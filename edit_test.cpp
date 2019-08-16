#include "edit.hpp"
#define CYBOZU_TEST_DISABLE_AUTO_RUN
#include <cybozu/test.hpp>
#include <mutex>
#include <condition_variable>
#include <thread>
#include <cybozu/xorshift.hpp>
#include <cybozu/option.hpp>
#include <cybozu/benchmark.hpp>

/*
	len / sec
	  20   40    80   160
	1.15 4.32 17.02 67.40 ; single thread
	0.78 2.81 11.11 42.50 ; OpenMP
	diagonal parallel
	 32   64   128 256
	2.5 9.69 38.06     w/o OpenMP
	       3    11  42 w OpenMP
	edit_test
	len = 64
	7.31

	S = {-2, -1, 0, 1, 2}
	len =
	8    16   32   64
	0.2  0.57 1.98 7.57 ; diffMin = -2, diffMax = 2
	0.18 0.47 1.62 6.10 ; diffMin = -2, diffMax = 1
*/

int g_n;
SecretKey g_sec;

struct MemoryStream {
	std::string s;
	std::mutex m_;
	std::condition_variable cv_;
	void write(const void *p, size_t n)
	{
		std::lock_guard<std::mutex> lk(m_);
		bool isEmpty = s.empty();
		s.append((const char*)p, n);
		if (isEmpty) cv_.notify_one();
	}
	size_t readSome(void *buf, size_t n)
	{
		std::unique_lock<std::mutex> lk(m_);
		cv_.wait(lk, [&] { return !s.empty(); });
		if (n > s.size()) {
			n = s.size();
		}
		memcpy(buf, s.data(), n);
		s = s.substr(n);
		return n;
	}
};

struct Pipe {
	MemoryStream *in;
	MemoryStream *out;
	Pipe(MemoryStream *in, MemoryStream *out)
		: in(in)
		, out(out)
	{
	}
	void write(bool *pb, const void *p, size_t n)
	{
		out->write(p, n);
		*pb = true;
	}
	void write(const void *p, size_t n)
	{
		out->write(p, n);
	}
	size_t readSome(void *buf, size_t n)
	{
		return in->readSome(buf, n);
	}
	void read(void *buf, size_t n)
	{
		char *p = (char *)buf;
		while (n > 0) {
			size_t readSize = readSome(p, n);
			if (readSize == 0) throw cybozu::Exception("MemoryStream:read:readSize is zero");
			p += readSize;
			n -= readSize;
		}
	}
};

CYBOZU_TEST_AUTO(init)
{
}

void setRand(std::string& s, cybozu::XorShift& rg, size_t n)
{
	s.resize(n);
	for (size_t i = 0; i < n; i++) {
		s[i] = 'a' + (rg() % 27);
	}
}

CYBOZU_TEST_AUTO(test)
{
	std::string s1;
	std::string s2;
	IntVec v1, v2;
	int n = 127;

	cybozu::XorShift rg;
	int len = 8;
	for (int i = 0; i < g_n; i++) {
		printf("len=%d\n", len);
		setRand(s1, rg, len);
		setRand(s2, rg, len);
		len *= 2;
		convertStringToIntVec(v1, s1, n);
		convertStringToIntVec(v2, s2, n);

		MemoryStream ms1, ms2;
		Pipe p1(&ms1, &ms2);
		Pipe p2(&ms2, &ms1);
		Timer tm;
		tm.begin("start");
		bool ok;
		std::thread t1([&]{ ok = clientProcess(p1, g_sec, v1, n); });
		std::thread t2([&]{ serverProcess(p2, v2); });
		t1.join();
		t2.join();
		tm.end();
		CYBOZU_TEST_ASSERT(ok);
	}
}

CYBOZU_TEST_AUTO(bench)
{
	PublicKey pub;
	g_sec.getPublicKey(pub);
	PrecomputedPublicKey ppub;
	ppub.init(pub);

	Fr gamma;
	gamma.setByCSPRNG();
	CipherTextG1 c1;
	puts("G1");
	CYBOZU_BENCH_C("enc", 1000, ppub.enc, c1, 123);
	CYBOZU_BENCH_C("isZero", 1000, g_sec.isZero, c1);
	CYBOZU_BENCH_C("add", 1000, CipherTextG1::add, c1, c1, c1);
	CYBOZU_BENCH_C("mul", 1000, CipherTextG1::mul, c1, c1, gamma);
	char buf[256];
	size_t n = c1.serialize(buf, sizeof(buf));
	printf("n=%zd\n", n);
	CYBOZU_BENCH_C("seri", 1000, c1.serialize, buf, sizeof(buf));
	n = c1.deserialize(buf, n);
	CYBOZU_BENCH_C("deseri", 1000, c1.deserialize, buf, sizeof(buf));
	{
		CipherPack cp;
		CipherTextG1 c;
		CipherTextG1 cTbl[edit::diffNum];
		cybozu::RandomGenerator rg;
		for (int i = 0; i < edit::diffNum; i++) {
			ppub.enc(cTbl[i], i + edit::diffMin);
		}
		IntVec idxVec(4);
		ppub.enc(c, 3);
		CYBOZU_BENCH_C("mixEnc", 1000, mixEnc, &idxVec[0], cp, c, rg, cTbl);
		const size_t n = 4;
		mpz_class rVec[n];
		for (size_t i = 0; i < n; i++) {
			mcl::gmp::getRand(rVec[i], 256);
		}
		G1 vec[n];
		CYBOZU_BENCH_C("winP.mul", 1000, edit::g_winP.mul, vec[0], rVec[0]);
		CYBOZU_BENCH_C("multiMul4", 10000, (multiMul<n, 4>), vec, c.getS(), rVec);
		CYBOZU_BENCH_C("multiMul5", 10000, (multiMul<n, 5>), vec, c.getS(), rVec);
		CYBOZU_BENCH_C("multiMul6", 10000, (multiMul<n, 6>), vec, c.getS(), rVec);
		CYBOZU_BENCH_C("multiMul7", 10000, (multiMul<n, 7>), vec, c.getS(), rVec);
		CYBOZU_BENCH_C("multiMul8", 10000, (multiMul<n, 8>), vec, c.getS(), rVec);
	}
#if 0
	puts("G2");
	CipherTextG2 c2;
	CYBOZU_BENCH_C("enc", 1000, ppub.enc, c2, 123);
	CYBOZU_BENCH_C("isZero", 1000, g_sec.isZero, c2);
	CYBOZU_BENCH_C("add", 1000, CipherTextG2::add, c2, c2, c2);
	CYBOZU_BENCH_C("mul", 1000, CipherTextG2::mul, c2, c2, gamma);
	CYBOZU_BENCH_C("mul160", 1000, CipherTextG2::mul, c2, c2, r);
	CipherTextGT ct;
	CYBOZU_BENCH_C("mulG1xG2", 1000, CipherTextGT::mul, ct, c1, c2);
#endif
}

int main(int argc, char *argv[])
	try
{
	cybozu::Option opt;
	opt.appendOpt(&g_n, 6, "n", ": num of loop");
	opt.appendHelp("h", ": show this message");
	if (!opt.parse(argc, argv)) {
		opt.usage();
		return 1;
	}

	edit::init();
	g_sec.setByCSPRNG();

	return cybozu::test::autoRun.run(argc, argv);
} catch (std::exception& e) {
	printf("ERR %s\n", e.what());
	return 1;
}
