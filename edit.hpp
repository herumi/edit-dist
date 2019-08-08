#pragma once
//#define MCLSHE_WIN_SIZE 10
#define MCLBN_FP_UNIT_SIZE 4
#define MCLBN_FR_UNIT_SIZE 4
#include <mcl/she.hpp>
#include <mcl/ecparam.hpp>
#include <cybozu/serializer.hpp>
#include <cybozu/time.hpp>
#include <fstream>
#include <mutex>
#include <atomic>

#define AAA

using namespace mcl::she;
using namespace mcl::bn;

typedef std::vector<int> IntVec;
typedef std::vector<CipherTextG1> CipherTextG1Vec;
typedef std::vector<CipherTextG1Vec> CipherTextG1VecVec;

namespace edit {

// fixed parameters
#ifdef AAA
const int diffMin = -5;
const int diffMax = 10;
const int diffNum = diffMax - diffMin + 1;
#else
const int diffMin = -1;
const int diffMax = 2;
const int diffNum = diffMax - diffMin + 1;
#endif

} // edit

struct CipherPack {
	int id;
	CipherTextG1 c[edit::diffNum];
};

struct Timer {
	const char *msg_;
	double begin_;
	void begin(const char *msg)
	{
		msg_ = msg;
		begin_ = cybozu::GetCurrentTimeSec();
	}
	void end() const
	{
		printf("%s %.2f sec\n", msg_, cybozu::GetCurrentTimeSec() - begin_);
	}
};

template<class V>
int editDist(const V& a, const V& b)
{
	size_t n = a.size();
	size_t m = b.size();
	std::vector<IntVec> dp(n + 1);
	for (size_t i = 0; i < n + 1; i++) {
		dp[i].resize(m + 1);
		dp[i][0] = int(i);
	}
	for (size_t j = 0; j < m + 1; j++) {
		dp[0][j] = int(j);
	}
	for (size_t i = 1; i < n + 1; i++) {
		for (size_t j = 1; j < m + 1; j++) {
			int D = dp[i - 1][j - 1] + (a[i - 1] == b[j - 1] ? 0 : 1);
			int U = dp[i][j - 1] + 1;
			int L = dp[i - 1][j] + 1;
			int t = (std::min)(U, L);
			dp[i][j] = (std::min)(t, D);
		}
	}
	return dp[n][m];
}

/*
	assume
	n = max - min + 1
	idxVec.size() == n
	out.size() == n
	cTbl[i - min] = Enc(i) for i in [min, max]
*/
template<class INT_VEC, class CIPHER>
void mixEnc(INT_VEC *idxVec, CipherPack& cp, const CipherTextG1& in, cybozu::RandomGenerator& rg, const CIPHER *cTbl)
{
	const int min = edit::diffMin;
	const int max = edit::diffMax;
	const int n = edit::diffNum;
	for (int i = 0; i < n; i++) {
		idxVec[i] = i;
	}
	cybozu::shuffle(idxVec, n, rg);
	CipherTextG1 v[n];
	mpz_class gamma;
	for (int i = min; i <= max; i++) {
		CipherTextG1 c;
		sub(c, in, cTbl[i - min]);
		mcl::gmp::getRand(gamma, 256);
		mul(v[i - min], c, gamma);
	}
	for (int i = 0; i < n; i++) {
		cp.c[idxVec[i]] = v[i];
	}
}

/*
	out[i] = min(out[i], bv[i]) for i in [0, m)
*/
template<class Stream, class CIPHER>
void serverMinVec(Stream& soc, CIPHER *out, const CIPHER *bv, int m, cybozu::RandomGenerator& rg, const CIPHER *cTbl)
{
	const int n = edit::diffNum;
	IntVec idxVec(m * n);

	cybozu::save(soc, m);
	std::mutex mw;
#pragma omp parallel for
	for (int i = 0; i < m; i++) {
		CipherPack cp;
		CipherTextG1 c;
		sub(c, out[i], bv[i]);
		mixEnc(&idxVec[i * n], cp, c, rg, cTbl);
		cp.id = i;

		{
			std::lock_guard<std::mutex> lk(mw);
			soc.write(&cp, sizeof(cp));
		}
	}

	/*
		min0(i) = min(i, 0) for i = -1, 0, 1, 2
		min(Enc(x), Enc(0)) = sum_i min0(i) Enc(delta_ix) = Enc(delta_(-1)x)
	*/
	for (int i = 0; i < m; i++) {
		CipherPack cp;
		soc.read(&cp, sizeof(cp));
		if (0 <= cp.id && cp.id < m) {
			sub(out[cp.id], bv[cp.id], cp.c[idxVec[cp.id * n]]);
		}
	}
}

// xv[] = min(0, xv[], yv[])
template<class Stream, class CIPHER>
void serverMinVec2(Stream& soc, CIPHER *xv, const CIPHER *yv, int m, cybozu::RandomGenerator& rg, const CIPHER *cTbl)
{
	const int n = edit::diffNum;
	IntVec idxVec(m * n);

	cybozu::save(soc, m);
	std::mutex mw;
#pragma omp parallel for
	for (int i = 0; i < m; i++) {
		CipherPack cp;
		CipherTextG1 c;
		add(c, xv[i], xv[i]);
		add(c, c, c); // 4xv[i]
		add(c, c, yv[i]); // 4xv[i] + yv[i] in [-5, 10]
		mixEnc(&idxVec[i * n], cp, c, rg, cTbl);
		cp.id = i;

		{
			std::lock_guard<std::mutex> lk(mw);
			soc.write(&cp, sizeof(cp));
		}
	}

	for (int i = 0; i < m; i++) {
		CipherPack cp;
		soc.read(&cp, sizeof(cp));
		if (0 <= cp.id && cp.id < m) {
			int idx = cp.id * n;
			CipherTextG1 c;
			add(c, cp.c[idxVec[idx + 0]], cp.c[idxVec[idx + 1]]);
			add(c, c, cp.c[idxVec[idx + 2]]);
			add(c, c, cp.c[idxVec[idx + 4]]);
			add(xv[cp.id], c, cp.c[idxVec[idx + 8]]);
		}
	}
}


/*
	BitEnc(x) := (Enc(delta_xi)) for i = 0, ..., n-1
	cv = [BitEnc(v[i]]
*/
static CipherTextG1 g_c[2];

void encIntVec(CipherTextG1Vec& cv, const IntVec& v, int charN, const PrecomputedPublicKey& ppub)
{
	cv.resize(charN * v.size());
	for (size_t j = 0; j < v.size(); j++) {
		int x = v[j];
		for (int i = 0; i < charN; i++) {
			ppub.enc(cv[j * charN + i], x == i ? 1 : 0);
		}
	}
}

void clientReEncPack(CipherPack *cp, const SecretKey& sec, const PrecomputedPublicKey& ppub)
{
	for (size_t i = 0; i < CYBOZU_NUM_OF_ARRAY(cp->c); i++) {
		bool isZero = sec.isZero(cp->c[i]);
		int v = isZero ? 1 :0;
#if 0
		// emulate precomputed version
		(void)ppub;
		cp->c[i] = g_c[v];
#else
		ppub.enc(cp->c[i], v);
#endif
	}
}
/*
	{Enc(c_i)} -> {Enc(1) if decoded else Enc(0)}
*/
template<class Stream>
bool clientReEnc(Stream& soc, const SecretKey& sec, const PrecomputedPublicKey& ppub)
{
	CipherPack cp;
	soc.read(&cp, sizeof(cp));
	if (cp.id < 0) return false;
	clientReEncPack(&cp, sec, ppub);
	soc.write(&cp, sizeof(cp));
	return true;
}

/*
	client has secret key
*/
template<class Stream>
bool clientProcess(Stream& soc, const SecretKey& sec, const IntVec& v, int charN)
{
	const int clientN = (int)v.size();
	PublicKey pub;
	sec.getPublicKey(pub);
	pub.save(soc);

	PrecomputedPublicKey ppub;
	ppub.init(pub);

	// precomputed
	ppub.enc(g_c[0], 0);
	ppub.enc(g_c[1], 1);

	CipherTextG1Vec cv;
	puts("prepare");
	Timer t;
	t.begin("prepare");
	encIntVec(cv, v, charN, ppub);
	t.end();
	cybozu::save(soc, v.size());
	cybozu::save(soc, charN);
	int serverN;
	cybozu::load(serverN, soc);
	printf("client serverN=%d, clientN=%d, charN=%d\n", serverN, clientN, charN);

	// send Enc(v)
	soc.write(cv.data(), sizeof(cv[0]) * cv.size());

	std::mutex mr, mw;
	for (;;) {
		int m;
		cybozu::load(m, soc);
		if (m == 0) break;
#pragma omp parallel for
		for (int i = 0; i < m; i++) {
			CipherPack cp;
			{
				std::lock_guard<std::mutex> lk(mr);
				soc.read(&cp, sizeof(cp));
			}
			clientReEncPack(&cp, sec, ppub);
			{
				std::lock_guard<std::mutex> lk(mw);
				soc.write(&cp, sizeof(cp));
			}
		}
	}

	CipherTextG1 c;
	c.load(soc);
	int dist = (int)sec.dec(c);
	printf("client result dist=%d\n", dist);

	IntVec vv(serverN);
	soc.read(vv.data(), sizeof(vv[0]) * vv.size());
	int okDist = editDist(v, vv);
	if (dist != okDist) {
		printf("err!!! %d %d\n", dist, okDist);
	}
	return dist == okDist;
}

template<class Stream>
void serverProcess(Stream& soc, const IntVec& v)
{
	const int serverN = (int)v.size();
	PublicKey pub;
	pub.load(soc);

	PrecomputedPublicKey ppub;
	ppub.init(pub);

	cybozu::RandomGenerator rg;

	CipherTextG1 one;
	ppub.enc(one, 1);

	SecretKey sec;
	int clientN;
	int charN;
	cybozu::load(clientN, soc);
	cybozu::load(charN, soc);
	cybozu::save(soc, serverN);
	printf("server serverN=%d, clientN=%d, charN=%d\n", serverN, clientN, charN);

	/*
		csMat[i][j] = client[j] == server[i] for j in [0, clientN), i in [0, serverN)
	*/
	CipherTextG1VecVec csMat;
	{
		CipherTextG1Vec cv(charN * clientN);
		soc.read(cv.data(), sizeof(cv[0]) * cv.size());
		csMat.resize(serverN);
		for (size_t i = 0; i < csMat.size(); i++) {
			csMat[i].resize(clientN);
			int x = v[i];
			if (v[i] >= charN) throw cybozu::Exception("serverProcess:bad v") << i << v[i] << charN;
			for (int j = 0; j < clientN; j++) {
				csMat[i][j] = cv[charN * j + x];
			}
		}
	}

	// init dp[serverN + 1][clientN + 1]
	CipherTextG1VecVec dp(serverN + 1);
	for (int i = 0; i <= serverN; i++) {
		dp[i].resize(clientN + 1);
		ppub.enc(dp[i][0], int(i));
	}
	for (int j = 1; j <= clientN; j++) {
		ppub.enc(dp[0][j], int(j));
	}

	printf("server");

	/*

		D = dp[i-1][j-1] + 1 - (a[i] == b[j]) | U = dp[i-1][j] + 1
		L = dp[i][j - 1] + 1                  | dp[i][j] = min(D, U, L)

		|
		v

		D = dp[i-1][j-1] - (a[i] == b[j]) | U = dp[i-1][j]
		L = dp[i][j - 1]                  | dp[i][j] = min(D, U, L) + 1
	*/
	/*
		U - D =  -1 0 1 2
		min(D, U) - L =  -1 0 1 2
	*/

	CipherTextG1 cTbl[edit::diffNum];
	for (int i = 0; i < edit::diffNum; i++) {
		ppub.enc(cTbl[i], i + edit::diffMin);
	}
	const int minN = (std::min)(clientN, serverN);
	const int maxN = (std::max)(clientN, serverN);
	CipherTextG1Vec av(minN);
	CipherTextG1Vec bv(minN);
#ifdef AAA
	CipherTextG1Vec cv(minN);
#endif
	for (int k = 1; k < minN + maxN; k++) {
		printf("%d ", k);
		const int beginI = (std::min)(k, serverN);
		const int beginJ = 1 + k - beginI;
		const int endJ = (std::min)(k, clientN);
		const int m = endJ - beginJ + 1;
#ifdef AAA
		for (int s = 0; s < m; s++) {
			int i = beginI - s;
			int j = beginJ + s;
			sub(av[s], dp[i - 1][j - 1], csMat[i - 1][j - 1]); // D
			sub(bv[s], dp[i - 1][j], av[s]); // U - D
			sub(cv[s], dp[i][j - 1], av[s]); // L - D
		}
		serverMinVec2(soc, bv.data(), cv.data(), m, rg, cTbl);
		// bv[] = min(0, U - D, L - D)
		for (int s = 0; s < m; s++) {
			int i = beginI - s;
			int j = beginJ + s;
			add(dp[i][j], av[s], one);
			add(dp[i][j], dp[i][j], bv[s]); // 1 + min(D, U, L)
		}
#else
		for (int s = 0; s < m; s++) {
			int i = beginI - s;
			int j = beginJ + s;
			sub(av[s], dp[i - 1][j - 1], csMat[i - 1][j - 1]);
			bv[s] = dp[i - 1][j];
		}
		serverMinVec(soc, bv.data(), av.data(), m, rg, cTbl);
		for (int s = 0; s < m; s++) {
			int i = beginI - s;
			int j = beginJ + s;
			av[s] = dp[i][j - 1];
		}
		serverMinVec(soc, av.data(), bv.data(), m, rg, cTbl);
		for (int s = 0; s < m; s++) {
			int i = beginI - s;
			int j = beginJ + s;
			add(dp[i][j], av[s], one);
		}
#endif
	}
	printf("\n");
	cybozu::save(soc, 0); // finish
	// to client
	dp[serverN][clientN].save(soc);

	soc.write(v.data(), sizeof(v[0]) * v.size());
}

/*
	conver ascii code [a-z] to int
	n ; max characters
*/
void convertStringToIntVec(IntVec& v, const std::string& s, int n)
{
	v.resize(s.size());
	for (size_t i = 0; i < v.size(); i++) {
		int x = s[i];
		if (x < 0 || x >= n) throw cybozu::Exception("convertStringToIntVec:bad s") << s;
		v[i] = x;
	}
}
