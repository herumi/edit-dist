/*
	normal dp
	n =   10  20  30
	sec   12  50 115
	diagonal parallel
	n =   20  30    40    80   160  512 1024
	sec 5.52 8.7 12.61 25.86 41.65  222  720
*/
#include "edit.hpp"
#include <mcl/ecparam.hpp>
#include <cybozu/socket.hpp>
#include <cybozu/serializer.hpp>
#include <cybozu/option.hpp>

const char *g_secretKeyName = "secretkey.txt";

int main(int argc, char *argv[])
	try
{
	cybozu::Option opt;
	int port;
	bool saveSecretKey;
	bool newSecretKey;
	bool swapRole;
	std::string ip;
	std::string s;
	int slen;
	bool debug;
	int n;
	opt.appendOpt(&ip, "", "ip", ": ip address");
	opt.appendOpt(&port, 10000, "p", ": port");
	opt.appendOpt(&s, "", "s", ": string");
	opt.appendOpt(&n, 128, "n", ": # of kind of characters");
	opt.appendOpt(&slen, 0, "slen", ": string length");
	opt.appendBoolOpt(&saveSecretKey, "save-sec", ": save secretKey");
	opt.appendBoolOpt(&newSecretKey, "new", ": new secretKey");
	opt.appendBoolOpt(&debug, "d", ": debug");
	opt.appendBoolOpt(&swapRole, "swap", ": swap role");
	opt.appendHelp("h", "show this message");
	if (!opt.parse(argc, argv)) {
		opt.usage();
		return 1;
	}
	const size_t tryNum = 1;
	mcl::she::initG1only(mcl::ecparam::secp256k1, 2048, tryNum);

	if (saveSecretKey) {
		SecretKey sec;
		sec.setByCSPRNG();
		std::ofstream ofs(g_secretKeyName, std::ios::binary);
		sec.save(ofs);
		return 0;
	}
	if (slen) {
		s.resize(slen);
		cybozu::RandomGenerator rg;
		for (int i = 0; i < slen; i++) {
			s[i] = 'a' + (rg() % 27);
		}
	}

	IntVec v;
	convertStringToIntVec(v, s, n);
	for (size_t i = 0; i < v.size(); i++) {
		printf("%d ", v[i]);
	}
	printf("\n");

	SecretKey sec;
	if (newSecretKey) {
		sec.setByCSPRNG();
	} else {
		std::ifstream ifs(g_secretKeyName, std::ios::binary);
		sec.load(ifs);
	}

	if (ip.empty()) {
		printf("server port=%d\n", port);
		cybozu::Socket server;
		server.bind(uint16_t(port));
		for (;;) {
			while (!server.queryAccept()) {
			}
			cybozu::Socket client;
			server.accept(client);
			client.setSocketOption(TCP_NODELAY, 1, IPPROTO_TCP);
#if 0
			serverProcess(client, v);
#else
			char c;
			client.read(&c, 1);
			printf("c=%c\n", c);
			if (c == '0') {
				serverProcess(client, v);
			} else {
				clientProcess(client, sec, v, n);
			}
#endif
		}
	} else {
		printf("client ip=%s port=%d\n", ip.c_str(), port);
		cybozu::SocketAddr sa(ip, uint16_t(port));
		printf("addr=%s\n", sa.toStr().c_str());
		cybozu::Socket client;
		client.connect(sa);
		Timer t;
		t.begin("total");
		client.setSocketOption(TCP_NODELAY, 1, IPPROTO_TCP);
#if 0
		clientProcess(client, sec, v, n);
#else
		client.write(swapRole ? "1" : "0", 1);
		if (swapRole) {
			serverProcess(client, v);
		} else {
			clientProcess(client, sec, v, n);
		}
#endif
		t.end();
	}
} catch (std::exception& e) {
	printf("err %s\n", e.what());
	return 1;
}
