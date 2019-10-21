# Editing Distance by Re-Encryption Protocol

# Support architecture

* x86-64 Windows + Visual Studio
* x86, x86-64 Linux/macOS + gcc/clang

# Installation Requirements

* [GMP](https://gmplib.org/) and OpenSSL
```
apt install libgmp-dev libssl-dev
```

Create a working directory (e.g., work) and clone the following repositories.
```
mkdir work
cd work
git clone git://github.com/herumi/mcl
git clone git://github.com/herumi/edit-dist
cd edit-dist
make
```
# How to use

```
usage:edit [opt]
  -h show this message
  -ip : ip address
  -p : port
  -s : string
  -n : # of kind of characters
  -slen : string length
  -save-sec : save secretKey
  -new : new secretKey
  -d : debug
  -delay : tcp delay
  -swap : swap role
```

sample test
- server
  - ./edit -slen 128
- client
  - ./edit -ip <server> -slen 128

# How to use VTune
```
cd mcl
make clean
cd ../edit-dist
make MCL_USE_PROF=2
./edit ...
```

# Benchmark
Xeon SP(Platinam 8280)
- env OMP_NUM_THREADS=56 numactl --cpunodebind=0 --membind=0 ./edit -slen 1024
- env OMP_NUM_THREADS=40 numactl --cpunodebind=1 --membind=1 ./edit -ip <server> -slen 1024

total 43 sec

# Reference
[Arbitrary Univariate Function Evaluation and Re-Encryption Protocols over Lifted-ElGamal Type Ciphertexts](https://eprint.iacr.org/2019/1233), Koji Nuida and Satsuya Ohata and Shigeo Mitsunari and Nuttapong Attrapadung


# Acknowledgements
Thanks to Masashi Horikoshi at Intel to provide Xeon SP environments to evaluate this protocol.

# Author

光成滋生 MITSUNARI Shigeo(herumi@nifty.com)
