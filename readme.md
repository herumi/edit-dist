# editing distance by reenc protocol

# Support architecture

* x86-64 Windows + Visual Studio
* x86, x86-64 Linux + gcc/clang

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
  - ./edit -ip localhost -slen 128

# Benchmark
Xeon SP(Platinam 8280)
- env OMP_NUM_THREADS=56 ./edit -slen 1024
- env OMP_NUM_THREADS=40 ./edit -ip localhost -slen 1024

total 44.72 sec

# Author

光成滋生 MITSUNARI Shigeo(herumi@nifty.com)
