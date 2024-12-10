# Discrete-Log
Deterministic Calculation Of Discrete Logarithms.


dlog.c
======

This program calculates discrete logarithms deterministically. 

Given a, b, n, it returns the smallest e such that a^e = b mod n, or an error message if it is impossible. 

0 <= a,b < n < 2^32

Usage:- ./dlog.bin a b n

E.g. ./dlog.bin 3 19 2800000051

3^2055010318 = 19 mod 2800000051

Time = 4502 us


./dlog.bin 3 19 4000000001

3^e = 19 mod 4000000001

No Solutions!

Time = 6042 us


./dlog.bin 35 6171875 10000000

35^17 = 6171875 mod 10000000

Time = 1 us

[CPU i7-6700 @3.4GHz]
