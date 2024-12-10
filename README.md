# Discrete-Log
Deterministic Calculation Of Discrete Logarithms.

This program calculates discrete logarithms deterministically. 

Given a, b, n, it returns an e such that a^e = b mod n, or an error message if it is impossible. 

There may be a smaller positive e which also works. 

The smallest non-negative solution is e mod o where o is the multiplicative order of a mod n.

0 < a,b,n < 2^32

Usage:- ./dlog.bin a b n

E.g. ./dlog.bin 3 19 2800000051

3^2055010318 = 19 mod 2800000051

Time = 4502 us


./dlog.bin 3 19 4000000001

3^e = 19 mod 4000000001

No Solutions!

Time = 6042 us

[CPU i7-6700 @3.4GHz]
