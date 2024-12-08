// Flags to use with gcc are -lm -std=c11
#include <math.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>

#ifndef alignof
  #define alignof _Alignof
#endif

#define DLOGU32MMDL_APOWERS_SIZE 131101
#define DLOGU32MMDL_NEXT(i) i = ((i + 1) < DLOGU32MMDL_APOWERS_SIZE ? i + 1 : 0)

void gcdu32(uint32_t A, uint32_t B, uint32_t *C) {
  uint32_t temp;
  if ((A == 0) || (B == 0)) {
    *C = 0;
    return;
  }
  while (1) {
    if (A < B) {
      temp = A;
      A = B;
      B = temp;
    }
    A = A % B;
    if (A == 0) {
      *C = B;
      return;
    }
  }
}

void egcdu32(uint32_t A, uint32_t B, int64_t *S, int64_t *T, uint32_t *C) {
  if ((A == 0) || (B == 0)) {
    *C = 0;
    return;
  }
  int64_t Ai, Bi, Qi, Siprev, Tiprev, Si, Ti, temp;  
  _Bool swapAB = false;
  Ai = A;
  Bi = B;
  if (Ai < Bi) {
    temp = Ai;
    Ai = Bi;
    Bi = temp;
    swapAB = true;
  }
  Si = 0;
  Ti = 1;
  Siprev = 1;
  Tiprev = 0;
  while (1) {
    if (Ai < Bi) {
      temp = Ai;
      Ai = Bi;
      Bi = temp;
    }
    Qi = Ai / Bi;
    Ai = Ai - (Qi * Bi);
    temp = Siprev - Qi*Si;
    Siprev = Si;
    Si = temp;
    temp = Tiprev - Qi*Ti;
    Tiprev = Ti;
    Ti = temp;
    if (Ai == 0) {
      *C = Bi;
      if (swapAB) {
        *T = Siprev;
        *S = Tiprev;
      } else {
        *S = Siprev;
        *T = Tiprev;
      }
      return;
    }
  }
}

_Bool modinvu32(uint32_t A, uint32_t N, uint32_t *C, uint32_t *gcd) {
  // C = A**(-1) mod N
  // A and N must be both +ve and relatively prime.
  if (N <= 1) return false;
  A %= N;
  if (A == 0) return false;
  int64_t S, T;
  egcdu32(A, N, &S, &T, gcd);
  if (*gcd != 1) return false;
  if (S >= 0) {
    *C = (S % N);
  } else {
    S *= -1;
    *C = N - (S % N);
  }
  return true;
}

uint32_t modpowu32(uint32_t a, uint32_t e, uint32_t n) {
// Returns a^e mod n
  if (n == 0) return 0;
  if (a < 2) return a;
  uint32_t res = 1;
  uint32_t sq = a % n;
  while (e) {
    if (e & 1U) res = ((uint64_t)res * sq) % n;
    sq = ((uint64_t)sq*sq) % n;
    e >>= 1;
  }
  return res;
}

uint32_t isqrt(uint64_t n) {
  if (n < 2) return n;
  uint64_t ai = sqrt(n);
  while (!((ai <= n/ai) && ((ai+1) > n/(ai+1)))) {    
    ai = (ai + n/ai)/2;
  }
  return ai;
}


_Bool dlogu32naive(uint32_t a, uint32_t b, uint32_t n, uint32_t count, uint32_t *e) {
  uint32_t i = 0;
  uint32_t apow = 1;
  for (; i<count; i++) {
    if (b == apow) {
      *e = i;
      return true;
    }
    apow = ((uint64_t)apow*a) % n;
  }
  return false;
}

typedef struct {
  uint32_t e, atoemodn;
} dlogu32mmdl_t;

_Bool dlogu32mmdl_insert(dlogu32mmdl_t el, dlogu32mmdl_t *apowers) {  
  uint32_t ix = el.atoemodn % DLOGU32MMDL_APOWERS_SIZE;
  while (apowers[ix].atoemodn != 0) DLOGU32MMDL_NEXT(ix);
  apowers[ix] = el;
  return true;  
}

_Bool dlogu32mmdl_select(uint32_t atoemodn, uint32_t *ix, dlogu32mmdl_t *apowers) {
  uint32_t i = atoemodn % DLOGU32MMDL_APOWERS_SIZE;
  while (apowers[i].atoemodn && (apowers[i].atoemodn != atoemodn)) DLOGU32MMDL_NEXT(i);
  if (apowers[i].atoemodn) {
    *ix = i;
    return true;  
  }
  return false;
}

_Bool dlogu32mmdl(uint32_t a, uint32_t b, uint32_t n, uint32_t *e) {
  // Returns an e if one exists, such that a^e = b mod n using a meet-in-the-middle algorithm. (n and a must be co-prime)  
  // The returned e is not necessarily the smallest. 
  // Take its residue modulo the multiplicative order of a mod n for the smallest possible e.
  // Requires ~1MiB RAM.
  if (n == 0) return false;
  if (a >= n) a %= n;
  if (a == 0) return false;
  if (b >= n) b %= n;
  if (b == 0) return false;
  if (n < 100000u) return dlogu32naive(a, b, n, n, e);
  if (dlogu32naive(a, b, n, 100, e)) return true;
  //uint32_t ainv = modpowu32(a, n-2, n); // Assumes n is prime.
  uint32_t ainv, gcd;
  if (!modinvu32(a, n, &ainv, &gcd)) return false;
  uint32_t m = isqrt(n);
  m += (n > m*m);
  uint32_t atominusm = modpowu32(ainv, m, n);
  uint32_t ix;
  dlogu32mmdl_t *apowers = aligned_alloc(alignof(dlogu32mmdl_t), DLOGU32MMDL_APOWERS_SIZE * sizeof(dlogu32mmdl_t));
  if (apowers == NULL) return false;   
  memset(apowers, 0, DLOGU32MMDL_APOWERS_SIZE * sizeof(dlogu32mmdl_t));
  dlogu32mmdl_t newel;
  newel.e = 0;
  newel.atoemodn = 1;
  do {
    dlogu32mmdl_insert(newel, apowers);
    newel.e++;
    newel.atoemodn = ((uint64_t)newel.atoemodn * a) % n;
  } while ((newel.atoemodn != 1) && (newel.e < m));
  uint32_t i = 0;
  uint32_t bi = b;
  uint32_t atominusmpowers = 1;
  while (i < m) {
    if (dlogu32mmdl_select(bi, &ix, apowers)) break;
    atominusmpowers = ((uint64_t)atominusmpowers*atominusm) % n;
    bi = ((uint64_t)b*atominusmpowers) % n;
    i++;
  }
  if (i >= m) {
    free(apowers);
    return false;
  }
  *e = i*m + apowers[ix].e;
  free(apowers);
  return true;
}

_Bool dlogu32(uint32_t a, uint32_t b, uint32_t n, uint32_t *e) {
  // Sets e such that a^e = b mod n if one exists and returns true, returns false otherwise.
  // If o is the multiplicative order of a mod n, the smallest possible solution is e mod o.
  if (n == 0) return false;
  uint32_t gcd;
  gcdu32(a, n, &gcd);
  if (gcd > 1) {
    if ((b % gcd) != 0) return false;
    uint32_t k = gcd;
    // (sk)^e = tk mod uk
    // (sk)^(e-1) = t/s mod u
    uint32_t invs = 1;
    uint32_t u = n/k;
    if(!modinvu32(a/k, u, &invs, &gcd)) return false;
    uint32_t newb = ((uint64_t)(b/k)*invs) % u;
    if (!dlogu32(a,newb,u,e)) return false;
    (*e)++;
    return b == modpowu32(a,*e,n);
  } else {
    return dlogu32mmdl(a,b,n,e);
  }
}
