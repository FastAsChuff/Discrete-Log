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
  if (n == 0) return false;
  if (a >= n) a %= n;
  if (b >= n) b %= n;
  *e = 1u;
  if ((a == 0) && (b == 0)) return true;
  *e = 0;
  if (a < 2u) {
    if (b == a) return true;
    return false;
  }
  uint32_t gcd;
  if (b == 1) {
    gcdu32(a, n, &gcd);
    if (gcd != 1) return false;
    return true;
  }
  uint32_t i = 1;
  uint32_t apow = a;
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
  // The returned e is the smallest non-negative solution. 
  // Requires ~1MiB RAM.
  if (n < 100000u) return dlogu32naive(a,b,n,n,e);
  if (n == 0) return false;
  if (a >= n) a %= n;
  if (b >= n) b %= n;
  *e = 0;
  if (a < 2u) {
    if (b == a) return true;
    return false;
  }
  uint32_t m = isqrt(n);
  m += (n > m*m);
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
  } while ((newel.atoemodn != b) && (newel.atoemodn != 1) && (newel.e < m));
  if (newel.atoemodn == b) {
    *e = newel.e;
    free(apowers);
    return true;
  }
  if (newel.atoemodn == 1) {
    free(apowers);
    return false;
  }
  // a(im+j) = b mod n => b(a^(-m))^i = a^j mod n
  uint32_t ainv, gcd;
  if (!modinvu32(a, n, &ainv, &gcd)) {
    free(apowers);
    return false;
  }
  uint32_t i = 1;
  uint32_t atominusm = modpowu32(ainv, m, n);
  uint32_t bi = ((uint64_t)b*atominusm) % n;
  uint32_t ix;
  while (i < m) {
    if (dlogu32mmdl_select(bi, &ix, apowers)) break;
    bi = ((uint64_t)bi*atominusm) % n;
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

uint32_t largestcoprimefactoru32(uint32_t n, uint32_t k) {
  // Returns the largest factor of n that is co-prime to k.
  uint32_t gcd;
  do {
    gcdu32(k, n, &gcd);
    n /= gcd;
  } while (gcd > 1);
  return n;
}

_Bool dlogu32(uint32_t a, uint32_t b, uint32_t n, uint32_t *e) {
  // Sets e to smallest non-negative solution such that a^e = b mod n and returns true if a solution exists.
  // Returns false otherwise.
  if (n == 0) return false;
  if (a >= n) a %= n;
  if (b >= n) b %= n;
  *e = 1u;
  if ((a == 0) && (b == 0)) return true;
  *e = 0;
  if (a < 2u) {
    if (b == a) return true;
    return false;
  }
  uint32_t gcd;
  gcdu32(a, n, &gcd);
  if (b == 1) {
    if (gcd != 1) return false;
    return true;
  }
  if (gcd > 1) {
    uint32_t v = largestcoprimefactoru32(n, gcd);
    uint32_t u = n/v;
    uint32_t bmodu = b % u;
    uint32_t mine = 1u;
    uint32_t apow = a % u;
    while ((apow != 0) && (apow != bmodu)) {
      mine++;
      apow = ((uint64_t)apow*a) % u;
    }
    *e = mine;
    if (v == 1) return b == modpowu32(a,*e,n);
    if (bmodu != 0) return b == modpowu32(a,*e,n);
    if(!dlogu32mmdl(a,b,v,e)) return false;
    if (mine > *e) {
      uint32_t invamodv, ov;
      if(!modinvu32(a, v, &invamodv, &gcd)) return false;
      if(!dlogu32mmdl(a,invamodv,v,&ov)) return false;
      ov++; // ov is now the multiplicative order of a mod v
      *e += ov*(((mine + (uint64_t)ov) - (1ULL + *e))/ov); // make >= mine by adding multiples of ov.
    }
    return b == modpowu32(a,*e,n);
  } else {
    return dlogu32mmdl(a,b,n,e);
  }
}
