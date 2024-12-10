#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <sys/time.h>
#include </home/simon/dlogu32.c>

// gcc dlog.c -o dlog.bin -lm -O3 -march=native -Wall -std=c11

uint64_t gettimeus(void) {
  struct timeval thistime;
  gettimeofday(&thistime, NULL);
  return thistime.tv_sec*1000000ULL + thistime.tv_usec;
}

int main(int argc, char* argv[]) {
  uint32_t a = 101;
  uint32_t b = 654321;
  uint32_t n = 100000007;
  uint32_t e = 0;
  _Bool success;
  if (argc >= 4) {
    a = atol(argv[1]);
    b = atol(argv[2]);
    n = atol(argv[3]);
  } 
  if ((argc < 4) || (n == 0)) {
    printf("This program calculates discrete logarithms deterministically. Given a, b, n, it returns the smallest non-negative e such that a^e = b mod n, or an error message if it is impossible. \n0 <= a,b < n < 2^32\nUsage:- %s a b n\n", argv[0]);
    exit(0);
  }
  a %= n;
  b %= n;
  uint64_t starttime, endtime;
  starttime = gettimeus();
  success = dlogu32(a, b, n, &e);
  endtime = gettimeus();
  if (!success) {
    printf("%u^e = %u mod %u\n", a,b,n); 
    printf("No Solutions!\n");
  } else {
    printf("%u^%u = %u mod %u\n", a,e,b,n); 
  }
  printf("Time = %lu us\n", endtime - starttime);
}
