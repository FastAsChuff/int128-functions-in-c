#define INT128_LOW 0
#define INT128_HIGH 1
#define INT128_PI_HIGH 0x3243f6a8885a308
#define INT128_PI_LOW 0xd313198a2e037073

typedef union {
  __int128 i128;
  uint64_t hilo[2];
} myint128_t;

typedef struct {
  myint128_t modulus; // modulus < 0 is invalid.
  myint128_t phi;
  uint8_t primecount; // primecount = 0 => (modulus < 2 || modulus not factorized)
                      // primecount > 0 => (modulus is factorized)
  uint32_t SQRT1count;
  uint8_t *S; // greatest exponent s.t. 2^S[i] | (P[i] - 1)
  uint8_t *primepower;  // greatest exponent s.t. P[i]^primepower[i] | modulus
  myint128_t X;  // general result mod modulus.
  myint128_t *Ppow;  // prime power Ppow[i] = P[i]^primepower[i]
  myint128_t *P;  // prime. If 2|modulus, P[0] = 2.
  myint128_t *Q;  // Q[i] = (P[i] - 1)/(2^S[i])
  myint128_t *z;  // some non-zero non-square mod P[i]
  myint128_t *SQRT1;  // SQRT1count square roots of 1 mod modulus 
} modulus128_t;

uint32_t lzcnt128(myint128_t n) {
  // WARNING!!! Don't call __builtin_clzll(0)
  if (n.hilo[INT128_HIGH] == 0) {
    if (n.hilo[INT128_LOW] == 0) return 128;
    return 64 + __builtin_clzll(n.hilo[INT128_LOW]);
  }
  return __builtin_clzll(n.hilo[INT128_HIGH]);
}

myint128_t nchoosek_simple(uint64_t n, uint64_t k) {
  // OK for all k if n <= 122
  // Returns -1 if overflow protection triggered.
  myint128_t res;
  uint64_t i, clz;
  uint64_t temp;
  res.i128 = 0;
  if (k > n) return res;
  if (n-k < k) k = n-k;
  res.i128 = 1;
  if (k == 0) return res;
  res.i128 = n;
  for (i=2; i<=k; i++) {    
    temp = n-i+1;
    clz = __builtin_clzll(temp);
    if (res.i128 >> (62 + clz)) {
      res.i128 = -1;
      return res;
    } else {
      res.i128 = (res.i128 * temp)/i;
    }
  }
  return res;
}
uint64_t nchoosek(uint64_t n, uint64_t k) {
  // No overflow for n <= 66.
  // Initialize with  nchoosek(0, 0);
  myint128_t res;
  static _Bool isinerror;
  if (n == 0) {
    isinerror = false;
    return 0;
  }
  if (k > n) return 0;
  if (k > n - k) k = n - k;
  if (k == 0) return 1;
  uint64_t nm1ckm1 = nchoosek(n-1,k-1);
  uint64_t nm1ckm1_lz = 64;
  if (nm1ckm1 > 0) {
    nm1ckm1_lz = __builtin_clzll(nm1ckm1);
    //nm1ckm1_lz = lzcnt(nm1ckm1);
  } else {
    return 0;
  }
  uint64_t n_lz = 64;
  if (n > 0) n_lz = __builtin_clzll(n);
  //if (n > 0) n_lz = lzcnt(n);
  if ((128 - nm1ckm1_lz - n_lz) > 127) {
    if (!isinerror) {
      fprintf(stderr, "Overflow in nchoosek().\n");
      isinerror = true;
    }    
    return 0;
  }
  res.i128 = (nm1ckm1 * (__int128)n) / k;
  if (lzcnt128(res) >= 65) {
    return res.i128;
  } else {
    return 0;
  }
}


myint128_t nchoosek128(uint64_t n, int64_t k) {
  // Ok for all n, k where n <= 128.
  // Returns -1 if overflow protection triggered.
  myint128_t res3, res2, res;
  uint32_t choosemax;
  int32_t i, kstart;
  choosemax = 120;
  if ((k < 0) || (k > n))   {
    res.i128 = 0;
    return res;
  }
  if (k > n - k) k = n - k;
  if (k == 0) {
    res.i128 = 1;
    return res;
  }
  if (k == 1) {
    res.i128 = n;
    return res;
  }
  if (n <= choosemax) {
    res2 = nchoosek128(n-1,k-1);
    res.i128 = (res2.i128 * n) / k;
    return res;
  }
  
  kstart = k - (n - choosemax);
  res.i128 = 0;
  for (i = kstart; i <= k; i++) {
    if ((128 - lzcnt128(res)) > 126) {
      res.i128 = -1;
      return res;
    }
    res2 = nchoosek128(choosemax, i);
    if (res2.i128 < 0) {
      res.i128 = -1;
      return res;
    }
    res3 = nchoosek128(n-choosemax, i-kstart);
    if (res3.i128 < 0) {
      res.i128 = -1;
      return res;
    }
    if ((256 - lzcnt128(res2) - lzcnt128(res3)) > 126) {
      res.i128 = -1;
      return res;
    }
    res.i128 += res2.i128*res3.i128;
  }
  return res;
}

myint128_t nchoosek128x(uint64_t n, int64_t k) {
  myint128_t res;
  uint64_t res64 = nchoosek(n,k);
  if (res64 > 0) {
    res.i128 = res64;
    return res;
  }
  return nchoosek128(n,k);
}

void printint128hex(__int128 A) {
  myint128_t Aunion;
  Aunion.i128 = A;
  printf("0x%.16lx%.16lx\n", Aunion.hilo[INT128_HIGH], Aunion.hilo[INT128_LOW]);
}

void strtomyint128(unsigned char* input, myint128_t *res) {
  // Converts a string representation of signed decimal integer to myint128_t.
  // The string can have any number of leading spaces or zeroes.
  // No formatted input such as '+' or ',' is supported.
  // Up to and including 38 decimal digits is safe.
  uint8_t i;
  res->i128 = 0;
  i = 0;
  _Bool ispositive = true;
  while (input[i] == 0x20) {
    i++;
  }
  if (input[i] == 0x2d){
    ispositive = false;
    i++;
  }
  if (input[i] == 0x22) {
    i++;
    if ((input[i] == 0x10) || (input[i] == 0x12)) {
      ispositive = false;
      i++;
    } else {
      return;
    }
  }
  if (input[i] == 0xe2) {
    i++;
    if (input[i] == 0x88) {
      i++;
      if (input[i] == 0x92) {
        ispositive = false;
        i++;
      } else {
        return;
      }
    } else {
      return;
    }
  }
  while (input[i] != 0) {
    if ((input[i] >= '0') && (input[i] <= '9')) {
      res->i128 = res->i128*10 + (input[i] - 48);
    } else {
      break;
    }
    i++;
  }
  if (!ispositive) res->i128 *= -1;
}

void myint128todec(char* res, int8_t reslen, myint128_t A, int8_t *firstcharres) {
  //Puts decimal signed value of A in res starting from index firstcharres. 
  //e.g.
  //  char decres[41];
  //  uint8_t firstchar;
  //  myint128_t C_i128; 
  //  myint128todec(decres, sizeof(decres), C_i128, &firstchar);
  //  printf("0x%.16lx%.16lx\n", C_i128.hilo[INT128_HIGH], C_i128.hilo[INT128_LOW]);
  //  printf("%s\n", &decres[firstchar]);
  _Bool ispositive = true;
  if (A.i128 < 0) {
    ispositive = false;
    A.i128 *= -1;
  }
  *firstcharres = 0;
  if (reslen == 1) res[0] = 0;
  if (reslen <= 1) return;
  res[reslen - 1] = 0;
  int8_t i = reslen - 2;
  int8_t firstchar = i;
  uint32_t j;
  for (; i >= 0; i--) {
    j = (A.i128 % 10);
    if (A.i128 == 0) {
      if (ispositive) {
        if (i == reslen - 2) {
          res[i] = 0x30;
        } else {
          res[i] = 0x20;
        }
      } else {
        ispositive = true;
        res[i] = 0x2d;
        firstchar = i;
      }
    } else {
      res[i] = 0x30 + j;
      firstchar = i;
    }
    A.i128 = (A.i128 - j) / 10;
  }
  *firstcharres = firstchar;
};


_Bool isqrt128(const __int128 x, __int128 *isqrtx) {
  // NOTE: isqrtx always highest int with square <= x.
  // Supports 0 <= x < (1 << 126) 
  if (x == 0) {
    isqrtx[0] = 0;
    return 1;
  }
  if (x >= ((__int128)1 << 126)) return 0;
  if (x < 0) return 0;
  uint32_t i;
  __int128 y, ai;
  y = x;
  i = 0;
  while (!(y & ((__int128)3 << 124))) {
    y = y << 2;
    i++;
  }
  y = y >> 122; // x approx. = y << (122 - 2i)
  ai = ((__int128)(2 | (y >= 9))) << (61 - i);   
  ai = (ai + x/ai)/2;
  ai = (ai + x/ai)/2;
  ai = (ai + x/ai)/2;
  ai = (ai + x/ai)/2;
  while (!((ai*ai <= x) && ((ai+1)*(ai+1) > x))) {    
    ai = (ai + x/ai)/2;
  }
  if (isqrtx != NULL) isqrtx[0] = ai;
  return 1;
}

void fxp128_div(myint128_t A, myint128_t B, myint128_t *C) {
  __int128 temp, sign, Q, R;
  sign = 1;
  if (A.i128 < 0) {
    sign = -1;
    A.i128 *= -1;
  }
  if (B.i128 < 0) {
    sign *= -1;
    B.i128 *= -1;
  }  
  C->i128 = 0;
  if (A.i128 == 0) return;
  while (!((B.i128 | A.i128) >> 126)) {
    B.i128 = B.i128 << 1;
    A.i128 = A.i128 << 1;
  }
  if (B.hilo[INT128_HIGH] != 0) {
    Q = A.hilo[INT128_HIGH];
    Q = Q << 64;
    temp = Q;
    Q /= (__int128)B.hilo[INT128_HIGH];
    R = temp - Q*(__int128)B.hilo[INT128_HIGH];   
    C->i128 = Q << 56;
    temp = R + (__int128)A.hilo[INT128_LOW];
    temp = temp << 56;
    temp /= (__int128)B.hilo[INT128_HIGH];
    C->i128 += temp;
    temp = Q;
    temp *= (__int128)B.hilo[INT128_LOW] >> 8;
    temp /= (__int128)B.hilo[INT128_HIGH];
    C->i128 -= temp;
  }
  C->i128 *= sign;
}

void print_myint128(myint128_t N) {
  char decres[41];
  int8_t firstchar;
  myint128todec(decres, sizeof(decres), N, &firstchar);
  printf("%s", &decres[firstchar]);    
}

void fxp128_mul(myint128_t A, myint128_t B, myint128_t *C) {
/*
Fixed point 128bit arithmetic
=============================

Integers represent [-128, 128) range.
*/
  __int128 temp, sign;
  sign = 1;
  if (A.i128 < 0) {
    sign = -1;
    A.i128 *= -1;
  }
  if (B.i128 < 0) {
    sign *= -1;
    B.i128 *= -1;
  }
  C->i128 = A.hilo[INT128_HIGH];
  C->i128 *= B.hilo[INT128_HIGH];
  C->i128 = C->i128 << 8;
  temp = A.hilo[INT128_HIGH];
  temp *= B.hilo[INT128_LOW];
  temp = temp >> 56;
  C->i128 += temp;
  temp = A.hilo[INT128_LOW];
  temp *= B.hilo[INT128_HIGH];
  temp = temp >> 56;
  C->i128 += temp;
  temp = A.hilo[INT128_LOW];
  temp *= B.hilo[INT128_LOW];
  temp = temp >> 120;
  C->i128 += temp;
  C->i128 *= sign;
}

void egcdint128(myint128_t A, myint128_t B, myint128_t *S, myint128_t *T, myint128_t *C) {
  if ((A.i128 <= 0) || (B.i128 <= 0)) {
    C->i128 = 0;
    return;
  }
  __int128 Ai, Bi, Qi, Siprev, Tiprev, Si, Ti, temp;  
  _Bool swapAB = false;
  Ai = A.i128;
  Bi = B.i128;
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
      C->i128 = Bi;
      if (swapAB) {
        T->i128 = Siprev;
        S->i128 = Tiprev;
      } else {
        S->i128 = Siprev;
        T->i128 = Tiprev;
      }
      return;
    }
  }
}

void gcdint128(myint128_t A, myint128_t B, myint128_t *C) {
  __int128 temp;
  if ((A.i128 <= 0) || (B.i128 <= 0)) {
    C->i128 = 0;
    return;
  }
  while (1) {
    if (A.i128 < B.i128) {
      temp = A.i128;
      A.i128 = B.i128;
      B.i128 = temp;
    }
    A.i128 = A.i128 % B.i128;
    if (A.i128 == 0) {
      C->i128 = B.i128;
      return;
    }
  }
}


int32_t findsmallprimefactor(myint128_t A, myint128_t *C, uint32_t L) {
  if (A.i128 <= 1) {
    C->i128 = 0;
    return -1;
  }
  if ((A.i128 % 2) == 0) {
    C->i128 = 2;
    return 0;
  }
  if ((A.i128 % 3) == 0) {
    C->i128 = 3;
    return 0;
  }
  uint64_t Atemp;
  int64_t i, step;
  //__int128 i;
  step = 2;
  i = 6*((__int128)L / 6) - 1;
  if (i < 5) i = 5;
  if (A.hilo[INT128_HIGH] == 0) {
    Atemp = A.i128;
    while(i < 0x10000) {
      if (i*i > Atemp) return -1;
      if ((Atemp % i) == 0) {
        C->i128 = i;
        return 0;
      }
      i += step;
      step ^= 6;
    }  
    return -1;
  }
  while(i < 0x10000) {
    if (i*i > A.i128) return -1;
    if ((A.i128 % i) == 0) {
      C->i128 = i;
      return 0;
    }
    i += step;
    step ^= 6;
  }  
  return -1;
}

void modprod(myint128_t A, myint128_t B, myint128_t N, myint128_t *C) {
  // C = A*B mod N
  // Won't overflow if N < 2^95.
  A.i128 = A.i128 % N.i128;
  B.i128 = B.i128 % N.i128;
  if ((N.i128 >> 63) != 0) {
    __int128 a, b, c, d, e, f, x, x2, x3, x4;
    // (ax2 + bx + c)(dx2 + ex + f) 
    // = (a*d)x4 + (a*e + b*d)x3 + (a*f + c*d + b*e)x2 + (b*f + e*c)x + c*f
    x = (__int128)1 << 32;
    x2 = (x << 32) % N.i128;
    x3 = (x2 << 32) % N.i128;
    x4 = (x3 << 32) % N.i128;
    a = A.i128 >> 64;
    d = B.i128 >> 64;
    b = (A.i128 >> 32) & 0xffffffff;
    e = (B.i128 >> 32) & 0xffffffff;
    c = A.i128 & 0xffffffff;
    f = B.i128 & 0xffffffff;
    if ((N.i128 >> 91) != 0) {
      C->i128 = (a*((d*x4) % N.i128)) % N.i128;
      C->i128 = (C->i128 + (a*((e*x3) % N.i128))) % N.i128;
      C->i128 = (C->i128 + (b*((d*x3) % N.i128))) % N.i128;
      C->i128 = (C->i128 + (a*((f*x2) % N.i128))) % N.i128;
      C->i128 = (C->i128 + (c*((d*x2) % N.i128))) % N.i128;
      C->i128 = (C->i128 + (b*((e*x2) % N.i128))) % N.i128;
      C->i128 = (C->i128 + (b*((f*x) % N.i128))) % N.i128;
      C->i128 = (C->i128 + (c*((e*x) % N.i128))) % N.i128;
      C->i128 = (C->i128 + (c*f)) % N.i128;
    } else {
      C->i128 = ((a*((d*x4) % N.i128)) 
            + (a*((e*x3) % N.i128))
            + (b*((d*x3) % N.i128))
            + (a*((f*x2) % N.i128))
            + (c*((d*x2) % N.i128))
            + (b*((e*x2) % N.i128))
            + (b*((f*x) % N.i128)) 
            + (c*((e*x) % N.i128))
            + (c*f)) % N.i128;
    }
  } else {
    C->i128 = (A.i128 * B.i128) % N.i128;
  }
}


void mod0toNminus1(myint128_t A, myint128_t N, myint128_t *C) {
  if (A.i128 >= (__int128)0) {
    C->i128 = A.i128 % N.i128;
  } else {
    C->i128 = A.i128 + N.i128*((N.i128 - A.i128)/N.i128);
  }
}

int32_t crt(myint128_t A, myint128_t M, myint128_t B, myint128_t N, myint128_t *C, myint128_t *D) {
  // Chinese Remainder Theorem. 
  // Result = C mod MN where
  // (C = A mod M) && (C = B mod N) 
  myint128_t temp;
  if (M.i128 <= 0) return -1;
  if (N.i128 <= 0) return -1;
  if (A.i128 < 0) return -1;
  if (B.i128 < 0) return -1;
  if (M.i128 < N.i128) {
    temp = M;
    M = N;
    N = temp;
    temp = A;
    A = B;
    B = temp;
  }
  if ((M.i128 >> 95) != 0) {
    fprintf(stderr, "Possible Overflow In crt().\n");
  }
  A.i128 = A.i128 % M.i128;
  B.i128 = B.i128 % N.i128;
  if (A.i128 == 0) return -1;
  if (B.i128 == 0) return -1;
  if (M.i128 == 1) {
    if (N.i128 == 1) {
      D->i128 = 1;
      C->i128 = 0;
      return 0;
    } else {
      D->i128 = N.i128;
      C->i128 = B.i128;
      return 0;
    }
  } else {
    if (N.i128 == 1) {
      D->i128 = M.i128;
      C->i128 = A.i128;
      return 0;
    }
  }
  myint128_t S, T, BmA;
  egcdint128(M, N, &S, &T, &temp);
  if (temp.i128 != 1) return -1;
  BmA.i128 = (B.i128 - A.i128);  
  D->i128 = M.i128 * N.i128;
  mod0toNminus1(BmA, N, &BmA);
  mod0toNminus1(S, N, &S);
  modprod(S, BmA, N, &temp);
  C->i128 = A.i128 + M.i128*temp.i128;
  if (C->i128 >= D->i128) {
    C->i128 -= D->i128;
  }
  return 0;
  // C = A + XM = B - YN
  // XM + YN = B - A
  // SM + TN = 1
  // X = S(B - A)
  // Y = T(B - A)
}


int32_t modinv(myint128_t A, myint128_t N, myint128_t *C) {
  // C = A**(-1) mod N
  // A and N must be both +ve and relatively prime.
  if (N.i128 <= 0) return -1;
  if (A.i128 <= 0) return -1;
  myint128_t S, T, gcd;
  egcdint128(A, N, &S, &T, &gcd);
  if (gcd.i128 != 1) return -1;
  if (S.i128 >= 0) {
    C->i128 = (S.i128 % N.i128);
  } else {
    S.i128 *= -1;
    C->i128 = N.i128 - (S.i128 % N.i128);
  }
  return 0;
}



void modpow(myint128_t A, myint128_t B, myint128_t N, myint128_t *C) {
  // C = A**B mod N
  // Won't overflow if N < 2^95.
  if (A.i128 <= 0) {
    C->i128 = 0;
    return;
  }
  A.i128 = A.i128 % N.i128;
  B.i128 = B.i128 % N.i128;
  myint128_t res, temp;
  res.i128 = 1;
  temp.i128 = B.i128;
  while (temp.i128 != 0) {
    if (temp.i128 & 1) modprod(res, A, N, &res);
    //printf("%li\n", (uint64_t)res.i128);
    modprod(A, A, N, &A);
    temp.i128 = temp.i128 >> 1;
  }
  C->i128 = res.i128;
}


int32_t primetest(myint128_t A, uint32_t strength) {
  // strength 0-6 for more/less certainty.
  // Returns 0 for (probable) prime, -1 for definite composite;
  // Won't overflow if A < 2^95.
  uint16_t s[] = {2, 3, 5, 7, 11, 13, 17, 19, 
23, 29, 31, 37, 41, 43, 47, 53, 
59, 61, 67, 71, 73, 79, 83, 89, 
97, 101, 103, 107, 109, 113, 127, 131, 
137, 139, 149, 151, 157, 163, 167, 173, 
179, 181, 191, 193, 197, 199, 211, 223, 
227, 229, 233, 239, 241, 251, 257, 263, 
269, 271, 277, 281, 283, 293, 307, 311, 
313, 317, 331, 337};
  uint32_t i;
  for (i=0; i<sizeof(s)/sizeof(s[0]); i++) {
    if (A.i128 == (__int128)s[i]) return 0;
    if ((A.i128 % (__int128)s[i]) == 0) return -1;
  }
  if (strength > 6) strength = 6;
  strength = (1 << strength) + 3;
  myint128_t x, pow, res;
  pow.i128 = (A.i128 - 1)/2;
  for (i=0; i<=strength; i++) {
    x.i128 = (__int128)s[i];
    modpow(x, pow, A, &res);
    if ((res.i128 != 1) && (res.i128 != (A.i128 - 1))) return -1;
  }
  return 0;
}

void findlargeprimefactor(myint128_t A, myint128_t *C) {
  // test 19807040628568054723222962221
  // Test for compositeness before calling.
  int32_t res;
  int64_t iter, itermax;
  if (A.i128 <= 4294836225) return;  
  if ((A.i128 >> 95) != 0) {
    fprintf(stderr, "Argument out of range.\n");
    C->i128 = 0;
    return;
  }
  /*
  printf("Factorizing ");
  print_myint128(A);
  printf(" With Pollard/Rho.\n");
  */
  myint128_t c, xbig, xlittle, xbig2, xlittle2;
  myint128_t test, gcd;
  isqrt128(A.i128, &test.i128);
  isqrt128(test.i128, &test.i128);
  itermax = 5*test.hilo[INT128_LOW];
  c.i128 = 1;
  xbig.i128 = 2;
  xlittle.i128 = xbig.i128;
  gcd.i128 = 1;
  iter = 0;
  //printf("Countmax = %li\n", (uint64_t)countmax);
  while (gcd.i128 == 1) {
    modprod(xbig, xbig, A, &xbig2);
    xbig.i128 = (xbig2.i128 + c.i128);
    modprod(xbig, xbig, A, &xbig2);
    xbig.i128 = (xbig2.i128 + c.i128) % A.i128;
    modprod(xlittle, xlittle, A, &xlittle2);
    xlittle.i128 = (xlittle2.i128 + c.i128) % A.i128;
    test.i128 = ((A.i128 - xlittle.i128) + xbig.i128) % A.i128;
    if (test.i128 != 0) gcdint128(A, test, &gcd);
    iter++;
    if (iter > itermax) {
      xbig.i128 = 2;
      xlittle.i128 = xbig.i128;
      gcd.i128 = 1;
      c.i128++;
      iter = 0;
    }
  }  
  res = primetest(gcd, 6);
  if (res == 0) {
    C->i128 = gcd.i128;
  } else {
    findlargeprimefactor(gcd, C);
  }
}

/*
void findlargeprimefactor(myint128_t A, myint128_t *C) {
  // test 10720264024669830996237285949
  // Test for compositeness before calling.
  int32_t res;
  __int128 count;
  if (A.i128 <= 4294836225) return;  
  if ((A.i128 >> 95) != 0) {
    fprintf(stderr, "Argument out of range.\n");
    C->i128 = 0;
    return;
  }
  myint128_t c, xbig, xlittle, xbig2, xlittle2;
  myint128_t test, gcd;
  c.i128 = 1;
  xbig.i128 = 2;
  xlittle.i128 = xbig.i128;
  gcd.i128 = 1;
  count = 0;
  //printf("Countmax = %li\n", (uint64_t)countmax);
  while (gcd.i128 == 1) {
    modprod(xbig, xbig, A, &xbig2);
    xbig.i128 = (xbig2.i128 + c.i128);
    modprod(xbig, xbig, A, &xbig2);
    xbig.i128 = (xbig2.i128 + c.i128) % A.i128;
    modprod(xlittle, xlittle, A, &xlittle2);
    xlittle.i128 = (xlittle2.i128 + c.i128) % A.i128;
    test.i128 = ((A.i128 - xlittle.i128) + xbig.i128) % A.i128;
    if (test.i128 != 0) gcdint128(A, test, &gcd);
    count++;
  }  
  res = primetest(gcd, 6);
  if (res == 0) {
    C->i128 = gcd.i128;
  } else {
    findlargeprimefactor(gcd, C);
  }
}
*/
int32_t findlargeprimefactorwithtotient(myint128_t A, myint128_t T, myint128_t *C) {
  // Test for compositeness before calling.  
  // Assumes all small primes have already been factored out.
  int32_t i, Tpowersof2, res, resfnd;
  myint128_t Ttemp, start, temp, tempm1, temp2;
  start.i128 = 2;
  if (T.i128 < 1) {
    findlargeprimefactor(A, C);
    return 1;
  }
  modpow(start, T, A, &temp); 
  if (temp.i128 != 1) {
    findlargeprimefactor(A, C);
    return 1;
  }
  gcdint128(T, A, &temp);
  if (temp.i128 != 1) {
    res = primetest(temp, 6);
    if (res == 0) {
      C->i128 = temp.i128;
      return 0;
    } 
    findlargeprimefactor(temp, C);
    return 1;
  }
  Ttemp = T;
  Tpowersof2 = 0;
  while((Ttemp.i128 & 1) == 0) {
    Ttemp.i128 = Ttemp.i128 >> 1;
    Tpowersof2++;
  }
  if (Tpowersof2 < 2) {    
    findlargeprimefactor(A, C);
    return 1;
  }
  while(1) {
    modpow(start, Ttemp, A, &temp); 
    while (temp.i128 == 1) {
      if (start.i128 < 0x1000) {
        start.i128++;   
        modpow(start, Ttemp, A, &temp); 
      } else {
        findlargeprimefactor(A, C);
        return 1;
      }
    }
    resfnd = 0;
    modpow(start, Ttemp, A, &temp); 
    for(i=0; (i<Tpowersof2) && (resfnd == 0); i++) {
      tempm1.i128 = temp.i128 - 1;
      if (tempm1.i128 > 0) {
        gcdint128(tempm1, A, &temp2);
        if (temp2.i128 > 1) {
          resfnd = 1;
        } 
        modprod(temp, temp, A, &temp);
      }
    }
    if (resfnd == 1) {
      res = primetest(temp2, 6);
      if (res == 0) {
        C->i128 = temp2.i128;
        return 0;
      } else { 
        temp.i128 = A.i128/temp2.i128;
        res = primetest(temp, 6);
        if (res == 0) {
          C->i128 = temp.i128;
          return 0;
        }
        if (temp2.i128 < temp.i128) {
          findlargeprimefactor(temp2, C);
          return 1;
        } else {
          findlargeprimefactor(temp, C);
          return 1;
        }
      }
    } else {
      if (start.i128 < 0x1000) {
        start.i128++;   
      } else {
        findlargeprimefactor(A, C);
        return 1;
      }
    }
  }
}

void modulus128_init(myint128_t A, modulus128_t *N) {
  if (A.i128 >= 0) {
    N->modulus.i128 = A.i128;
  } else {
    N->modulus.i128 = 0;
  }
  N->X.i128 = 0;
  N->phi.i128 = 0;
  N->Q = NULL;
  N->S = NULL;
  N->z = NULL;
  N->P = NULL;
  N->Ppow = NULL;
  N->SQRT1 = NULL;
  N->primepower = NULL;
  N->primecount = 0;
  N->SQRT1count = 0;
}


void modulus128_nullify(modulus128_t *N) {
  N->modulus.i128 = 0;
  N->X.i128 = 0;
  N->phi.i128 = 0;
  N->primecount = 0;
  N->SQRT1count = 0;
  if (N->Q != NULL) {
    free(N->Q);
    N->Q = NULL;
  }
  if (N->S != NULL) {
    free(N->S);
    N->S = NULL;
  }
  if (N->z != NULL) {
    free(N->z);
    N->z = NULL;
  }
  if (N->P != NULL) {
    free(N->P);
    N->P = NULL;
  }
  if (N->Ppow != NULL) {
    free(N->Ppow);
    N->Ppow = NULL;
  }
  if (N->SQRT1 != NULL) {
    free(N->SQRT1);
    N->SQRT1 = NULL;
  }
  if (N->primepower != NULL) {
    free(N->primepower);
    N->primepower = NULL;
  }
}

void printfactors(modulus128_t *N) {
  char decres[41];
  uint8_t i;
  int8_t firstchar;
  for (i=0; i<N->primecount; i++) {
    myint128todec(decres, sizeof(decres), N->P[i], &firstchar);
    printf("(%s ^ %i) ", &decres[firstchar], N->primepower[i]);
  }
  printf("\n");
}

int32_t factorize(modulus128_t *N, myint128_t *T) {
  if (N->primecount != 0) return 0;
  if (N->modulus.i128 < 2) return -1;
  myint128_t L, Tlocal, Tsmallphi;
  if (T != NULL) Tlocal = *T;
  L.i128 = N->modulus.i128;
  Tsmallphi.i128 = 1;
  myint128_t primefactor[30];
  uint32_t primepower[30], primefactorcount = 0;
  int32_t i, j, res = findsmallprimefactor(L, primefactor + primefactorcount, 0);
  while (res == 0) {
    Tsmallphi.i128 *= primefactor[primefactorcount].i128 - 1;
    L.i128 /= primefactor[primefactorcount].i128;
    for (primepower[primefactorcount]=1; (L.i128 % primefactor[primefactorcount].i128) == 0; primepower[primefactorcount]++) {      
      L.i128 /= primefactor[primefactorcount].i128;
      Tsmallphi.i128 *= primefactor[primefactorcount].i128;
    }
    primefactorcount++;
    res = findsmallprimefactor(L, primefactor + primefactorcount, (uint32_t)primefactor[primefactorcount - 1].i128);
  }
  if ((L.i128 >> 95) != 0) {
    fprintf(stderr, "Argument out of range in factorize().\n");
    return -1;
  }
  if (T != NULL) {
    Tlocal.i128 /= Tsmallphi.i128;
  }
  if (L.i128 > (__int128)4294836225) {
   res = primetest(L, 6);
   while (res != 0) { 
    if (T == NULL) {
      findlargeprimefactor(L, primefactor + primefactorcount);
    } else {
      if (T->i128 > 2) {
        res = findlargeprimefactorwithtotient(L, Tlocal, primefactor + primefactorcount);
        if (res == 0) {
          Tlocal.i128 /= (primefactor[primefactorcount].i128 - 1);
        } else {
          T = NULL;
        }
      } else {
        findlargeprimefactor(L, primefactor + primefactorcount);
      }
    }
    L.i128 /= primefactor[primefactorcount].i128;
    for (primepower[primefactorcount]=1; (L.i128 % primefactor[primefactorcount].i128) == 0; primepower[primefactorcount]++) {      
      L.i128 /= primefactor[primefactorcount].i128;
    }
    primefactorcount++;
    res = primetest(L, 6);
   }
  }
  if (L.i128 > (__int128)1) {    
    primepower[primefactorcount] = 1;
    primefactor[primefactorcount].i128 = L.i128;
    primefactorcount++;
  }
  N->primecount = primefactorcount;
  N->P = aligned_alloc(16, primefactorcount*sizeof(myint128_t));
  N->Ppow = aligned_alloc(16, primefactorcount*sizeof(myint128_t));
  N->primepower = malloc(primefactorcount*sizeof(uint8_t));
  for (i=0; i<primefactorcount; i++) {
    N->P[i] = primefactor[i];
    N->Ppow[i] = primefactor[i];
    N->primepower[i] = primepower[i];
    for (j = 1; j < primepower[i]; j++) {
      N->Ppow[i].i128 *= N->P[i].i128;
    }
  }  
  //printfactors(N);
  return 0;
}

int32_t modulus128_getphi(modulus128_t *N) {
  if (N->modulus.i128 <= 1) return -1;
  if (N->phi.i128 != 0) return 0;
  int32_t res, pix;
  if (N->primecount == 0) {
    res = factorize(N, NULL);
    if (res != 0) return -1;
  }
  N->phi.i128 = 1;
  for (pix = 0; pix < N->primecount; pix++) {
    N->phi.i128 *= N->Ppow[pix].i128 - N->Ppow[pix].i128/N->P[pix].i128;
  }
  return 0;
}


void print_SQRT1(modulus128_t *N) {
  printf("Square Roots Of Unity Mod ");   
  print_myint128(N->modulus);
  printf(" = \n");   
  uint32_t SQRT1no;
  for (SQRT1no = 0; SQRT1no < N->SQRT1count; SQRT1no++) { 
    print_myint128(N->SQRT1[SQRT1no]);
    printf("  ");
  }
  printf("\n");  
}

void print_primes(modulus128_t *N) {
  uint8_t i;
  print_myint128(N->modulus);
  printf(" = ");    
  for (i = 0; i < N->primecount; i++) {
    printf(" (");    
    print_myint128(N->P[i]);
    printf(" ^ %i)", N->primepower[i]);    
  }
  printf("\n");  
}

int32_t modulus128_SQRT1(modulus128_t *N) {
  if (N->modulus.i128 <= 1) return -1;
  int32_t res;
  if (N->primecount == 0) {
    res = factorize(N, NULL);
    if (res != 0) return -1;
  }
  if (N->modulus.i128 > ((__int128)1 << 95)) {
    return -1;
  }
  myint128_t SQRT1powof2s[4], SQRT1, SQRT1ppowprod;
  uint32_t SQRT1powof2count, i, temp, primeno, SQRT1no, SQRT1ix;
  if (N->P[0].i128 == 2) {
    if (N->primepower[0] == 1) {
      N->SQRT1count = (uint32_t)1 << (N->primecount - 1);
      SQRT1powof2s[0].i128 = 1;
      SQRT1powof2count = 1;
    }
    if (N->primepower[0] == 2) {
      N->SQRT1count = (uint32_t)1 << N->primecount;
      SQRT1powof2s[0].i128 = 1;
      SQRT1powof2s[1].i128 = 3;
      SQRT1powof2count = 2;
    }
    if (N->primepower[0] > 2) {
      N->SQRT1count = (uint32_t)1 << (N->primecount + 1);
      SQRT1powof2s[0].i128 = 1;
      SQRT1powof2s[1].i128 = N->Ppow[0].i128 - 1;
      SQRT1powof2s[2].i128 = N->Ppow[0].i128/2 - 1;
      SQRT1powof2s[3].i128 = N->Ppow[0].i128/2 + 1;
      SQRT1powof2count = 4;
    }
  } else {
    //printf("primecount = %i\n", N->primecount);
    N->SQRT1count = (uint32_t)1 << N->primecount;
    SQRT1powof2count = 0;
  }
  N->SQRT1 = aligned_alloc(16, N->SQRT1count*sizeof(myint128_t));
  if (SQRT1powof2count == 0) {
   for (SQRT1no = 0; SQRT1no < N->SQRT1count; SQRT1no++) {
    temp = SQRT1no;
    for (primeno = 0; primeno < N->primecount; primeno++) {
      if (temp & 1) {
        SQRT1.i128 = 1;
      } else {
        SQRT1.i128 = N->Ppow[primeno].i128 - 1;
      }
      if (primeno > 0) {
        crt(SQRT1, N->Ppow[primeno], N->SQRT1[SQRT1no], SQRT1ppowprod, &N->SQRT1[SQRT1no], &SQRT1ppowprod);
      } else {
        N->SQRT1[SQRT1no] = SQRT1;
        SQRT1ppowprod = N->Ppow[0];
      }
      temp = temp >> 1;      
    }
   }
   return 0;
  } else {
   if (N->primecount == 1) {
     memcpy(N->SQRT1, SQRT1powof2s, N->SQRT1count*sizeof(myint128_t));
     return 0;
   }
   for (i = 0; i < SQRT1powof2count; i++) {
     for (SQRT1no = 0; SQRT1no < N->SQRT1count/SQRT1powof2count; SQRT1no++) {
       temp = SQRT1no;
       SQRT1ix = SQRT1no*SQRT1powof2count + i;
       N->SQRT1[SQRT1ix] = SQRT1powof2s[i];
       SQRT1ppowprod = N->Ppow[0];
       for (primeno = 1; primeno < N->primecount; primeno++) {
         if (temp & 1) {
           SQRT1.i128 = 1;
         } else {
           SQRT1.i128 = N->Ppow[primeno].i128 - 1;
         }
         crt(SQRT1, N->Ppow[primeno], N->SQRT1[SQRT1ix], SQRT1ppowprod, &N->SQRT1[SQRT1ix], &SQRT1ppowprod);
         temp = temp >> 1;      
       }
     }
   }
  }
  return 0;
}

int32_t perfectdivisor(modulus128_t *N, myint128_t *D, modulus128_t *M) {
  // M = N/D
  if (D->i128 <= 0) {
    return -1;
  }
  if ((N->modulus.i128) <= 1) {
    return -1;
  }
  if ((N->modulus.i128 % D->i128) != 0) {
    return -1;
  }
  modulus128_nullify(M);
  M->modulus.i128 = N->modulus.i128 / D->i128;
  if (N->primecount == 0) return 0;
  uint32_t iN, iM;
  myint128_t tempM, temp;
  tempM = M->modulus;
  M->primecount = 0;
  iM = 0;
  for (iN=0; iN < N->primecount; iN++) {
    if ((M->modulus.i128 % N->P[iN].i128) == 0) M->primecount++;
  }
  M->P = aligned_alloc(16, M->primecount*sizeof(myint128_t));
  M->Ppow = aligned_alloc(16, M->primecount*sizeof(myint128_t));
  M->primepower = malloc(M->primecount*sizeof(uint8_t));
  for (iN=0; iN < N->primecount; iN++) {
    temp.i128 = tempM.i128 % N->P[iN].i128;
    if (temp.i128 == 0) {    
      M->primepower[iM] = 1;
      M->P[iM] = N->P[iN];
      M->Ppow[iM] = N->P[iN];
      tempM.i128 /= N->P[iN].i128;
      temp.i128 = tempM.i128 % N->P[iN].i128;
      while (temp.i128 == 0) {
        M->primepower[iM]++;
        M->Ppow[iM].i128 *= N->P[iN].i128;
        tempM.i128 /= N->P[iN].i128;
        temp.i128 = tempM.i128 % N->P[iN].i128;
      }
      iM++;
    }
  }
  return 0;
}

void modsqrt0(modulus128_t *N, myint128_t *M) {
  // (x = 0 mod M) <=> (x^2 = 0 mod N->modulus)
  uint32_t res, i, j;
  if (N->modulus.i128 <= 1) {
    M->i128 = 0;
    return;
  }
  if (N->primecount == 0) {
    res = factorize(N, NULL);
    if (res != 0) {
      M->i128 = 0;
      return;
    }
  }
  M->i128 = 1;
  for (i=0; i<N->primecount; i++) {
    for (j=0; j<(N->primepower[i] + 1)/2; j++) {
      M->i128 *= N->P[i].i128;
    }
  }
}

int32_t modulusprod(modulus128_t *A, modulus128_t *B, modulus128_t *C) {
  // C=A*B
  uint32_t sorted, i, res, primelistix, Cmaxprimecount;
  modulus128_nullify(C);
  C->modulus.i128 = A->modulus.i128*B->modulus.i128;
  if (B->modulus.i128 != C->modulus.i128 / A->modulus.i128) return -1;
  if (A->modulus.i128 != C->modulus.i128 / B->modulus.i128) return -1;
  if (C->modulus.i128 % A->modulus.i128 != 0) return -1;
  if (C->modulus.i128 % B->modulus.i128 != 0) return -1;
  if (C->modulus.i128 <= 1) return -1;
  if ((A->primecount == 0) && (B->primecount == 0)) return 0;
  if ((A->primecount == 0) && (B->primecount > 0)) {
    res = factorize(A, NULL);
    if (res != 0) return -1;
  }
  if ((A->primecount > 0) && (B->primecount == 0)) {
    res = factorize(B, NULL);
    if (res != 0) return -1;
  }
  Cmaxprimecount = A->primecount + B->primecount;
  myint128_t tempP, tempPpow;
  myint128_t *Plist = aligned_alloc(16, Cmaxprimecount*sizeof(myint128_t));
  myint128_t *Ppowlist = aligned_alloc(16, Cmaxprimecount*sizeof(myint128_t));
  uint8_t temppower;
  uint8_t *primepowerlist = malloc(Cmaxprimecount*sizeof(uint8_t));
  primelistix = 0;
  for (i=0; i<A->primecount; i++) {
    Plist[primelistix] = A->P[i];
    Ppowlist[primelistix] = A->Ppow[i];
    primepowerlist[primelistix] = A->primepower[i];
    primelistix++;
  }
  for (i=0; i<B->primecount; i++) {
    Plist[primelistix] = B->P[i];
    Ppowlist[primelistix] = B->Ppow[i];
    primepowerlist[primelistix] = B->primepower[i];
    primelistix++;
  }
  sorted = 0;
  while(!sorted) {
    sorted = 1;
    for (primelistix=0; primelistix<Cmaxprimecount-1; primelistix++) {
      if (Plist[primelistix].i128 > Plist[primelistix+1].i128) {
        sorted = 0;
        tempP = Plist[primelistix];
        tempPpow = Ppowlist[primelistix];
        temppower = primepowerlist[primelistix];
        Plist[primelistix] = Plist[primelistix+1];
        Ppowlist[primelistix] = Ppowlist[primelistix+1];
        primepowerlist[primelistix] = primepowerlist[primelistix+1];
        Ppowlist[primelistix+1] = tempPpow;
        Plist[primelistix+1] = tempP;
        primepowerlist[primelistix+1] = temppower;
      }
    }
  }
  C->primecount = 1;
  for (primelistix=1; primelistix<Cmaxprimecount; primelistix++) {
    if (Plist[primelistix-1].i128 != Plist[primelistix].i128) C->primecount++;
  }
  C->P = aligned_alloc(16, C->primecount*sizeof(myint128_t));
  C->Ppow = aligned_alloc(16, C->primecount*sizeof(myint128_t));
  C->primepower = malloc(C->primecount*sizeof(uint8_t));
  C->P[0] = Plist[0];
  C->Ppow[0] = Ppowlist[0];
  C->primepower[0] = primepowerlist[0];
  i = 0;
  for (primelistix=1; primelistix<Cmaxprimecount; primelistix++) {
    if (Plist[primelistix-1].i128 != Plist[primelistix].i128) {
      i++;
      C->Ppow[i] = Ppowlist[primelistix];
      C->P[i] = Plist[primelistix];
      C->primepower[i] = primepowerlist[primelistix];
    } else {
      C->primepower[i] += primepowerlist[primelistix];
      C->Ppow[i].i128 *= Ppowlist[primelistix].i128;
    }
  }
  free(Plist);
  free(Ppowlist);
  free(primepowerlist);
  return 0;
}

int32_t gcoprimed(modulus128_t *N, myint128_t *D, modulus128_t *M, modulus128_t *W) {
  // Greatest Co-Prime Divisor.
  // For all prime factors P of N divisor D, 
  //   copy all of N's P prime power factors to W
  // Copy remaining prime powers to M.
  // M is the largest divisor of N than is also co-prime to D.
  // (D|W) && (MW = N) && (gcd(D, M) = 1) && (E|W => gcd(E, D) > 1)
  if (D->i128 <= 1) {
    return -1;
  }
  if ((N->modulus.i128 % D->i128) != 0) {
    return -1;
  }
  uint32_t res, i;  
  modulus128_t Dtemp, NoverDtemp;
  myint128_t temp, tempM, tempW;
  if (N->primecount == 0) {
    modulus128_init(*D, &Dtemp);
    res = factorize(&Dtemp, NULL);
    if (res != 0) {
      modulus128_nullify(&Dtemp);
      return -1;
    }
    temp.i128 = N->modulus.i128 / D->i128;
    modulus128_init(temp, &NoverDtemp);
    res = factorize(&NoverDtemp, NULL);
    if (res != 0) {
      modulus128_nullify(&Dtemp);
      modulus128_nullify(&NoverDtemp);
      return -1;
    }
    res = modulusprod(&NoverDtemp, &Dtemp, N);
    if (res != 0) {
      modulus128_nullify(&Dtemp);
      modulus128_nullify(&NoverDtemp);
      return -1;
    }
    modulus128_nullify(&Dtemp);
    modulus128_nullify(&NoverDtemp);
  }
  modulus128_nullify(M);
  modulus128_nullify(W);
  tempM.i128 = 1;
  tempW.i128 = 1;
  for (i=0; i<N->primecount; i++) {
    if ((D->i128 % N->P[i].i128) == 0) {
      tempW.i128 *= N->Ppow[i].i128;
    } else {
      tempM.i128 *= N->Ppow[i].i128;
    }
  }
  res = perfectdivisor(N, &tempM, W);
  if (res != 0) {      
    return -1;
  }
  res = perfectdivisor(N, &tempW, M);
  if (res != 0) {      
    return -1;
  }
  return 0;
}

int32_t modpsqrt(myint128_t A, modulus128_t *N, uint8_t primeno, myint128_t *C) {
  // C^2 = A mod P
  int32_t tsaM;
  uint32_t i;
  myint128_t tsab, tsac, tsaR, tsat, pow, temp;
  if (N->primecount == 0) return -1;
  if (primeno >= N->primecount) return -1;
  if (N->P[primeno].i128 < 2) {
    return -1;
  }
  A.i128 = A.i128 % N->P[primeno].i128;
  if (A.i128 == 0) {
    C->i128 = 0;
    return 0;
  }
  if (A.i128 == 1) {
    C->i128 = 1;
    return 0;
  }
  if ((N->P[primeno].i128 & (__int128)3) == 3) {
    pow.i128 = (N->P[primeno].i128 + 1)/4;
    modpow(A, pow, N->P[primeno], C);
    modprod(*C, *C, N->P[primeno], &temp);
    if (temp.i128 != A.i128) return -1;
    return 0;
  }
  if ((N->P[primeno].i128 & 7) == 5) {
    pow.i128 = (N->P[primeno].i128 + 3)/8;
    modpow(A, pow, N->P[primeno], C);
    modprod(*C, *C, N->P[primeno], &temp);
    if (temp.i128 == A.i128) return 0;
    temp.i128 = 2;
    pow.i128 = (N->P[primeno].i128 - 1)/4;
    modpow(temp, pow, N->P[primeno], &temp);
    modprod(temp, *C, N->P[primeno], C);
    modprod(*C, *C, N->P[primeno], &temp);
    if (temp.i128 == A.i128) return 0;
    return -1;
  }
  pow.i128 = (N->P[primeno].i128 - 1)/2;
  modpow(A, pow, N->P[primeno], &temp);
  if (temp.i128 != 1) return -1;
  if (N->P[primeno].i128 < 0x100) {
    temp.i128 = 2;
    while (temp.i128 < N->P[primeno].i128) {
      if (((temp.i128 * temp.i128) % N->P[primeno].i128) == A.i128) {
        C->i128 = temp.i128;
        return 0;
      }
      temp.i128++;
    }
    return -1;
  }
  if (N->S == NULL) {
    N->S = calloc(N->primecount, sizeof(uint8_t));
    N->Q = aligned_alloc(16, N->primecount*sizeof(myint128_t));
    memset(N->Q, 0, N->primecount*sizeof(myint128_t));
    N->z = aligned_alloc(16, N->primecount*sizeof(myint128_t));
    memset(N->z, 0, N->primecount*sizeof(myint128_t));
  }
  if (N->S[primeno] == 0) {
    N->z[primeno].i128 = 2;
    modpow(N->z[primeno], pow, N->P[primeno], &temp);
    while (temp.i128 != N->P[primeno].i128-1) {
      N->z[primeno].i128++;
      modpow(N->z[primeno], pow, N->P[primeno], &temp);
    }
    N->S[primeno] = 1;
    while ((pow.i128 & 1) == 0) {
      N->S[primeno]++;
      pow.i128 = pow.i128 >> 1;
    }
    N->Q[primeno] = pow;    
  }
  modpow(A, N->Q[primeno], N->P[primeno], &tsat);
  modpow(N->z[primeno], N->Q[primeno], N->P[primeno], &tsac);
  pow.i128 = (N->Q[primeno].i128 + 1)/2;
  modpow(A, pow, N->P[primeno], &tsaR); 
  tsaM = N->S[primeno];
  while (tsat.i128 != 1) {
    temp = tsat;
    i = 0;
    while (temp.i128 != 1) {
      modprod(temp, temp, N->P[primeno], &temp);
      i++;
    }
    pow.i128 = ((__int128)1) << (tsaM - i - 1);
    modpow(tsac, pow, N->P[primeno], &tsab);
    modprod(tsaR, tsab, N->P[primeno], &tsaR);
    modprod(tsab, tsab, N->P[primeno], &tsac);
    modprod(tsat, tsac, N->P[primeno], &tsat);
    tsaM = i;
  }
  C->i128 = tsaR.i128;
  modprod(tsaR, tsaR, N->P[primeno], &temp);
  if (temp.i128 == A.i128) return 0;  
  print_myint128(tsaR);
  printf("  -- Error occurred in modpsqrt().\n");
  return -1;
}

int32_t mod2powsqrt(myint128_t A, uint8_t exp, myint128_t *root) {  
  // A = 1 mod 8
  // exp >= 3
  if (exp == 3) {
    root->i128 = 1;
    return 0;
  }  
  int32_t i, e, newrootsix;
  myint128_t temp, localA, twopow, roots[4], newroots[4];
  roots[0].i128 = 1;
  roots[1].i128 = 3;
  roots[2].i128 = 5;
  roots[3].i128 = 7;
  e = 4;
  while (e < exp) {
    twopow.i128 = (__int128)1 << e;
    localA.i128 = (A.i128 % twopow.i128);
    newrootsix = 0;
    for (i = 0; i < 4; i++) {
      modprod(roots[i], roots[i], twopow, &temp);
      if (temp.i128 == localA.i128) {
        newroots[newrootsix] = roots[i];
        newrootsix++;
      }
    }
    for (i = 0; i < 4; i++) {
      temp.i128 = roots[i].i128 + (twopow.i128 >> 1);
      modprod(temp, temp, twopow, &temp);
      if (temp.i128 == localA.i128) {
        newroots[newrootsix].i128 = roots[i].i128 + (twopow.i128 >> 1);
        newrootsix++;
      }
    }
    if (newrootsix != 4) return -1;
    for (i = 0; i < 4; i++) {
      roots[i] = newroots[i];
    }
    e++;
  }
  *root = newroots[0];
  return 0;
}
/*
int32_t mod2powsqrt(myint128_t A, uint8_t exp, myint128_t *roots) {  
  // A = 1 mod 8
  // exp >= 3
  if (exp == 3) {
    roots[0].i128 = 1;
    roots[1].i128 = 3;
    roots[2].i128 = 5;
    roots[3].i128 = 7;
    return 0;
  }  
  int32_t i, res, newrootsix;
  myint128_t temp, localA, twopow, newroots[4];
  twopow.i128 = (__int128)1 << exp;
  localA.i128 = (A.i128 % twopow.i128);
  res = mod2powsqrt(A, exp-1, roots);
  if (res != 0) return -1;
  newrootsix = 0;
  for (i = 0; i < 4; i++) {
    modprod(roots[i], roots[i], twopow, &temp);
    if (temp.i128 == localA.i128) {
      newroots[newrootsix] = roots[i];
      newrootsix++;
    }
  }
  for (i = 0; i < 4; i++) {
    temp.i128 = roots[i].i128 + (twopow.i128 >> 1);
    modprod(temp, temp, twopow, &temp);
    if (temp.i128 == localA.i128) {
      newroots[newrootsix].i128 = roots[i].i128 + (twopow.i128 >> 1);
      newrootsix++;
    }
  }
  for (i = 0; i < 4; i++) {
    roots[i] = newroots[i];
  }
  if (newrootsix == 4) return 0;
  return -1;
}
*/
int32_t modppowsqrt(myint128_t A, modulus128_t *N, uint8_t primeno, myint128_t *C) {
  // If gcd(A,N->Ppow[primeno]) = 1, then C^2 = A mod N->Ppow[primeno] if soln exists.
  // (x + bP^i)^2 = x^2 + 2bxP^i + b^2P^(2i) = A mod P^n
  // bP^i = ((x + bP^i)^2 - x^2)/2x mod P^(2i)
  myint128_t bPpower, inv2x, temp;
  A.i128 = A.i128 % N->Ppow[primeno].i128;
  if ((A.i128 % N->P[primeno].i128) == 0) return -1;
  int32_t i, res;
  if (N->P[primeno].i128 == 2) {
    if (N->primepower[primeno] == 1) {
      C->i128 = 1;
      return 0;
    }
    if (N->primepower[primeno] == 2) {
      if (A.i128 == 1) {
        C->i128 = 1;
        return 0;
      } else {
        return -1;
      }
    }
    res = mod2powsqrt(A, N->primepower[primeno], C);
    /*
    newroots = aligned_alloc(16, 4*sizeof(myint128_t));
    res = mod2powsqrt(A, N->primepower[primeno], newroots);
    C->i128 = newroots[0].i128;
    free(newroots);
    */
    if (res != 0) return -1;
    return 0;
  }
  res = modpsqrt(A, N, primeno, C);
  if (res != 0) return -1;
  if (N->primepower[primeno] == 1) return 0;
  i = 1;
  while (1) {
    if (i <= 2) {
      temp.i128 = C->i128 << 1;
      res = modinv(temp, N->Ppow[primeno], &inv2x);
      if (res != 0) return -1;
    }
    temp.i128 = C->i128;
    modprod(temp, temp, N->Ppow[primeno], &temp);
    if (A.i128 >= temp.i128) {
      temp.i128 = A.i128 - temp.i128;
    } else {
      temp.i128 = N->Ppow[primeno].i128 - temp.i128 + A.i128;
    }
    modprod(inv2x, temp, N->Ppow[primeno], &bPpower);
    C->i128 += bPpower.i128;
    modprod(*C, *C, N->Ppow[primeno], &temp);
    if (temp.i128 == A.i128) {
      C->i128 %= N->Ppow[primeno].i128;
      return 0;
    }
    i++;
  }
}

int32_t modnsqrt(myint128_t A, modulus128_t *N, modulus128_t *M, modulus128_t *W) {
/*==========================================================================
    N = prod_{i = 0 to n-1} p_i^e_i
    M = 1
    W = 1
    For each p_i
      If (p_i | A) 
        K = A mod p_i^e_i
        If (K == 0)
          W *= p_i^((e_i + 1)/2)
        Else          
          j = number of times p_i divides K
          If (j odd) return -1 (no solutions)
          W *= p_i^(j/2)
          If ((K / p^j) not a square mod p) return -1 (no solutions)
        EndIf 
      Else
        X = modppowsqrt(A, p_i^e_i)
        If (M == 1) 
          M = p_i^e_i
          MX = X
        Else
          MX = crt(MX, M, X, p_i^e_i)
          M *= p_i^e_i
        EndIf  
      EndIf
    EndFor
    invW = 1
    If (W > 1)
      invW = modinv(W, M)
    EndIf
    WX = invW
    
    A successful call to this function (returned 0), only guarantees
    that square roots exist mod MW where MW | N.
    
    If M = 1 then all solutions are 0 mod W. For M > 1...
    
    All square roots of A mod MW are
    { W*((invW*M->SQRT1[i]*M->X) mod M) } 
      for 0 <= i < M->SQRT1count.
      
    All possible square roots of A mod N are contained in
    { W*((invW*M->SQRT1[i]*M->X) mod M) + jMW} 
      for 0 <= i < N->SQRT1count
          0 <= j < N/(MW)
      
    The square roots of unity are not created by this function. 
============================================================================
*/
  if (N->modulus.i128 <= 1) {
    return -1;
  }
  A.i128 = A.i128 % N->modulus.i128;
  myint128_t Atemp, temp, K;
  modulus128_t Mtemp;
  int32_t res;
  uint32_t j, i, primeno;
  if (N->primecount == 0) {
    res = factorize(N, NULL);
    if (res != 0) return -1;
  }
  Mtemp.modulus.i128 = 1;
  W->modulus.i128 = 1;
  for (primeno = 0; primeno < N->primecount; primeno++) {
    if ((A.i128 % N->P[primeno].i128) == 0) {
      K.i128 = A.i128 % N->Ppow[primeno].i128;
      if (K.i128 == 0) {        
        temp.i128 = N->P[primeno].i128;
        for (i = 1; i < (N->primepower[primeno] + 1)/2; i++) {
          temp.i128 *= N->P[primeno].i128; 
        }       
        W->modulus.i128 *= temp.i128;
      } else {
        temp = A;
        j = 0;
        while ((temp.i128 % N->P[primeno].i128) == 0) {
          temp.i128 /= N->P[primeno].i128;
          j++;
        }
        if ((j & 1) != 0) return -1;
        temp.i128 = N->P[primeno].i128;
        for (i = 1; i < j/2; i++) {
          temp.i128 *= N->P[primeno].i128; 
        }       
        W->modulus.i128 *= temp.i128;
        Atemp.i128 = A.i128 / temp.i128;
        res = modpsqrt(Atemp, N, primeno, &temp);
        if (res != 0) return -1;
      }
    }
  }
  for (primeno = 0; primeno < N->primecount; primeno++) {
    if ((A.i128 % N->P[primeno].i128) != 0) {
      res = modppowsqrt(A, N, primeno, &temp);
      if (res != 0) return -1;
      if (Mtemp.modulus.i128 == 1) {
        Mtemp.modulus.i128 = N->Ppow[primeno].i128;
        Mtemp.X = temp;
      } else {
        res = crt(Mtemp.X, Mtemp.modulus, temp, N->Ppow[primeno], &Mtemp.X, &temp);
        if (res != 0) {
          fprintf(stderr, "Error in modnsqrt() calling crt().\n");
          return -1;
        }
        Mtemp.modulus.i128 *= N->Ppow[primeno].i128;
      }
    }
  }
  if ((W->modulus.i128 > 1) && (Mtemp.modulus.i128 > 1)) {
    res = modinv(W->modulus, Mtemp.modulus, &W->X);
    if (res != 0) {
      fprintf(stderr, "Error in modnsqrt() calling modinv().\n");
      return -1;
    }
  } else {
    W->X.i128 = 1;
  }
  temp.i128 = N->modulus.i128 / Mtemp.modulus.i128;
  res = perfectdivisor(N, &temp, M);
  if (res != 0) {      
    return -1;
  }
  M->X = Mtemp.X;
  return 0;
}

void double_to_fxp128(double input, myint128_t *A) {
  if (input == 0) {
    A->i128 = 0;
    return;
  }
  if (input >= 128) {    
    A->hilo[INT128_HIGH] = 0x7fffffffffffffff;
    A->hilo[INT128_LOW] = 0xffffffffffffffff;
    return;
  }
  if (input <= -128) {    
    A->hilo[INT128_HIGH] = 0x8000000000000000;
    A->hilo[INT128_LOW] = 0;
    return;
  }
  __int128 sign = 1;  
  if (input < 0) {
    sign = -1;
    input *= -1;
  }
  int32_t pows = 0;
  while (input < 64) {
    pows++;
    input *= 2;
  }
  A->hilo[INT128_HIGH] = floor(input * pow(2,56));
  A->hilo[INT128_LOW] = 0;
  A->i128 = (A->i128 >> pows);
  A->i128 *= sign;
}

double fxp128_to_double(myint128_t A) {
  myint128_t temp;
  temp.hilo[INT128_HIGH] = 1UL << 56;
  temp.hilo[INT128_LOW] = 0;
  return (((double)A.i128) / temp.i128);
}

void fiba_slow(long i, myint128_t * res) {
  long j;
  myint128_t xj, xjp1, temp;
  if (i <= 0) {
    res->i128 = 0;
    return;
  }
  if (i > 186) {
    res->i128 = 0;
    return;
  }
  xj.i128 = 0;
  xjp1.i128 = 1;
  for (j = 0; j<i-2; j++) {
    temp = xj;
    xj = xjp1;
    xjp1.i128 += temp.i128;
  }
  res->i128 = xj.i128+xjp1.i128;
}

