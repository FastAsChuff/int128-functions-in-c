#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include </home/simon/int128fns.c>

int main(int32_t argc, char* argv[]) {
  int32_t i, printsqrtsmax, res;
  printsqrtsmax = 100000000;
  myint128_t A_i128;
  myint128_t B_i128;
  myint128_t T_i128;
  myint128_t C_i128;
  myint128_t temp;
  if (argc < 3) {
    printf("Usage: %s Quadraticresidue Modulus [Totient]\n", argv[0]);
    printf("This program finds square roots of quadratic residues in arbitrary Modulus >= 2.\nFor faster solutions, provide the Euler Totient if known.\nModuli up to 10^27 are supported.\nCopyright Simon Goater April 2023.\n");
    exit(0);
  }
  strtomyint128((unsigned char*)argv[1], &A_i128);
  strtomyint128((unsigned char*)argv[2], &B_i128);
  modulus128_t BMOD, MMOD, WMOD;
  modulus128_init(B_i128, &BMOD);
  modulus128_init(B_i128, &MMOD);
  modulus128_init(B_i128, &WMOD);
  if (B_i128.i128 < 2) exit(1);
  if (A_i128.i128 < 0) {
    mod0toNminus1(A_i128, B_i128, &A_i128);
  } else {
    if (A_i128.i128 >= B_i128.i128) A_i128.i128 %= B_i128.i128;
  }
  if (argc > 3) {
    strtomyint128((unsigned char*)argv[3], &T_i128);
    res = factorize(&BMOD, &T_i128);
  } else {
    res = factorize(&BMOD, NULL);
  }
  if (res != 0) {    
    printf("factorize() failed.\n");
    exit(1);
  }
  print_primes(&BMOD);
  modulus128_getphi(&BMOD);
  printf("Totient = ");
  print_myint128(BMOD.phi);
  printf("\n");
  modsqrt0(&BMOD, &C_i128);
  printf("x^2 = 0 mod %s iff x = 0 mod ", argv[2]);  
  print_myint128(C_i128);
  printf("\n");
  if (A_i128.i128 == 0) exit(0);
  res = modnsqrt(A_i128, &BMOD, &MMOD, &WMOD);
  if (res == 0) {    
    res = modulus128_SQRT1(&MMOD);
    if (res != 0) {    
      if (WMOD.modulus.i128 == 1) {
        printf("One possible solution is ");
        print_myint128(MMOD.X);
        printf(" mod ");
        print_myint128(MMOD.modulus);
        printf("\n");
      }
      printf("modulus128_SQRT1() failed.\n");
      exit(1);
    }
    C_i128.i128 = MMOD.modulus.i128*WMOD.modulus.i128;
    if (WMOD.modulus.i128 > 1) {
      printf("Square roots must be zero mod ");  
      print_myint128(WMOD.modulus);
      printf("\n");    
    }
    printf("Found %i square roots of %s mod ", MMOD.SQRT1count, argv[1]);  
    print_myint128(C_i128);
    printf("\n");    
    if (MMOD.SQRT1count > printsqrtsmax) printf("Truncating after %i results.\n", printsqrtsmax);
    for (i=0; i<MMOD.SQRT1count && i<printsqrtsmax; i++) {
      printf("    ");  
      modprod(MMOD.SQRT1[i], MMOD.X, MMOD.modulus, &temp);
      modprod(temp, WMOD.X, MMOD.modulus, &temp);
      temp.i128 *= WMOD.modulus.i128; 
      print_myint128(temp);
      /*
      printf("    ");   
      modprod(temp, temp, C_i128, &temp);
      print_myint128(temp);
      */
      printf("\n");   
    }
  } else {
    printf("Failed to find square roots of %s mod %s\n", argv[1], argv[2]);  
  }
  exit(0);
}
