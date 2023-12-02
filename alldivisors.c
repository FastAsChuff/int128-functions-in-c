#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include </home/simon/int128fns.c>

// gcc alldivisors.c -o alldivisors.bin -O3 -march=native -lm -Wall

int main(int argc, char **argv) {
  myint128_t A_i128, S_i128;
  modulus128_t A, Divisor;
  strtomyint128((unsigned char*)argv[1], &A_i128);
  if (A_i128.i128 < 2) exit(0);
  modulus128_init(A_i128, &A);
  if (factorize(&A, NULL) != 0) exit(1);
  all_divisors_init(A, &Divisor);
  printf("1 ");
  S_i128.i128 = Divisor.modulus.i128 + 1;
  print_myint128(Divisor.modulus);
  while (all_divisors_next(A, &Divisor) == 0) {
    printf(" ");
    print_myint128(Divisor.modulus);
    S_i128.i128 += Divisor.modulus.i128;
  }
  printf("\n");
  printf("Sum of all positive divisors = ");    
  print_myint128(S_i128);
  printf("\n");  
  modulus128_nullify(&Divisor);
  modulus128_nullify(&A);
}
