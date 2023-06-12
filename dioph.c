#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include </home/simon/int128fns.c>

int main(int32_t argc, char* argv[]) {
  myint128_t A_i128;
  myint128_t B_i128;
  myint128_t S_i128;
  myint128_t T_i128;
  myint128_t C_i128;
  myint128_t X_i128;
  myint128_t Y_i128;
  myint128_t temp;
  if (argc < 4) {
    printf("Usage: %s A B C\n", argv[0]);
    printf("This program finds x and y, if they exist, such that Ax + By = C.\n");
    exit(0);
  }
  strtomyint128((unsigned char*)argv[1], &A_i128);
  strtomyint128((unsigned char*)argv[2], &B_i128);
  strtomyint128((unsigned char*)argv[3], &C_i128);
  if ((A_i128.i128 <= 0) ||(B_i128.i128 <= 0) ||(C_i128.i128 <= 0)) {
    printf("Arguments must all be positive and less than 2^63.\n");
    exit(0);
  }
  if ((A_i128.i128 >= ((__int128)1 << 63)) ||(B_i128.i128 >= ((__int128)1 << 63)) ||(C_i128.i128 >= ((__int128)1 << 63))) {
    printf("Arguments must all be positive and less than 2^63.\n");
    exit(0);
  }
  printf("%sx + %sy = %s\n", argv[1], argv[2], argv[3]);
  egcdint128(A_i128, B_i128, &S_i128, &T_i128, &temp);
  if ((C_i128.i128 % temp.i128) != 0) {
    printf("No Solutions.\n");
    exit(0);
  }
  if (temp.i128 > 1) {
    A_i128.i128 /= temp.i128;
    B_i128.i128 /= temp.i128;
    C_i128.i128 /= temp.i128;
  }
  // if x0,y0 is a solution,
  // x = x0 + nB, y = y0 - nA are all solutions.
  X_i128.i128 = S_i128.i128*C_i128.i128;
  Y_i128.i128 = T_i128.i128*C_i128.i128;
  if (X_i128.i128 < 0) {
    temp.i128 = Y_i128.i128/A_i128.i128;
    X_i128.i128 += temp.i128*B_i128.i128;  
    Y_i128.i128 -= temp.i128*A_i128.i128;  
  } else {
    temp.i128 = X_i128.i128/B_i128.i128;
    X_i128.i128 -= temp.i128*B_i128.i128;
    Y_i128.i128 += temp.i128*A_i128.i128;
  }
  if (X_i128.i128 < 0) {
    temp.i128 = (B_i128.i128 - X_i128.i128)/B_i128.i128; 
    printf("x +ve for n >= ");
    print_myint128(temp);
    printf("\n"); 
  } else {
    temp.i128 = -1*(X_i128.i128/B_i128.i128); 
    printf("x +ve for n >= ");
    print_myint128(temp);
    printf("\n"); 
  }
  if (Y_i128.i128 < 0) {
    temp.i128 = (Y_i128.i128 - A_i128.i128)/A_i128.i128; 
    printf("y +ve for n <= ");
    print_myint128(temp);
    printf("\n"); 
  } else {
    temp.i128 = Y_i128.i128/A_i128.i128; 
    printf("y +ve for n <= ");
    print_myint128(temp);
    printf("\n"); 
  }
  printf("x = "); 
  print_myint128(X_i128);
  printf(" + "); 
  print_myint128(B_i128);
  printf("n\n"); 
  printf("y = "); 
  print_myint128(Y_i128);
  printf(" - "); 
  print_myint128(A_i128);
  printf("n\n"); 
  exit(0);
}
