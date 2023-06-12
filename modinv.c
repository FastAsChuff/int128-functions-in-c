#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include </home/simon/int128fns.c>

int main(int32_t argc, char* argv[]) {
  int32_t res;
  myint128_t A, C, N;
  A.i128 = 1;
  N.i128 = 1;
  if (argc > 1) strtomyint128((unsigned char*)argv[1], &A);
  if (argc > 2) strtomyint128((unsigned char*)argv[2], &N);
  res = modinv(A, N, &C);
  if (res == 0) {
    print_myint128(C);
    printf("\n");
    printint128hex(C.i128);
    printf("\n");
  } else {
    printf("No inverse\n");
  }
  exit(0);
}
