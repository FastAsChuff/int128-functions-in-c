#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include </home/simon/int128fns.c>


int main(int32_t argc, char* argv[]) {
  nchoosek(0, 0);
  uint64_t k = 1, n = 321;
  myint128_t res128;
  if (argc > 1) n = atol(argv[1]);
  if (argc > 2) k = atol(argv[2]);
  res128 = nchoosek_simple(n, k);
  printf("nchoosek_simple(%lu, %lu) = ", n, k);
  print_myint128(res128);
  printf("\n");
  printint128hex(res128.i128);
  printf("\n");
  exit(0);
}
