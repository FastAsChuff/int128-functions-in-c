#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>

#define INT128_LOW 0
#define INT128_HIGH 1

typedef union {
  __int128 i128;
  uint64_t hilo[2];
} myint128_t;

void strtomyint128(char* input, myint128_t *res) {
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
  if (input[i] == 0x2d) {
    ispositive = false;
    i++;
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

void main(int32_t argc, char* argv[]) {
  /*
    char decres[42];
    uint8_t firstchar;
    myint128_t C_i128; 
    myint128todec(decres, 42, C_i128, &firstchar);
    printf("0x%.16lx%.16lx\n", C_i128.hilo[INT128_HIGH], C_i128.hilo[INT128_LOW]);
    printf("%s\n", &decres[firstchar]);
  */
  myint128_t C_i128;
  char decres[41];
  uint64_t i;
  uint8_t firstchar;
  for (i=0; i<187; i++) {    
    fiba_slow(i, &C_i128);
    myint128todec(decres, sizeof(decres), C_i128, &firstchar);
    printf("fiba(%li) = 0x%.16lx%.16lx %s\n", i, C_i128.hilo[INT128_HIGH], C_i128.hilo[INT128_LOW], &decres[firstchar]);
  }
  char *somenum = "31415926535897932384626433832795028842";
  strtomyint128(somenum, &C_i128);
  myint128todec(decres, sizeof(decres), C_i128, &firstchar);
  printf("%s\n", &decres[firstchar]);
}
