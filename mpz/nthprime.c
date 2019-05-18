/* mpz_nthprime_ui(p, n) - compute the nth prime and store it in p.

Copyright 2018 Free Software Foundation, Inc.

Contributed to the GNU project by Seth Troisi

This file is part of the GNU MP Library.

The GNU MP Library is free software; you can redistribute it and/or modify
it under the terms of either:

  * the GNU Lesser General Public License as published by the Free
    Software Foundation; either version 3 of the License, or (at your
    option) any later version.

or

  * the GNU General Public License as published by the Free Software
    Foundation; either version 2 of the License, or (at your option) any
    later version.

or both in parallel, as here.

The GNU MP Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received copies of the GNU General Public License and the
GNU Lesser General Public License along with the GNU MP Library.  If not,
see https://www.gnu.org/licenses/.  */

#include "gmp-impl.h"
#include "longlong.h"
#include "math.h"
#include "stdio.h"

/*************************************************************/
/* Section macros: common macros, for swing/fac/bin (&sieve) */
/*************************************************************/

#define LOOP_ON_SIEVE_CONTINUE(prime,end,sieve)			\
    __max_i = (end);						\
								\
    do {							\
      ++__i;							\
      if (((sieve)[__index] & __mask) == 0)			\
	{							\
          mp_limb_t prime;					\
	  prime = id_to_n(__i)

#define LOOP_ON_SIEVE_BEGIN(prime,start,end,off,sieve)		\
  do {								\
    mp_limb_t __mask, __index, __max_i, __i;			\
								\
    __i = (start)-(off);					\
    __index = __i / GMP_LIMB_BITS;				\
    __mask = CNST_LIMB(1) << (__i % GMP_LIMB_BITS);		\
    __i += (off);						\
								\
    LOOP_ON_SIEVE_CONTINUE(prime,end,sieve)

#define LOOP_ON_SIEVE_STOP					\
	}							\
      __mask = __mask << 1 | __mask >> (GMP_LIMB_BITS-1);	\
      __index += __mask & 1;					\
    }  while (__i <= __max_i)

#define LOOP_ON_SIEVE_END					\
    LOOP_ON_SIEVE_STOP;						\
  } while (0)


/*********************************************************/
/* Section sieve: sieving functions and tools for primes */
/*********************************************************/

/* id_to_n (x) = bit_to_n (x-1) = (id*3+1)|1*/
static mp_limb_t
id_to_n  (mp_limb_t id)  { return id*3+1+(id&1); }

/* n_to_bit (n) = ((n-1)&(-CNST_LIMB(2)))/3U-1 */
static mp_limb_t
n_to_bit (mp_limb_t n) { return ((n-5)|1)/3U; }

static mp_size_t
primesieve_size (mp_limb_t n) { return n_to_bit(n) / GMP_LIMB_BITS + 1; }



/* This threshold determines when to stop using gmp_nextprime and switch to
   primepi method.  */
#define NTH_PRIME_SIEVE_THRESHOLD 1000000

static unsigned long
mpz_primepi_isqrt(unsigned long x)
{
  // Enhancement: faster sqrt, call mpn_sqrtrem1 directly?
  mp_limb_t res;
  const mp_limb_t src = x;
  mpn_sqrtrem(&res, NULL, &src, 1);
  return (unsigned long) res;
}

/* Enhancements:

   - Use sieve for smallish n
   - Use precomputed lookup values as started point
   - Implement a more modern algorithm */

/* P1(x, z) = count of numbers <= x, <s>not</s> divisible by any prime <= Pend
 */
static unsigned long
mpz_primepi_p1(mp_ptr primes, unsigned long x, mp_size_t end)
{
/* Enhancements:
   - can be precomputed when z is small
   - cache results
   - loop unroll when primepi_p1_ui(x, z) = 1 */

  if (end == 0) { return 0; }
  if (x <= 0) { return 0; }

  unsigned long count = 0;

  if (x <= primes[end-1]) {
    return x - 1;
  }

  for (mp_size_t pi = 0; pi < end; pi++)
    {
      unsigned long p = primes[pi];
      unsigned long temp = x / p;

      count += temp;
      count -= mpz_primepi_p1(primes, temp, pi);
    }

/*
  mp_size_t pi = end - 1;
  for (; pi >= 0; pi--)
    {
      unsigned long p = primes[pi];
      unsigned long temp = x / p;

      if (temp > p);
        break;

      count += temp;
    }
  for(; pi >= 0 ; pi--)
    {
      unsigned long p = primes[pi];
      unsigned long temp = x / p;
      count += temp;
      count -= mpz_primepi_p1(primes, temp, pi);
    }
*/
/*
  unsigned long count = 0;
  mp_size_t pi = start;
  for (; pi < end; pi++)
    {
      unsigned long p = primes[pi];
      unsigned long temp = x / p;
      count += temp;

      if (temp < p)
        break;
      count -= mpz_primepi_p1(primes, temp, pi + 1, end);
    }

  pi++;
  for(; pi < end; pi++)
    {
      unsigned long p = primes[pi];
      unsigned long temp = x / p;
      count += temp;
    }
*/

//  printf("%lu, %lu = %lu\n", x, start, count);
  return count;
}

/* Legendre's Formula
 * let x_2 = x^(1/2)
 * primepi(x) = primepi(x_2) + P1(x, x_2) - 1
 * P1(x, z) = count of numbers <= x, not divisible by any prime <= Pz
 */
static unsigned long
mpz_primepi_ui(mp_ptr primes, mp_size_t primepi, unsigned long x)
{
  // Primes should have been sieved

  // Check if last prime >= x
  //    if so binary search for answer
  mp_limb_t last_prime = primes[primepi - 1];
  if (last_prime > x)
    {
      // Binary search
      mp_size_t low = 0;
      mp_size_t high = primepi -  1;
      while (low < high)
        {
          mp_size_t test = (low + high + 1) / 2;
          if (primes[test] <= x)
            low = test;
          else
            high = test - 1;
        }
      ASSERT ( low == high );
      ASSERT ( primes[low] <= x );
      ASSERT ( primes[low + 1] > x );
      return low + 1;
    }

  unsigned long x_2 = mpz_primepi_isqrt(x);
  unsigned long primepi_2 = mpz_primepi_ui(primes, primepi, x_2);

  // requires oversizing by prime_gap
  ASSERT (last_prime >= x_2);
  unsigned long p1 = x - mpz_primepi_p1(primes, x, primepi_2);
//  unsigned long p1 = x - mpz_primepi_p1(primes, x, 0, primepi_2);

  printf("x: %lu, x_2: %lu => %lu, %lu to %lu\n",
      x, x_2, primepi_2, primes[0], primes[primepi_2-1]);
  printf("primepi(%lu) = %lu + %lu = %lu\n", x, primepi_2, p1, primepi_2 + p1 - 1);

  return primepi_2 + p1 - 1;
}

void
mpz_nthprime_ui (mpz_ptr p, unsigned long n)
{
  // UNDERSTAND: this condition appears for mpz_primorial_ui
  ASSERT (n <= GMP_NUMB_MAX);

  if (n == 0)
    {
      // FIXME: is 0 reasonable here?
      MPZ_REALLOC (p, 0);
      return;
    }

  // Allocate 1 LIMBS for now.
  MPZ_NEWALLOC (p, 1)[0] = 1;
  SIZ (p) = 1;

#if 0
  while (n-- > 0)
    {
      mpz_nextprime(p, p);
    }
  return;
#endif

  if (n <=  NTH_PRIME_SIEVE_THRESHOLD)
    {
      gmp_primesieve_t ps;
      gmp_init_primesieve (&ps);

      while (n-- > 1)
        {
          gmp_nextprime (&ps);
        }

      PTR (p)[0] = gmp_nextprime (&ps);
      return;
    }

  // Choose X such that primepi(X) is slightly below n
  // if primepi(X) = n_x is very close to n
  //    call nextprime(X,X) n - n_x times
  // Refine a better exstimate of X and re run primepi(x)

  // P_n-est ~= n * (Log[n] + (Log[Log[n]] - 1))
  // P_n - P_n-est > 0 (for n < 10^18)
  // P_n-est / P_n ~= 0.999 for large n
  // Enhancement figure out real log.
  double logN = log(n);
  unsigned long X = n * ( logN + (log(logN) - 1));

  // Abuse knowledge of implementation needing (x+e)^(1/2) primes.
  // +e is case we readjust x because we guessed too low.
  unsigned long high_x_2 = mpz_primepi_isqrt(1.01 * X);

  mp_ptr primes;
  mp_size_t primepi;

  TMP_DECL;
  TMP_MARK;
  {
    mp_size_t end = high_x_2;
    mp_size_t sieve_size = primesieve_size(end);

    mp_ptr sieve = TMP_ALLOC_LIMBS (sieve_size);
    primepi = 2 + gmp_primesieve (sieve, end);

    primes = TMP_ALLOC_LIMBS(primepi);

    /* Store primes from 5 to n */
    primes[0] = 2;
    primes[1] = 3;

    mp_size_t j = 2;

    LOOP_ON_SIEVE_BEGIN (prime, n_to_bit(5), n_to_bit (end), 0, sieve);
    primes[j++] =  prime;
    LOOP_ON_SIEVE_END;
    // how to unallocate sieve here?
  }

//  printf("Primes <= %lu\n", high_x_2);
  for (mp_size_t i = 0; i < primepi; i += 75)
    {
//      printf("\tprime(%lu) = %lu\n", i, primes[i]);
    }

  unsigned long n_est = mpz_primepi_ui(primes, primepi, X);
  long delta = n - n_est;
//  printf("Looking for %lu\n", n);
//  printf("prime(%lu) = %lu, delta=%ld\n", X, n_est, delta);
  ASSERT (delta > 0);

  // Do a refined estimate
  unsigned long Y = X + (n - n_est) * 0.99 * log(X);
  // Make sure we sieved far enough
  ASSERT ( mpz_primepi_isqrt(Y) <= high_x_2);

  unsigned long n_est2 = mpz_primepi_ui(primes, primepi, Y);
  long delta2 = n - n_est2;
//  printf("prime(%lu) = %lu, delta=%ld\n", Y, n_est2, delta2);
  ASSERT (delta2 > 0);

  TMP_FREE;

  // TODO make this support prev prime.
  SIZ(p) = 1;
  PTR (p)[0] = Y;
  for (int i = 0; i < delta2; i++)
    {
      mpz_nextprime(p, p);
    }

  return;
}
