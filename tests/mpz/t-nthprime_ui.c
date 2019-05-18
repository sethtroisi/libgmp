/* Test mpz_nthprime_ui.

Copyright 2019 Free Software Foundation, Inc.

This file is part of the GNU MP Library test suite.

The GNU MP Library test suite is free software; you can redistribute it
and/or modify it under the terms of the GNU General Public License as
published by the Free Software Foundation; either version 3 of the License,
or (at your option) any later version.

The GNU MP Library test suite is distributed in the hope that it will be
useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
Public License for more details.

You should have received a copy of the GNU General Public License along with
the GNU MP Library test suite.  If not, see https://www.gnu.org/licenses/.  */

#include <stdio.h>
#include <stdlib.h>
#include "gmp-impl.h"
#include "tests.h"


/* Verify nth_prime(n) = want and isprime(want). */
void
check_one (long n, long want) {
  mpz_t got;
  mpz_init(got);

  mpz_nthprime_ui(got, n);
  MPZ_CHECK_FORMAT (got);

  if (mpz_cmp_si(got, want) != 0)
    {
      printf ("mpz_nthprime wrong\n");
      printf ("  n=%lu\n", n);
      printf ("  got="); mpz_out_str (stdout, 10, got); printf ("\n");
      printf ("  want=%lu\n", want);
      abort();
    }

  if (!mpz_probab_prime_p(got, 25))
    {
      printf ("mpz_nthprime_ui not prime\n");
      printf ("  n=%lu\n", n);
      printf ("  p="); mpz_out_str (stdout, 10, got); printf ("\n");
      abort();
    }

  mpz_clear(got);
}

void
check_data (void)
{
  static const struct {
    int n;
    long want;
  } data[] = {
    { 1,      2 },
    { 2,      3 },
    { 3,      5 },
    { 4,      7 },
    { 10,     29 },
    { 25,     97 },
    { 100,    541 },
    { 1000,   7919 },
    { 3141,   28843 },
  };

  for (int i = 0; i < numberof (data); i++)
    {
      check_one(data[i].n, data[i].want);
    }
}

/* check nthprime(0) is correct */
void
check_zero (void)
{
  mpz_t test;
  mpz_init(test);

  mpz_nthprime_ui(test, 0);
  if (mpz_cmp_si(test, 2))
    {
      printf ("mpz_nthprime_ui(0) != 2\n");
      printf ("  got  "); mpz_out_str (stdout, 10, test); printf("\n");
      abort ();
    }

  MPZ_CHECK_FORMAT (test);

  mpz_clear(test);
}

/* check small n's. */
void
check_small (void)
{
  mpz_t test;
  mpz_init_set_ui(test, 1);

  for (int i = 1; i < 1000; i++)
    {
      mpz_nextprime(test, test);

      check_one(i, mpz_get_ui(test));
    }

  mpz_clear(test);
}

int
main (int argc, char **argv)
{
  tests_start ();

  check_zero ();
  check_small ();
  check_data ();

  tests_end ();
  exit (0);
}
