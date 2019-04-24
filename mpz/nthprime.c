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


/* Enhancements:

   - Use sieve for smallish n
   - Use precomputed lookup values as started point
   - Implement a more modern algorithm */


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

  gmp_primesieve_t ps;
  gmp_init_primesieve (&ps);

  while (n-- > 1)
    {
      gmp_nextprime (&ps);
    }

  PTR (p)[0] = gmp_nextprime (&ps);

  return;
}
