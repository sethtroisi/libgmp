/* mpz_nthprime_ui(p, n) - compute the nth prime and store it in p.

Copyright 2019 Free Software Foundation, Inc.

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
  /* UNDERSTAND: this condition appears for mpz_primorial_ui */
  ASSERT (n <= GMP_NUMB_MAX);

  /* Allocate 1 LIMBS */
  MPZ_NEWALLOC (p, 1);

  PTR (p)[0] = 1;
  SIZ (p) = 1;

  if (n <= 1)
    {
      /* nth_prime(p, 0) = 2 */
      PTR (p)[0] = 2;
      return;
    }

  /* Simple proof of concept implementation, soon to be replaced. */
  while (n-- > 0)
    {
      mpz_nextprime(p, p);
    }
}
