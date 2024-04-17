/* rng/gsl_rng.h
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004, 2007 James Theiler, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */
#include <stdio.h>
#include <stddef.h>
#include <rancommon.h>
__import rancommon;

/* These structs also need to appear in default.c so you can select
   them via the environment variable GSL_RNG_TYPE */

GSL_VAR const gsl_rng_type *gsl_rng_borosh13;
GSL_VAR const gsl_rng_type *gsl_rng_coveyou;
GSL_VAR const gsl_rng_type *gsl_rng_cmrg;
GSL_VAR const gsl_rng_type *gsl_rng_fishman18;
GSL_VAR const gsl_rng_type *gsl_rng_fishman20;
GSL_VAR const gsl_rng_type *gsl_rng_fishman2x;
GSL_VAR const gsl_rng_type *gsl_rng_gfsr4;
GSL_VAR const gsl_rng_type *gsl_rng_knuthran;
GSL_VAR const gsl_rng_type *gsl_rng_knuthran2;
GSL_VAR const gsl_rng_type *gsl_rng_knuthran2002;
GSL_VAR const gsl_rng_type *gsl_rng_lecuyer21;
GSL_VAR const gsl_rng_type *gsl_rng_minstd;
GSL_VAR const gsl_rng_type *gsl_rng_mrg;
GSL_VAR const gsl_rng_type *gsl_rng_mt19937;
GSL_VAR const gsl_rng_type *gsl_rng_mt19937_1999;
GSL_VAR const gsl_rng_type *gsl_rng_mt19937_1998;
GSL_VAR const gsl_rng_type *gsl_rng_r250;
GSL_VAR const gsl_rng_type *gsl_rng_ran0;
GSL_VAR const gsl_rng_type *gsl_rng_ran1;
GSL_VAR const gsl_rng_type *gsl_rng_ran2;
GSL_VAR const gsl_rng_type *gsl_rng_ran3;
GSL_VAR const gsl_rng_type *gsl_rng_rand;
GSL_VAR const gsl_rng_type *gsl_rng_rand48;
GSL_VAR const gsl_rng_type *gsl_rng_random128_bsd;
GSL_VAR const gsl_rng_type *gsl_rng_random128_glibc2;
GSL_VAR const gsl_rng_type *gsl_rng_random128_libc5;
GSL_VAR const gsl_rng_type *gsl_rng_random256_bsd;
GSL_VAR const gsl_rng_type *gsl_rng_random256_glibc2;
GSL_VAR const gsl_rng_type *gsl_rng_random256_libc5;
GSL_VAR const gsl_rng_type *gsl_rng_random32_bsd;
GSL_VAR const gsl_rng_type *gsl_rng_random32_glibc2;
GSL_VAR const gsl_rng_type *gsl_rng_random32_libc5;
GSL_VAR const gsl_rng_type *gsl_rng_random64_bsd;
GSL_VAR const gsl_rng_type *gsl_rng_random64_glibc2;
GSL_VAR const gsl_rng_type *gsl_rng_random64_libc5;
GSL_VAR const gsl_rng_type *gsl_rng_random8_bsd;
GSL_VAR const gsl_rng_type *gsl_rng_random8_glibc2;
GSL_VAR const gsl_rng_type *gsl_rng_random8_libc5;
GSL_VAR const gsl_rng_type *gsl_rng_random_bsd;
GSL_VAR const gsl_rng_type *gsl_rng_random_glibc2;
GSL_VAR const gsl_rng_type *gsl_rng_random_libc5;
GSL_VAR const gsl_rng_type *gsl_rng_randu;
GSL_VAR const gsl_rng_type *gsl_rng_ranf;
GSL_VAR const gsl_rng_type *gsl_rng_ranlux;
GSL_VAR const gsl_rng_type *gsl_rng_ranlux389;
GSL_VAR const gsl_rng_type *gsl_rng_ranlxd1;
GSL_VAR const gsl_rng_type *gsl_rng_ranlxd2;
GSL_VAR const gsl_rng_type *gsl_rng_ranlxs0;
GSL_VAR const gsl_rng_type *gsl_rng_ranlxs1;
GSL_VAR const gsl_rng_type *gsl_rng_ranlxs2;
GSL_VAR const gsl_rng_type *gsl_rng_ranmar;
GSL_VAR const gsl_rng_type *gsl_rng_slatec;
GSL_VAR const gsl_rng_type *gsl_rng_taus;
GSL_VAR const gsl_rng_type *gsl_rng_taus2;
GSL_VAR const gsl_rng_type *gsl_rng_taus113;
GSL_VAR const gsl_rng_type *gsl_rng_transputer;
GSL_VAR const gsl_rng_type *gsl_rng_tt800;
GSL_VAR const gsl_rng_type *gsl_rng_uni;
GSL_VAR const gsl_rng_type *gsl_rng_uni32;
GSL_VAR const gsl_rng_type *gsl_rng_vax;
GSL_VAR const gsl_rng_type *gsl_rng_waterman14;
GSL_VAR const gsl_rng_type *gsl_rng_zuf;

const gsl_rng_type ** gsl_rng_types_setup(void);

GSL_VAR const gsl_rng_type *gsl_rng_default;
GSL_VAR unsigned long int gsl_rng_default_seed;

gsl_rng *gsl_rng_alloc (const gsl_rng_type * T);
int gsl_rng_memcpy (gsl_rng * dest, const gsl_rng * src);
gsl_rng *gsl_rng_clone (const gsl_rng * r);

void gsl_rng_free (gsl_rng * r);

void gsl_rng_set (const gsl_rng * r, unsigned long int seed);
unsigned long int gsl_rng_max (const gsl_rng * r);
unsigned long int gsl_rng_min (const gsl_rng * r);
const char *gsl_rng_name (const gsl_rng * r);

int gsl_rng_fread (FILE * stream, gsl_rng * r);
int gsl_rng_fwrite (FILE * stream, const gsl_rng * r);

size_t gsl_rng_size (const gsl_rng * r);
void * gsl_rng_state (const gsl_rng * r);

void gsl_rng_print_state (const gsl_rng * r);

const gsl_rng_type * gsl_rng_env_setup (void);

INLINE_DECL unsigned long int gsl_rng_get (const gsl_rng * r);
INLINE_DECL double gsl_rng_uniform (const gsl_rng * r);
INLINE_DECL double gsl_rng_uniform_pos (const gsl_rng * r);
INLINE_DECL unsigned long int gsl_rng_uniform_int (const gsl_rng * r, unsigned long int n);

//~ #ifdef HAVE_INLINE

INLINE_FUN unsigned long int
gsl_rng_get (const gsl_rng * r)
{
  return (r->type->get) (r->state);
}

INLINE_FUN double
gsl_rng_uniform (const gsl_rng * r)
{
  return (r->type->get_double) (r->state);
}

INLINE_FUN double
gsl_rng_uniform_pos (const gsl_rng * r)
{
  double x ;
  do
    {
      x = (r->type->get_double) (r->state) ;
    }
  while (x == 0) ;

  return x ;
}

/* Note: to avoid integer overflow in (range+1) we work with scale =
   range/n = (max-min)/n rather than scale=(max-min+1)/n, this reduces
   efficiency slightly but avoids having to check for the out of range
   value.  Note that range is typically O(2^32) so the addition of 1
   is negligible in most usage. */

INLINE_FUN unsigned long int
gsl_rng_uniform_int (const gsl_rng * r, unsigned long int n)
{
  unsigned long int offset = r->type->min;
  unsigned long int range = r->type->max - offset;
  unsigned long int scale;
  unsigned long int k;

  if (n > range || n == 0) 
    {
      GSL_ERROR_VAL ("invalid n, either 0 or exceeds maximum value of generator",
                     GSL_EINVAL, 0) ;
    }

  scale = range / n;

  do
    {
      k = (((r->type->get) (r->state)) - offset) / scale;
    }
  while (k >= n);

  return k;
}
//~ #endif /* HAVE_INLINE */

gsl_rng *
gsl_rng_alloc (const gsl_rng_type * T)
{

  gsl_rng *r = (gsl_rng *) malloc (sizeof (gsl_rng));

  if (r == 0)
    {
      GSL_ERROR_VAL ("failed to allocate space for rng struct",
                        GSL_ENOMEM, 0);
    };

  r->state = calloc (1, T->size);

  if (r->state == 0)
    {
      free (r);         /* exception in constructor, avoid memory leak */

      GSL_ERROR_VAL ("failed to allocate space for rng state",
                        GSL_ENOMEM, 0);
    };

  r->type = T;

  gsl_rng_set (r, gsl_rng_default_seed);        /* seed the generator */

  return r;
}

int
gsl_rng_memcpy (gsl_rng * dest, const gsl_rng * src)
{
  if (dest->type != src->type)
    {
      GSL_ERROR ("generators must be of the same type", GSL_EINVAL);
    }

  memcpy (dest->state, src->state, src->type->size);

  return GSL_SUCCESS;
}

gsl_rng *
gsl_rng_clone (const gsl_rng * q)
{
  gsl_rng *r = (gsl_rng *) malloc (sizeof (gsl_rng));

  if (r == 0)
    {
      GSL_ERROR_VAL ("failed to allocate space for rng struct",
                        GSL_ENOMEM, 0);
    };

  r->state = malloc (q->type->size);

  if (r->state == 0)
    {
      free (r);         /* exception in constructor, avoid memory leak */

      GSL_ERROR_VAL ("failed to allocate space for rng state",
                        GSL_ENOMEM, 0);
    };

  r->type = q->type;

  memcpy (r->state, q->state, q->type->size);

  return r;
}

void
gsl_rng_set (const gsl_rng * r, unsigned long int seed)
{
  (r->type->set) (r->state, seed);
}

unsigned long int
gsl_rng_max (const gsl_rng * r)
{
  return r->type->max;
}

unsigned long int
gsl_rng_min (const gsl_rng * r)
{
  return r->type->min;
}

const char *
gsl_rng_name (const gsl_rng * r)
{
  return r->type->name;
}

size_t
gsl_rng_size (const gsl_rng * r)
{
  return r->type->size;
}

void *
gsl_rng_state (const gsl_rng * r)
{
  return r->state;
}

void
gsl_rng_print_state (const gsl_rng * r)
{
  size_t i;
  unsigned char *p = (unsigned char *) (r->state);
  const size_t n = r->type->size;

  for (i = 0; i < n; i++)
    {
      /* FIXME: we're assuming that a char is 8 bits */
      printf ("%.2x", *(p + i));
    }

}

void
gsl_rng_free (gsl_rng * r)
{
  RETURN_IF_NULL (r);
  free (r->state);
  free (r);
}
