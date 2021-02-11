/*
 *	Copyright 2021 MTG AG
 *
 *	Licensed under the Apache License, Version 2.0 (the "License");
 *	you may not use this file except in compliance with the License.
 *	You may obtain a copy of the License at
 *
 *		http://www.apache.org/licenses/LICENSE-2.0
 *
 *	Unless required by applicable law or agreed to in writing, software
 *	distributed under the License is distributed on an "AS IS" BASIS,
 *	WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *	See the License for the specific language governing permissions and
 *	limitations under the License.
*/

#include <stdint.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

#include "controlbits.h"
#include "pk_gen.h"
#include "params.h"
#include "root.h"
#include "util.h"
#include "lu_decomp.h"

/*
 * like the original root() function but only for the first SYS_T*GFBITS (instead of SYS_N) elements
 */
static void root_reduced(gf *out, gf *f, gf *L)
{
  int i;

  for (i = 0; i < SYS_T*GFBITS; i++)
    out[i] = eval(f, L[i]);
}

/*
 * multiply mat (with dimension rows x cols) by vec
 */
static void mat_vec_mul(unsigned char * out, unsigned const char * mat, size_t rows, size_t cols, unsigned const char * vec)
{
  unsigned char b;
  const unsigned char *mat_ptr = mat;

  int i, j;
  for (i = 0; i < cols/8; i++)
    out[i] = 0;

  for (i = 0; i < rows; i++)
  {
    b = 0;
    for (j = 0; j < cols/8; j++)
      b ^= mat_ptr[j] & vec[j];

    b ^= b >> 4;
    b ^= b >> 2;
    b ^= b >> 1;
    b &= 1;

    out[ i/8 ] |= (b << (i%8));

    mat_ptr += cols/8;
  }
}

/*
 * generate a column of \hat{H} from g and L
 */
static int gen_hat_H_col(unsigned char * pk_col, size_t col, gf * g, gf * L)
{
  gf inv;
  int i,k;
  unsigned char b;

  // inv[7] is for first column
  col = (col/8 * 8) + 7 - (col % 8);

  inv = eval(g, L[col]);
  inv = gf_inv(inv);

  for (i = 0; i < SYS_T*GFBITS/8; i++)
    pk_col[i] = 0;

  size_t cur_index = 0;
  b = 0;
  size_t left = 8;

  for (i = 0; i < SYS_T; i++)
  {
    for (k = 0; k < GFBITS ; k++)
    {
      if(left == 0)
      {
        pk_col[ cur_index ] = b;

        cur_index++;
        b = 0;
        left = 8;
      }

      b >>= 1;
      b  |= ((inv >> k) & 1) << 7;

      left--;
    }
    inv = gf_mul(inv, L[col]);
  }

  if(left == 0)
  {
    pk_col[ cur_index ] = b;
  }

  return 0;
}


/*
 * Compute the inverse matrix inv_mat in-place (small constant overhead) and constant time with LU decomposition
 */
int sk_inv_gen(unsigned char * sk, unsigned char *inv_mat, gf *g, gf *L)
{
  int i, j, k;
//  int row, c;

  //uint64_t buf[ 1 << GFBITS ];

//	unsigned char mat[ GFBITS * SYS_T ][ SYS_N/8 ];
  unsigned char b;

  gf inv[ SYS_T*GFBITS ];

  //

  g[ SYS_T ] = 1;

  for (i = 0; i < SYS_T; i++) { g[i] = load2(sk); g[i] &= GFMASK; sk += 2; }

  // filling the matrix

  root_reduced(inv, g, L);


  for (i = 0; i < SYS_T*GFBITS; i++)
    inv[i] = gf_inv(inv[i]);

  for (i = 0; i < SYS_T*GFBITS; i++)
    for (j = 0; j < SYS_T*GFBITS/8; j++)
      inv_mat[i*SYS_T*GFBITS/8 + j] = 0;

  for (i = 0; i < SYS_T; i++)
  {
    //    for (j = 0; j < SYS_N; j+=8)
    for (j = 0; j < SYS_T*GFBITS; j+=8)
      for (k = 0; k < GFBITS;  k++)
      {
//			b  = (inv[j+7] >> k) & 1; b <<= 1;
//			b |= (inv[j+6] >> k) & 1; b <<= 1;
//			b |= (inv[j+5] >> k) & 1; b <<= 1;
//			b |= (inv[j+4] >> k) & 1; b <<= 1;
//			b |= (inv[j+3] >> k) & 1; b <<= 1;
//			b |= (inv[j+2] >> k) & 1; b <<= 1;
//			b |= (inv[j+1] >> k) & 1; b <<= 1;
//			b |= (inv[j+0] >> k) & 1;
// different bit-order for LU code
        b  = (inv[j+0] >> k) & 1; b <<= 1;
        b |= (inv[j+1] >> k) & 1; b <<= 1;
        b |= (inv[j+2] >> k) & 1; b <<= 1;
        b |= (inv[j+3] >> k) & 1; b <<= 1;
        b |= (inv[j+4] >> k) & 1; b <<= 1;
        b |= (inv[j+5] >> k) & 1; b <<= 1;
        b |= (inv[j+6] >> k) & 1; b <<= 1;
        b |= (inv[j+7] >> k) & 1;

        inv_mat[ (i*GFBITS + k)*SYS_T*GFBITS/8 + j/8 ] = b;
      }

    for (j = 0; j < SYS_T*GFBITS; j++)
      inv[j] = gf_mul(inv[j], L[j]);
  }


  int P[SYS_T*GFBITS];

  // compute PA = LU
  if(Doolittle32(inv_mat, P, SYS_T*GFBITS))
  {
    // singular, try again
    return -1;
  }

  // compute U^-1 * L^-1
  invert_U8(inv_mat, SYS_T*GFBITS);
  invert_L8(inv_mat, SYS_T*GFBITS);
  multiply_UL_inplace8(inv_mat, SYS_T*GFBITS);

  // permute columns to compute A^-1 = U^-1 L^-1 P
  transpose8(inv_mat, SYS_T*GFBITS);
  undo_perm32(inv_mat, P, SYS_T*GFBITS);
  transpose8(inv_mat, SYS_T*GFBITS);

  // reverse all bits since Classic McEliece uses a different bit order...
  // TODO: unify
  for(i = 0; i < SYS_T*GFBITS; i++)
    for(j = 0; j < SYS_T*GFBITS/8; j++)
    {
      b = inv_mat[i*SYS_T*GFBITS/8 + j];
      b = (b & 0xF0) >> 4 | (b & 0x0F) << 4;
      b = (b & 0xCC) >> 2 | (b & 0x33) << 2;
      b = (b & 0xAA) >> 1 | (b & 0x55) << 1;
      inv_mat[i*SYS_T*GFBITS/8 + j] = b;
    }

  return 0;
}

/*
 * retrieve 8 columns of the public key (col, col+1, ... col+7)
 * the SYS_T*GFBITS bytes are written to cols_out.
 *
 * 8 columns at a time were chosen in order to prevent costly accesses to single columns (single bits in a byte)
 */
void pk_retrieve_8cols(unsigned char *cols_out, unsigned int col, unsigned char *inv, gf *g, gf *L)
{
  unsigned char b;
  unsigned char tmp[SYS_T * GFBITS];
  unsigned char cols_hat_H[SYS_T * GFBITS/8];

  // compute 8 columns of \hat{H} and multiply by the inverse to obtain the pubkey columns
  for(int j = 0; j < 8; j++)
  {
    gen_hat_H_col(cols_hat_H, col+j, g, L);
    mat_vec_mul(tmp + (SYS_T*GFBITS/8) * j, inv, SYS_T * GFBITS, SYS_T * GFBITS, cols_hat_H);
  }

  // reassamble the 8 bit columns into the corresponding bytes
  for(int i = 0; i < SYS_T*GFBITS/8; i++)
    for(int j = 0; j < 8; j++)
    {
      b  = (tmp[SYS_T*GFBITS/8 * 0 + i] >> j) & 1; b <<= 1;
      b |= (tmp[SYS_T*GFBITS/8 * 1 + i] >> j) & 1; b <<= 1;
      b |= (tmp[SYS_T*GFBITS/8 * 2 + i] >> j) & 1; b <<= 1;
      b |= (tmp[SYS_T*GFBITS/8 * 3 + i] >> j) & 1; b <<= 1;
      b |= (tmp[SYS_T*GFBITS/8 * 4 + i] >> j) & 1; b <<= 1;
      b |= (tmp[SYS_T*GFBITS/8 * 5 + i] >> j) & 1; b <<= 1;
      b |= (tmp[SYS_T*GFBITS/8 * 6 + i] >> j) & 1; b <<= 1;
      b |= (tmp[SYS_T*GFBITS/8 * 7 + i] >> j) & 1;

      cols_out[i*8 + j] = b;
    }
}


/*
 * Example usage that shows how the public key can be retrieved from the extended private key.
 * Writing to the pk_out buffer can be replaced by writing to any kind of stream or storage.
 */
int retrieve_pk(unsigned char *pk_out, unsigned char *inv, gf *g, gf *L)
{
  unsigned char col[SYS_T*GFBITS];

  // skip producing the identity matrix, i.e. the first SYS_T*GFBITS columns
  for(int k = SYS_T*GFBITS; k < SYS_N; k += 8)
  {
    pk_retrieve_8cols(col, k, inv, g, L);

    for(int i = 0; i < SYS_T*GFBITS; i++)
    {
      pk_out[i*PK_ROW_BYTES + (k - SYS_T*GFBITS)/8] = col[i];
    }
  }

  return 0;
}


