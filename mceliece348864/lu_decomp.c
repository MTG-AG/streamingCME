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

#include <stdio.h>
#include <stdint.h>

#include "lu_decomp.h"


/*
 * undo permutation P in constant time
 */
void undo_perm8(unsigned char *mat, int *P, int n)
{
  unsigned char mask_cpy;
  unsigned char tmp;

  for(int i = n-1; i >= 0; i--)
  {
    for(int j = 0; j < n; j++)
    {
      mask_cpy = (j == P[i]);
      mask_cpy = -mask_cpy;

      for(int c = 0; c < n/8; c++)
      {
        tmp = mat[(n/8)*i + c];
        mat[(n/8)*i + c] = (mat[(n/8)*i + c] & ~mask_cpy) | (mat[(n/8)*j + c] & mask_cpy);
        mat[(n/8)*j + c] = (mat[(n/8)*j + c] & ~mask_cpy) | (tmp & mask_cpy);
      }
    }
  }
}


/*
 * undo permutation P in constant time
 * assumes (n%32) == 0
 */
void undo_perm32(unsigned char *mat, int *P, int n)
{
  uint32_t mask_cpy;
  uint32_t tmp;
  uint32_t *mat_32_ptr = (uint32_t*)mat;

  for(int i = n-1; i >= 0; i--)
  {
    for(int j = 0; j < n; j++)
    {
      mask_cpy = (j == P[i]);
      mask_cpy = -mask_cpy;

      for(int c = 0; c < n/32; c++)
      {
        tmp = mat_32_ptr[(n/32)*i + c];
        mat_32_ptr[(n/32)*i + c] = (mat_32_ptr[(n/32)*i + c] & ~mask_cpy) | (mat_32_ptr[(n/32)*j + c] & mask_cpy);
        mat_32_ptr[(n/32)*j + c] = (mat_32_ptr[(n/32)*j + c] & ~mask_cpy) | (tmp & mask_cpy);
      }
    }
  }
}


/*
 * undo permutation P in constant time
 * assumes (n%64) == 0
 */
void undo_perm64(unsigned char *mat, int *P, int n)
{
  uint64_t mask_cpy;
  uint64_t tmp;
  uint64_t *mat_64_ptr = (uint64_t*)mat;

  for(int i = n-1; i >= 0; i--)
  {
    for(int j = 0; j < n; j++)
    {
      mask_cpy = (j == P[i]);
      mask_cpy = -mask_cpy;

      for(int c = 0; c < n/64; c++)
      {
        tmp = mat_64_ptr[(n/64)*i + c];
        mat_64_ptr[(n/64)*i + c] = (mat_64_ptr[(n/64)*i + c] & ~mask_cpy) | (mat_64_ptr[(n/64)*j + c] & mask_cpy);
        mat_64_ptr[(n/64)*j + c] = (mat_64_ptr[(n/64)*j + c] & ~mask_cpy) | (tmp & mask_cpy);
      }
    }
  }
}

/*
 * transposes the matrix mat
 * assumes that n is divisible by 8 and operates on 8x8 chunks
 *
 * in-place (constant overhead)
 */
void transpose8(unsigned char *mat, int n)
{
  // naive solution ...

  // temporarily store a 8x8 chunk
  unsigned char tmp_rows[8];

  // 8x8 block at once
  for(int i = 0; i < n; i += 8)
  {

    // handle i == j case (on diagonal)
    for(int k = 0; k < 8; k++)
    {
      tmp_rows[k] = mat[(n/8)*(i+k) + i/8];
      mat[(n/8)*(i+k) + i/8] = 0;
    }
    for(int k = 0; k <= 7; k++)
    {
      for(int l = k; l < 8; l++)
      {
        mat[(n/8)*(i+k) + i/8] |= (((0x80 >> k) & tmp_rows[l]) << k) >> l;
      }
    }
    for(int k = 7; k > 0; k--)
    {
      for(int l = 0; l <= k-1; l++)
      {
        mat[(n/8)*(i+k) + i/8] |= (((0x80 >> k) & tmp_rows[l]) << k) >> l;
      }
    }

    // i != j
    for (int j = i+8; j < n; j+= 8)
    {
      for(int k = 0; k < 8; k++)
      {
        tmp_rows[k] = mat[(n/8)*(i+k) + j/8];
        mat[(n/8)*(i+k) + j/8] = 0;
      }

      // write from columns to rows
      for(int k = 0; k < 8; k++)
      {
        for(int l = 0; l < 8; l++)
        {
          mat[(n/8)*(i+k) + j/8] |= (((0x80 >> k) & mat[(n/8)*(j+l) + i/8]) << k) >> l;
        }
      }

      // write from rows to columns
      for(int k = 0; k < 8; k++)
      {
        mat[(n/8)*(j+k) + i/8] = 0;
        for(int l = 0; l < 8; l++)
        {
          mat[(n/8)*(j+k) + i/8] |= (((0x80 >> k) & tmp_rows[l]) << k) >> l;
        }
      }
    }
  }
}

/*
 * Computes the LU Decomposition of a binary matrix in-place in constant time.
 *
 * mat: the matrix that shall be decomposed into L and U (in-place).
 *      one unsigned char represents 8 bits of the binary matrix, i.e. mat[0] equals the first 8 elements of the first row.
 *      L will be in the lower triangle matrix of mat and U in the upper triangle matrix.
 * pivot: int array that determines how the rows are switched (length N)
 *      on 32 bit system this will probably result in N * 4 bytes additional memory
 * n: size of the square matrix mat.
 *
 * assumes dimensions (in bit) of mat are divisible by 8.
 */
int Doolittle8(unsigned char *mat, int *pivot, int n)
{

  /* LU pseudo code (inplace, without pivoting)
      // n-1 Iterationsschritte
      for k := 1 to n-1
       // Zeilen der Restmatrix werden durchlaufen
       for i := k+1 to n
         // Berechnung von L
         A(i,k) := A(i,k) / A(k,k) // Achtung: vorher Prüfung auf Nullwerte notwendig
         // Spalten der Restmatrix werden durchlaufen
         for j := k+1 to n
           // Berechnung von R
           A(i,j) := A(i,j) - A(i,k) * A(k,j)
   */

  for(int k = 0; k < n-1; k++)
  {

    /*
     * swap rows in constant time ...
     */

    pivot[k] = k;

    unsigned char tmp;
    for(unsigned int j = k+1; j < n; j++)
    {
      unsigned int mask_cpy = (mat[(n/8)*j + k/8] >> (7 - (k % 8))) & 1; // note: type has to be large enough to use as mask with j in pivot computing
      mask_cpy = mask_cpy & !((mat[(n/8)*k + k/8] >> (7 - (k % 8))) & 1);  // add this in order to only actually switch with the first suitable row => consistent with "ref" non-constant-time impl
      mask_cpy = -mask_cpy;
      for(int c = 0; c < n/8; c++)
      {
        tmp = mat[(n/8)*j + c];
        mat[(n/8)*j + c] = (mat[(n/8)*j + c] & ~mask_cpy) | (mat[(n/8)*k + c] & mask_cpy);
        mat[(n/8)*k + c] = (mat[(n/8)*k + c] & ~mask_cpy) | (tmp & mask_cpy);
      }

      // update pivot
      pivot[k] = (j & mask_cpy) | (pivot[k] & ~mask_cpy);
    }

    if (!((mat[k*n/8 + (k/8)] >> (7 - (k % 8))) & 1))
    {
      return -1;
    }

    // for i := k+1 to n
    //      for j := k+1 to n
    //          A(i,j) := A(i,j) - A(i,k) * A(k,j)
    //
    // Bitwise: A(i,j) := A(i,j) ^ A(i,k) & A(k,j)
    // idea: take element A(i,k). It's either 0 or 1. Expand to 8 bits. Multiply with A(k,j), ... A(k,j+7)
    // basically, bit A(i,k) decides whether row A(k) (starting from the k+1'th column) is added to the i+1'th, ..., n'th rows below
    // use a mask for constant time instead of if(bit==1)


    for(int i = k+1; i < n; i++)
    {
      // 0x00 iff A(i,k) is zero. Otherweise 0xFF
      unsigned char mask_ik = 0x00 - ((mat[i*n/8 + (k/8)] >> (7 - (k % 8))) & 1);

      /*
       * First, handle the case that j = k+1 is not divisible by 8; i.e. only handle the last j%8 bits of this byte
       */
      unsigned char mask_prepending = (unsigned char) (0x100 >> ((k+1) % 8)) - 1;

      unsigned char tmp = mat[n/8*i + (k+1)/8] & ~mask_prepending; // store the left side of the splitted byte in tmp
      mat[i*(n/8) + (k+1)/8]  ^= (mask_ik & mat[k*(n/8) + (k+1)/8]);

      // restore the left bits: first set bits 0, then 'or' tmp
      mat[n/8*i + (k+1)/8] &= mask_prepending;
      mat[n/8*i + (k+1)/8] |= tmp;


      /*
       * now handle the rest of the row efficiently in byte-aligned byte-chunks
       */
      for(int j = k+1 + (8-((k+1)%8)); j < n; j += 8)
      {
        // A(i,j) = A(i,j) + mask_ik * A(k,j)
        // => XOR Bits A(i,j) ... A(i,j+7) with (mask_ik & (A(k,j), ... A(k,j+7))
        mat[i*(n/8) + j/8] ^= mask_ik & mat[k*(n/8) + j/8];
      }
    }
  }
  pivot[n-1] = n-1;

  if ((mat[(n/8)*(n-1) + n/8-1] & 0x01) == 0)
  {
    return -1;
  }
  return 0;
}


/* assumes (n % 32) == 0 */
int Doolittle32(unsigned char *mat, int *pivot, int n)
{
  uint32_t *mat_32_ptr = (uint32_t*)mat;

  for(int k = 0; k < n-1; k++)
  {

    /*
     * swap rows in constant time ...
     */

    pivot[k] = k;

    uint32_t tmp;
    for(uint32_t j = k+1; j < n; j++)
    {
      uint32_t mask_cpy = (mat[(n/8)*j + k/8] >> (7 - (k % 8))) & 1; // note: type has to be large enough to use as mask with j in pivot computing
      mask_cpy = mask_cpy & !((mat[(n/8)*k + k/8] >> (7 - (k % 8))) & 1);  // add this in order to only actually switch with the first suitable row => consistent with "ref" non-constant-time impl
      mask_cpy = -mask_cpy;
      for(int c = 0; c < n/32; c++)
      {
        // copy in 32 bit chunks
        tmp = mat_32_ptr[(n/32)*j + c];
        mat_32_ptr[(n/32)*j + c] = (mat_32_ptr[(n/32)*j + c] & ~mask_cpy) | (mat_32_ptr[(n/32)*k + c] & mask_cpy);
        mat_32_ptr[(n/32)*k + c] = (mat_32_ptr[(n/32)*k + c] & ~mask_cpy) | (tmp & mask_cpy);
      }

      // update pivot
      pivot[k] = (j & mask_cpy) | (pivot[k] & ~mask_cpy);
    }

    if (!((mat[k*n/8 + (k/8)] >> (7 - (k % 8))) & 1))
    {
      return -1;
    }

    for(int i = k+1; i < n; i++)
    {
      // 0x00 iff A(i,k) is zero. Otherweise 0xFF
      uint32_t mask_ik = 0x00000000 - ((mat[i*n/8 + (k/8)] >> (7 - (k % 8))) & 1);

      /*
       * First, handle the case that j = k+1 is not divisible by 8; i.e. only handle the last j%8 bits of this byte
       */
      unsigned char mask_prepending = (unsigned char) (0x100 >> ((k+1) % 8)) - 1;

      unsigned char tmp = mat[n/8*i + (k+1)/8] & ~mask_prepending; // store the left side of the splitted byte in tmp
      mat[i*(n/8) + (k+1)/8]  ^= (mask_ik & mat[k*(n/8) + (k+1)/8]);

      // restore the left bits: first set bits 0, then 'or' tmp
      mat[n/8*i + (k+1)/8] &= mask_prepending;
      mat[n/8*i + (k+1)/8] |= tmp;

      // handle byte-aligned chunks until we are 4-byte aligned
      int k_1_aligned = k+1 + (8-((k+1)%8));
      int k_4_aligned = k_1_aligned;
      for(int j = k_1_aligned; (j < n) && ((j%32) != 0); j += 8)
      {
        mat[i*(n/8) + j/8] ^= mask_ik & mat[k*(n/8) + j/8];
        k_4_aligned = j+8;
      }

      /*
       * now handle the rest of the row efficiently in 4-byte-aligned byte-chunks
       */
      for(int j = k_4_aligned; j < n; j += 32)
      {
        mat_32_ptr[i*(n/32) + j/32] ^= mask_ik & mat_32_ptr[k*(n/32) + j/32];
      }
    }
  }
  pivot[n-1] = n-1;

  if ((mat[(n/8)*(n-1) + n/8-1] & 0x01) == 0)
  {
    return -1;
  }
  return 0;
}



/* assumes (n % 64) == 0 */
int Doolittle64(unsigned char *mat, int *pivot, int n)
{
  uint64_t *mat_64_ptr = (uint64_t*)mat;

  for(int k = 0; k < n-1; k++)
  {

    /*
     * swap rows in constant time ...
     */

    pivot[k] = k;

    uint64_t tmp;
    for(uint64_t j = k+1; j < n; j++)
    {
      uint64_t mask_cpy = (mat[(n/8)*j + k/8] >> (7 - (k % 8))) & 1; // note: type has to be large enough to use as mask with j in pivot computing
      mask_cpy = mask_cpy & !((mat[(n/8)*k + k/8] >> (7 - (k % 8))) & 1);  // add this in order to only actually switch with the first suitable row => consistent with "ref" non-constant-time impl
      mask_cpy = -mask_cpy;
      for(int c = 0; c < n/64; c++)
      {
        // copy in 32 bit chunks
        tmp = mat_64_ptr[(n/64)*j + c];
        mat_64_ptr[(n/64)*j + c] = (mat_64_ptr[(n/64)*j + c] & ~mask_cpy) | (mat_64_ptr[(n/64)*k + c] & mask_cpy);
        mat_64_ptr[(n/64)*k + c] = (mat_64_ptr[(n/64)*k + c] & ~mask_cpy) | (tmp & mask_cpy);
      }

      // update pivot
      pivot[k] = (j & mask_cpy) | (pivot[k] & ~mask_cpy);
    }

    if (!((mat[k*n/8 + (k/8)] >> (7 - (k % 8))) & 1))
    {
      return -1;
    }

    for(int i = k+1; i < n; i++)
    {
      // 0x00 iff A(i,k) is zero. Otherweise 0xFF
      uint64_t mask_ik = 0x0000000000000000 - ((mat[i*n/8 + (k/8)] >> (7 - (k % 8))) & 1);

      /*
       * First, handle the case that j = k+1 is not divisible by 8; i.e. only handle the last j%8 bits of this byte
       */
      unsigned char mask_prepending = (unsigned char) (0x100 >> ((k+1) % 8)) - 1;

      unsigned char tmp = mat[n/8*i + (k+1)/8] & ~mask_prepending; // store the left side of the splitted byte in tmp
      mat[i*(n/8) + (k+1)/8]  ^= (mask_ik & mat[k*(n/8) + (k+1)/8]);

      // restore the left bits: first set bits 0, then 'or' tmp
      mat[n/8*i + (k+1)/8] &= mask_prepending;
      mat[n/8*i + (k+1)/8] |= tmp;

      // handle byte-aligned chunks until we are 8-byte aligned
      int k_1_aligned = k+1 + (8-((k+1)%8));
      int k_8_aligned = k_1_aligned;
      for(int j = k_1_aligned; (j < n) && ((j%64) != 0); j += 8)
      {
        mat[i*(n/8) + j/8] ^= mask_ik & mat[k*(n/8) + j/8];
        k_8_aligned = j+8;
      }

      /*
       * now handle the rest of the row efficiently in 4-byte-aligned byte-chunks
       */
      for(int j = k_8_aligned; j < n; j += 64)
      {
        mat_64_ptr[i*(n/64) + j/64] ^= mask_ik & mat_64_ptr[k*(n/64) + j/64];
      }
    }
  }
  pivot[n-1] = n-1;

  if ((mat[(n/8)*(n-1) + n/8-1] & 0x01) == 0)
  {
    return -1;
  }
  return 0;
}


/*
 * invert an upper triangle matrix with diagonal only 1s
 * Assumes n%8 == 0
 * in-place, constant time
 */
void invert_U8(unsigned char *U, int n)
{
  /*
   * inversion == backwards substitution
   *
   * For all columns c = n-1 ... 0
   *   For all rows c-1 ... 0
   *     add (row c) * (leading (c'th) element of row r) to row r.
   *       only do it for the k elements k = c+1 ... n-1
   */
  for(int c = n-2; c >= 0; c--)
  {
    /*
     * Basically:
     *      copy elements c+1 to n-1 from c'th row to r'th row iff leading element of row r is 1
     */

    for(int r = c-1; r >= 0; r--)
    {
      /*
       * first copy all "full bytes" at the right
       */
      unsigned char mask_rc = 0x00 - ((U[r*n/8 + (c/8)] >> (7 - (c % 8))) & 1); // mask for the element U(r,c)
      for(int k = n-1; k >= c+1 + 8; k -= 8)
      {
        U[r*(n/8) + k/8] ^= mask_rc & U[c*(n/8) + k/8];
      }

      /*
       * then handle the left-most "partial byte"
       */
      int k = c+1;
      unsigned char mask_prepending = (unsigned char) (0x100 >> (k % 8)) - 1;

      unsigned char tmp_rk = U[r*(n/8) + k/8] & ~mask_prepending; // store the left side of the splitted byte in tmp
      U[r*(n/8) + k/8]  ^= (mask_rc & U[c*(n/8) + k/8]);

      // restore the left bits: first set bits 0, then 'or' tmp
      U[r*(n/8) + k/8] &= mask_prepending;
      U[r*(n/8) + k/8] |= tmp_rk;
    }
  }
}


/*
 * invert a lower triangle matrix with diagonal only 1s
 * Assumes n%8 == 0
 * in-place, constant time
 */
void invert_L8(unsigned char *L, int n)
{
  /*
   * inversion == forward substitution
   *
   * For all columns c = 1 ... n-1
   *   For all rows c+1 ... n-1
   *     add (row c) * (leading (c'th) element of row r) to row r.
   *       only do it for the k elements k = 0 ... c-1
   */

  for(int c = 1; c <= n-1; c++)
  {
    /*
     * Basically:
     *      copy elements 0 to c-1  from c'th row to r'th row iff leading element of row r is 1
     */
    for(int r = c+1; r <= n-1; r++)
    {
      /*
       * first copy all "full bytes" at the left
       */
      unsigned char mask_rc = 0x00 - ((L[r*n/8 + (c/8)] >> (7 - (c % 8))) & 1); // mask for the element U(r,c)
      for(int k = 0; k <= c-1 - 8; k += 8)
      {
        L[r*(n/8) + k/8] ^= mask_rc & L[c*(n/8) + k/8];
      }

      /*
       * then handle the right-most "partial byte"
       */
      int k = c-1;
      unsigned char mask_prepending = ((unsigned char) (0x80 >> (k % 8)) - 1);

      unsigned char tmp_rk = L[r*(n/8) + k/8] & mask_prepending; // store the right side of the splitted byte in tmp
      L[r*(n/8) + k/8]  ^= (mask_rc & L[c*(n/8) + k/8]);

      // restore the right bits: first set bits 0, then 'or' tmp
      L[r*(n/8) + k/8] &= ~mask_prepending;
      L[r*(n/8) + k/8] |= tmp_rk;
    }
  }
}

/*
 * multiplies the matrices L and U that are inside mat (after LU  decomposition)
 * in-place, constant time
 *
 * uses additional space of n/8 bytes (one row / col)
 */
void multiply_UL_inplace8(unsigned char *mat, int n)
{
  // compute A_11 first, then the surrounding elements, and so on
  // => the upper-left 1x1 matrix
  // => the upper-left 2x2 matrix
  // => ...
  // for the i x i submatrix, the new elements are computed by the product of a (n-i) dimensional row and column vector each.
  // In the i'th step, we compute the column A(j, i) for 0 <= j <= i
  // In the i'th step, we compute the row A(i, j) for 0 <= j < i
  // A column element A(j, i) is computed by A(j, i) ^= A(j, k) & A(k, i) for k = i+1, ..., n
  // A row element A(i, j) is computed by A(i, j) ^= A(k, j) & A(i, k) for k = i+1, ..., n

  // store one colum intermediately
  unsigned char col[n/8];

  // compute upper-left i x i matrix
  for(int i = 0; i < n-1; i++)
  {
    // compute i'th column that is used i times for the computation of the column elements A(i, j) j <= i
    col[(i+1)/8] = 0;
    for(int c = i+1; c < n; c++)
    {
      if(c % 8 == 0)
      {
        col[c/8] = 0;
      }

      col[c/8] |= (((0x80U >> i%8) & mat[(n/8)*c + i/8]) << i%8) >> c%8;
    }

    unsigned char col_mul_mask = (0x100 >> ((i+1) % 8)) - 1;

    // decide to view the diagonal element i*i as part of the column (and not of the row)
    // compute column elements
    for(int j = 0; j <= i; j++)
    {
      // compute the j'th column element A(j, i)
      // first handle the possible partial byte at the beginning
      unsigned char parity_bit = col[(i+1)/8] & mat[(n/8)*j + (i+1)/8];
      parity_bit &= col_mul_mask;

      // handle the rest of the row
      for(int k = (i+1) + (8-(i+1)%8); k < n; k+=8)
      {
        // ith_col times the kth row: A(j, i) ^= A(j, k) & A(k, i)
        parity_bit ^= mat[(n/8)*j + k/8] & col[k/8];
      }
      // now actually compute and set the bit A(j, i)
      parity_bit ^= parity_bit >> 4;
      parity_bit ^= parity_bit >> 2;
      parity_bit ^= parity_bit >> 1;
      // TODO: simplify....?
      mat[(n/8)*j + i/8] = (mat[(n/8)*j + i/8] & (unsigned char)~(0x80U >> (i%8))) | (((parity_bit << (7-i%8)) ^ mat[(n/8)*j + i/8]) & (1 << (7-(i%8))));
    }

    // compute row elements
    for(int j = 0; j < i; j+=8)
    {
      // compute the j'th row element A(i, j)
      for(int k = i+1; k < n; k++)
      {
        // multiply the 8 elements A(k, j) ... A(k, j+7) by A(i, k)
        unsigned char mask_ik = 0x00 - ((mat[(n/8)*i + k/8] >> (7 - (k % 8))) & 1); // 0x00 if bit is zero, 0xFF otherwise

        if(j + 8 > i)
        {
          // handle the case that i is not divisible by 8 and we only partially set the first i % 8 bits
          unsigned char mask_prepending = ~((unsigned char) (0x100 >> (i % 8)) - 1);
          unsigned char tmp_ij = mat[(n/8)*i + j/8] & ~mask_prepending;
          mat[(n/8)*i + j/8] ^= mat[(n/8)*k + j/8] & mask_ik;

          mat[(n/8)*i + j/8] &= mask_prepending;
          mat[(n/8)*i + j/8] |= tmp_ij;
        }
        else
        {
          mat[(n/8)*i + j/8] ^= mat[(n/8)*k + j/8] & mask_ik;
        }
      }
    }
  }
}
