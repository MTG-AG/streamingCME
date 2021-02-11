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

#include "params.h"
#include "gf.h"

#include "encrypt.h"
#include "encrypt_streaming.h"


void encrypt_init(unsigned char *acc_synd, unsigned char *e)
{
  gen_e(e);

  // first handle the identity matrix on the left
  for(int i = 0; i < SYS_T*GFBITS/8; i++)
  {
    acc_synd[i] = e[i];
  }
}

/*
 * col_nr: number of the column. we operate on 8 columns at a time, i.e. col_nr means we
 * compute the partial sum for the syndrome for columns (col_nr, col_nr+1, ..., col_nr+7)
 */
void encrypt_8cols(unsigned char *acc_synd, unsigned char *col, unsigned int col_nr, unsigned char *e)
{
  unsigned char b;
  unsigned char cur_byte;

  /*
   * compute syndrome s = He, s \in F_2^{n-k}, e \in F_2^n
   */

  for(int i = 0; i < SYS_T*GFBITS/8; i++)
  {
    cur_byte = 0;
    for(int j = 0; j < 8; j++)
    {
      b = col[i*8+j] & e[col_nr/8];

      b ^= b >> 4;
      b ^= b >> 2;
      b ^= b >> 1;
      b &= 1;

      cur_byte |= b << j;
    }

    acc_synd[i] = acc_synd[i] ^ cur_byte;
  }
}