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

#ifndef LU_DECOMP_H
#define LU_DECOMP_H


int Doolittle8(unsigned char *mat, int *pivot, int n);
int Doolittle32(unsigned char *mat, int *pivot, int n);
int Doolittle64(unsigned char *mat, int *pivot, int n);
void invert_U8(unsigned char *U, int n);
void invert_L8(unsigned char *U, int n);
void multiply_UL_inplace8(unsigned char *mat, int n);
void transpose8(unsigned char *mat, int n);
void undo_perm8(unsigned char *mat, int *P, int n);
void undo_perm32(unsigned char *mat, int *P, int n);

#endif //LU_DECOMP_H
