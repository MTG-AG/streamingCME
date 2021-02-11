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

#ifndef PK_GEN_H
#define PK_GEN_H

#include "gf.h"

int sk_inv_gen(unsigned char * sk, unsigned char *inv_mat, gf *g, gf *L);
int retrieve_pk(unsigned char *pk_out, unsigned char *inv, gf *g, gf *L);
void pk_retrieve_8cols(unsigned char *cols_out, unsigned int col, unsigned char *inv, gf *g, gf *L);


#endif