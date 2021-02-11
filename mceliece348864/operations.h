#ifndef OPERATIONS_H
#define OPERATIONS_H

#include "crypto_kem.h"

int crypto_kem_enc(
       unsigned char *c,
       unsigned char *key,
       const unsigned char *pk
);

int crypto_kem_dec(
       unsigned char *key,
       const unsigned char *c,
       const unsigned char *sk,
       gf *L
);

int crypto_kem_keypair
(
       unsigned char *pk,
       unsigned char *sk 
);


int crypto_kem_keypair_sk_only(unsigned char *sk, unsigned char *inv_mat, gf *g, gf *L);
int crypto_kem_enc_col_start(unsigned char *c, unsigned char *e);
int crypto_kem_enc_col_update(unsigned char *c, unsigned char *e, unsigned char *col, unsigned int col_nr);
int crypto_kem_enc_col_finish(unsigned char *c, unsigned char *key, unsigned char *e);

#endif

