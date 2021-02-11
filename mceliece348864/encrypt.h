/*
  This file is for Niederreiter encryption
*/

#ifndef ENCRYPT_H
#define ENCRYPT_H

void encrypt(unsigned char *, const unsigned char *, unsigned char *);

void encrypt_init(unsigned char *acc_synd, unsigned char *e);
void encrypt_8cols(unsigned char *acc_synd, unsigned char *col, unsigned int col_nr, unsigned char *e);

void gen_e(unsigned char *e);

#endif

