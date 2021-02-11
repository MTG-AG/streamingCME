#ifndef SHAKE256_H
#define SHAKE256_H

#include <string.h>

extern void hash_shake256_variable_output(unsigned char *output, size_t outputLen, const unsigned char *input, size_t inputLen);

extern void hash_shake256(unsigned char *output, const unsigned char *input, size_t inputLen);

#endif
