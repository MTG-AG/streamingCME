#include "shake256.h"
#include "shake256_impl.h"

void hash_shake256_variable_output(unsigned char *output, size_t outputLen, const unsigned char *input, size_t inputLen)
{
    shake256(output, outputLen, input, inputLen);
}

void hash_shake256(unsigned char *output, const unsigned char *input, size_t inputLen)
{
    hash_shake256_variable_output(output, 32, input, inputLen);
}
