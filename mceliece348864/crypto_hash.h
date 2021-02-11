#include <shake256.h>

#define crypto_hash_32b(out,in,inlen) \
  hash_shake256(out,in,inlen)
