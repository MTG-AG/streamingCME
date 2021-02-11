#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "mceliece348864/api.h"
#include "mceliece348864/operations.h"
#include "mceliece348864/params.h"
#include "mceliece348864/pk_gen.h"


int test_LUKeyGen_StreamingEncaps();

int main()
{
  return test_LUKeyGen_StreamingEncaps();
}

int test_LUKeyGen_StreamingEncaps()
{
  int ret_val;
  unsigned char *ct = malloc(CRYPTO_CIPHERTEXTBYTES);
  unsigned char *ss = malloc(CRYPTO_BYTES);
  unsigned char *ss1 = malloc(CRYPTO_BYTES);
  unsigned char *pk = malloc(CRYPTO_PUBLICKEYBYTES);
  unsigned char *sk = malloc(CRYPTO_SECRETKEYBYTES);
  unsigned char *inverseMat = malloc(SYS_T*GFBITS * SYS_T*GFBITS/8);
  gf g[SYS_T + 1];
  gf L[SYS_N];
  clock_t begin;
  clock_t end;
  double time_spent;

  printf("PUBLIC KEY SIZE = %u\n", CRYPTO_PUBLICKEYBYTES);
  printf("PRIVATE KEY SIZE = %u\n", CRYPTO_SECRETKEYBYTES);
  printf("CIPHERTEXT SIZE = %u\n", CRYPTO_CIPHERTEXTBYTES);
  printf("SHARED SECRET SIZE = %u\n", CRYPTO_BYTES);

  /*
   * compute (extended) private key
   */

  printf("\ncompute extended private key\n");

  begin = clock();
  if ((ret_val = crypto_kem_keypair_sk_only(sk, inverseMat, g, L)) != 0)
  //if ((ret_val = crypto_kem_keypair(pk, sk)) != 0)
  {
      fprintf(stderr, "crypto_kem_keypair_sk_only returned <%d>\n", ret_val);
      return -1;
  }
  end = clock();
  time_spent = (double) (end - begin) / CLOCKS_PER_SEC;
  printf("KeyGen Duration: %f\n", time_spent);


  /*
   * compute public key from (extended) private key
   */

  printf("\nreconstruct public key\n");

  begin = clock();
  retrieve_pk(pk, inverseMat, g, L);
  end = clock();
  time_spent = (double) (end - begin) / CLOCKS_PER_SEC;
  printf("Pubkey Reconstruction Duration: %f\n", time_spent);


  /*
   * test encapsulation and decapsulation
   */

  printf("\ntest encapsulation with reconstructed public key\n");

  if ((ret_val = crypto_kem_enc(ct, ss, pk)) != 0)
  {
      fprintf(stderr, "crypto_kem_enc returned <%d>\n", ret_val);
      return -1;
  }
  if ((ret_val = crypto_kem_dec(ss1, ct, sk, L)) != 0)
  {
      fprintf(stderr, "crypto_kem_dec returned <%d>\n", ret_val);
      return -1;
  }
  if (memcmp(ss, ss1, CRYPTO_BYTES) != 0)
  {
      fprintf(stderr, "decrypt returned bad 'ss' value\n");
      return -1;
  }

  /*
   * test streaming encapsulation
   */

  printf("\ntest streaming encapsulation\n");

  unsigned char tmp_8cols[SYS_T*GFBITS];
  unsigned char e[SYS_N];
  crypto_kem_enc_col_start(ct, e);

  // takes 8 columns of the public key at a time
  // skip the identity matrix (SYS_T*GFBITS columns) because we handle this separately in `crypto_kem_enc_col_start()`
  for(int col = SYS_T*GFBITS; col < SYS_N; col += 8)
  {
    pk_retrieve_8cols(tmp_8cols, col, inverseMat, g, L);
    crypto_kem_enc_col_update(ct, e, tmp_8cols, col);
  }
  crypto_kem_enc_col_finish(ct, ss, e);


  if ((ret_val = crypto_kem_dec(ss1, ct, sk, L)) != 0)
  {
    fprintf(stderr, "crypto_kem_dec returned <%d>\n", ret_val);
    return -1;
  }
  if (memcmp(ss, ss1, CRYPTO_BYTES) != 0)
  {
    fprintf(stderr, "decrypt returned bad 'ss' value\n");
    return -1;
  }


  printf("\nSuccess!\n");
  return 0;
}