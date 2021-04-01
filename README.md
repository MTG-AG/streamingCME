This project contains the code base that was used for the paper Classic McEliece Implementation with Low Memory Footprint ([Springer](https://link.springer.com/chapter/10.1007/978-3-030-68487-7_3) | [ePrint](https://eprint.iacr.org/2021/138)).

The code can be built with cmake. Has been tested with GCC (WSL).


# Difference to the Reference Code

The code is based on the [Classic McEliece round2 reference code](https://classic.mceliece.org/nist.html). Some of the changes are outlined:

* A freestanding shake256 and AES implementation have been added (removes OpenSSL dependency).
* The key generation process has been changed according to the algorithms that are detailled in the paper.
  * Code for in-place matrix inversion (using the LU Decomposition in lu_decomp.c/h).
  * Since the Benes network is not used in the reference implementation and the existing implementation has a huge memory overhead, it has been removed from the code. Instead, the support elements are simply stored as an array.
* The streaming encapsulation has been implemented in operations.c/h and encrypt_streaming.c/h.
* The RNG is a simple linear congruent RNG that can (and should) be replaced.


## Other Parameter Sets

This code is for the mceliece348864 parameter set which produces the smallest keys out of the [proposed parameter sets](https://classic.mceliece.org/nist/mceliece-20190331.pdf#section.3).
It can be used for other parameter sets as well. However, mind the small differences in the reference code for the parameter sets.
Also, the mceliece6960119 parameter set produces matrices with dimensions that are not divisible by 8. 
This implementation does not currently take this into account and therefore does not work for mceliece6960119.


# Additional Notes

This code does not contain all tweaks that were used on the ARM board. This is for two reasons. First, the code is more clear and readable this way. Second, some parts are device specific.
The most notable tweak will be outlined next.

The Classic McEliece reference code makes some large memory allocations in the key generation process. Almost all of them can be avoided if using the memory that already is allocated by the matrix S which is part of the extended private key. 
Fortunately, generating S is the last step of the key generation process. By using this technique, the "Peak Memory Usage" can be kept low. Some places to look out for are:

* crypto_kem_keypair_sk_only(): 
  * unsigned char r[ SYS_T*2 + (1 << GFBITS)*sizeof(uint32_t) + SYS_N/8 + 32 ];
  * uint32_t perm[ 1 << GFBITS ]; 
* gen_L(): uint64_t buf[ 1 << GFBITS ];
* genpoly_gen(): gf mat[ SYS_T+1 ][ SYS_T ];
* perm_check(): uint64_t list[1 << GFBITS];


# License

See LICENSE.txt file. Parts of the code are licensed under the CC0 / MIT / Apache 2.0 license.