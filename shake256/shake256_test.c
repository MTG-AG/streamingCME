//
// Created by sdeligeorgopoulos on 21.09.2018.
//

#include <stdlib.h>
#include <stdio.h>
#include "shake256.h"

void printBstr(char *S, unsigned char *A, unsigned long long L)
{
    unsigned long long i;

    printf("%s", S);

    for (i = 0; i < L; i++)
        printf("%02X", A[i]);

    if (L == 0)
        printf("00");

    printf("\n");
}

int main()
{
    const unsigned char msg[] = "test";

    unsigned char *hash = calloc(32, sizeof(unsigned char));

    hash_shake256(hash, msg, strlen(msg));

    printBstr("hash = ", hash, 32);

    return 0;
}

