//
// Created by jroth on 13.06.2019.
//

#include "rng.h"
#include <stdlib.h>


static unsigned long int next = 100;

int rand2(void) // RAND_MAX assumed to be 32767
{
    next = next * 1103515245 + 12345;
    return (unsigned int)(next/65536) % 32768;
}

void srand2(unsigned int seed)
{
    next = seed;
}

int randombytes(unsigned char *x, unsigned long long xlen)
{
    for(size_t pos = 0; pos < xlen; pos++)
    {
        x[pos] = rand2() % 256;
    }

    return 0;
}