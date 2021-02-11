cmake_minimum_required(VERSION 3.10)
project(cme348864 C)

set(CMAKE_C_STANDARD 99)


# McEliece reference code
set(MCELIECE_DIR mceliece348864)
set(MCELIECE_SRC
        ${MCELIECE_DIR}/operations.c
        ${MCELIECE_DIR}/encrypt.c
        ${MCELIECE_DIR}/encrypt_streaming.c
        ${MCELIECE_DIR}/decrypt.c
        ${MCELIECE_DIR}/sk_gen.c
        ${MCELIECE_DIR}/pk_gen.c
        ${MCELIECE_DIR}/controlbits.c
        ${MCELIECE_DIR}/gf.c
        ${MCELIECE_DIR}/root.c
        ${MCELIECE_DIR}/synd.c
        ${MCELIECE_DIR}/util.c
        ${MCELIECE_DIR}/bm.c
        ${MCELIECE_DIR}/aes.c
        ${MCELIECE_DIR}/lu_decomp.c
        )

# shake256 code
set(SHAKE_DIR shake256)
set(SHAKE_SRC
        ${SHAKE_DIR}/shake256.c
        ${SHAKE_DIR}/shake256_impl.c
        )


set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")

add_executable(cme348864 main.c rng.c ${MCELIECE_SRC} ${SHAKE_SRC})

target_include_directories(cme348864 PRIVATE ${SHAKE_DIR})


add_executable(shake256_test shake256/shake256.c shake256/shake256.h shake256/shake256_impl.c shake256/shake256_impl.h
        shake256/shake256_test.c)