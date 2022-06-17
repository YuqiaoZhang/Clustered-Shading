#!/bin/bash

MY_DIR="$(cd "$(dirname "$0")" 1>/dev/null 2>/dev/null && pwd)"  

# https://github.com/ispc/ispc
ispc --emit-llvm -o "${MY_DIR}/kernels_x86.bc" -h "${MY_DIR}/kernels_x86.h" -O0 --arch=x86 --target=sse2 "${MY_DIR}/kernels.ispc"
ispc --emit-llvm -o "${MY_DIR}/kernels_x64.bc" -h "${MY_DIR}/kernels_x64.h" -O0 --arch=x86-64 --target=sse2 "${MY_DIR}/kernels.ispc"

clang++ -c -o "${MY_DIR}/kernels_x86.obj" -m32 "${MY_DIR}/kernels_x86.bc" 
clang++ -c -o "${MY_DIR}/kernels_x64.obj" "${MY_DIR}/kernels_x64.bc" 