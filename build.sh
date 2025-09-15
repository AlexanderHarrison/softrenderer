#!/bin/bash

set -e

FLAGS="-fmax-errors=1 -mavx2 -mavx -ggdb"
WARN_FLAGS="-Wall -Wextra -Wuninitialized -Wcast-qual -Wdisabled-optimization -Winit-self -Wlogical-op -Wmissing-include-dirs -Wredundant-decls -Wshadow -Wundef -Wstrict-prototypes -Wpointer-to-int-cast -Wint-to-pointer-cast -Wconversion -Wduplicated-cond -Wduplicated-branches -Wformat=2 -Wshift-overflow=2 -Wint-in-bool-context -Wvector-operation-performance -Wvla -Wdisabled-optimization -Wredundant-decls -Wmissing-parameter-type -Wold-style-declaration -Wlogical-not-parentheses -Waddress -Wmemset-transposed-args -Wmemset-elt-size -Wsizeof-pointer-memaccess -Wwrite-strings -Wtrampolines -Werror=implicit-function-declaration"
if [ "$1" = 'release' ]; then
    BASE_FLAGS="-O2"
else
    BASE_FLAGS=""
fi
PATH_FLAGS="-I."
LINK_FLAGS="-lm libSDL3.so.0 -Wl,-rpath=."

mkdir -p build

# build emu
/usr/bin/gcc ${FLAGS} ${WARN_FLAGS} ${PATH_FLAGS} ${BASE_FLAGS} src/window.c ${LINK_FLAGS} -o build/softrenderer
