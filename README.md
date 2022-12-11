# A soft rasteriser

```txt

Added:

1. a check box to enable and disable parallel

2. time consuming output on each step, which can be seen in terminal

```

## build on linux

```bash
mkdir build && cd build
module load legacy-eng
module add qt/5.13.0
qmake ../LeedsGLRenderWindow/
make -j16
```