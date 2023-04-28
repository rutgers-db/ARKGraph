# ARKGraph: All-Range Approximate K-Nearest-Neighbor Graph

ARKGraph is a efficient library for all-range k-nearest-neighbor graph (KNNG) construction. It implements the algorithms orf our paper: [ARKGraph: All-Range Approximate K-Nearest-Neighbor Graph](https://github.com/rutgers-db/ARKGraph).

## Complie and Run

```bash
mkdir build && cd build
cmake ..
make exp
```

Running example experiment: `./tests/exp -N 10000 -K 10`

## Data Reading Wrapper

We predefined some data reader wrapper in function `ReadDataWrapper()`, it supports reading vector data from several common datasets.

## Index Algorithms

Our core ARKGraph and its variants algorithms are under the directory `src/CompactKNNG/`. They support the same `buildIndex()` method with different compact strategies.
