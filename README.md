# ARKGraph: All-Range Approximate K-Nearest-Neighbor Graph

ARKGraph is a efficient library for all-range k-nearest-neighbor graph (KNNG) indexing, which can process range KGraph queries in almost real-time. For more details please refer to our paper: [ARKGraph: All-Range Approximate K-Nearest-Neighbor Graph](https://dl.acm.org/doi/10.14778/3603581.3603601).


## Usage

### Complie and Run

```bash
mkdir build && cd build
cmake ..
make exp
```

Running example experiment: `./tests/exp -N 10000 -K 10`

### Data Reading Wrapper

We predefined some data reader wrapper in function `ReadDataWrapper()`, it supports reading vector data from several common datasets.

### Index Algorithms

Our core ARKGraph and its variants algorithms are under the directory `src/CompactKNNG/`. They support the same `buildIndex()` method with different compact strategies.

## Reference

ARKGraph paper:

```
@article{10.14778/3603581.3603601,
author = {Zuo, Chaoji and Deng, Dong},
title = {ARKGraph: All-Range Approximate K-Nearest-Neighbor Graph},
year = {2023},
issue_date = {June 2023},
publisher = {VLDB Endowment},
volume = {16},
number = {10},
issn = {2150-8097},
url = {https://doi.org/10.14778/3603581.3603601},
doi = {10.14778/3603581.3603601},
abstract = {Given a collection of vectors, the approximate K-nearest-neighbor graph (KGraph for short) connects every vector to its approximate K-nearest-neighbors (KNN for short). KGraph plays an important role in high dimensional data visualization, semantic search, manifold learning, and machine learning. The vectors are typically vector representations of real-world objects (e.g., images and documents), which often come with a few structured attributes, such as times-tamps and locations. In this paper, we study the all-range approximate K-nearest-neighbor graph (ARKGraph) problem. Specifically, given a collection of vectors, each associated with a numerical search key (e.g., a timestamp), we aim to build an index that takes a search key range as the query and returns the KGraph of vectors whose search keys are within the query range. ARKGraph can facilitate interactive high dimensional data visualization, data mining, etc. A key challenge of this problem is the huge index size. This is because, given n vectors, a brute-force index stores a KGraph for every search key range, which results in O(Kn3) index size as there are O(n2) search key ranges and each KGraph takes O(Kn) space. We observe that the KNN of a vector in nearby ranges are often the same, which can be grouped together to save space. Based on this observation, we propose a series of novel techniques that reduce the index size significantly to just O(Kn log n) in the average case. Furthermore, we develop an efficient indexing algorithm that constructs the optimized ARKGraph index directly without exhaustively calculating the distance between every pair of vectors. To process a query, for each vector in the query range, we only need O(log log n + K log K) to restore its KNN in the query range from the optimized ARKGraph index. We conducted extensive experiments on real-world datasets. Experimental results show that our optimized ARKGraph index achieved a small index size, low query latency, and good scalability. Specifically, our approach was 1000x faster than the baseline method that builds a KGraph for all the vectors in the query range on-the-fly.},
journal = {Proc. VLDB Endow.},
month = {jun},
pages = {2645–2658},
numpages = {14}
}
```
