Verifying the Bilu-Linial conjecture
====================================

This program verifies the Bilu-Linial conjecture. [0] The Bilu-Linial conjecture states that every d-regular (d > 1) graph has a signing with spectral radius at most 2 sqrt(d-1).

A _signing_ of a graph G is a mapping from EG to {-1, +1}. It is naturally identified with a real symmetric n x n matrix by taking the entries which are not edges as zero.
The _spectral radius_ of a matrix is the absolutely largest eigenvalue.

Luckily, nauty [1] can generate all regular graphs up to isomorphism for us. This program expects such a file as input. For example, verifying the Bilu-Linial conjecture for n = 10 can be done using:

    F=`mktemp`
    for i in `seq 1 9`; do nauty24r2/geng -d$i -D$i 10 >> "$F"; done
    ./bilu-linial "$F"
    rm "$F"


This program is super naive and many improvements are possible. Using nauty to generate all signings of a graph would be an obvious one.

In order to compile it, you will need to download nauty 2.4r2 from http://cs.anu.edu.au/~bdm/nauty/nauty24r2.tar.gz
and unpack it. Relative to bilu-linial.c, there should be a file nauty24r2/nauty.c. After that has been completed, you can compile this program:

    pushd nauty24r2
    ./configure
    make
    popd
    cmake .
    make


[0] Yonatan Bilu and Nathan Linial, Lifts, discrepancy and nearly optimal spectral gap, Combinatorica 26 (5) (2006) 495-519.

[1] http://cs.anu.edu.au/~bdm/nauty/


