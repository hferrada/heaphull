# heaphull
This code computes the convex hull in 2D

# Authors
Hector Ferrada, Cristobal Navarro and Nancy Hitschfeld {hferrada@inf.uach.cl, cnavarro@inf.uach.cl, nancy@dcc.uchile.cl}

# Description
This is a small library with a method to compute the hull. 
This work presents an optimization technique that reduces the computational cost for building the Convex Hull from a set of points. The proposed method pre-processes the input set, filtering all points inside an eight-vertex polygon in O(n) time and returns a reduced set of candidate points, ordered and distributed across four priority queues. Experimental results show that for a normal distribution of points in two-dimensional space, the filtering approach in conjunction with the Graham scan is up to 10× faster than the qhull library, and between 1.7× to 10× faster than the Convex Hull methods available in the CGAL library.
In terms of memory efficiency, the proposed implementation manages to use from 3× to 6× less memory
than the other methods. The reason behind this memory improvement is because the proposed method stores indices of the input arrays, avoiding duplicates of the original floating points.
This library works only for 64 bits and supporting long sequences.

# Parameters 
For the testing code, the unique parameter is the value for n, the points in this method are generated randomly.
For the method HeapCHull(), which compute the convex hull, the parameter are the number o points and the two float arrays X and Y with the points (see the example in the file file test.cpp).

# Input & Output
Input: Two arrays of floating points X[0..n-1] and Y [0..n-1].
Ouput: An array of nCH integer with the indexes of the nCH points which form the convex hull for the input.
(see the example in the file file test.cpp).

# Make
To make the library just give the command 'make', this will create the lib: 'heapchull.a'.

# Compile & Linking
To use the library you must compile your program linking 'heapchull.a' and include the the header "HeapCHull.h". For example, compiling the file test.cpp included here: "g++ -std=c++11 test.cpp -o testHeaphull -O3 heapchull.a" or simply run the command 'make test'. it will create the binary 'testHeaphull'. The parameters are documented in the code

# Requisites
A compiler C++11 like than g++ version 4.9 or higher . 

# References
Please, if you want to include this tool as part of your experiments, in your references include the paper [1]

[1]. H. Ferrada, C. Navarro and N. Hitschfeld. A Filtering Technique for Fast Convex Hull Construction in R^2. To appear in Journal of Computational and Applied Mathematics, 2019.
