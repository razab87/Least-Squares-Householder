# Least-Squares-Householder
A simple C program solving bounded linear least square problems of the form

|Ax-b|=min s.t. Cx=d

where A,C are given matrices, b,d are given vectors and we search for the unknown vector x. The algorithm is based on householder transformations and implemented in C. Example problems to which it is applied can be found at the end of the code in the main method. It outputs the entries of the corresponding vectors x (written in a row).

Focus of the programming activity was not user friendliness but implementation of a non-trivial algorithm from numerical analysis.

Comments in the code and output remarks are in German.

Wrote this code a few years ago as a project in a numerical analysis course.
