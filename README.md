# cozyhome-distance-gjk-js
This repository houses my implementation of the Gilbert Johnson Keerthi algorithm for convex sets in two dimensions. The output of this algorithm returns the simplex used to minimize the norm of the configuration space obstacle as well as the minimized point in configuration space. 

Using this simplex, it is trivial to construct a nearest pair of vertices via barycentric coordinates. I've designed this algorithm to return as much information as possible to ensure ease of usability. As well as this, it is optionally able to exploit frame coherency by reinserting a simplex from a prior iteration.
