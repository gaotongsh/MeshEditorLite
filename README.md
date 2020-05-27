# Mesh Editor Lite

This is a simple mesh editor. It can load a `obj` mesh file, and perform two operations: mesh simplification and subdivision Loop scheme.

## How to compile

This project is compiled with CMake. On macOS, run

```bash
mkdir build
cd build
cmake ..
make
./mesh_editor <path_to_obj_model>
```

On Windows, use CMake GUI.

## Usage

Use `H` to switch among the four display modes: wireframe mode, vertex mode, face mode, face and edge mode.

Use `W`, `S`, `A`, `D`, `J`, and `K` to translate the model along the three axes.

Use `U`, `M`, `I`, `,`, `O`, and `.` to rotate the model about the three axes.

Use `C` to change the color of points and the wireframe. 

**Use Right Arrow key to carry out mesh subdivision. **

**Use Left Arrow key to carry out mesh simplification.**

Two sample models are attached in the `model` directory, `eight` and `cube`.

## Framework

This project is based on [Scotty3D](https://github.com/cmu462/Scotty3D), the homework framework of the Computer Graphics course at Carnegie Mellon University. The borrowed part only includes the *definition* of the halfedge data structure in `halfEdgeMesh.h/.cpp`. 

## Simplification

The mesh simplification algorithm is described [here](http://mgarland.org/research/quadrics.html). A detailed tutorial can be found in [Scotty3D's wiki](https://github.com/cmu462/Scotty3D/wiki/Simplification). In short, in every iteration, we need to find the best edge to collapse according to a specific quadric error. In this project, we create a priority queue with all the edges and their quadrics (the priority queue is simulated with a sorted set). Our algorithm is:

1. Take out the best edge from the queue.
2. Compute the position of the new vertex.
3. Also remove any edge adjacent to the best edge (we will re-insert them later).
4. Collapse the edge into a vertex.
5. Re-calculate quadrics of edges previously removed.
6. Re-insert these edges to the priority queue. 

Each simplification will reduce the number of faces to 1/4 of the original mesh.

The algorithm described above would sometimes fail, because the local mesh will sometimes form a tetrahedron, and collapsing one edge in it will break the mesh. We add one extra step in the loop to remove any potential tetrahedron before we go on to the next loop.

## Subdivision

The mesh subdivision Loop scheme is also described in [Scotty3D's wiki](https://github.com/cmu462/Scotty3D/wiki/Loop-Subdivision). The connectivity change is rather simple: we split each and every existing edge, and connect the new vertices of each of the original face. We also calculate the positions of the new vertices with a specific weight. Our algorithm is:

1. Iterate over all vertices, calculate their new positions, and mark them as *old* vertices.
2. Iterate over all edges and calculate the positions of the vertices to be created on them.
3. Iterate over all *old* edges, and split them all. The order is not important. Mark all created vertices and edges not along the original edges as *new*.
4. Iterate over all *new* edges, and flip any edge that touches an *old* and a *new* vertex.
5. Finally, move all the old and new vertices to their new positions.

After each subdivision, the mesh will have 4 times the number of faces of the original mesh. So after one simplification and one subdivision, the number of faces in the mesh should remain the same.

