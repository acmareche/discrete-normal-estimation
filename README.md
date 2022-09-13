# discrete-normal-estimation
Normal vector estimation on a digital surface with a local sectorial approach.

Using DGtal 1.3 http://dgtal.org (in beta at the time of writing) with libQGLViewer and Qt5.

Usage example :
```shell
./estimation goursat025.vol 0 0 33 12
```

 ```goursat025.vol``` is your digital surface in .vol format.

 ```0 0 33``` are the (x, y, z) coordinates of the point on the digital surface for which you want the normal estimation. Inputing the coordinates of a point that is NOT on the surface might result in undefined behavior and sad error messages.

 ```12``` is the maximum radius you want the program to try for the size of the neighborhood of the point.

The program will compute the smallest radius for which the point's circular neighborhood on the digital surface is no longer planar and split it in planar angular sectors, from which it will estimate normals using COBA and Chord's algorithm. It will also compute normal estimations with the Integral Invariants and Plane Probing methods (see [DGtal's documentation](https://dgtal-team.github.io/doc-nightly/) for more details).

