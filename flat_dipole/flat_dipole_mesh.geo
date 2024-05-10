mesh_sizing = 0.002;
radius = 0.002;
length = 0.055;


Point(1) = {radius, 0, -length/2, mesh_sizing};
Point(2) = {-radius, 0, -length/2, mesh_sizing};
Point(3) = {-radius, 0, 0, mesh_sizing};
Point(4) = {radius, 0, 0, mesh_sizing};

Point(7) = {-radius, 0, length/2, mesh_sizing};
Point(8) = {radius, 0, length/2, mesh_sizing};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line(5) = {4,3};
Line(6) = {3,7};
Line(7) = {7,8};
Line(8) = {8,4};

Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Line Loop(2) = {5,6,7,8};
Plane Surface(2) = {2};
Physical Surface(1) = {{1,2}};

Mesh.Algorithm = 6;
//+
Coherence;
