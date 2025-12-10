SetFactory("OpenCASCADE");

// --------------------------------------
// Parameters
// --------------------------------------
a = 1.0;       // outside semi-axis X
b = 1.0;       // outside semi-axis Y
c = 3.0;       // outside semi-axis Z
t = 0.3;       // wall thickness
lc = 0.2;

// --------------------------------------
// OUTER ELLIPSOID
// --------------------------------------
Sphere(1) = {0,0,0,1};
Dilate {{0,0,0}, {a, b, c}} { Volume{1}; }

// --------------------------------------
// INNER ELLIPSOID (shrunken version)
// --------------------------------------
Sphere(2) = {0,0,0,1};
Dilate {{0,0,0}, {a - t, b - t, c - t}} { Volume{2}; }

// --------------------------------------
// Create half-space box: z < 0
// --------------------------------------
Box(10) = {-10, -10, -10, 20, 20, 10};

// --------------------------------------
// Cut OUTER ellipsoid with half-space
// --------------------------------------
BooleanIntersection{ Volume{1}; Delete; } { Volume{10}; }

// --------------------------------------
// Cut INNER ellipsoid with half-space
// --------------------------------------
BooleanIntersection{ Volume{2}; Delete; } { Volume{10}; Delete; }

// --------------------------------------
// Subtract inner half from outer half 
// to create a concave shell
// --------------------------------------
BooleanDifference{ Volume{1}; Delete; } { Volume{2}; Delete; }


// Get all boundary surfaces of the concave shell
shell_surfaces[] = Boundary{ Volume{1}; };

Printf("Number of boundary surfaces: %g", #shell_surfaces[]);

// IMPORTANT: After generating, you need to check in GMSH GUI which surface
// is which, then adjust the indices below. Typically:
// - The flat top (z=0) will be one surface
// - The outer curved surface will be one or more surfaces
// - The inner curved surface will be one or more surfaces

// For a simple case with 3 surfaces total:
// Usually: surface[0] = top rim, surface[1] = outside, surface[2] = inside
// But you should verify this in GMSH GUI!
Physical Surface("outside", 1) = {shell_surfaces[0]};
Physical Surface("top_rim", 2) = {shell_surfaces[1]};
Physical Surface("inside", 3) = {shell_surfaces[2]};

// --------------------------------------
// Volume physical ID
// --------------------------------------
Physical Volume("shell", 10) = {1};

// --------------------------------------
Mesh.CharacteristicLengthMin = lc;
Mesh.CharacteristicLengthMax = lc;

// Generate 3D mesh
Mesh 3;