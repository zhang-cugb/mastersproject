
// Define points
p0 = newp; Point(p0) = {0.0, -2.0, 0.0, 0.3 };
p1 = newp; Point(p1) = {5.0, -2.0, 0.0, 0.3 };
p2 = newp; Point(p2) = {5.0, 3.0, 0.0, 0.3 };
p3 = newp; Point(p3) = {0.0, 3.0, 0.0, 0.3 };
p4 = newp; Point(p4) = {1.0, 0.0, 0.0, 0.2 };
p5 = newp; Point(p5) = {4.0, 0.0, 0.0, 0.2 };
// End of point specification

// Start of specification of domain// Define lines that make up the domain boundary
bound_line_0 = newl;
Line(bound_line_0) ={p0, p1};
Physical Line("DOMAIN_BOUNDARY_0") = { bound_line_0 };
bound_line_1 = newl;
Line(bound_line_1) ={p1, p2};
Physical Line("DOMAIN_BOUNDARY_1") = { bound_line_1 };
bound_line_2 = newl;
Line(bound_line_2) ={p2, p3};
Physical Line("DOMAIN_BOUNDARY_2") = { bound_line_2 };
bound_line_3 = newl;
Line(bound_line_3) ={p3, p0};
Physical Line("DOMAIN_BOUNDARY_3") = { bound_line_3 };

// Line loop that makes the domain boundary
Domain_loop = newll;
Line Loop(Domain_loop) = {bound_line_0, bound_line_1, bound_line_2, bound_line_3};
domain_surf = news;
Plane Surface(domain_surf) = {Domain_loop};
Physical Surface("DOMAIN") = {domain_surf};
// End of domain specification

// Start specification of fractures/compartment boundary/auxiliary elements
frac_line_0 = newl; Line(frac_line_0) = {p4, p5};
Line{frac_line_0} In Surface{domain_surf};
Physical Line("FRACTURE_4") = { frac_line_0 };

// End of /compartment boundary/auxiliary elements specification

// Start physical point specification
// End of physical point specification

