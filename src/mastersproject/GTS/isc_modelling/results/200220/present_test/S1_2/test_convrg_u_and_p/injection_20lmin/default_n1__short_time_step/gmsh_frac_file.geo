Geometry.Tolerance = 1e-10;
// Define points
p0 = newp; Point(p0) = {10.869386500225858, 95.69772585839029, 34.410362762094856, 10.0 };
p1 = newp; Point(p1) = {75.47721099067371, 133.76440049383834, 34.25338188367685, 10.0 };
p2 = newp; Point(p2) = {53.895135783950764, 114.4752293110289, 8.594096929400413, 8.721907427899884 };
p3 = newp; Point(p3) = {48.42258592286497, 109.95228536778566, 3.5279499432713717, 8.721907427899884 };
p4 = newp; Point(p4) = {-20.0, 50.0, 75.0, 10.0 };
p5 = newp; Point(p5) = {-20.0, 50.0, -25.0, 10.0 };
p6 = newp; Point(p6) = {-20.0, 150.0, -25.0, 10.0 };
p7 = newp; Point(p7) = {-20.0, 150.0, 75.0, 10.0 };
p8 = newp; Point(p8) = {80.0, 50.0, 75.0, 10.0 };
p9 = newp; Point(p9) = {80.0, 50.0, -25.0, 10.0 };
p10 = newp; Point(p10) = {80.0, 150.0, -25.0, 10.0 };
p11 = newp; Point(p11) = {80.0, 150.0, 75.0, 10.0 };
// End of point specification

// Define lines 
frac_line_0= newl; Line(frac_line_0) = {p0, p1};
Physical Line("FRACTURE_TIP_0") = {frac_line_0};

frac_line_1= newl; Line(frac_line_1) = {p0, p3};
Physical Line("FRACTURE_TIP_1") = {frac_line_1};

frac_line_2= newl; Line(frac_line_2) = {p1, p2};
Physical Line("FRACTURE_TIP_2") = {frac_line_2};

frac_line_3= newl; Line(frac_line_3) = {p2, p3};
Physical Line("FRACTURE_TIP_3") = {frac_line_3};

frac_line_4= newl; Line(frac_line_4) = {p4, p5};
Physical Line("AUXILIARY_LINE_4") = {frac_line_4};

frac_line_5= newl; Line(frac_line_5) = {p4, p7};
Physical Line("AUXILIARY_LINE_5") = {frac_line_5};

frac_line_6= newl; Line(frac_line_6) = {p4, p8};
Physical Line("AUXILIARY_LINE_6") = {frac_line_6};

frac_line_7= newl; Line(frac_line_7) = {p5, p6};
Physical Line("AUXILIARY_LINE_7") = {frac_line_7};

frac_line_8= newl; Line(frac_line_8) = {p5, p9};
Physical Line("AUXILIARY_LINE_8") = {frac_line_8};

frac_line_9= newl; Line(frac_line_9) = {p6, p7};
Physical Line("AUXILIARY_LINE_9") = {frac_line_9};

frac_line_10= newl; Line(frac_line_10) = {p6, p10};
Physical Line("AUXILIARY_LINE_10") = {frac_line_10};

frac_line_11= newl; Line(frac_line_11) = {p7, p11};
Physical Line("AUXILIARY_LINE_11") = {frac_line_11};

frac_line_12= newl; Line(frac_line_12) = {p8, p9};
Physical Line("AUXILIARY_LINE_12") = {frac_line_12};

frac_line_13= newl; Line(frac_line_13) = {p8, p11};
Physical Line("AUXILIARY_LINE_13") = {frac_line_13};

frac_line_14= newl; Line(frac_line_14) = {p9, p10};
Physical Line("AUXILIARY_LINE_14") = {frac_line_14};

frac_line_15= newl; Line(frac_line_15) = {p10, p11};
Physical Line("AUXILIARY_LINE_15") = {frac_line_15};

// End of line specification 

// Start domain specification
frac_loop_1 = newll; 
Line Loop(frac_loop_1) = { frac_line_4, frac_line_7, frac_line_9, -frac_line_5};
auxiliary_1 = news; Plane Surface(auxiliary_1) = {frac_loop_1};
Physical Surface("AUXILIARY_1") = {auxiliary_1};

frac_loop_2 = newll; 
Line Loop(frac_loop_2) = { frac_line_12, frac_line_14, frac_line_15, -frac_line_13};
auxiliary_2 = news; Plane Surface(auxiliary_2) = {frac_loop_2};
Physical Surface("AUXILIARY_2") = {auxiliary_2};

frac_loop_3 = newll; 
Line Loop(frac_loop_3) = { frac_line_4, frac_line_8, -frac_line_12, -frac_line_6};
auxiliary_3 = news; Plane Surface(auxiliary_3) = {frac_loop_3};
Physical Surface("AUXILIARY_3") = {auxiliary_3};

frac_loop_4 = newll; 
Line Loop(frac_loop_4) = { frac_line_9, frac_line_11, -frac_line_15, -frac_line_10};
auxiliary_4 = news; Plane Surface(auxiliary_4) = {frac_loop_4};
Physical Surface("AUXILIARY_4") = {auxiliary_4};

frac_loop_5 = newll; 
Line Loop(frac_loop_5) = { frac_line_7, frac_line_10, -frac_line_14, -frac_line_8};
auxiliary_5 = news; Plane Surface(auxiliary_5) = {frac_loop_5};
Physical Surface("AUXILIARY_5") = {auxiliary_5};

frac_loop_6 = newll; 
Line Loop(frac_loop_6) = { frac_line_5, frac_line_11, -frac_line_13, -frac_line_6};
auxiliary_6 = news; Plane Surface(auxiliary_6) = {frac_loop_6};
Physical Surface("AUXILIARY_6") = {auxiliary_6};

domain_loop = newsl;
Surface Loop(domain_loop) = {auxiliary_1,auxiliary_2,auxiliary_3,auxiliary_4,auxiliary_5,auxiliary_6};
Volume(1) = {domain_loop};
Physical Volume("DOMAIN") = {1};
// End of domain specification

// Start fracture specification
frac_loop_0 = newll; 
Line Loop(frac_loop_0) = { frac_line_0, frac_line_2, frac_line_3, -frac_line_1};
fracture_0 = news; Plane Surface(fracture_0) = {frac_loop_0};
Physical Surface("FRACTURE_0") = {fracture_0};
Surface{fracture_0} In Volume{1};


// End of fracture specification

// Start physical point specification
// End of physical point specification

