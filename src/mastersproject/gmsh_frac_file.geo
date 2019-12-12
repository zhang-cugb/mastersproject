Geometry.Tolerance = 1.0526315789473684e-10;
// Define points
p0 = newp; Point(p0) = {9.735, 88.36, 35.42, 10.0 };
p1 = newp; Point(p1) = {25.898207534469616, 95.81544745655496, 25.39682476739433, 10.0 };
p2 = newp; Point(p2) = {51.230000000000004, 107.49999999999999, 9.687999999999995, 10.0 };
p3 = newp; Point(p3) = {57.41893851709343, 111.57966494114152, 17.06669612262349, 10.0 };
p4 = newp; Point(p4) = {72.62, 121.6, 35.19, 10.0 };
p5 = newp; Point(p5) = {23.752748050532844, 103.62131116664274, 34.31506181035821, 10.0 };
p6 = newp; Point(p6) = {10.920000000000005, 95.62, 34.43, 10.0 };
p7 = newp; Point(p7) = {41.1, 109.2, 11.219999999999997, 10.0 };
p8 = newp; Point(p8) = {74.56, 135.3, 33.86, 10.0 };
p9 = newp; Point(p9) = {18.560000000000002, 107.69999999999999, 34.92, 10.0 };
p10 = newp; Point(p10) = {44.650000000000006, 117.29999999999998, 10.940000000000003, 10.0 };
p11 = newp; Point(p11) = {74.84, 143.3, 33.61, 10.0 };
p12 = newp; Point(p12) = {49.86, 100.6, 6.763999999999996, 10.0 };
p13 = newp; Point(p13) = {72.09, 106.6, 32.81, 10.0 };
p14 = newp; Point(p14) = {22.42, 113.2, 33.61, 10.0 };
p15 = newp; Point(p15) = {52.96, 105.5, 12.279999999999998, 10.0 };
p16 = newp; Point(p16) = {72.19, 110.0, 33.64, 10.0 };
p17 = newp; Point(p17) = {25.13, 118.2, 33.76, 10.0 };
p18 = newp; Point(p18) = {-6.0, 55.0, 50.0, 10.0 };
p19 = newp; Point(p19) = {-6.0, 55.0, 0.0, 10.0 };
p20 = newp; Point(p20) = {-6.0, 150.0, 0.0, 10.0 };
p21 = newp; Point(p21) = {-6.0, 150.0, 50.0, 10.0 };
p22 = newp; Point(p22) = {80.0, 55.0, 50.0, 10.0 };
p23 = newp; Point(p23) = {80.0, 55.0, 0.0, 10.0 };
p24 = newp; Point(p24) = {80.0, 150.0, 0.0, 10.0 };
p25 = newp; Point(p25) = {80.0, 150.0, 50.0, 10.0 };
p26 = newp; Point(p26) = {43.616395750791206, 103.988145913246, 14.409370395002789, 10.0 };
p27 = newp; Point(p27) = {50.184695990470104, 109.51071082067439, 33.16281343280902, 10.0 };
p28 = newp; Point(p28) = {56.252515058537654, 112.77703732511668, 33.68063957061146, 10.0 };
p29 = newp; Point(p29) = {49.41051368350976, 107.11977995757908, 15.019596337700687, 10.0 };
p30 = newp; Point(p30) = {36.42906488116287, 111.33851765219096, 33.384365776023145, 10.0 };
p31 = newp; Point(p31) = {34.695139413689105, 107.56345639167337, 21.600583356417722, 10.0 };
p32 = newp; Point(p32) = {42.45182660558299, 115.18174717029791, 33.715830446394605, 10.0 };
p33 = newp; Point(p33) = {40.415697558854106, 111.224493029197, 21.962055926547386, 10.0 };
p34 = newp; Point(p34) = {26.874103352109337, 112.60815216178938, 33.53826086809568, 10.0 };
p35 = newp; Point(p35) = {26.642490825852114, 111.26110115139444, 29.478914405582145, 10.0 };
p36 = newp; Point(p36) = {33.30004223657774, 116.77640573013308, 33.739166913123896, 10.0 };
p37 = newp; Point(p37) = {32.783317840667614, 114.70746904144885, 27.852947638607965, 10.0 };
// End of point specification

// Define lines 
frac_line_0= newl; Line(frac_line_0) = {p0, p1};
Physical Line("FRACTURE_TIP_0") = {frac_line_0};

frac_line_1= newl; Line(frac_line_1) = {p0, p4};
Physical Line("FRACTURE_TIP_1") = {frac_line_1};

frac_line_2= newl; Line(frac_line_2) = {p1, p26};
Physical Line("FRACTURE_TIP_2") = {frac_line_2};

frac_line_3= newl; Line(frac_line_3) = {p2, p3};
Physical Line("FRACTURE_TIP_3") = {frac_line_3};

frac_line_4= newl; Line(frac_line_4) = {p2, p26};
Physical Line("FRACTURE_TIP_4") = {frac_line_4};

frac_line_5= newl; Line(frac_line_5) = {p3, p4};
Physical Line("FRACTURE_TIP_5") = {frac_line_5};

frac_line_6= newl; Line(frac_line_6) = {p5, p6};
Physical Line("FRACTURE_TIP_6") = {frac_line_6};

frac_line_7= newl; Line(frac_line_7) = {p5, p8};
Physical Line("FRACTURE_TIP_7") = {frac_line_7};

frac_line_8= newl; Line(frac_line_8) = {p6, p7};
Physical Line("FRACTURE_TIP_8") = {frac_line_8};

frac_line_9= newl; Line(frac_line_9) = {p7, p8};
Physical Line("FRACTURE_TIP_9") = {frac_line_9};

frac_line_10= newl; Line(frac_line_10) = {p9, p10};
Physical Line("FRACTURE_TIP_10") = {frac_line_10};

frac_line_11= newl; Line(frac_line_11) = {p9, p11};
Physical Line("FRACTURE_TIP_11") = {frac_line_11};

frac_line_12= newl; Line(frac_line_12) = {p10, p11};
Physical Line("FRACTURE_TIP_12") = {frac_line_12};

frac_line_13= newl; Line(frac_line_13) = {p12, p13};
Physical Line("FRACTURE_TIP_13") = {frac_line_13};

frac_line_14= newl; Line(frac_line_14) = {p12, p31};
Physical Line("FRACTURE_TIP_14") = {frac_line_14};

frac_line_15= newl; Line(frac_line_15) = {p13, p27};
Physical Line("FRACTURE_TIP_15") = {frac_line_15};

frac_line_16= newl; Line(frac_line_16) = {p14, p34};
Physical Line("FRACTURE_TIP_16") = {frac_line_16};

frac_line_17= newl; Line(frac_line_17) = {p14, p35};
Physical Line("FRACTURE_TIP_17") = {frac_line_17};

frac_line_18= newl; Line(frac_line_18) = {p15, p16};
Physical Line("FRACTURE_TIP_18") = {frac_line_18};

frac_line_19= newl; Line(frac_line_19) = {p15, p29};
Physical Line("FRACTURE_TIP_19") = {frac_line_19};

frac_line_20= newl; Line(frac_line_20) = {p16, p28};
Physical Line("FRACTURE_TIP_20") = {frac_line_20};

frac_line_21= newl; Line(frac_line_21) = {p17, p36};
Physical Line("FRACTURE_TIP_21") = {frac_line_21};

frac_line_22= newl; Line(frac_line_22) = {p17, p37};
Physical Line("FRACTURE_TIP_22") = {frac_line_22};

frac_line_23= newl; Line(frac_line_23) = {p18, p19};
Physical Line("AUXILIARY_LINE_23") = {frac_line_23};

frac_line_24= newl; Line(frac_line_24) = {p18, p21};
Physical Line("AUXILIARY_LINE_24") = {frac_line_24};

frac_line_25= newl; Line(frac_line_25) = {p18, p22};
Physical Line("AUXILIARY_LINE_25") = {frac_line_25};

frac_line_26= newl; Line(frac_line_26) = {p19, p20};
Physical Line("AUXILIARY_LINE_26") = {frac_line_26};

frac_line_27= newl; Line(frac_line_27) = {p19, p23};
Physical Line("AUXILIARY_LINE_27") = {frac_line_27};

frac_line_28= newl; Line(frac_line_28) = {p20, p21};
Physical Line("AUXILIARY_LINE_28") = {frac_line_28};

frac_line_29= newl; Line(frac_line_29) = {p20, p24};
Physical Line("AUXILIARY_LINE_29") = {frac_line_29};

frac_line_30= newl; Line(frac_line_30) = {p21, p25};
Physical Line("AUXILIARY_LINE_30") = {frac_line_30};

frac_line_31= newl; Line(frac_line_31) = {p22, p23};
Physical Line("AUXILIARY_LINE_31") = {frac_line_31};

frac_line_32= newl; Line(frac_line_32) = {p22, p25};
Physical Line("AUXILIARY_LINE_32") = {frac_line_32};

frac_line_33= newl; Line(frac_line_33) = {p23, p24};
Physical Line("AUXILIARY_LINE_33") = {frac_line_33};

frac_line_34= newl; Line(frac_line_34) = {p24, p25};
Physical Line("AUXILIARY_LINE_34") = {frac_line_34};

frac_line_35= newl; Line(frac_line_35) = {p26, p27};
Physical Line("FRACTURE_LINE_35") = {frac_line_35};

frac_line_36= newl; Line(frac_line_36) = {p27, p30};
Physical Line("FRACTURE_TIP_36") = {frac_line_36};

frac_line_37= newl; Line(frac_line_37) = {p28, p29};
Physical Line("FRACTURE_LINE_37") = {frac_line_37};

frac_line_38= newl; Line(frac_line_38) = {p28, p32};
Physical Line("FRACTURE_TIP_38") = {frac_line_38};

frac_line_39= newl; Line(frac_line_39) = {p29, p33};
Physical Line("FRACTURE_TIP_39") = {frac_line_39};

frac_line_40= newl; Line(frac_line_40) = {p30, p31};
Physical Line("FRACTURE_LINE_40") = {frac_line_40};

frac_line_41= newl; Line(frac_line_41) = {p30, p34};
Physical Line("FRACTURE_TIP_41") = {frac_line_41};

frac_line_42= newl; Line(frac_line_42) = {p31, p35};
Physical Line("FRACTURE_TIP_42") = {frac_line_42};

frac_line_43= newl; Line(frac_line_43) = {p32, p33};
Physical Line("FRACTURE_LINE_43") = {frac_line_43};

frac_line_44= newl; Line(frac_line_44) = {p32, p36};
Physical Line("FRACTURE_TIP_44") = {frac_line_44};

frac_line_45= newl; Line(frac_line_45) = {p33, p37};
Physical Line("FRACTURE_TIP_45") = {frac_line_45};

frac_line_46= newl; Line(frac_line_46) = {p34, p35};
Physical Line("FRACTURE_LINE_46") = {frac_line_46};

frac_line_47= newl; Line(frac_line_47) = {p36, p37};
Physical Line("FRACTURE_LINE_47") = {frac_line_47};

// End of line specification 

// Start domain specification
frac_loop_5 = newll; 
Line Loop(frac_loop_5) = { frac_line_23, frac_line_26, frac_line_28, -frac_line_24};
auxiliary_5 = news; Plane Surface(auxiliary_5) = {frac_loop_5};
Physical Surface("AUXILIARY_5") = {auxiliary_5};

frac_loop_6 = newll; 
Line Loop(frac_loop_6) = { frac_line_31, frac_line_33, frac_line_34, -frac_line_32};
auxiliary_6 = news; Plane Surface(auxiliary_6) = {frac_loop_6};
Physical Surface("AUXILIARY_6") = {auxiliary_6};

frac_loop_7 = newll; 
Line Loop(frac_loop_7) = { frac_line_23, frac_line_27, -frac_line_31, -frac_line_25};
auxiliary_7 = news; Plane Surface(auxiliary_7) = {frac_loop_7};
Physical Surface("AUXILIARY_7") = {auxiliary_7};

frac_loop_8 = newll; 
Line Loop(frac_loop_8) = { frac_line_28, frac_line_30, -frac_line_34, -frac_line_29};
auxiliary_8 = news; Plane Surface(auxiliary_8) = {frac_loop_8};
Physical Surface("AUXILIARY_8") = {auxiliary_8};

frac_loop_9 = newll; 
Line Loop(frac_loop_9) = { frac_line_26, frac_line_29, -frac_line_33, -frac_line_27};
auxiliary_9 = news; Plane Surface(auxiliary_9) = {frac_loop_9};
Physical Surface("AUXILIARY_9") = {auxiliary_9};

frac_loop_10 = newll; 
Line Loop(frac_loop_10) = { frac_line_24, frac_line_30, -frac_line_32, -frac_line_25};
auxiliary_10 = news; Plane Surface(auxiliary_10) = {frac_loop_10};
Physical Surface("AUXILIARY_10") = {auxiliary_10};

domain_loop = newsl;
Surface Loop(domain_loop) = {auxiliary_5,auxiliary_6,auxiliary_7,auxiliary_8,auxiliary_9,auxiliary_10};
Volume(1) = {domain_loop};
Physical Volume("DOMAIN") = {1};
// End of domain specification

// Start fracture specification
frac_loop_0 = newll; 
Line Loop(frac_loop_0) = { frac_line_0, frac_line_2, -frac_line_4, frac_line_3, frac_line_5, -frac_line_1};
fracture_0 = news; Plane Surface(fracture_0) = {frac_loop_0};
Physical Surface("FRACTURE_0") = {fracture_0};
Surface{fracture_0} In Volume{1};

Line{frac_line_35} In Surface{fracture_0};
Line{frac_line_37} In Surface{fracture_0};

frac_loop_1 = newll; 
Line Loop(frac_loop_1) = { frac_line_6, frac_line_8, frac_line_9, -frac_line_7};
fracture_1 = news; Plane Surface(fracture_1) = {frac_loop_1};
Physical Surface("FRACTURE_1") = {fracture_1};
Surface{fracture_1} In Volume{1};

Line{frac_line_40} In Surface{fracture_1};
Line{frac_line_43} In Surface{fracture_1};

frac_loop_2 = newll; 
Line Loop(frac_loop_2) = { frac_line_10, frac_line_12, -frac_line_11};
fracture_2 = news; Plane Surface(fracture_2) = {frac_loop_2};
Physical Surface("FRACTURE_2") = {fracture_2};
Surface{fracture_2} In Volume{1};

Line{frac_line_46} In Surface{fracture_2};
Line{frac_line_47} In Surface{fracture_2};

frac_loop_3 = newll; 
Line Loop(frac_loop_3) = { frac_line_13, frac_line_15, frac_line_36, frac_line_41, -frac_line_16, frac_line_17, -frac_line_42, -frac_line_14};
fracture_3 = news; Plane Surface(fracture_3) = {frac_loop_3};
Physical Surface("FRACTURE_3") = {fracture_3};
Surface{fracture_3} In Volume{1};

Line{frac_line_35} In Surface{fracture_3};
Line{frac_line_40} In Surface{fracture_3};
Line{frac_line_46} In Surface{fracture_3};

frac_loop_4 = newll; 
Line Loop(frac_loop_4) = { frac_line_18, frac_line_20, frac_line_38, frac_line_44, -frac_line_21, frac_line_22, -frac_line_45, -frac_line_39, -frac_line_19};
fracture_4 = news; Plane Surface(fracture_4) = {frac_loop_4};
Physical Surface("FRACTURE_4") = {fracture_4};
Surface{fracture_4} In Volume{1};

Line{frac_line_37} In Surface{fracture_4};
Line{frac_line_43} In Surface{fracture_4};
Line{frac_line_47} In Surface{fracture_4};

// End of fracture specification

// Start physical point specification
// End of physical point specification

