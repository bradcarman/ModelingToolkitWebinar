model Workshop
	parameter Real res1_C_0 = 2.7;
	parameter Real res1_A_0 = 0.00094;
	parameter Real res1_rho_0 = 1000;
	parameter Real res2_C_0 = 2.7;
	parameter Real res2_A_0 = 0.00094;
	parameter Real res2_rho_0 = 1000;
	parameter Real act_A = 0.1;
	parameter Real act_vol1_A = 0.1;
	parameter Real act_vol1_rho_0 = 1000;
	parameter Real act_vol1_beta = 2.0e9;
	parameter Real act_vol1_direction = -1;
	parameter Real act_vol1_L = 0.5;
	parameter Real act_vol2_A = 0.1;
	parameter Real act_vol2_rho_0 = 1000;
	parameter Real act_vol2_beta = 2.0e9;
	parameter Real act_vol2_direction = 1;
	parameter Real act_vol2_L = 0.5;
	parameter Real act_mass_m = 100;
	parameter Real src_p = 3.0e7;
	parameter Real snk_p = 0;
	parameter Real dmp_c = 1000;
	Real act_mass_x(start = 0.0);
	Real act_vol1_x(start = 0.0);
	Real act_vol1_m(start = 50.74999999999999);
	Real act_vol2_x(start = 0.0);
	Real act_vol2_m(start = 50.0);
	Real act_mass_dx(start = 0.0);
	Real act_vol1_r_t(start = 0.0);
	Real act_vol2_r_t(start = 0.0);
	Real act_vol1_r(start = 1014.9999999999999);
	Real act_vol2_r(start = 1000.0);
	
initial equation
	act_mass_x = 0;
    der(act_mass_x) = 0;
    
    act_vol1_x = 0;
    der(act_vol1_m) = 0;
    
    act_vol2_x = 0;
    der(act_vol2_m) = 0;

equation
	der(act_mass_x) = act_mass_dx;
	der(act_vol1_x) = act_vol1_direction*act_mass_dx;
	der(act_vol1_m) = (2*res1_A_0*(src_p + (act_vol1_beta*(act_vol1_rho_0 - act_vol1_r)) / act_vol1_rho_0)) / (res1_C_0*((0.0001 + (4*((src_p + (act_vol1_beta*(act_vol1_rho_0 - act_vol1_r)) / act_vol1_rho_0)^2)) / ((res1_C_0^2)*(res1_rho_0^2)))^0.25));
	der(act_vol2_x) = act_vol2_direction*act_mass_dx;
	der(act_vol2_m) = (-2*res2_A_0*(-snk_p + (-act_vol2_beta*(act_vol2_rho_0 - act_vol2_r)) / act_vol2_rho_0)) / (res2_C_0*((0.0001 + (4*((-snk_p + (-act_vol2_beta*(act_vol2_rho_0 - act_vol2_r)) / act_vol2_rho_0)^2)) / ((res2_C_0^2)*(res2_rho_0^2)))^0.25));
	der(act_mass_dx) = ((-act_vol2_A*act_vol2_beta*(act_vol2_rho_0 - act_vol2_r)) / (act_vol2_direction*act_vol2_rho_0) + (-act_vol1_A*act_vol1_beta*(act_vol1_rho_0 - act_vol1_r)) / (act_vol1_direction*act_vol1_rho_0) - dmp_c*act_mass_dx) / act_mass_m;
	0 = -act_vol1_m + act_vol1_A*(act_vol1_L + act_vol1_x)*act_vol1_r;
	0 = -act_vol2_m + act_vol2_A*(act_vol2_L + act_vol2_x)*act_vol2_r;
	0 = (-2*res1_A_0*(src_p + (act_vol1_beta*(act_vol1_rho_0 - act_vol1_r)) / act_vol1_rho_0)) / (res1_C_0*((0.0001 + (4*((src_p + (act_vol1_beta*(act_vol1_rho_0 - act_vol1_r)) / act_vol1_rho_0)^2)) / ((res1_C_0^2)*(res1_rho_0^2)))^0.25)) + act_vol1_A*(act_vol1_L + act_vol1_x)*act_vol1_r_t + act_vol1_A*act_vol1_direction*act_mass_dx*act_vol1_r;
	0 = (2*res2_A_0*(-snk_p + (-act_vol2_beta*(act_vol2_rho_0 - act_vol2_r)) / act_vol2_rho_0)) / (res2_C_0*((0.0001 + (4*((-snk_p + (-act_vol2_beta*(act_vol2_rho_0 - act_vol2_r)) / act_vol2_rho_0)^2)) / ((res2_C_0^2)*(res2_rho_0^2)))^0.25)) + act_vol2_A*(act_vol2_L + act_vol2_x)*act_vol2_r_t + act_vol2_A*act_vol2_direction*act_mass_dx*act_vol2_r;
end Workshop;
