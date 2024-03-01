model MTK
	parameter Real rho_0 = 1000;
	parameter Real beta = 2.0e9;
	parameter Real A = 0.1;
	parameter Real m = 100;
	parameter Real L = 1;
	parameter Real p_s = 1.0e7;
	parameter Real p_r = 1.0e6;
	parameter Real C = 1.35;
	parameter Real c = 1000;
	parameter Real A_p = 0.00094;
	Real x(start = 0);
	Real dx(start = 0);
	Real rho_1(start = 1004.9999999999999);
	Real rho_2(start = 1000.5);
	Real drho_1(start = 0);
	Real drho_2(start = 0);
	Real dm_1(start = 0);
	Real dm_2(start = 0);
equation
	der(x) = dx;
	der(dx) = (A*((-beta*(-rho_0 + rho_2)) / rho_0 + (beta*(-rho_0 + rho_1)) / rho_0) - c*dx) / m;
	der(rho_1) = drho_1;
	der(rho_2) = drho_2;
	0 = -dm_1 + (L + x)*drho_1 + dx*rho_1;
	0 = dm_2 + (L - x)*drho_2 - dx*rho_2;
	0 = -p_s + (beta*(-rho_0 + rho_1)) / rho_0 + C*rho_0*((dm_1 / (A_p*rho_0))^2);
	0 = p_r + (-beta*(-rho_0 + rho_2)) / rho_0 + C*rho_0*((dm_2 / (A_p*rho_0))^2);
end MTK;
