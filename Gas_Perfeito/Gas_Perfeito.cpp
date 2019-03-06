#include <stdio.h>
#include <math.h>

int main() {
	double R, Rmol, Mmol, Tr, Pr, Tc, Pc;
	double T, P, v, vr, vnew1, vnew2, vold;
	double Ao, a, Bo, b, c,Pbb, Rbb;
	double Z ,w, f_acent;
	double Bref, Cref, Dref, Zref, b1ref, b2ref, b3ref, b4ref, c1ref;
	double c2ref, c3ref, c4ref, d1ref, d2ref, betaref, gamaref;
	double Bsim, Csim, Dsim, Zsim, b1sim, b2sim, b3sim, b4sim, c1sim;
	double c2sim, c3sim, c4sim, d1sim, d2sim, betasim, gamasim;
	double check, precisao, count;

ini:
	// Dados 
	//Universais
	precisao = 10E-6;
	check = precisao * 100;
	Rmol = 8314.5; //J/kmol*K
	Rbb = 0.08206; // atm*l/mol*K
	count = 0;
	
	//Constantes de Lee-Kesler fluido simples
	b1sim = 0.1181193;
	b2sim = 0.265728;
	b3sim = 0.154790;
	b4sim = 0.030323;
	c1sim = 0.0236744;
	c2sim = 0.0186984;
	c3sim = 0.0;
	c4sim = 0.042724;
	d1sim = 0.155488E-4;
	d2sim = 0.623689E-4;
	betasim = 0.65392;
	gamasim = 0.060167;

	//Constantes de Lee-Kesler de referência
	b1ref = 0.2026579;
	b2ref = 0.331511;
	b3ref = 0.027655;
	b4ref = 0.203488;
	c1ref = 0.0313385;
	c2ref = 0.0503618;
	c3ref = 0.016901;
	c4ref = 0.041577;
	d1ref = 0.48736E-4;
	d2ref = 0.07040336E-4;
	betaref = 1.226;
	gamaref = 0.03754;
	
	printf("\nEscolha o gas \n1 - Amonia \n2 - Propano \n:");
	scanf("%lf", &T);
	
	if (T==1) {
		// Amônia
		Tc = 405.5; // K
		Pc = 11.35E6;  // Pa
		Mmol = 17.031; // g/mol
		Ao = 2.3930;
		a = 0.17031;
		Bo = 0.03415;
		b = 0.19112;
		c = 476.87E4;
		f_acent = 0.250;
	}
	
	else {
		// Propano
		Tc = 369.8; // K
		Pc = 4.257E6;  // Pa
		Mmol = 44.097; // g/mol
		Ao = 11.9200;
		a = 0.07321;
		Bo = 0.18100;
		b = 0.04293;
		c = 120E4;
		f_acent = 0.153;
	}

	printf("\nEntre com o valor de T em K: ");
	scanf("%lf", &T);
	printf("Entre com o valor de P em bar: ");
	scanf("%lf", &P);

	P = P *1E5; //Pa
	R = Rmol / Mmol; // J/kgK
	Pbb = P/101325; //atm
	Tr = T / Tc;
	Pr = P / Pc;

	// Gás Perfeito 

	v = R*T / P;
	printf("\nGas Perfeito: v = %lf m3/kg\n", v);

	// Beattie-Bridgeman
	vold = v*Mmol;
	while (check > precisao) {
		
		vnew1 = (Rbb*T*(vold + Bo*(1 - b / vold))*(1 - c / (vold*pow(T, 3))) - Ao*(1 - a / vold)) / (Pbb*vold);
		vnew2 = 0.1*vnew1 + 0.9*vold;

		if (vold > vnew2) {
			check = vold - vnew2;
		}
		else {
			check = vnew2 - vold;
		}

		vold = vnew2;
	}
	v = vold / Mmol;
	printf("Beattie-Bridgeman: v = %lf m3/kg\n", v);

	// Lee-Kesler
	
	Bref=b1ref-b2ref/Tr-b3ref/pow(Tr,2)-b4ref/pow(Tr,3);
	Cref=c1ref-c2ref/Tr+c3ref/pow(Tr,2);
	Dref=d1ref+d2ref/Tr;
	
	Bsim=b1sim-b2sim/Tr-b3sim/pow(Tr,2)-b4sim/pow(Tr,3);
	Csim=c1sim-c2sim/Tr+c3sim/pow(Tr,2);
	Dsim=d1sim+d2sim/Tr;
	
	check = precisao*100;
	while (check > precisao) {
		
		Zref = 1 + Bref / vold + Cref / pow(vold, 2) + Dref / pow(vold, 3) + c4ref*(betaref + gamaref / pow(vold, 2))*exp(-gamaref / pow(vold, 2)) / (pow(Tr, 3)*pow(vold, 2));
		Zsim = 1 + Bsim / vold + Csim / pow(vold, 2) + Dsim / pow(vold, 3) + c4sim*(betasim + gamasim / pow(vold, 2))*exp(-gamasim / pow(vold, 2)) / (pow(Tr, 3)*pow(vold, 2));
		Z=Zsim+f_acent*(Zref-Zsim)/0.3978;
		vnew1 = Z*Tr / Pr;
		vnew2 = 0.1*vnew1 + 0.9*vold;

		if (vold > vnew2) {
			check = vold - vnew2;
		}
		else {
			check = vnew2 - vold;
		}
		vold = vnew2;
	}
	v = vold*R*Tc/Pc;
	printf("Lee-Kesler: v = %lf m3/kg\n", v);

	printf("\n\nEntre com 1 para parar, ou para 0 continuar: "); 
	scanf("%lf", &P);
	if (P == 1) return 0;
	goto ini;
	return 0;
}
