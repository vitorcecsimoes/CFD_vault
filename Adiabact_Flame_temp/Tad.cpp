# include <stdio.h>
# include <math.h>

//Declaração de constantes
const double Rcal=1.987; //kcal/(mol*K)
const double Rmol=8.314; //kJ/(kmol*K)

const double hfCO=-110530; //kj/kmol [nist.gov]
const double hfCO2=-393522.4; //kj/kmol [nist.gov]
const double hfH2O=-241826.4; //kj/kmol [nist.gov]

const double hfCH4=-74873.1; //kJ/kmol [nist.gov]
const double hfC2H2=226731.4; //kJ/kmol [nist.gov]
const double hfC2H4=-52467; //kJ/kmol [nist.gov]
const double hfC2H6=-83846; //kJ/kmol [NASA]
const double hfC3H8=-104700; //kJ/kmol [nist.gov]
const double hfC4H10=-125782; //kJ/kmol [NASA]
const double hfC5H12=-146753; //kJ/kmol [NASA]

//cabeçalhos das funções
int CalcComp(double phi, double nC, double nH, double *nO2real, double *nCO2real, double *nH2Oreal, double *nexcomb, double *nexO2);
int DeltahComb(int select, double T,double  *dhcomb);
int PropertyO2 (char prop,double T);
int PropertyN2(char prop,double T);
int PropertyCO(char prop,double T);
int PropertyCO2(char prop,double T);
int PropertyH2(char prop,double T);
int PropertyH2O(char prop,double T);
int KpH2OCO2 (double T, double *Kp);


//inicialização da variável de arquivo
FILE *file;


int main (){
//======Declaração de variáveis=======//
	double precision, dif, TadP, TadV, TadPdiss, TadVdiss, Tadnew, fn;
	double A, B, D, lower, upper;
	double phi, Ti, Tad, dG, Kp;
	double hfreact, Hreact, Hprod;
	double dhCO, dhCO2, dhH2, dhH2O, dhO2, dhN2, dhreact;
	double nC, nCO, nCO2, nH, nH2, nH2O, nN2, nO2, nexO2, ncomb, nexcomb, Nreact, Nprod;
	double sCO, sCO2, sH2, sH2O, sN2, sO2;
	
	int select,loop;


//loop principal, para cada valor de select um combustível é utilizado
	for (select=1; select<=7; select++){
		
		precision = 1E-3;
		Ti=298; //K
		ncomb=1; //mol
		
		//fecha o arquivo anteriormente aberto
		if (select>1){
			fclose (file);
		}
		
		//inicialização dos arquivos, guardando cada combustível no seu respectivo arquivo
		//os arquivos ficarão localizados na mesma pasta que o programa
		if (select==1){
			file = fopen("CH4.txt","w");
		}
		else if (select==2){
			file = fopen("C2H2.txt","w");
		}
		else if (select==3){
			file = fopen("C2H4.txt","w");
		}
		else if (select==4){
			file = fopen("C2H6.txt","w");
		}
		else if (select==5){
			file = fopen("C3H8.txt","w");
		}
		else if (select==6){
			file = fopen("C4H10.txt","w");
		}
		else if (select==7){
			file = fopen("C5H12.txt","w");
		}
		
		fprintf(file,"phi	TadV	TadVdiss	TadP	TadPdiss\n");
		
		for (phi = 0.25; phi<2.2; phi=phi+0.01){

			if (select==1){
			//combustível = metano = CH4
				nC=ncomb;
				nH=4*ncomb;
				CalcComp(phi, nC, nH, &nO2, &nCO2, &nH2O, &nexcomb, &nexO2); //cálculo da composição dos gases de entrada e saída
				hfreact=hfCH4;
				DeltahComb(select,Ti,&dhreact);
			}
			if (select==2){
			//combustível = etino = C2H2
				nC=2*ncomb;
				nH=2*ncomb;
				CalcComp(phi, nC, nH, &nO2, &nCO2, &nH2O, &nexcomb, &nexO2);
				hfreact=hfC2H2;
				DeltahComb(select,Ti,&dhreact);
			}
			
			if (select==3){
			//combustível = eteno = C2H4
				nC=2*ncomb;
				nH=4*ncomb;
				CalcComp(phi, nC, nH, &nO2, &nCO2, &nH2O, &nexcomb, &nexO2);
				hfreact=hfC2H4;
				DeltahComb(select,Ti,&dhreact);
			}
			if (select==4){
			//combustível = etano = C2H6
				nC=2*ncomb;
				nH=6*ncomb;
				CalcComp(phi, nC, nH, &nO2, &nCO2, &nH2O, &nexcomb, &nexO2);
				hfreact=hfC2H6;
				DeltahComb(select,Ti,&dhreact);
			}
			
			if (select==5){
			//combustível = propano = C3H8
				nC=3*ncomb;
				nH=8*ncomb;
				CalcComp(phi, nC, nH, &nO2, &nCO2, &nH2O, &nexcomb, &nexO2);
				hfreact=hfC3H8;
				DeltahComb(select,Ti,&dhreact);
			}
			if (select==6){
			//combustível = butano = C4H10
				nC=4*ncomb;
				nH=10*ncomb;
				CalcComp(phi, nC, nH, &nO2, &nCO2, &nH2O, &nexcomb, &nexO2);
				hfreact=hfC4H10;
				DeltahComb(select,Ti,&dhreact);
			}
			if (select==7){
			//combustível = pentano = C5H12
				nC=5*ncomb;
				nH=12*ncomb;
				CalcComp(phi, nC, nH, &nO2, &nCO2, &nH2O, &nexcomb, &nexO2);
				hfreact=hfC5H12;
				DeltahComb(select,Ti,&dhreact);
			}
			
			
			//cálculos preliminares
					
			nN2=3.76*nO2;
			
			Nreact=ncomb+nO2+nN2;
						
			dhO2=PropertyO2('h', Ti);
			dhN2=PropertyN2('h', Ti);
			
			Hreact = ncomb*(hfreact+dhreact)+nO2*dhO2+nN2*dhN2; //kJ

//======V constante=======//
			
			//parâmetros de entrada
			dif=10;
			lower=500+Ti;
			upper=5000;
			Tad=(lower+upper)/2;
			
			while (dif>precision){
				
				//cálculo das entalpias
				dhO2=PropertyO2('h', Tad);
				dhN2=PropertyN2('h', Tad);
				dhCO2=PropertyCO2('h', Tad);
				dhH2O=PropertyH2O('h', Tad);
				DeltahComb(select,Tad,&dhreact);
				
				Hprod=nCO2*(hfCO2+dhCO2)+nH2O*(hfH2O+dhH2O)+nexO2*dhO2+nN2*dhN2+nexcomb*dhreact; //kJ
				
				//Avaliação da temperatura
				if(Hreact>Hprod){
					lower=Tad;
					Tadnew=(lower+upper)/2;
				}
				else{
					upper=Tad;
					Tadnew=(lower+upper)/2;
				}
				
				//Avaliação da variação da temperatura, se a temperatura varia muito pouco, então a temperatura alvo já foi atingida
				dif=fabs(Tad-Tadnew);
				
				Tad=Tadnew;
			
			}
			
			TadV = Tad;
			
			
//======P constante=======//
			
			dif=10;
			lower=500+Ti;
			upper=5000;
			Tad=(lower+upper)/2;
			
			while (dif>precision){
				
				dhO2=PropertyO2('h', Tad);
				dhN2=PropertyN2('h', Tad);
				dhCO2=PropertyCO2('h', Tad);
				dhH2O=PropertyH2O('h', Tad);
				DeltahComb(select,Tad,&dhreact);
				
				
				Hprod=nCO2*(hfCO2+dhCO2)+nH2O*(hfH2O+dhH2O)+nexO2*dhO2+nN2*dhN2+nexcomb*dhreact; //kJ
				
				Nprod=nCO2+nH2O+nN2+nexO2+nexcomb;
				fn=Hprod+Rmol*(Nreact*Ti-Nprod*Tad);
					
					if(Hreact>fn){
						lower=Tad;
						Tadnew=(lower+upper)/2;
					}
					else{
						upper=Tad;
						Tadnew=(lower+upper)/2;
					}
					dif=fabs(Tad-Tadnew);
					Tad=Tadnew;
			}
			TadP = Tad;
			
			TadPdiss=0;
			TadVdiss=0;

//======Com dissociação=======//
			
			if (phi>1){

//======V constante com dissociação=======//
				dif=10;
				lower=500+Ti;
				upper=TadV+1000;
				Tad=(lower+upper)/2;
				
				while (dif>precision){
					
					dhN2=PropertyN2('h', Tad);
					dhCO=PropertyCO('h', Tad);
					dhCO2=PropertyCO2('h', Tad);
					dhH2=PropertyH2('h', Tad);
					dhH2O=PropertyH2O('h', Tad);
					
					//Cálculo das entropias
					sCO=PropertyCO('s', Tad);
					sCO2=PropertyCO2('s', Tad);
					sH2=PropertyH2('s', Tad);
					sH2O=PropertyH2O('s', Tad);
					
					//Cálculo da dissociação, o sistema foi resolvido algebricamente
					dG=-(hfCO+dhCO)-(hfH2O+dhH2O)+(hfCO2+dhCO2)+dhH2-Tad*(-sCO-sH2O+sCO2+sH2);
					Kp=exp(-dG/(Rmol*Tad));
					
					A=1-Kp;
					B=nC*Kp+Kp*(2*nO2-nC)-2*nO2+nC+nH/2;
					D=-nC*Kp*(2*nO2-nC);
					
					nCO2=(-B+sqrt(pow(B,2)-4*A*D))/(2*A);
					nCO=nC-nCO2;
					nH2O=2*nO2-nC-nCO2;
					nH2=-2*nO2+nC+nH/2+nCO2;
					
					
					Hprod=nCO2*(hfCO2+dhCO2)+nCO*(hfCO+dhCO)+nH2O*(hfH2O+dhH2O)+nH2*dhH2+nN2*dhN2; //kJ
					
					if(Hreact>Hprod){
						lower=Tad;
						Tadnew=(lower+upper)/2;
					}
					else{
						upper=Tad;
						Tadnew=(lower+upper)/2;
					}
					dif=fabs(Tad-Tadnew);
					Tad=Tadnew;
				}
				
				
				TadVdiss = Tad;
				
//======P constante com dissociação=======//
				dif=10;
				lower=500+Ti;
				upper=TadP+1000;
				Tad=(lower+upper)/2;
				
				while (dif>precision){

					dhN2=PropertyN2('h', Tad);
					dhCO=PropertyCO('h', Tad);
					dhCO2=PropertyCO2('h', Tad);
					dhH2=PropertyH2('h', Tad);
					dhH2O=PropertyH2O('h', Tad);
					
					sCO=PropertyCO('s', Tad);
					sCO2=PropertyCO2('s', Tad);
					sH2=PropertyH2('s', Tad);
					sH2O=PropertyH2O('s', Tad);
					
					dG=-(hfCO+dhCO)-(hfH2O+dhH2O)+(hfCO2+dhCO2)+dhH2-Tad*(-sCO-sH2O+sCO2+sH2);
					Kp=exp(-dG/(Rmol*Tad));
					
					A=1-Kp;
					B=nC*Kp+Kp*(2*nO2-nC)-2*nO2+nC+nH/2;
					D=-nC*Kp*(2*nO2-nC);
					
					nCO2=(-B+sqrt(pow(B,2)-4*A*D))/(2*A);
					nCO=nC-nCO2;
					nH2O=2*nO2-nC-nCO2;
					nH2=-2*nO2+nC+nH/2+nCO2;
					
					
					Hprod=nCO2*(hfCO2+dhCO2)+nCO*(hfCO+dhCO)+nH2O*(hfH2O+dhH2O)+nH2*dhH2+nN2*dhN2; //kJ
					
					Nprod=nCO+nCO2+nH2+nH2O+nN2;
					fn=Hprod+Rmol*(Nreact*Ti-Nprod*Tad);
					
					if(Hreact>fn){
						lower=Tad;
						Tadnew=(lower+upper)/2;
					}
					else{
						upper=Tad;
						Tadnew=(lower+upper)/2;
					}
					
					dif=fabs(Tad-Tadnew);
					Tad=Tadnew;
					
				}
				TadPdiss = Tad;
			}
				
			fprintf(file,"%.2lf	%lf	%lf	%lf	%lf\n",phi, TadV, TadVdiss, TadP, TadPdiss);
		}
	}
	//fecha o último arquivo
	fclose (file);
	return 0;
}


int CalcComp(double phi, double nC, double nH, double *nO2real, double *nCO2real, double *nH2Oreal, double *nexcomb, double *nexO2){
	double nO2est;
	
	nO2est=nC+nH/4;
	
	*nO2real=nO2est/phi;
	
	if (phi>1){
		*nCO2real=nC/phi; 
		*nH2Oreal=nH/(2*phi);
		*nexcomb=1-1/phi;
		*nexO2=0.0;
	}
			
	else {
		*nCO2real=nC;
		*nH2Oreal=nH/2;
		*nexO2=nO2est-nC-nH/4;
		*nexcomb=0;
	}
}

int DeltahComb(int select, double T,double  *dhcomb) {
    double A,B,C,D,E,F,G,H,theta, theta2, theta3,theta4;
	double a1, a2, a3, a4, a5, b1, b2, T2, T3, T4, T5;    
    
    theta = T/1000;
    theta2= pow(theta,2);
    theta3= pow(theta,3);
    theta4= pow(theta,4);
    
    T2=pow(T,2);
    T3=pow(T,3);
    T4=pow(T,4);
    T5=pow(T,5);
    
    if (select==1){
		    
	    //kJ/kmol [nist.gov]
	    //CH4
	   
	    if (T>=1300&&T<=6000){
	    	A =	85.81217;
			B =	11.26467;
			C =	-2.114146;
			D =	0.138190;
			E =	-26.42221;
			F =	-153.5327;
			G =	224.4143;
			H =	-74.87310;
			}
	    else if (T>=298){
	    	A =	-0.703029;
			B =	108.4773;
			C =	-42.52157;
			D =	5.862788;
			E =	0.678565;
			F =	-76.84376;
			G =	158.7163;
			H =	-74.87310;
		}
		else{
			printf("CH4 - temperatura imprópria (T= %lf)\n", T);
		}
		*dhcomb =(A*theta+B*theta2/2+C*theta3/3+D*theta4/4-E/theta+F-H)*1000;
	}
	else if (select==2){	
	    //C2H2
	    //kJ/kmol [nist.gov]
	
	    if (T>=1100&&T<=6000){
			A = 67.47244;
			B = 11.7511;
			C = -2.02147;
			D = 0.136195;
			E = -9.806418;
			F = 185.455;
			G = 253.5337;
			H = 226.7314;
	
			}
	    else if (T>=298){
			A = 40.68697;
			B = 40.73279;
			C = -16.1784;
			D = 3.669741;
			E = -0.658411;
			F = 210.7067;
			G = 235.0052;
			H = 226.7314;
	
		}
		else{
			printf("C2H2 - temperatura imprópria (T= %lf)\n", T);
		}
		*dhcomb =(A*theta+B*theta2/2+C*theta3/3+D*theta4/4-E/theta+F-H)*1000;
	}
	else if (select==3){	
	    //C2H4
	    //kJ/kmol [nist.gov]
	
	    if (T>=1200&&T<=6000){
			A = 106.5104;
			B = 13.7326;
			C = -2.628481;
			D = 0.174595;
			E = -26.14469;
			F = -35.36237;
			G = 275.0424;
			H = 52.46694;
	
			}
	    else if (T>=298){
			A = -6.38788;
			B = 184.4019;
			C = -112.9718;
			D = 28.49593;
			E = 0.31554;
			F = 48.17332;
			G = 163.1568;
			H = 52.46694;
	
		}
		else{
			printf("C2H4 - temperatura imprópria (T= %lf)\n", T);
		}
		*dhcomb =(A*theta+B*theta2/2+C*theta3/3+D*theta4/4-E/theta+F-H)*1000;
	}
	else if (select==4){
		//C2H6
		//kJ/kmol [NASA]
		if (T>=1000&&T<=6000){
			a1 = 4.04666674;
			a2 = 0.0153538766;
			a3 = -5.47039321E-6;
			a4 = 8.77826228E-10;
			a5 = -6.23167306E-14;
			b1 = -12447.3512;
			b2 = -0.968683607;
		}
		else if (T>=200){
			a1 = 4.29142492;
			a2 = -0.0055015427;
			a3 = 5.99438288E-5;
			a4 = -7.08466286E-8;
			a5 = 2.68685771E-11;
			b1 = -11522.2055;
			b2 = 2.66682316;
		}
		else{
			printf("C2H6 - temperatura imprópria (T= %lf)\n", T);
		}
		*dhcomb =Rmol*(a1*T+a2*T2/2+a3*T3/3+a4*T4/4+a5*T5/5+b1);
	}
	else if (select==5){
    //C3H8
    //kJ/kmol [Turns]
    
	    A = -1.4867;
		B = 74.339;
		C = -39.065;
		D = 8.0543;
		E = 0.01219;
		F = -27.313;
		G = 8.852;
	
		*dhcomb =4184*(A*theta+B*theta2/2+C*theta3/3+D*theta4/4-E/theta+F);
	}
	else if (select==6){
		//C4H10 - butane
		//kJ/kmol [NASA]
		
		if (T>=1000&&T<=6000){
			a1 = 9.44535834;
			a2 = 0.0257858073;
			a3 = -9.23619122E-6;
			a4 = 1.48632755E-9;
			a5 = -8.87897158E-14;
			b1 = -20138.2165;
			b2 = -26.3470076;
		}
		else if (T>=200){
			a1 = 6.14746806;
			a2 = 0.000155947389;
			a3 = 9.67913517E-5;
			a4 = -1.2548391E-7;
			a5 = 4.97816555E-11;
			b1 = -17599.4402;
			b2 = -1.09409879;
		}
		else{
			printf("C4H10 - temperatura imprópria (T= %lf)\n", T);
		}
		*dhcomb =Rmol*(a1*T+a2*T2/2+a3*T3/3+a4*T4/4+a5*T5/5+b1);
	}
	else if (select==7){
		//C5H12 - pentane
		//kJ/kmol [NASA]
		
		if (T>=1000&&T<=5000){
			a1 = 13.546998;
			a2 = 0.028421786;
			a3 = -9.4174648E-6;
			a4 = 1.3893589E-9;
			a5 = -7.4212609E-14;
			b1 = -24581.68;
			b2 = -47.021185;
		}
		else if (T>=298){
			a1 = 1.8983679;
			a2 = 0.041203037;
			a3 = 0.000012312175;
			a4 = -3.6589501E-8;
			a5 = 1.5042509E-11;
			b1 = -20091.5;
			b2 = 18.679072;
		}
		else{
			printf("C5H12 - temperatura imprópria (T= %lf)\n", T);
		}
		*dhcomb =Rmol*(a1*T+a2*T2/2+a3*T3/3+a4*T4/4+a5*T5/5+b1);
	}
}

int PropertyO2(char prop,double T) {
    double A,B,C,D,E,F,G,H,theta, theta2, theta3,theta4, propO2;
      
    //[nist.gov]
    
    theta = T/1000;
    theta2= pow(theta,2);
    theta3= pow(theta,3);
    theta4= pow(theta,4);

   
    if (T>=2000&&T<=6000){
    	A =	20.91111;
		B =	10.72071;
		C =	-2.020498;
		D =	0.146449;
		E =	9.245722;
		F =	5.337651;
		G =	237.6185;
		H =	0.0;
	}
    else if (T>=700){
    	A =	30.03235;
		B =	8.772972;
		C =	-3.988133;
		D =	0.788313;
		E =	-0.741599;
		F =	-11.32468;
		G =	236.1663;
		H =	0.0;
	}
	else if (T>=100){
    	A =	31.32234;
		B =	-20.23531;
		C =	57.86644;
		D =	-36.50624;
		E =	-0.007374;
		F =	-8.903471;
		G =	246.7945;
		H =	0.0;
	}
	else{
		printf("O2 - temperatura imprópria (T= %lf)\n", T);
	}
	
	if (prop=='h'){
		//dH [kJ/kmol]
		propO2 =(A*theta+B*theta2/2+C*theta3/3+D*theta4/4-E/theta+F-H)*1000;
	}
	if (prop=='s'){
		//S0 [kJ/kmol]
		propO2 =A*log(theta)+B*theta+C*theta2/2+D*theta3/3-E/(2*theta2)+G;	
	}
	return propO2;
}

int PropertyN2(char prop,double T) {
	
    double A,B,C,D,E,F,G,H,theta, theta2, theta3,theta4, propN2;
    
    //[nist.gov]
    
    theta = T/1000;
    theta2= pow(theta,2);
    theta3= pow(theta,3);
    theta4= pow(theta,4);

   
    if (T>=2000&&T<=6000){
    	A = 35.51872;
		B = 1.128728;
		C = -0.196103;
		D = 0.014662;
		E = -4.55376;
		F = -18.97091;
		G = 224.981;
		H = 0;
	}
    else if (T>=700){
    	A = 19.50583;
		B = 19.88705;
		C = -8.598535;
		D = 1.369784;
		E = 0.527601;
		F = -4.935202;
		G = 212.39;
		H = 0;
	}
	else if (T>=100){
    	A = 28.98641;
		B = 1.853978;
		C = -9.647459;
		D = 16.63537;
		E = 0.000117;
		F = -8.671914;
		G = 226.4168;
		H = 0;
	}
	else{
		printf("N2 - temperatura imprópria (T= %lf)\n", T);
	}
	if (prop=='h'){
		//dH [kJ/kmol]
		propN2 =(A*theta+B*theta2/2+C*theta3/3+D*theta4/4-E/theta+F-H)*1000;
	}
	if (prop=='s'){
		//S0 [kJ/kmol]
		propN2 =A*log(theta)+B*theta+C*theta2/2+D*theta3/3-E/(2*theta2)+G;	
	}
	return propN2;
}

int PropertyCO(char prop,double T) {
	
    double A,B,C,D,E,F,G,H,theta, theta2, theta3,theta4, propCO;
    
    //kJ/kmol [nist.gov]
    
    theta = T/1000;
    theta2= pow(theta,2);
    theta3= pow(theta,3);
    theta4= pow(theta,4);

   
    if (T>=1300&&T<=6000){
		A = 35.1507;
		B = 1.300095;
		C = -0.205921;
		D = 0.01355;
		E = -3.28278;
		F = -127.8375;
		G = 231.712;
		H = -110.5271;
	}
	else if (T>=298){
		A = 25.56759;
		B = 6.09613;
		C = 4.054656;
		D = -2.671301;
		E = 0.131021;
		F = -118.0089;
		G = 227.3665;
		H = -110.5271;
	}
	else{
		printf("CO - temperatura imprópria (T= %lf)\n", T);
	}
	if (prop=='h'){
		//dH [kJ/kmol]
		propCO =(A*theta+B*theta2/2+C*theta3/3+D*theta4/4-E/theta+F-H)*1000;
	}
	if (prop=='s'){
		//S0 [kJ/kmol]
		propCO =A*log(theta)+B*theta+C*theta2/2+D*theta3/3-E/(2*theta2)+G;	
	}
	return propCO;
}

int PropertyCO2(char prop,double T) {
	
    double A,B,C,D,E,F,G,H,theta, theta2, theta3,theta4, propCO2;
    
    //kJ/kmol [nist.gov]
    
    theta = T/1000;
    theta2= pow(theta,2);
    theta3= pow(theta,3);
    theta4= pow(theta,4);

   
    if (T>=1200&&T<=6000){
		A = 58.16639;
		B = 2.720074;
		C = -0.492289;
		D = 0.038844;
		E = -6.447293;
		F = -425.9186;
		G = 263.6125;
		H = -393.5224;
	}
	else if (T>=298){
    	A = 24.99735;
		B = 55.18696;
		C = -33.69137;
		D = 7.948387;
		E = -0.136638;
		F = -403.6075;
		G = 228.2431;
		H = -393.5224;
	}
	else{
		printf("CO2 - temperatura imprópria (T= %lf)\n", T);
	}
	if (prop=='h'){
		//dH [kJ/kmol]
		propCO2 =(A*theta+B*theta2/2+C*theta3/3+D*theta4/4-E/theta+F-H)*1000;
	}
	if (prop=='s'){
		//S0 [kJ/kmol]
		propCO2 =A*log(theta)+B*theta+C*theta2/2+D*theta3/3-E/(2*theta2)+G;	
	}
	return propCO2;
}

int PropertyH2(char prop,double T) {
	
    double A,B,C,D,E,F,G,H,theta, theta2, theta3, theta4, propH2;
    
    //kJ/kmol [nist.gov]
    
    theta = T/1000;
    theta2= pow(theta,2);
    theta3= pow(theta,3);
    theta4= pow(theta,4);

   
    if (T>=2500&&T<=6000){
		A = 43.41356;
		B = -4.293079;
		C = 1.272428;
		D = -0.096876;
		E = -20.533862;
		F = -38.515158;
		G = 162.081354;
		H = 0;

	}
	else if (T>=1000){
		A = 18.563083;
		B = 12.257357;
		C = -2.859786;
		D = 0.268238;
		E = 1.97799;
		F = -1.147438;
		G = 156.288133;
		H = 0;

	}
	else if (T>=298){
		A = 33.066178;
		B = -11.363417;
		C = 11.432816;
		D = -2.772874;
		E = -0.158558;
		F = -9.980797;
		G = 172.707974;
		H = 0;
	}
	else{
		printf("H2 - temperatura imprópria (T= %lf)\n", T);
	}
	if (prop=='h'){
		//dH [kJ/kmol]
		propH2 =(A*theta+B*theta2/2+C*theta3/3+D*theta4/4-E/theta+F-H)*1000;
	}
	if (prop=='s'){
		//S0 [kJ/kmol]
		propH2 =A*log(theta)+B*theta+C*theta2/2+D*theta3/3-E/(2*theta2)+G;	
	}
	return propH2;
}

int PropertyH2O(char prop,double T) {
	
    double A,B,C,D,E,F,G,H,theta, theta2, theta3, theta4, propH2O;
    
    //kJ/kmol [nist.gov]
    
    theta = T/1000;
    theta2= pow(theta,2);
    theta3= pow(theta,3);
    theta4= pow(theta,4);

   
    if (T>=1700&&T<=6000){
		A = 41.96426;
		B = 8.622053;
		C = -1.49978;
		D = 0.098119;
		E = -11.15764;
		F = -272.1797;
		G = 219.7809;
		H = -241.8264;
	}
	else if (T>=500){
		A = 30.092;
		B = 6.832514;
		C = 6.793435;
		D = -2.53448;
		E = 0.082139;
		F = -250.881;
		G = 223.3967;
		H = -241.8264;
	}
	else{
		printf("H2O - temperatura imprópria (T= %lf)\n", T);
	}
	if (prop=='h'){
		//dH [kJ/kmol]
		propH2O =(A*theta+B*theta2/2+C*theta3/3+D*theta4/4-E/theta+F-H)*1000;
	}
	if (prop=='s'){
		//S0 [kJ/kmol]
		propH2O =A*log(theta)+B*theta+C*theta2/2+D*theta3/3-E/(2*theta2)+G;	
	}
	return propH2O;
}
