#include <stdio.h>
#include <math.h>
/*   HIPÓTESYS
1) Dry air
2) Free falling
3) Straight trajectory */


//enthalpy evaluation (valid for T between 283 to 299 K)
int entalp (double T, double *hl, double *hv){
    *hl=-0.0006*T*T+4.5406*T-1194.9464;
    *hv=219.999*T*T-126501.777*T+18184080.665;
}

//specific volume evaluation (valid for T between 283 to 299 K)
int prop (double T, double *nil, double *nig){
    *nil=3.180957E-21*T*T+2E-7*T+0.0009434;
    *nig=0.13911*T*T-85.00449*T+13021.31973;
}

int main (){

    double d, A, Tinf, T, Dab, pinf,v,M, Vini, V, mi,na,Sc,Re,visc,Sh,h, nil, nig;
    double g, dt, t, D, roar,dV, T1, Nu, Pr, k, hl, hv, dM, s;
    int i, a, c;
	
    //File handling
    FILE *Vol, *dist, *temp, *vel;
    Vol = fopen("volume.txt","w");
    vel = fopen("velocidade.txt","w");
    dist = fopen("distancia.txt","w");
    temp = fopen("temperatura.txt","w");

    //Initial inputs
    Tinf=295; //K
    T=Tinf;
    d=1E-3; //m (diâmetro)
    roar=1.18473; //kg/m^3
    nil=1.002E-3; //m^3/kg
    nig=51.94; //m^3/kg
    Dab=0.26E-4; //m/s^2
    pinf=1.013E5; //Pa
    visc= 15.445E-6; //m/s^2
    mi=182.1E-7; //Ns/m^2
    v=0; //m/s
    g=9.806; //m/s^2
    Pr=0.70803; //ar a 295K
    k=0.606; //W/mK
    dt=1E-8; //s (passo)
    t=0;
    i=0;
    s=0;

    // Preliminary evaluation
    A=M_PI*pow(d,2); //m^2
    Vini=M_PI*pow(d,3)/6; //m^3
    M=Vini/nil; //kg
    V=Vini;


    while(V>0.001*Vini) {

        //Mass convection
        prop(T, &nil, &nig);
        s=v*dt+s; //Distance
        D=3*M_PI*mi*v*d+9*M_PI*roar*pow(v,2)*pow(d,2)/16; //Drag
        v=(g-D/M)*dt+v; //Velocity derived from: mass*acc = Weight - Drag
        Re=v*d/visc;
        Sc=nig/Dab;
        Sh=2+(0.4*pow(Re,0.5)+0.06*pow(Re,2/3))*pow(Sc,0.4);
        h=Sh*Dab/d;
        na=h*A/nig; //Mass loss

        //Final values
        dV=na*dt*nil;
        V=V-dV;
        M=V/nil;
        d=pow(6*V/M_PI,1/3);
        A=M_PI*pow(d,2);
        dM=dV/nil;

        //Heat excange
        Re=v*d/visc;
        Nu=2+(0.4*pow(Re,0.5)+0.06*pow(Re,2/3))*pow(Pr,0.4);
        h=Nu*k/d;
        a=0;
        T1=T;
        while (a==0) {
            entalp (T, &hl, &hv);
            T=dM*(hl-hv)/(dt*h*A)+Tinf;
            if (T<=T1) a=1;
            else T1=T;
        }

        //Saving values to the .txt file
        c=i%100;
        if (c==0) {
            fprintf (Vol, "%.8f %.15f\n", t, V);
            fprintf (vel, "%.8f %.15f\n", t, v);
            fprintf (dist, "%.8f %.15f\n", t, s);
            fprintf (temp, "%.8f %.15f\n", t, T);
        }

        t=t+dt;
        i++;
    }

    fprintf (Vol, "%.8f %.15f", t, V);
    fprintf (vel, "%.8f %.15f\n", t, v);
    fprintf (dist, "%.8f %.15f\n", t, s);
    fprintf (temp, "%.8f %.15f\n", t, T);

    close (Vol);
    close (dist);
    close (vel);
    close (temp);

    return 0;
}
