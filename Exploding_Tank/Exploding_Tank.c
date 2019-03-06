#include <stdio.h>
#include <math.h>

int deltah (double T, double *dhcomb,double  *dhCO2,double  *dhCO,double  *dhO2,double  *dhN2,double  *dhH20)
{
    double T2, T3, T4, T5,A,B,C;
    T2=pow(T,2);
    T3=pow(T,3);
    T4=pow(T,4);
    T5=pow(T,5);

    *dhcomb=20.056*T+1.4077E-1*T2-4.7143E-6*T3-2.3643E-8*T4+6.8398E-12*T5-1.8198E4;
    if (T>=1000){
        *dhCO=(0.03025078E2*T +0.14426885E-2/2*T2 -0.05630827E-05/3*T3+ 0.10185813E-09/4*T4 -0.06910951E-13/5*T5-0.14268350E05)*8.31446;
        *dhCO2=(0.04453623E+02*T+0.03140168E-01/2*T2-0.12784105E-05/3*T3+0.02393996E-08/4*T4-0.16690333E-13/5*T5-0.04896696E+06)*8.31446;
        *dhH20=(0.02672145E+02*T+0.03056293E-01/2*T2-0.08730260E-05/3*T3+0.12009964E-09/4*T4-0.06391618E-13/5*T5-0.02989921E+06)*8.31446;
        *dhN2=(0.02926640E+02*T+0.14879768E-02/2*T2-0.05684760E-05/3*T3+ 0.10097038E-09/4*T4-0.06753351E-13/5*T5-0.09227977E+04)*8.31446;
        *dhO2=(0.03697578E+02*T+0.06135197E-02/2*T2-0.12588420E-06/3*T3+ 0.01775281E-09/4*T4-0.11364354E-14/5*T5-0.12339301E+04)*8.31446;
    }
    else if (T>=350){
        *dhCO=(0.03262451E+02*T+0.15119409E-02/2*T2 -0.03881755E-04/3*T3+ 0.05581944E-07/4*T4 -0.02474951E-10/5*T5 -0.14310539E+05)*8.31446;
        *dhCO2=(0.02275724E+02*T+0.09922072E-01/2*T2 -0.10409113E-04/3*T3+0.06866686E-07/4*T4 -0.02117280E-10/5*T5 -0.04837314E+06)*8.31446;
        *dhH20=(0.02672145E+02*T+0.03056293E-01/2*T2-0.08730260E-05/3*T3+ 0.12009964E-09/4*T4 -0.06391618E-13/5*T5 -0.02989921E+06)*8.31446;
        *dhN2=(0.03298677E+02*T+0.14082404E-02/2*T2 -0.03963222E-04/3*T3+ 0.05641515E-07/4*T4 -0.02444854E-10/5*T5 -0.10208999E+04)*8.31446;
        *dhO2=(0.03212936E+02*T+0.11274864E-02/2*T2 -0.05756150E-05/3*T3+ 0.13138773E-08/4*T4 -0.08768554E-11/5*T5 -0.10052490E+04)*8.31446;
    }
    else printf("deu merda");

}

int main()
{
    //variáveis
    double Vint, pini, Tini, ycomb, Rcal, Rmol, nN2Ar, passo, R, m, t, v;
    double hfcomb, hfCO2, hfH20, hfCO, dhcomb, dhCO2, dhCO, dhO2, dhN2, dhH20, MO2, MN2, Mcomb, Mmist, MCO, Mag, MCO2;
    double ncomb, nAr, nO2, nN2, nCO2, nag, ccombini, cO2ini, cN2ini, ccomb, cO2, cN2, temp, temp1, U2, p;
    double cCO, cCO2, cag, ecCO, ecCO2, ecag, ecO2, eccomb, fcCO, fcCO2, fcag, fcO2, fccomb, ro,  U1;
    double dif, dif1, difO2, difCO,difCO2, difag, T, ftemp, ftemp1, fdif, fdif1, fdifO2, fdifCO, fdifCO2, fdifag;
    int i,a,c;

    //valores iniciais
    Vint=2; //m^3
    pini=0.7; //bar
    Tini=750; //K
    ycomb=2.835E-2;
    Rcal=1.987; //kcal/(mol*K)
    Rmol=8.31446; //kJ/(kmol*K)
    nN2Ar=3.76; //mol
    passo=1E-9;

    //propriedades constantes
    MO2=31.999;
    MN2=28.013;
    Mcomb=58.123;
    Mag=18.016;
    MCO=28.010;
    MCO2=44.011;
    hfcomb=-126200; //kj/kmol [Yaws]
    hfCO2=-393546; //kj/kmol [Turns]
    hfH20=-241854; //kj/kmol [Turns]
    hfCO=-110541; //kj/kmol [Turns]

    //calculos preliminares
    ncomb=1;
    nAr=(1/ycomb-1)/4.76;
    nO2=nAr;
    nN2=nN2Ar*nAr;
    ro=pini/(Rmol*Tini)*1E-1;
    ccombini=ycomb*ro;
    cO2ini=nAr*ccombini;
    cN2ini=nN2Ar*cO2ini;
    deltah(Tini, &dhcomb, &dhCO2, &dhCO, &dhO2, &dhN2, &dhH20);
    U1=(ncomb*(hfcomb+dhcomb)+nN2*(dhN2)+nO2*(dhO2)-(ncomb+nO2+nN2)*(Rmol*Tini))*ccombini;
    printf("U1= %f\n", U1);
    ccomb=ccombini;
    cO2=cO2ini;
    cCO=0;
    cCO2=0;
    cag=0;
    cN2=cN2ini;
    T=Tini;
    p=pini;

    i=0;
    t=0;
    while(ccomb>0.02*ccombini){

        t=t+passo;
        i++;

        v=(ccomb+cO2+cCO+cCO2+cag+cN2);
        p=v*Rmol*T*10;



        //predictor
        ftemp=-8.8E11*pow(p,1.2)*exp(-30000/(Rcal*T))*pow(ccomb,0.15)*pow(cO2,1.6);
        ftemp1=-pow(10, 14.6)*exp(-40000/(Rcal*T))*cCO*pow(cO2,0.25)*pow(cag,0.5)+5E8*exp(-40000/(Rcal*T))*cCO2;

        fdif=ftemp*passo;
        fdifO2=4.5*fdif;
        fdifCO=-4*fdif;
        fdifag=-5*fdif;

        fdif1=ftemp1*passo;
        fdifO2=difO2+0.5*fdif1;
        fdifCO=difCO+fdif1;
        fdifCO2=-fdif1;

        fccomb=ccomb+fdif;
        fcO2=cO2+fdifO2;
        fcCO=cCO+fdifCO;
        fcCO2=cCO2+fdifCO2;
        fcag=cag+fdifag;

        //corrector
        eccomb=(ccomb+fccomb)/2;
        ecO2=(cO2+fcO2)/2;
        ecCO=(cCO+fcCO)/2;
        ecCO2=(cCO2+fcCO2)/2;
        ecag=(cag+fcag)/2;

        a=0;
        if (i==1) a=1;
        while (a==0){
            deltah(T, &dhcomb, &dhCO2, &dhCO, &dhO2, &dhN2, &dhH20);
            U2=(eccomb*(hfcomb+dhcomb)+cN2*(dhN2)+ecO2*(dhO2)+ecCO*(dhCO)+ecCO2*(dhCO2)+ecag*(dhH20)-(eccomb+ecO2+cN2+ecCO+ecCO2+ecag)*(Rmol*T));
            if (U2<U1) T=T+0.000001;
            else a=1;
        }

        temp=-8.8E11*pow(p,1.2)*exp(-30000/(Rcal*T))*pow(eccomb,0.15)*pow(ecO2,1.6);
        temp1=-pow(10, 14.6)*exp(-40000/(Rcal*T))*ecCO*pow(ecO2,0.25)*pow(ecag,0.5)+5E8*exp(-40000/(Rcal*T))*ecCO2;

        dif=temp*passo;
        difO2=4.5*dif;
        difCO=-4*dif;
        difag=-5*dif;

        dif1=temp1*passo;
        difO2=difO2+0.5*dif1;
        difCO=difCO+dif1;
        difCO2=-dif1;

        ccomb=ccomb+dif;
        cO2=cO2+difO2;
        cCO=cCO+difCO;
        cCO2=cCO2+difCO2;
        cag=cag+difag;

        if(i==1) printf("t=%.9f\n ccomb=%.30f T=%f p=%.15f\n", t, ccomb, T, p);
        if(i==5) printf("t=%.9f\n ccomb=%.30f T=%f p=%.15f\n", t,ccomb, T, p);
        if(i==10) printf("t=%.9f\n ccomb=%.30f T=%f p=%.15f\n", t,ccomb, T, p);
        if(i==50) printf("t=%.9f\n ccomb=%.30f T=%f p=%.15f\n",t, ccomb, T, p);
        if(i==100) printf("t=%.9f\n ccomb=%.30f T=%f p=%.15f\n",t, ccomb, T, p);
        if(i==500) printf("t=%.9f\n ccomb=%.30f T=%f p=%.15f\n",t, ccomb, T, p);
        if(i==1000) printf("t=%.9f\n ccomb=%.30f T=%f p=%.15f\n",t, ccomb, T, p);
        if(i==5000) printf("t=%.9f\n ccomb=%.30f T=%f p=%.15f\n",t, ccomb, T, p);
        if(i==10000) printf("t=%.9f\n ccomb=%.30f T=%f p=%.15f\n",t, ccomb, T, p);
        if(i==50000) printf("t=%.9f\n ccomb=%.30f T=%f p=%.15f\n",t, ccomb, T, p);
        if(i==100000) printf("t=%.9f\n ccomb=%.30f T=%f p=%.15f\n",t, ccomb, T, p);
        if(i==500000) printf("t=%.9f\n ccomb=%.30f T=%f p=%.15f\n",t, ccomb, T, p);
        if(i==1000000) printf("t=%.9f\n ccomb=%.30f T=%f p=%.15f\n",t, ccomb, T, p);
        if(i==5000000) printf("t=%.9f\n ccomb=%.30f T=%f p=%.15f\n",t, ccomb, T, p);
        if(i==10000000) printf("t=%.9f\n ccomb=%.30f T=%f p=%.15f\n",t, ccomb, T, p);
    }

    printf("tempo=%1.30f\n", t);
    printf("temperatura=%1.30f\n", T);
    v=(ccomb+cO2+cCO+cCO2+cag+cN2);
    p=v*Rmol*T*10;
    printf("pressao=%.30f\n", p);

    scanf ("%d", c);




    return 0;

}
