#include <iostream> 
#include <fstream> 
#include <TROOT.h>
void jpsi_t() {
// https://pdg.lbl.gov/2019/reviews/rpp2019-rev-kinematics.pdf
double E[9]={6,7,8,9,10,12,14,16,18};
for(int i=0;i<9;i++)
{
double Ebeam=E[i];
// double m1=0,m2=0.938,m3=3.1,m4=0.938; //proton target
double m1=0,m2=1.876,m3=3.1,m4=1.876; //deuteron target
double s=(Ebeam+m2)*(Ebeam+m2)-Ebeam*Ebeam;
double E1cm=(s+m1*m1-m2*m2)/(2.*sqrt(s));
double E3cm=(s+m3*m3-m4*m4)/(2.*sqrt(s));
double p1cm=sqrt(E1cm*E1cm-m1*m1);
double p3cm=sqrt(E3cm*E3cm-m3*m3);
double t0=pow((m1*m1-m3*m3-m2*m2+m4*m4)/(2.*sqrt(s)),2)-pow(p1cm-p3cm,2);
double t1=pow((m1*m1-m3*m3-m2*m2+m4*m4)/(2.*sqrt(s)),2)-pow(p1cm+p3cm,2);
double t=t0-4*p1cm*p3cm*sin(28./2./180.*3.1416); //solid jpsi setup max acceptance angle at 28deg
cout << Ebeam << "\t" << t0 << "\t" << t1 << "\t" << t << endl;
}
}