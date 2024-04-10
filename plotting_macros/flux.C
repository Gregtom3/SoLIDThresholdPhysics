using namespace std;

Double_t brem(Double_t *x, Double_t *p) //p{d,X0,Eb}
{
   return (0.5*p[0]/p[1])*(1/x[0])*((4./3.) - (4./3.)*(x[0]/p[2]) + x[0]*x[0]/(p[2]*p[2]));
    //    (0.5*d/X0)*(1/Eg)*((4./3.) - (4./3.)*(Eg/Eb) + Eg*Eg/(Eb*Eb));
}

Double_t epa(Double_t *x, Double_t *p) //p{Eb,Q2_max}
{
  const double alpha = 1./137.;
  const double PI = 3.14159265358979312;
  double y = x[0]/p[0]; // Eg/Eb
  double me = 0.00051;
  double Mp = 0.9383;
  double Q2_min = me*me*y*y/(1 - y);
  return (1/p[0])*alpha/(PI*y)*( (1 - y + y*y/2)*log(p[1]/Q2_min) - (1 - y));
//            (1/Eb)*alpha/(PI*x)*( (1 - x + x*x/2)*log(Q2_max/Q2_min) - (1 - x));
}

void flux(double Eb=8.8)
{
// double Emin=5,Emax=12;
// TH1F *flux_brem=new TH1F("flux_brem",";E_{#gamma} (GeV); photon flux",35,Emin,Emax);
// TH1F *flux_epa=new TH1F("flux_epa",";E_{#gamma} (GeV); photon flux",35,Emin,Emax);

// double Eb=8.8;
TF1 *fbrem  = new TF1("flux",brem,5,Eb,3);
TF1 *fepa  = new TF1("flux",epa,5,Eb,2);
fbrem->SetNpx(500);
fepa->SetNpx(500);

fbrem->SetParameters(15,769,Eb); //15cm LD2
fepa->SetParameters(Eb,0.5); //Q2_max=0.5 and 0.1 differ by a few percent only for EPA
// cout << fepa->Eval(8) << endl;

gStyle->SetOptStat(0);

gStyle->SetOptStat(0);
gStyle->SetOptFit(1);

  gStyle->SetPadColor(0);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.1);    
  gStyle->SetLabelSize(0.04,"xyz"); // size of axis values
  gStyle->SetTitleSize(0.06,"xyz");  
  gStyle->SetTitleOffset(1.3,"y");
  gStyle->SetTitleOffset(1,"x");    
  gStyle->SetTitleSize(0.06,"t"); 
  
TCanvas *c_flux = new TCanvas("c_flux","c_flux",1000,600);
TH1F *hflux=new TH1F("hflux",Form("e^{-} beam energy %.01f GeV;E_{#gamma} (GeV);photo flux dN/dE_{#gamma} (GeV^{-1})",Eb),60,5,Eb);
hflux->SetMinimum(0);
hflux->SetMaximum(0.005);
hflux->Draw();
fepa->SetLineColor(kRed);
fepa->Draw("same");
fbrem->SetLineColor(kBlue);
fbrem->Draw("same");

    TLegend* leg1 = new TLegend(0.7, 0.75, 0.95, 0.9);  
    leg1->AddEntry(fepa,"flux EPA","l");
    leg1->AddEntry(fbrem,"flux Brem","l");    
    leg1->Draw();    

}

// def Bremmstrahlung(kmin,kmax,beamE,target_length,X0):
//     ymin = kmin/beamE
//     ymax = kmax/beamE
//     y = random.uniform(ymin,ymax)
//     # We use "0.5*d" to simulate half the target as the radiator
//     flux = 0.5*target_length/X0/(y*beamE)*(4/3 - 4/3 * y + y**2)
//     return  flux, y*beamE
// 
// def N_EquivalentPhotonApproximation(gammaE,beamE,q2max):
//     # (https://lss.fnal.gov/archive/other/lpc-94-35.pdf)
//     # Factors of mE^2 are needed atop the fractions 1/q2min and 1/q2max for units
//     # I believe this is a mistake in the derivation
//     y = gammaE/beamE
//     q2min = mE**2 * y**2 / (1-y)
//     return (alpha_em*y)/(4*np.pi)*((2*q2min+4*mE**2)/(q2min)*np.log(q2max/q2min)-4*mE**2/q2min+4*mE**2/q2max)/beamE
// }
// 
// double KinFuncs::N_EPA(double Eb, double Eg, double Q2_max)
// {
//   const double alpha = 1./137.;
//   const double PI = 3.14159265358979312;
// 
//   double x = Eg/Eb;
//   double me = 0.00051;
//   double Mp = 0.9383;
//   double Q2_min = me*me*x*x/(1 - x);
//   return (1/Eb)*alpha/(PI*x)*( (1 - x + x*x/2)*log(Q2_max/Q2_min) - (1 - x));
// }
// 
// double KinFuncs::N_Brem(double Eg, double Eb, double d, double X0)
// {
//   // The factor 0.5 is because when one integrates over (l - x)*dx, then you get l^2/2
//   return (0.5*d/X0)*(1/Eg)*((4./3.) - (4./3.)*(Eg/Eb) + Eg*Eg/(Eb*Eb));
// }
