#include <iostream> 
#include <fstream> 
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TF1.h>
#include <TLegend.h>
#include <TLorentzVector.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TMinuit.h>
#include <TRandom3.h>
#include <TString.h>
#include <TGraph.h>
#include <TLegend.h>

using namespace std;

const double DEG=180./3.1416;

void plot(string input_filename){
//  set_style();
gROOT->Reset();
//   gStyle->SetPalette(1);
// // gStyle->SetPalette(kBird);
// gStyle->SetOptStat(11111111);
gStyle->SetOptStat(0);
gStyle->SetOptFit(1);

  gStyle->SetPadColor(0);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.1);    
  gStyle->SetLabelSize(0.05,"xyz"); // size of axis values
  gStyle->SetTitleSize(0.05,"xyz");  
  gStyle->SetTitleOffset(1.3,"y");
  gStyle->SetTitleOffset(1,"x");    
  gStyle->SetTitleSize(0.06,"t"); 

  double targetM,Emin,Ymax,Ymin;
  if (input_filename.find("_ed_",0) != string::npos)  {targetM=1.8761358; Emin=5.5; Ymax=1e-4; Ymin=1e-7;}
  else if(input_filename.find("_ep_",0) != string::npos)  {targetM=0.938272; Emin=8; Ymax=1e2; Ymin=1e-1;}
  else {cout << " wrong target " << endl; return;}
  
// TH2F *hacceptance_ThetaP_forwardangle,*hacceptance_ThetaP_largeangle;
// TFile *acceptancefile=new TFile("tables/acceptance_solid_JPsi_electron_target315_output.root");  
// hacceptance_ThetaP_forwardangle=(TH2F*) acceptancefile->Get("acceptance_ThetaP_forwardangle");  
// hacceptance_ThetaP_largeangle=(TH2F*) acceptancefile->Get("acceptance_ThetaP_largeangle");
// TCanvas *c_acceptance = new TCanvas("acceptance","acceptance",1200,900);
// c_acceptance->Divide(2,1);
// c_acceptance->cd(1);
// hacceptance_ThetaP_forwardangle->Draw("colz");
// c_acceptance->cd(2);
// hacceptance_ThetaP_largeangle->Draw("colz");   

TH1F *hcount=new TH1F("hcount","SoLID J/#psi-d, 1.25uA on 15cm LD2 for 50 days with eff 0.8;E_{#gamma} (GeV);counts /0.2 GeV",25,5,10);
TH1F *hcount_acc2=new TH1F("hcount_acc2",";E_{#gamma} (GeV);counts /0.2 GeV",25,5,10);
TH1F *hcount_acc3=new TH1F("hcount_acc3",";E_{#gamma} (GeV);counts /0.2 GeV",25,5,10);

TH1F *hsigma=new TH1F("hsigma",";E_{#gamma} (GeV);#sigma (nb)",100,5,10);

const int n=3;
TH1F *hdsigmadt[n];
for(int k=0;k<n;k++){
  char hstname[100];
  sprintf(hstname,"hdsigmadt_%i",k);
  hdsigmadt[k]=new TH1F(hstname,"compare generator with Harry's model;-(t-t0) (GeV);d#sigma/dt (nb/GeV^{2})",100,0,20);
}

TH2F *hgen_ePlus= new TH2F("hgen_ePlus",";#theta (deg);P (GeV)",180,0,180,100,0,10);
TH2F *hgen_eMinus= new TH2F("hgen_eMinus",";#theta (deg);P (GeV)",180,0,180,100,0,10);
TH2F *hgen_hOut= new TH2F("hgen_hOut",";#theta (deg);P (GeV)",180,0,180,100,0,10);
TH2F *hgenCM_ePlus= new TH2F("hgenCM_ePlus",";#theta (deg);P (GeV)",180,0,180,100,0,10);
TH2F *hgenCM_eMinus= new TH2F("hgenCM_eMinus",";#theta (deg);P (GeV)",180,0,180,100,0,10);
TH2F *hgenCM_hOut= new TH2F("hgenCM_hOut",";#theta (deg);P (GeV)",180,0,180,100,0,10);
TH2F *hacc1_ePlus= new TH2F("hacc1_ePlus",";#theta (deg);P (GeV)",40,0,40,100,0,10);
TH2F *hacc1_eMinus= new TH2F("hacc1_eMinus",";#theta (deg);P (GeV)",40,0,40,100,0,10);
TH2F *hacc1_hOut= new TH2F("hacc1_hOut",";#theta (deg);P (GeV)",40,0,40,100,0,10);
TH2F *hacc2_ePlus= new TH2F("hacc2_ePlus",";#theta (deg);P (GeV)",40,0,40,100,0,10);
TH2F *hacc2_eMinus= new TH2F("hacc2_eMinus",";#theta (deg);P (GeV)",40,0,40,100,0,10);
TH2F *hacc2_hOut= new TH2F("hacc2_hOut",";#theta (deg);P (GeV)",40,0,40,100,0,10);
TH2F *hacc3_ePlus= new TH2F("hacc3_ePlus",";#theta (deg);P (GeV)",40,0,40,100,0,10);
TH2F *hacc3_eMinus= new TH2F("hacc3_eMinus",";#theta (deg);P (GeV)",40,0,40,100,0,10);
TH2F *hacc3_hOut= new TH2F("hacc3_hOut",";#theta (deg);P (GeV)",40,0,40,100,0,10);

TH1F *hflux=new TH1F("hflux",";E_{#gamma} (GeV);flux",70,5,12);
TH1F *hflux_brem=new TH1F("hflux_brem",";E_{#gamma} (GeV);flux",70,5,12);
TH1F *hflux_epa=new TH1F("hflux_epa",";E_{#gamma} (GeV);flux",70,5,12);

TH2F *hangle[3][3];
for(int i=0;i<3;i++){
for(int k=0;k<3;k++){
  char hstname[100];
  sprintf(hstname,"hangle_%i_%i",i,k);
  hangle[i][k]=new TH2F(hstname,";deg;deg",180,0,180,180,0,180);
}}

TFile *inputfile=new TFile(input_filename.c_str());
TTree *tree = (TTree*) inputfile->Get("tree");

// have to initialize those pointers, otherwise weird result and segmentation error can happen
TLorentzVector *eIn=0, *hIn=0, *eOut=0, *hOut=0, *ePlus=0, *eMinus=0, *q=0;
TLorentzVector *smear_eOut=0, *smear_hOut=0, *smear_ePlus=0, *smear_eMinus=0, *smear_q=0, *smear_VM=0;
// TLorentzVector eIn, hIn, eOut, hOut, ePlus, eMinus, q;
// TLorentzVector smear_eOut, smear_hOut, smear_ePlus, smear_eMinus, smear_q, smear_VM;

// double  *dsigma=0,*weight=0,*decay_weight=0,*psf=0,*flux=0,*flux_brem=0,*flux_epa=0,*acc_eOut=0,*acc_hOut=0,*acc_ePlus=0,*acc_eMinus=0,*Q2=0,*smear_Q2=0,*gammaE=0,*smear_gammaE=0,*t=0,*smear_t=0,*W=0,*smear_W=0,*m_vm=0,*smear_m_vm=0,*J=0,*smear_J=0;
double  dsigma,weight,decay_weight,psf,flux,flux_brem,flux_epa,acc_eOut,acc_hOut,acc_ePlus,acc_eMinus,Q2,smear_Q2,gammaE,smear_gammaE,t,smear_t,W,smear_W,m_vm,smear_m_vm,J,smear_J;

// tree->SetBranchStatus("*",0);
tree->SetBranchAddress("eIn",&eIn);
tree->SetBranchAddress("hIn",&hIn);
tree->SetBranchAddress("eOut",&eOut);
tree->SetBranchAddress("hOut",&hOut);
tree->SetBranchAddress("ePlus",&ePlus);
tree->SetBranchAddress("eMinus",&eMinus);
tree->SetBranchAddress("smear_eOut",&smear_eOut);
tree->SetBranchAddress("smear_hOut",&smear_hOut);
tree->SetBranchAddress("smear_ePlus",&smear_ePlus);
tree->SetBranchAddress("smear_eMinus",&smear_eMinus);
tree->SetBranchAddress("q",&q);
tree->SetBranchAddress("smear_q",&smear_q);
tree->SetBranchAddress("smear_VM",&smear_VM);
tree->SetBranchAddress("dsigma",&dsigma);
tree->SetBranchAddress("weight",&weight);
tree->SetBranchAddress("decay_weight",&decay_weight);
tree->SetBranchAddress("psf",&psf);
tree->SetBranchAddress("flux",&flux);
tree->SetBranchAddress("flux_brem",&flux_brem);
tree->SetBranchAddress("flux_epa",&flux_epa);
tree->SetBranchAddress("acc_eOut",&acc_eOut);
tree->SetBranchAddress("acc_hOut",&acc_hOut);
tree->SetBranchAddress("acc_ePlus",&acc_ePlus);
tree->SetBranchAddress("acc_eMinus",&acc_eMinus);
tree->SetBranchAddress("Q2",&Q2);
tree->SetBranchAddress("smear_Q2",&smear_Q2);
tree->SetBranchAddress("gammaE",&gammaE);
tree->SetBranchAddress("smear_gammaE",&smear_gammaE);
tree->SetBranchAddress("t",&t);
tree->SetBranchAddress("smear_t",&smear_t);
tree->SetBranchAddress("W",&W);
tree->SetBranchAddress("smear_W",&smear_W);
tree->SetBranchAddress("m_vm",&m_vm);
tree->SetBranchAddress("smear_m_vm",&smear_m_vm);
tree->SetBranchAddress("J",&J);
tree->SetBranchAddress("smear_J",&smear_J);

int nevent = tree->GetEntries();
cout << "nevent " << nevent << endl;

// nevent=nevent/10;

for(Int_t i=0; i < nevent; i++){
//     Double_t id = i;
//     Double_t neventd = nevent;
//     if(floor(id/neventd*20) == countd){
//    
//       cout << "=" << flush;
//       countd++;
//       if(countd == 20){
// 	cout << "|" << endl;
//       }
//     }

//     cout << i << "\r";
    tree->GetEvent(i);
//         cout << q.E() << " " << ePlus.M() << " " << gammaE << endl;    
//         cout << (*ePlus+*eMinus).M() << " " << m_vm << " " << gammaE << endl;
//     cout << hOut->Pz() << " " << q->Pz()+hIn->Pz()-ePlus->Pz()-eMinus->Pz() << endl;
    
    // 3e-6/1.6e-19*15*0.071*6.02e23=1.2e37 neucleon lumi for 3uA 15cm LH2
    // 1.25e-6/1.6e-19*15*0.169*6.02e23=1.2e37  nucleon lumi for 1.25uA 15cm LD2
    double lumi=0.6e37; //nucleus lumi for 1.25uA 15cm LD2
    double rate_convert = 1e-33*lumi*0.8;  //nb = 1e-33 cm2 and 0.8 eff
    double count_convert = rate_convert*3600*24*50; //50 days

    //to get dcount/dE curve
//     hcount->Fill(gammaE,weight/nevent/hcount->GetXaxis()->GetBinWidth(1)*count_convert);
//     hcount_acc2->Fill(gammaE,weight/nevent/hcount_acc2->GetXaxis()->GetBinWidth(1)*count_convert*acc_ePlus*acc_eMinus);
//     hcount_acc3->Fill(gammaE,weight/nevent/hcount_acc3->GetXaxis()->GetBinWidth(1)*count_convert*acc_ePlus*acc_eMinus*acc_hOut);

    //to get count histogram
    hcount->Fill(gammaE,weight/nevent*count_convert);
    hcount_acc2->Fill(gammaE,weight/nevent*count_convert*acc_ePlus*acc_eMinus);    
    hcount_acc3->Fill(gammaE,weight/nevent*count_convert*acc_ePlus*acc_eMinus*acc_hOut);  

    //check acceptance
//     double acc=0;      
//     double acc_prot=0,acc_em=0,acc_ep=0;
//     double acc_em_FA=0,acc_em_LA=0,acc_ep_FA=0,acc_ep_LA=0;
//     double acc_prot_FA=0,acc_prot_LA=0;
//     acc_prot_FA=hacceptance_ThetaP_forwardangle->GetBinContent(hacceptance_ThetaP_forwardangle->FindBin(hOut->Theta()*DEG,hOut->P()));
//     acc_ep_FA=hacceptance_ThetaP_forwardangle->GetBinContent(hacceptance_ThetaP_forwardangle->FindBin(ePlus->Theta()*DEG,ePlus->P()));	  
//     acc_em_FA=hacceptance_ThetaP_forwardangle->GetBinContent(hacceptance_ThetaP_forwardangle->FindBin(eMinus->Theta()*DEG,eMinus->P()));	  
//     
//     acc_prot_LA = hacceptance_ThetaP_largeangle->GetBinContent(hacceptance_ThetaP_largeangle->FindBin(hOut->Theta()*DEG,hOut->P()));
//     acc_ep_LA = hacceptance_ThetaP_largeangle->GetBinContent(hacceptance_ThetaP_largeangle->FindBin(ePlus->Theta()*DEG,ePlus->P()));  
//     acc_em_LA = hacceptance_ThetaP_largeangle->GetBinContent(hacceptance_ThetaP_largeangle->FindBin(eMinus->Theta()*DEG,eMinus->P()));	
// 
//     acc_prot = acc_prot_FA+acc_prot_LA;	  	  
//     acc_ep= acc_ep_LA+acc_ep_FA;
//     acc_em= acc_em_LA+acc_em_FA;
//     acc=acc_prot*acc_ep*acc_em;
//     hcount_acc->Fill(gammaE,weight/nevent/hcount_acc->GetXaxis()->GetBinWidth(1)*count_convert*acc);

//     cout << gammaE << " " << flux  << " " << flux_brem << " " << flux_epa << endl;
    double Erange_flux=(8.8-5.67);
    hflux->Fill(gammaE,flux*Erange_flux/nevent/hflux->GetXaxis()->GetBinWidth(1));
    hflux_brem->Fill(gammaE,flux_brem*Erange_flux/nevent/hflux_brem->GetXaxis()->GetBinWidth(1));
    hflux_epa->Fill(gammaE,flux_epa*Erange_flux/nevent/hflux_epa->GetXaxis()->GetBinWidth(1));
    
//     cout << q->Px() << "\t" << q->Py() << "\t" << q->Pz() << "\t" << q->M() << "\t" << endl;
//     cout << hIn->Px() << "\t" << hIn->Py() << "\t" << hIn->Pz() << "\t" << hIn->M() << "\t" << endl;
//     cout << hOut->Px() << "\t" << hOut->Py() << "\t" << hOut->Pz() << "\t" << hOut->M() << "\t" << endl;
//     cout << ePlus->Px() << "\t" << ePlus->Py() << "\t" << ePlus->Pz() << "\t" << ePlus->M() << "\t" << endl;
//     cout << eMinus->Px() << "\t" << eMinus->Py() << "\t" << eMinus->Pz() << "\t" << eMinus->M() << "\t" << endl;

//     TLorentzVector *CM=0;  
//     *CM=*q+*hIn;
//     TLorentzVector *ePlusCM=0,*eMinusCM=0,*hOutCM=0;
//     *ePlusCM=*ePlus; *eMinusCM=*eMinus; *hOutCM=*hOut;
//     ePlusCM->Boost(CM->BoostVector());
//     eMinusCM->Boost(CM->BoostVector());
//     hOutCM->Boost(CM->BoostVector());

    TLorentzVector CM=*q+*hIn;
    TLorentzVector ePlusCM=*ePlus,eMinusCM=*eMinus,hOutCM=*hOut;
    ePlusCM.Boost(-CM.BoostVector());
    eMinusCM.Boost(-CM.BoostVector());
    hOutCM.Boost(-CM.BoostVector());

//     q->Boost(-CM.BoostVector());
//     hIn->Boost(-CM.BoostVector());    
//     cout << "CM " << q->Px() << "\t" << q->Py() << "\t" << q->Pz() << "\t" << q->M() << "\t" << endl;
//     cout << "CM " << hIn->Px() << "\t" << hIn->Py() << "\t" << hIn->Pz() << "\t" << hIn->M() << "\t" << endl;
    
    hgenCM_ePlus->Fill(ePlusCM.Theta()*DEG,ePlusCM.P());
    hgenCM_eMinus->Fill(eMinusCM.Theta()*DEG,eMinusCM.P());
    hgenCM_hOut->Fill(hOutCM.Theta()*DEG,hOutCM.P());        
    hgen_ePlus->Fill(ePlus->Theta()*DEG,ePlus->P());
    hgen_eMinus->Fill(eMinus->Theta()*DEG,eMinus->P());
    hgen_hOut->Fill(hOut->Theta()*DEG,hOut->P());    
    hacc1_ePlus->Fill(ePlus->Theta()*DEG,ePlus->P(),weight/nevent*count_convert*acc_ePlus);
    hacc1_eMinus->Fill(eMinus->Theta()*DEG,eMinus->P(),weight/nevent*count_convert*acc_eMinus);
    hacc1_hOut->Fill(hOut->Theta()*DEG,hOut->P(),weight/nevent*count_convert*acc_hOut);    
    hacc2_ePlus->Fill(ePlus->Theta()*DEG,ePlus->P(),weight/nevent*count_convert*acc_ePlus*acc_eMinus);
    hacc2_eMinus->Fill(eMinus->Theta()*DEG,eMinus->P(),weight/nevent*count_convert*acc_ePlus*acc_eMinus);
    hacc2_hOut->Fill(hOut->Theta()*DEG,hOut->P(),weight/nevent*count_convert*acc_ePlus*acc_eMinus);    
    hacc3_ePlus->Fill(ePlus->Theta()*DEG,ePlus->P(),weight/nevent*count_convert*acc_ePlus*acc_eMinus*acc_hOut);
    hacc3_eMinus->Fill(eMinus->Theta()*DEG,eMinus->P(),weight/nevent*count_convert*acc_ePlus*acc_eMinus*acc_hOut);
    hacc3_hOut->Fill(hOut->Theta()*DEG,hOut->P(),weight/nevent*count_convert*acc_ePlus*acc_eMinus*acc_hOut);    

    hangle[0][0]->Fill(ePlus->Theta()*DEG,eMinus->Theta()*DEG,weight/nevent*count_convert);
    hangle[0][1]->Fill(ePlus->Theta()*DEG,eMinus->Theta()*DEG,weight/nevent*count_convert*acc_ePlus*acc_eMinus);
    hangle[0][2]->Fill(ePlus->Theta()*DEG,eMinus->Theta()*DEG,weight/nevent*count_convert*acc_ePlus*acc_eMinus*acc_hOut);
    hangle[1][0]->Fill(ePlus->Theta()*DEG,hOut->Theta()*DEG,weight/nevent*count_convert);
    hangle[1][1]->Fill(ePlus->Theta()*DEG,hOut->Theta()*DEG,weight/nevent*count_convert*acc_ePlus*acc_eMinus);
    hangle[1][2]->Fill(ePlus->Theta()*DEG,hOut->Theta()*DEG,weight/nevent*count_convert*acc_ePlus*acc_eMinus*acc_hOut);
    hangle[2][0]->Fill(eMinus->Theta()*DEG,hOut->Theta()*DEG,weight/nevent*count_convert);
    hangle[2][1]->Fill(eMinus->Theta()*DEG,hOut->Theta()*DEG,weight/nevent*count_convert*acc_ePlus*acc_eMinus);
    hangle[2][2]->Fill(eMinus->Theta()*DEG,hOut->Theta()*DEG,weight/nevent*count_convert*acc_ePlus*acc_eMinus*acc_hOut);    
    
    hsigma->Fill(gammaE,weight/0.05971/flux/nevent/hsigma->GetXaxis()->GetBinWidth(1));

    int indexE=int(gammaE-Emin); //E range 5.5,6.5,7.5,8.5
    double Erange_dsigmadt=1.;    //because how indexE defined
    //     cout << gammaE << " " << indexE << endl;

    if (0<=indexE && indexE<n){
	double m1=0,m2=targetM,m3=3.0969,m4=targetM;
	double s=(gammaE+m2)*(gammaE+m2)-gammaE*gammaE;
	double E1cm=(s+m1*m1-m2*m2)/(2.*sqrt(s));
	double E3cm=(s+m3*m3-m4*m4)/(2.*sqrt(s));
	double p1cm=sqrt(E1cm*E1cm-m1*m1);
	double p3cm=sqrt(E3cm*E3cm-m3*m3);
	double t0=pow((m1*m1-m3*m3-m2*m2+m4*m4)/(2.*sqrt(s)),2)-pow(p1cm-p3cm,2);
	double t1=pow((m1*m1-m3*m3-m2*m2+m4*m4)/(2.*sqrt(s)),2)-pow(p1cm+p3cm,2);
	
	hdsigmadt[indexE]->Fill(-t+t0,weight/0.05971/Erange_dsigmadt/flux/nevent/hdsigmadt[indexE]->GetXaxis()->GetBinWidth(1)); 
    }
    
    //note about normalization: weight is total crosssection of photoproduction over a photon range, divide by binwidth is just standard way to convert histogram to curve

}
cout << "end of event" << endl;

  TCanvas *c_acc = new TCanvas("c_acc","c_acc",1900,1000);
  c_acc->Divide(3,5);  
  c_acc->cd(1);        
    hgenCM_ePlus->Draw("colz");
  c_acc->cd(2);        
    hgenCM_eMinus->Draw("colz");
  c_acc->cd(3);        
    hgenCM_hOut->Draw("colz");  
  c_acc->cd(4);        
    hgen_ePlus->Draw("colz");
  c_acc->cd(5);        
    hgen_eMinus->Draw("colz");
  c_acc->cd(6);        
    hgen_hOut->Draw("colz");
  c_acc->cd(7);        
    hacc1_ePlus->Draw("colz");
  c_acc->cd(8);        
    hacc1_eMinus->Draw("colz");
  c_acc->cd(9);        
    hacc1_hOut->Draw("colz");
  c_acc->cd(10);        
    hacc2_ePlus->Draw("colz");
  c_acc->cd(11);        
    hacc2_eMinus->Draw("colz");
  c_acc->cd(12);        
    hacc2_hOut->Draw("colz");
  c_acc->cd(13);        
    hacc3_ePlus->Draw("colz");
  c_acc->cd(14);        
    hacc3_eMinus->Draw("colz");
  c_acc->cd(15);        
    hacc3_hOut->Draw("colz");

  TCanvas *c_flux = new TCanvas("c_flux","c_flux",1900,800);
  hflux->SetLineColor(kBlack);
  hflux->Draw();
  hflux_brem->SetLineColor(kBlue);
  hflux_brem->Draw("same");
  hflux_epa->SetLineColor(kRed);
  hflux_epa->Draw("same");
    TLegend* leg_flux = new TLegend(0.7, 0.75, 0.95, 0.9);  
    leg_flux->AddEntry(hflux_epa,"flux EPA","l");
    leg_flux->AddEntry(hflux_brem,"flux Brem","l");    
    leg_flux->AddEntry(hflux,"flux Brem + EPA","l");    
    leg_flux->Draw();    
  
    
  TCanvas *c_angle = new TCanvas("c_angle","c_angle",1900,1000);
  c_angle->Divide(3,3);
    for(int i=0;i<3;i++){
    for(int k=0;k<3;k++){
        c_angle->cd(i*3+k+1);
        hangle[i][k]->Draw("colz");
    }}
  
  TCanvas *c_comp = new TCanvas("c_comp","c_comp",1900,800);
  c_comp->Divide(2,1);   
  c_comp->cd(1);      
  gPad->SetLogy(1);
  hdsigmadt[0]->SetMinimum(Ymin);  
  hdsigmadt[0]->SetMaximum(Ymax);
  hdsigmadt[0]->SetLineColor(kBlue);
  hdsigmadt[0]->Draw();
  hdsigmadt[1]->SetLineColor(kBlack);
  hdsigmadt[1]->Draw("same");
  hdsigmadt[2]->SetLineColor(kRed);
  hdsigmadt[2]->Draw("same");
    
    //Harry model 6gev
    const Int_t n0 = 16;
    Double_t x0[n0]  = {0.00000000E+00,0.22382910E-02,0.57264695E-02,0.96674839E-02,0.13158775E-01,0.15400605E-01,0.20984303E-01,0.41202025E-01,0.72799675E-01,0.10863564E+00,0.14050741E+00,0.16103693E+00,0.21755150E+00,0.42799070E+00,0.77177306E+00,0.11946204E+01};
    Double_t y0[n0]  = {0.10922545E-05,0.10924661E-05,0.10913836E-05,0.10897937E-05,0.10874356E-05,0.10856548E-05,0.10809539E-05,0.10652939E-05,0.10421790E-05,0.10176951E-05,0.99615133E-06,0.98208074E-06,0.94019130E-06,0.76727688E-06,0.53800241E-06,0.31991053E-06};
    TGraph *g0 = new TGraph(n0,x0,y0);
    g0->SetLineColor(kBlue);
    g0->Draw("same L");    
    // Harry model 7gev
    const Int_t n1 = 19;
    Double_t x1[n1]  = {0.00000000E+00,0.21404081E-02,0.54747490E-02,0.92400523E-02,0.12573995E-01,0.14713950E-01,0.20041058E-01,0.39295505E-01,0.69279391E-01,0.10312313E+00,0.13307617E+00,0.15229555E+00,0.20489792E+00,0.39664295E+00,0.69442498E+00,0.10296783E+01,0.13259881E+01,0.15160929E+01,0.20373169E+01};
    Double_t y1[n1]  = {0.69152152E-05,0.68924536E-05,0.68734699E-05,0.68631303E-05,0.68604722E-05,0.68612741E-05,0.68700797E-05,0.69554308E-05,0.71432671E-05,0.72423032E-05,0.70961239E-05,0.68764421E-05,0.59578511E-05,0.32131615E-05,0.17051351E-05,0.97224242E-06,0.62931268E-06,0.48696801E-06,0.25810302E-06};
    TGraph *g1 = new TGraph(n1,x1,y1);
    g1->SetLineColor(kBlack);
    g1->Draw("same L");  
    //Harry model 8gev
    const Int_t n2 = 19;
    Double_t x2[n2]  = {0.00000000E+00,0.25101144E-02,0.64201443E-02,0.10835190E-01,0.14744125E-01,0.17252996E-01,0.23497934E-01,0.46063521E-01,0.81183953E-01,0.12079647E+00,0.15582928E+00,0.17829525E+00,0.23973167E+00,0.46302514E+00,0.80768785E+00,0.11924419E+01,0.15293837E+01,0.17439027E+01,0.23247762E+01};
    Double_t y2[n2]  = {   0.30364978E-04,0.30080439E-04,0.29847900E-04,0.29682525E-04,0.29566810E-04,0.29494879E-04,0.29292661E-04,0.27924538E-04,0.24086627E-04,0.19634392E-04,0.16540368E-04,0.14942719E-04,0.11635170E-04,0.57423041E-05,0.25948127E-05,0.13157766E-05,0.78623046E-06,0.58203259E-06,0.28173396E-06};
    TGraph *g2 = new TGraph(n2,x2,y2);  
    g2->SetLineColor(kRed);
    g2->Draw("same L");  

    TLegend* leg1 = new TLegend(0.7, 0.75, 0.95, 0.9);  
    leg1->AddEntry(hdsigmadt[0],"E_{#gamma}=6 GeV","l");
    leg1->AddEntry(hdsigmadt[1],"E_{#gamma}=7 GeV","l");  
    leg1->AddEntry(hdsigmadt[2],"E_{#gamma}=8 GeV","l");
    leg1->Draw();    
    
  c_comp->cd(2);
  gPad->SetLogy(1);
  hsigma->SetMinimum(Ymin);  
  hsigma->SetMaximum(Ymax);
  hsigma->Draw();  

  TCanvas *c_count = new TCanvas("c_count","c_count",1200,800);  
  gPad->SetLogy(1);
  hcount->SetLineColor(kBlack);  
  hcount->SetMinimum(1e-1);
  hcount->SetMaximum(2e1);     
  hcount->Draw("HIST E1"); 
  hcount_acc2->SetLineColor(kBlue);  
  hcount_acc2->Draw("HIST E1 same");  
  hcount_acc3->SetLineColor(kRed);  
  hcount_acc3->Draw("HIST E1 same");
    TLegend* leg2 = new TLegend(0.2, 0.8, 0.6, 0.9);  
    leg2->AddEntry(hcount, Form("total counts %.0f, full phasespace",hcount->Integral()),"l");
    leg2->AddEntry(hcount_acc2, Form("total counts %.0f, e^{+}e^{-} detected",hcount_acc2->Integral()),"l");  
    leg2->AddEntry(hcount_acc3, Form("total counts %.0f, e^{+}e^{-}d detected",hcount_acc3->Integral()),"l");
    leg2->Draw();
  
  cout << "hcount " << hcount->Integral() << endl;
  cout << "hcount_acc2 " << hcount_acc2->Integral() << endl;    
  cout << "hcount_acc3 " << hcount_acc3->Integral() << endl;  

}