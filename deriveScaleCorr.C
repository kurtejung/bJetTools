#include "TF1.h"
#include "TH1.h"
#include "TRandom.h"
#include "TCanvas.h"

void deriveScaleCorr(long ntries=100000000){


  TF1 *dndpt = new TF1("dndpt","pow(x,-5.95)",40,250);
  TF1 *bscaleCorr = new TF1("bscaleCorr","[0]+[1]*log(x)+[2]*log(x)*log(x)",40,250);
  bscaleCorr->SetParameters(5.23359e-01,1.15226e-01,-7.96949e-03);

  TH1F *hBdndpt = new TH1F("hBdndpt","hBdndpt",210,40,200);
  TH1F *hBdndptNew = new TH1F("hBdndptNew","hBdndptNew",210,40,200);
  //bscaleCorr->Draw();


  TRandom *trand=new TRandom();
  for(int i=0;i<ntries;i++){
    double pt = (250.-40.)*trand->Rndm()+40.;
    double newpt = pt*bscaleCorr->Eval(pt);
    hBdndpt->Fill(pt,dndpt->Eval(pt));
    hBdndptNew->Fill(newpt,dndpt->Eval(pt));
  }

  hBdndpt->Draw();
  hBdndptNew->SetLineColor(2);
  hBdndptNew->Draw("same");

  TCanvas *c2=new TCanvas("c2","c2",1);
  TH1F *hratio = (TH1F*)hBdndpt->Clone("hratio");
  hratio->Divide(hBdndptNew);
  hratio->Draw();


  cout<<" pp bins "<<endl;
  cout<<" correction in 60,70 = "<<hBdndpt->Integral(60,70)/hBdndptNew->Integral(60,70)<<endl;
  cout<<" correction in 70,80 = "<<hBdndpt->Integral(70,80)/hBdndptNew->Integral(70,80)<<endl;
  cout<<" correction in 80,95 = "<<hBdndpt->Integral(80,95)/hBdndptNew->Integral(80,95)<<endl;
  cout<<" correction in 95,120 = "<<hBdndpt->Integral(95,120)/hBdndptNew->Integral(95,120)<<endl; 
  cout<<" correction in 120,200 = "<<hBdndpt->Integral(120,200)/hBdndptNew->Integral(120,200)<<endl;

  cout<<" PbPb bins "<<endl;
  cout<<" correction in 60,80 = "<<hBdndpt->Integral(60,80)/hBdndptNew->Integral(60,80)<<endl;
  cout<<" correction in 80,100 = "<<hBdndpt->Integral(80,100)/hBdndptNew->Integral(80,100)<<endl;
  cout<<" correction in 100,120 = "<<hBdndpt->Integral(100,120)/hBdndptNew->Integral(100,120)<<endl;
  cout<<" correction in 120,150 = "<<hBdndpt->Integral(120,150)/hBdndptNew->Integral(120,150)<<endl;
  cout<<" correction in 150,200 = "<<hBdndpt->Integral(150,200)/hBdndptNew->Integral(150,200)<<endl; 


  //hBdndpt->Scale(dndpt->Integral(60,250)/hBdndpt->Integral(21,250));

  //dndpt->Draw("same");
  //TF1 *outfit = new TF1("outfit","[0]*pow(x,[1])",60,200);
  //outfit->SetParameters(10.,-6.);
  //hBdndpt->Fit(outfit,"R");

}
