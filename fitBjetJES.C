#include <iostream>
#include "TProfile.h"
#include "TFile.h"
#include "TStyle.h"
#include "TNtuple.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TLegend.h"

void fitBjetJES(int ppPbPb=1, int cbinlo=12, int cbinhi=40){

  if(!ppPbPb){
    cbinlo=0;
    cbinhi=40;
  }

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  TFile *fL;
  
  if(!ppPbPb)fL=new TFile("histos/ppMC_hiReco_jetTrig_highPurity_JEC.root");
  else     fL=new TFile("histos/PbPbQCDMC_JEC.root");

  // these are dummy files for pp
  TFile *fB=new TFile("histos/PbPbBMC_JEC.root");
  TFile *fC=new TFile("histos/PbPbCMC_JEC.root");


  TNtuple *tL = (TNtuple*) fL->Get("nt");
  TNtuple *tB = (TNtuple*) fB->Get("nt");
  TNtuple *tC = (TNtuple*) fC->Get("nt");
  
  float jtptL, refptL, jtetaL, weightL, refparton_flavorForBL, binL;

  tL->SetBranchAddress("jtpt",&jtptL);
  tL->SetBranchAddress("jteta",&jtetaL);
  tL->SetBranchAddress("refpt",&refptL);
  tL->SetBranchAddress("weight",&weightL);
  if(ppPbPb)tL->SetBranchAddress("bin",&binL);
  tL->SetBranchAddress("refparton_flavorForB",&refparton_flavorForBL);

  float jtptB, refptB, jtetaB, weightB, refparton_flavorForBB, binB;

  tB->SetBranchAddress("jtpt",&jtptB);
  tB->SetBranchAddress("jteta",&jtetaB);
  tB->SetBranchAddress("refpt",&refptB);
  tB->SetBranchAddress("weight",&weightB);
  if(ppPbPb)tB->SetBranchAddress("bin",&binB);
  tB->SetBranchAddress("refparton_flavorForB",&refparton_flavorForBB);

  float jtptC, refptC, jtetaC, weightC, refparton_flavorForBC, binC;

  tC->SetBranchAddress("jtpt",&jtptC);
  tC->SetBranchAddress("jteta",&jtetaC);
  tC->SetBranchAddress("refpt",&refptC);
  tC->SetBranchAddress("weight",&weightC);
  if(ppPbPb)tC->SetBranchAddress("bin",&binC);
  tC->SetBranchAddress("refparton_flavorForB",&refparton_flavorForBC);

  TProfile  *hL = new TProfile("hL","hL",250,50,300,0,10);
  TProfile  *hB = new TProfile("hB","hB",250,50,300,0,10);
  TProfile  *hC = new TProfile("hC","hC",250,50,300,0,10);
  hL->Sumw2(),hB->Sumw2(),hC->Sumw2();


  for(int i=0;i<tL->GetEntries();i++){
    tL->GetEntry(i);
    if(!ppPbPb) binL=39;

    if(fabs(jtetaL)<2 && binL>=cbinlo && binL<cbinhi)
      hL->Fill(refptL,jtptL/refptL,weightL); 


    if(!ppPbPb){
      if(fabs(jtetaL)<2 && binL>=cbinlo && binL<cbinhi && abs(refparton_flavorForBL)==5)
	hB->Fill(refptL,jtptL/refptL,weightL);
      
      if(fabs(jtetaL)<2 && binL>=cbinlo && binL<cbinhi && abs(refparton_flavorForBL)==4)
	hC->Fill(refptL,jtptL/refptL,weightL);      
    }
  }

  if(ppPbPb){
    for(int i=0;i<tB->GetEntries();i++){
      tB->GetEntry(i);
      if(fabs(jtetaB)<2 && binB>=cbinlo && binB<cbinhi && abs(refparton_flavorForBB)==5)
	hB->Fill(refptB,jtptB/refptB,weightB);
    }
    for(int i=0;i<tC->GetEntries();i++){
      tC->GetEntry(i);
      if(fabs(jtetaC)<2 && binC>=cbinlo && binC<cbinhi && abs(refparton_flavorForBC)==4)
	hC->Fill(refptC,jtptC/refptC,weightC);
    }
  }

  
 
  hL->SetMinimum(0.);
  
  hL->SetLineColor(kBlue);
  hB->SetLineColor(kRed);
  hC->SetLineColor(kGreen);

  hL->SetMarkerColor(kBlue);
  hB->SetMarkerColor(kRed);
  hC->SetMarkerColor(kGreen);
  
  //hL->SetMarkerStyle(4);
  //hB->SetMarkerStyle(4);
  //hC->SetMarkerStyle(4);
  
  hL->SetXTitle("genJet p_{T} (GeV/c)");
  hL->SetYTitle("<reco p_{T} / gen p_{T} >");

  hL->GetXaxis()->SetRangeUser(50.,199.999);
  hL->GetYaxis()->SetRangeUser(0.5,1.05);
  
  TCanvas *c1=new TCanvas("c1","c1",800,600);
  c1->SetGridx(1);
  c1->SetGridy(1);

  hL->Draw("e1");
  hB->Draw("e1,same");
  hC->Draw("e1,same");

  TLegend *leg=new TLegend(0.4,0.15,0.9,0.45);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  if(ppPbPb&&cbinlo==0&&cbinhi==40)leg->SetHeader("Pythia+Hydjet, 0-100%");
  leg->AddEntry(hL,"Inclusive jets","pl");
  leg->AddEntry(hC,"c-jets","pl");
  leg->AddEntry(hB,"b-jets","pl");
  leg->Draw();

  TCanvas *c2=new TCanvas("c2","c2",1);
  /*
  TH1F *hL2 = (TH1F*)hL->Clone("hL2");
  TH1F *hB2 = (TH1F*)hB->Clone("hB2");
  hL2->Add(hB2,-1);
  hL2->Draw();
  */

  TH1F  *hcorr = new TH1F("hcorr","hcorr",250,50,300);
  hcorr->Sumw2();

  for(int i=0;i<hL->GetNbinsX();i++){
    cout<<" b resp "<<hB->GetBinContent(i+1)<<endl;
    cout<<" l resp "<<hL->GetBinContent(i+1)<<endl;
    cout<<" l offset "<<1.-hL->GetBinContent(i+1)<<endl;
    cout<<" corrected b resp "<<hB->GetBinContent(i+1)+1.-hL->GetBinContent(i+1)<<endl;
    float jesOffset = 1.-hL->GetBinContent(i+1);

    hcorr->SetBinContent(i+1,hB->GetBinContent(i+1)+jesOffset);

    hcorr->SetBinError(i+1,sqrt(hB->GetBinError(i+1)*hB->GetBinError(i+1)+hL->GetBinError(i+1)*hL->GetBinError(i+1)));


  }

  hcorr->SetMinimum(0.5);
  hcorr->SetMaximum(1.1);
      
  hcorr->SetLineColor(kRed);
  hcorr->SetMarkerColor(kRed);
  hcorr->SetMarkerStyle(4);
  hcorr->Draw();

  TF1 *fCorr = new TF1("fCorr","[0]+[1]*log(x)+[2]*log(x)*log(x)",50,300);
  fCorr->SetLineWidth(1);
  fCorr->SetLineColor(kBlue);
  hcorr->Fit(fCorr);

  TFile *fout;
  if(ppPbPb) fout =new TFile(Form("bJEShistos/bJetScale_PbPb_Cent_fineBin_%d_%d.root",cbinlo,cbinhi),"recreate");
  else fout =new TFile("bJEShistos/bJetScale_PP_fineBin.root","recreate");
  hcorr->Write();
  fCorr->Write();
  fout->Close();

}

