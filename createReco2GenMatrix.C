#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TH2.h"

void createReco2GenMatrix(){


  TFile *fin = new TFile("histos/PbPbBMC_pt30by3_ipHICalibCentWeight_noTrig.root");

  TTree *t=(TTree*) fin->Get("nt");
  
  double binBounds[7]  = {55,65,80,100,120,150,200};

  TH2F *hRecoVsGen = new TH2F("hRecoVsGen","hRecoVsGen;genJet p_{T} (GeV/c);recoJet p_{T} (GeV/c)",6,binBounds,6,binBounds);
  hRecoVsGen->Sumw2();

  double jtpt, jteta, refpt, discr_ssvHighEff, pthat, weight;
  int refparton_flavorForB, bin, jet65;
  
  t->SetBranchAddress("jtpt",&jtpt);
  t->SetBranchAddress("jteta",&jteta);
  t->SetBranchAddress("refpt",&refpt);
  t->SetBranchAddress("discr_ssvHighEff",&discr_ssvHighEff);
  t->SetBranchAddress("pthat",&pthat);
  t->SetBranchAddress("refparton_flavorForB",&refparton_flavorForB);
  t->SetBranchAddress("bin",&bin);
  t->SetBranchAddress("jet65",&jet65);
  t->SetBranchAddress("weight",&weight);

  for(int i=0;i<t->GetEntries();i++){
    t->GetEntry(i);

    if(!jet65) continue;
    if(abs(refparton_flavorForB)!=5) continue;
    if(discr_ssvHighEff<2)continue;
    if(abs(jteta)>2) continue;
    

    // cutoffs to keep stat errors under control -- still needed?
    if(jtpt>100&&pthat<50) continue;
    if(jtpt>120&&pthat<65) continue;
    if(jtpt>150&&pthat<80) continue;

    hRecoVsGen->Fill(refpt,jtpt,weight);


  }

  TFile *fout=new TFile("outputTowardsFinal/reco2GenMatrix.root","recreate");

  hRecoVsGen->Write();

  TH2F* hRecoVsGenNorm = (TH2F*) hRecoVsGen->Clone("hRecoVsGenNorm");
 
  

  for(int i=0;i<hRecoVsGenNorm->GetNbinsX();i++){

    float nGen = 0;    
    for(int j=0;j<hRecoVsGenNorm->GetNbinsY();j++){
      nGen+=hRecoVsGenNorm->GetBinContent(i+1,j+1);
    }
    for(int j=0;j<hRecoVsGenNorm->GetNbinsY();j++){
      hRecoVsGenNorm->SetBinContent(i+1,j+1,hRecoVsGenNorm->GetBinContent(i+1,j+1)/nGen);
      hRecoVsGenNorm->SetBinError(i+1,j+1,hRecoVsGenNorm->GetBinError(i+1,j+1)/nGen);
    }
  }

  hRecoVsGenNorm->Write();
  fout->Close();



}
