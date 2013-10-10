
#include "TTree.h"
#include "TH1D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"

void PUeffect(){
  
  TFile *f1 = new TFile("histos/ppMC_ppReco_ak3PF_QCDjetTrig_noIPupperCut.root");
  //TFile *f1 = new TFile("histos/ppdata_ppReco_ak3PF_jetTrig_noIPupperCut.root");
  TTree *nt = (TTree*)f1->Get("nt");

  TH1D *h_QCDpu = new TH1D("h_QCDpu","",50,0,600); h_QCDpu->Sumw2();
  TH1D *h_QCDnoPU = new TH1D("h_QCDnoPU","",50,0,600); h_QCDnoPU->Sumw2();
  TH1D *h_Bpu = new TH1D("h_Bpu","",50,0,600); h_Bpu->Sumw2();
  TH1D *h_BnoPU = new TH1D("h_BnoPU","",50,0,600); h_BnoPU->Sumw2();
  
  nt->Draw("jtpt>>h_QCDpu","weight*(pVertexFilterCutGplusUpsPP)");
  nt->Draw("jtpt>>h_QCDnoPU","weight*(1)");
  nt->Draw("jtpt>>h_Bpu","weight*(pVertexFilterCutGplusUpsPP && discr_ssvHighEff>2 )");
  nt->Draw("jtpt>>h_BnoPU","weight*(discr_ssvHighEff>2)");

  h_QCDpu->Divide(h_QCDpu,h_QCDnoPU,1,1,"B");
  h_Bpu->Divide(h_Bpu,h_BnoPU,1,1,"B");

  TCanvas *c1 = new TCanvas("c1","",800,600);
  c1->cd();
  h_QCDpu->SetXTitle("reco p_{T}");
  h_QCDpu->SetYTitle("PU Filter / No PU Filter [MC]");
  h_QCDpu->SetMaximum(1.3);
  h_QCDpu->SetMinimum(0.4);
  h_QCDpu->Draw();
  h_Bpu->SetMarkerColor(2);
  h_Bpu->SetLineColor(2);
  h_Bpu->Draw("same");

  TLegend *l1 = new TLegend(0.5,0.5,0.8,0.8);
  l1->AddEntry(h_QCDpu,"QCD Jets");
  l1->AddEntry(h_Bpu,"B Jets");
  l1->Draw();
}
