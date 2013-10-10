
#include <iostream>
#include <cstdio>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TChain.h"

using namespace std;

void JetTrgNorm(){

  gStyle->SetPadLeftMargin(0.2);

  TFile *fin = new TFile("histos/ppdata_ppReco_ak3PF_DinkoXchk_jetTrig_noIPupperCut.root");
  //TChain *nt = new TChain("hltanalysis/HltTree");
  //TChain *jt = new TChain("ak3PFJetAnalyzer/t");
  //nt->Add("/mnt/hadoop/store/user/kjung/pPb_BForest/hiForestTree_1*");
  //jt->Add("/mnt/hadoop/store/user/kjung/pPb_BForest/hiForestTree_1*");
  //nt->AddFriend(jt);
  TTree *nt = (TTree*)fin->Get("nt");
  if(!nt){ cout << "Can't get tree!" << endl; exit(0); }

  TH1D *jetHist[6];
  TH1D *jetNorm[6];

  jetHist[0] = new TH1D("jet20","",50,30,250);
  jetHist[1] = new TH1D("jet40","",50,30,250);
  jetHist[2] = new TH1D("jet60","",50,30,250);
  jetHist[3] = new TH1D("jet80","",50,30,250);
  jetHist[4] = new TH1D("jet100","",50,30,250);
  jetHist[5] = new TH1D("jet120","",50,30,250);
  
  jetNorm[0] = new TH1D("jet20Comb","",50,30,250);
  jetNorm[1] = new TH1D("jet40Comb","",50,30,250);
  jetNorm[2] = new TH1D("jet60Comb","",50,30,250);
  jetNorm[3] = new TH1D("jet80Comb","",50,30,250);
  jetNorm[4] = new TH1D("jet100Comb","",50,30,250); //Not used since ps(Jet80)=1.  Only here for completeness.
  jetNorm[5] = new TH1D("jet120Comb","",50,30,250);
  
  for(int i=0; i<6; i++){
    jetHist[i]->Sumw2();
    jetNorm[i]->Sumw2();
  }

  Int_t JetTrg20, JetTrg40, JetTrg60, JetTrg80, JetTrg100, JetTrg120;
  Double_t jtpt;
  Double_t jteta;

  nt->SetBranchAddress("HLT_Jet20_noJetID_v1",&JetTrg20);
  nt->SetBranchAddress("HLT_Jet40_noJetID_v1",&JetTrg40);
  nt->SetBranchAddress("HLT_Jet60_noJetID_v1",&JetTrg60);
  nt->SetBranchAddress("HLT_Jet80_noJetID_v1",&JetTrg80);
  nt->SetBranchAddress("HLT_Jet100_noJetID_v1",&JetTrg100);
  nt->SetBranchAddress("HLT_Jet120_noJetID_v1",&JetTrg120);
  nt->SetBranchAddress("jtpt",&jtpt);
  nt->SetBranchAddress("jteta",&jteta);

  double ov1 = nt->GetEntries("HLT_Jet20_noJetID_v1 && HLT_Jet80_noJetID_v1");
  cout <<"ov1: "<< ov1 << endl;
  double ov2 = nt->GetEntries("HLT_Jet40_noJetID_v1 && HLT_Jet80_noJetID_v1");
  cout <<"ov2: "<< ov2 << endl;
  double ov3 = nt->GetEntries("HLT_Jet60_noJetID_v1 && HLT_Jet80_noJetID_v1");
  cout <<"ov3: "<< ov3 << endl;
  double ov4 = nt->GetEntries("HLT_Jet80_noJetID_v1 && HLT_Jet80_noJetID_v1");
  cout <<"ov4: "<< ov4 << endl;
  double ov5 = nt->GetEntries("HLT_Jet100_noJetID_v1 && HLT_Jet80_noJetID_v1");
  cout <<"ov5: "<< ov5 << endl;
  double ov6 = nt->GetEntries("HLT_Jet80_noJetID_v1");
  cout <<"ov6: "<< ov6 << endl;

  double ps1 = ov6/ov1;
  double ps2 = ov6/ov2;
  double ps3 = ov6/ov3;
  double ps4 = ov6/ov4;
  double ps5 = ov6/ov5;

  for(int i=0; i<nt->GetEntries(); i++){
    nt->GetEntry(i);
    
    if(fabs(jteta)<3){
      //Fill Trigger Combination Spectra
      //if(JetTrg20 && !JetTrg40 && !JetTrg60 && !JetTrg80) jetNorm[0]->Fill(jtpt, 1./(1./ps1));
      if(JetTrg40 && !JetTrg60 && !JetTrg80 /*&& !JetTrg100 && !JetTrg120*/) jetNorm[1]->Fill(jtpt,(1./(1./ps2)));
      if(JetTrg60 && !JetTrg80 /*&& !JetTrg100 && !JetTrg120*/) jetNorm[2]->Fill(jtpt);//, 1./(1./ps2 + 1./ps3 - (1./(ps2*ps3))));
      if(JetTrg80 /*&& !JetTrg100 && !JetTrg120*/) jetNorm[3]->Fill(jtpt);//, 1./(1./ps2 + 1./ps3 + 1./ps4 - (1./(ps2*ps3)) - (1./(ps3*ps4)) - (1./(ps2*ps4)) + 1./(ps2*ps3*ps4)));
      //if(JetTrg100 && !JetTrg120) jetNorm[4]->Fill(jtpt, 1./(1./ps2 + 1./ps3 + 1./ps4 + 1./ps5 - (1./(ps2*ps3)) - (1./(ps2*ps4)) - 1./(ps2*ps5) - 1./(ps3*ps4) - 1./(ps3*ps5) - 1./(ps4*ps5) + 1./(ps2*ps3*ps4) + 1./(ps2*ps3*ps5) + 1./(ps2*ps4*ps5) + 1./(ps3*ps4*ps5) - 1./(ps2*ps3*ps4*ps5)));
      //if(JetTrg120) jetNorm[5]->Fill(jtpt);

      //Fill Trigger Normalization cross-check spectra
      if(JetTrg20) jetHist[0]->Fill(jtpt);
      if(JetTrg40) jetHist[1]->Fill(jtpt);
      if(JetTrg60) jetHist[2]->Fill(jtpt);
      if(JetTrg80) jetHist[3]->Fill(jtpt);
      if(JetTrg100) jetHist[4]->Fill(jtpt);
    }
  }

  TH1D *jetCln = (TH1D*)jetNorm[5]->Clone("jetCln");
  for(int i=1; i<4; i++){
    jetNorm[4]->Add(jetNorm[i]);
  }

  jetHist[0]->Scale(1/ov1);
  jetHist[1]->Scale(1/ov2);
  jetHist[2]->Scale(1/ov3);
  jetHist[3]->Scale(1/ov4);
  jetHist[4]->Scale(1/ov5);
  jetHist[5]->Scale(1/ov6);

  //Scale by bin width to get dN/dpt
  /* for(int i=1; i<5; i++){
    for(int j=1; j<=jetHist[i]->GetNbinsX(); j++){
      jetHist[i]->SetBinContent(j, jetHist[i]->GetBinContent(j)/jetHist[i]->GetBinWidth(j));
      jetHist[i]->SetBinError(j, jetHist[i]->GetBinError(j)/jetHist[i]->GetBinWidth(j));
    }
    }*/
      
  for(int i=1; i<6; i++){
    jetHist[i]->SetLineColor(pow(2,i));
    jetHist[i]->SetMarkerColor(pow(2,i));
  }
  jetHist[4]->SetMarkerColor(5);

  TCanvas *c1 = new TCanvas("c1","",800,600);
  c1->Divide(1,2);
  c1->cd(1);
  gPad->SetLogy();

  jetHist[0]->SetYTitle("1/N_{Jet80Trg} dN/dp_{T}");
  jetHist[0]->SetXTitle("Jet p_{T}");
  jetHist[0]->Draw();
  for(int i=1; i<4; i++){
    jetHist[i]->Draw("same");
  }
  TLegend *leg = new TLegend(0.6,0.5,0.9,0.9);
  leg->SetFillColor(0);
  //leg->AddEntry(jetHist[0],"Jet20 Trg");
  leg->AddEntry(jetHist[1],"Jet40 Trg");
  leg->AddEntry(jetHist[2],"Jet60 Trg");
  leg->AddEntry(jetHist[3],"Jet80 Trg");
  //leg->AddEntry(jetHist[4],"Jet100 Trg");
  //leg->AddEntry(jetHist[5],"Jet120 Trg");
  leg->Draw();

  c1->cd(2);
  TH1D *cln[5];
  char* histoname = new char[20];
  for(int i=1; i<6; i++){
    sprintf(histoname,"%s%d","clone_",i);
    cln[i] = (TH1D*)jetHist[i]->Clone(histoname);
    cln[i]->Divide(jetHist[i-1]);
  }

  cln[1]->SetXTitle("Jet p_{T}");
  cln[1]->SetYTitle("Efficiency");
  cln[1]->Draw();
  for(int i=1; i<5; i++){
    cln[i]->Draw("same");
  }

  TCanvas *c2 = new TCanvas("c2","",800,600);
  c2->cd();
  c2->SetLogy();
  //c2->SetLogx();
  jetNorm[4]->SetMinimum(10);
  jetNorm[4]->SetMaximum(1E8);
  jetNorm[4]->SetXTitle("Jet p_{T} (GeV/c)");
  jetNorm[4]->SetYTitle("dN/dp_{T} (GeV/c)^{-1}");
  //jetNorm[5]->Scale(1./12.);
  //double tempInt = jetNorm[4]->Integral();
  //jetNorm[4]->Scale(1./jetNorm[4]->Integral());
  jetNorm[4]->Draw();
  jetNorm[3]->SetMarkerColor(2);
  jetNorm[3]->SetLineColor(2);
  jetNorm[2]->SetMarkerColor(4);
  jetNorm[2]->SetLineColor(4);
  jetNorm[1]->SetMarkerColor(8);
  jetNorm[1]->SetLineColor(8);		     
  for(int i=1; i<4; i++){
    //jetNorm[i]->SetMarkerColor(2+(i*2));
    //jetNorm[i]->SetLineColor(2+(i*2));
    //jetNorm[i]->Scale(1./12.);
    //jetNorm[i]->Scale(1./tempInt);
    jetNorm[i]->Draw("same");
    
  }
  jetCln->SetMarkerColor(kGreen-2);
  jetCln->SetLineColor(kGreen-2);
  //jetCln->Draw("same");
  jetNorm[4]->Draw("same");
  //jetHist[1]->Draw("same");

  TLegend *leg1 = new TLegend(0.577,0.613,0.888,0.914);
  leg1->SetFillColor(0);
  leg1->AddEntry(jetNorm[4], "Combined");
  //leg1->AddEntry(jetNorm[0], "Jet20_NoID Trg");
  leg1->AddEntry(jetNorm[3], "ppJet80");
  leg1->AddEntry(jetNorm[2], "ppJet60");
  leg1->AddEntry(jetNorm[1], "ppJet40");
  //leg1->AddEntry(jetNorm[4], "Jet100 Trigger");
  //leg1->AddEntry(jetCln, "Jet120 Trigger");
  leg1->Draw();
}
