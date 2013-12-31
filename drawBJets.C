
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"

using namespace std;

void stackHistogram(TH1D **thist, int nHist){
  thist[0]->SetFillColor(6);
  thist[0]->SetLineColor(6);
  thist[1]->SetFillColor(2);
  thist[1]->SetLineColor(2);
  thist[2]->SetFillColor(3);
  thist[2]->SetLineColor(3);
  thist[3]->SetFillColor(8);
  thist[3]->SetLineColor(8);
  thist[4]->SetFillColor(4);
  thist[4]->SetLineColor(4);

  for(int i=0; i<nHist; i++){
    thist[i]->SetMarkerStyle(1);
    //thist[i]->Scale(1./thist[i]->Integral());
  }
  for(int i=0; i<nHist; i++){
    int j=i+1;
    while(j<5){
      thist[i]->Add(thist[j]);
      j++;
    }
  }
}

void drawBJets(int leadjet=65, int subleadjet=30){

  TFile *f1 = new TFile("histos/ppMC_ppReco_jetTrig_Dijet_noIPupperCut.root");
  TTree *nt = (TTree*)f1->Get("nt");
  TFile *f2 = new TFile("histos/ppdata_ppReco_jetTrig_Dijet_noIPupperCut.root");
  TTree *data = (TTree*)f2->Get("nt");

  char* histoname = new char[40];
  // [0] = b-bbar, [1] = b+X, [2] = c-cbar, [3] = c+X, [4] = udsg jets, [5] = data
  TH1D *X_all[6];
  TH1D *dphi_all[6];
  TH1D *X_leadB[6];
  TH1D *dphi_leadB[6];
  TH1D *X_subleadB[6];
  TH1D *dphi_subleadB[6];
  TH1D *X_doubleB[6];
  TH1D *dphi_doubleB[6];

  for(int i=0; i<6; i++){
    sprintf(histoname,"%s%d","X_all_",i);
    X_all[i] = new TH1D(histoname,"",10,0,1);
    X_all[i]->Sumw2();
    sprintf(histoname,"%s%d","dphi_all_",i);
    dphi_all[i] = new TH1D(histoname,"",30,0,3.14);
    dphi_all[i]->Sumw2();
    sprintf(histoname,"%s%d","X_leadB_",i);
    X_leadB[i] = new TH1D(histoname,"",10,0,1);
    X_leadB[i]->Sumw2();
    sprintf(histoname,"%s%d","dphi_leadB_",i);
    dphi_leadB[i] = new TH1D(histoname,"",30,0,3.14);
    dphi_leadB[i]->Sumw2();
    sprintf(histoname,"%s%d","X_subleadB_",i);
    X_subleadB[i] = new TH1D(histoname,"",10,0,1);
    X_subleadB[i]->Sumw2();
    sprintf(histoname,"%s%d","dphi_subleadB_",i);
    dphi_subleadB[i] = new TH1D(histoname,"",30,0,3.14);
    dphi_subleadB[i]->Sumw2();
    sprintf(histoname,"%s%d","X_doubleB_",i);
    X_doubleB[i] = new TH1D(histoname,"",10,0,1);
    X_doubleB[i]->Sumw2();
    sprintf(histoname,"%s%d","dphi_doubleB_",i);
    dphi_doubleB[i] = new TH1D(histoname,"",30,0,3.14);
    dphi_doubleB[i]->Sumw2();
  }

  double jtpt1, jtpt2, jtphi1, jtphi2, refparton_flavorForB1, refparton_flavorForB2, discr_prob1, discr_prob2, weight;
  double djtpt1, djtpt2, djtphi1, djtphi2, drefparton_flavorForB1, drefparton_flavorForB2, ddiscr_prob1, ddiscr_prob2;

  nt->SetBranchAddress("jtpt1",&jtpt1);
  nt->SetBranchAddress("jtpt2",&jtpt2);
  nt->SetBranchAddress("jtphi1",&jtphi1);
  nt->SetBranchAddress("jtphi2",&jtphi2);
  nt->SetBranchAddress("refparton_flavorForB1",&refparton_flavorForB1);
  nt->SetBranchAddress("refparton_flavorForB2",&refparton_flavorForB2);
  nt->SetBranchAddress("discr_prob1",&discr_prob1);
  nt->SetBranchAddress("discr_prob2",&discr_prob2);
  nt->SetBranchAddress("weight",&weight);

  data->SetBranchAddress("jtpt1",&djtpt1);
  data->SetBranchAddress("jtpt2",&djtpt2);
  data->SetBranchAddress("jtphi1",&djtphi1);
  data->SetBranchAddress("jtphi2",&djtphi2);
  data->SetBranchAddress("refparton_flavorForB1",&drefparton_flavorForB1);
  data->SetBranchAddress("refparton_flavorForB2",&drefparton_flavorForB2);
  data->SetBranchAddress("discr_prob1",&ddiscr_prob1);
  data->SetBranchAddress("discr_prob2",&ddiscr_prob2);

  cout << "Starting data..." << endl;
  for(int ientry=0; ientry<data->GetEntries(); ientry++){
    data->GetEntry(ientry);

    //Logic tree for Data
    if(djtpt1>leadjet && djtpt2>subleadjet){
      dphi_all[5]->Fill(acos(cos(djtphi1-djtphi2)));
      if(acos(cos(djtphi1-djtphi2))>2.0944) X_all[5]->Fill(djtpt2/djtpt1);
      
      if(ddiscr_prob1 > 0.6){
	dphi_leadB[5]->Fill(acos(cos(djtphi1-djtphi2)));
	if(acos(cos(djtphi1-djtphi2))>2.0944) X_leadB[5]->Fill(djtpt2/djtpt1);
      }
      if(ddiscr_prob2 > 0.6){
	dphi_subleadB[5]->Fill(acos(cos(djtphi1-djtphi2)));
	if(acos(cos(djtphi1-djtphi2))>2.0944) X_subleadB[5]->Fill(djtpt2/djtpt1);
      }
      if(ddiscr_prob1 > 0.6 && ddiscr_prob2 > 0.6){
	dphi_doubleB[5]->Fill(acos(cos(djtphi1-djtphi2)));
	if(acos(cos(djtphi1-djtphi2))>2.0944) X_doubleB[5]->Fill(djtpt2/djtpt1);
      }
    }
  }

  cout << "Starting MC..." << endl;
  for(int ientry=0; ientry<nt->GetEntries(); ientry++){
    nt->GetEntry(ientry);

    //Logic tree for All Jets
    if(jtpt1>leadjet && jtpt2>subleadjet){
      if(fabs(refparton_flavorForB1)==5){
	if(fabs(refparton_flavorForB2)==5){
	  dphi_all[0]->Fill(acos(cos(jtphi1-jtphi2)), weight);
	  if(acos(cos(jtphi1-jtphi2))>2.0944) X_all[0]->Fill(jtpt2/jtpt1, weight);
	}
	else{
	  dphi_all[1]->Fill(acos(cos(jtphi1-jtphi2)), weight);
	  if(acos(cos(jtphi1-jtphi2))>2.0944) X_all[1]->Fill(jtpt2/jtpt1, weight);
	}
      }
      else if(fabs(refparton_flavorForB1)==4){
	if(fabs(refparton_flavorForB2)==5){
	  dphi_all[1]->Fill(acos(cos(jtphi1-jtphi2)), weight);
	  if(acos(cos(jtphi1-jtphi2))>2.0944) X_all[1]->Fill(jtpt2/jtpt1, weight);
	}
	else if(fabs(refparton_flavorForB2)==4){
	  dphi_all[2]->Fill(acos(cos(jtphi1-jtphi2)), weight);
	  if(acos(cos(jtphi1-jtphi2))>2.0944) X_all[2]->Fill(jtpt2/jtpt1, weight);
	}
	else{
	  dphi_all[3]->Fill(acos(cos(jtphi1-jtphi2)), weight);
	  if(acos(cos(jtphi1-jtphi2))>2.0944) X_all[3]->Fill(jtpt2/jtpt1, weight);
	}
      }
      else{
	if(fabs(refparton_flavorForB1)==4){
	  dphi_all[3]->Fill(acos(cos(jtphi1-jtphi2)), weight);
	  if(acos(cos(jtphi1-jtphi2))>2.0944) X_all[3]->Fill(jtpt2/jtpt1, weight);
	}
	else{
	  dphi_all[4]->Fill(acos(cos(jtphi1-jtphi2)), weight);
	  if(acos(cos(jtphi1-jtphi2))>2.0944) X_all[4]->Fill(jtpt2/jtpt1, weight);
	}
      }
    }

      //Logic tree for Leading B-tagged Jet
    if(jtpt1>leadjet && jtpt2>subleadjet && discr_prob1 > 0.6){
      if(fabs(refparton_flavorForB1)==5){
	if(fabs(refparton_flavorForB2)==5){
	  dphi_leadB[0]->Fill(acos(cos(jtphi1-jtphi2)), weight);
	  if(acos(cos(jtphi1-jtphi2))>2.0944) X_leadB[0]->Fill(jtpt2/jtpt1, weight);
	}
	else{
	  dphi_leadB[1]->Fill(acos(cos(jtphi1-jtphi2)), weight);
	  if(acos(cos(jtphi1-jtphi2))>2.0944) X_leadB[1]->Fill(jtpt2/jtpt1, weight);
	}
      }
      else if(fabs(refparton_flavorForB1)==4){
	if(fabs(refparton_flavorForB2)==5){
	  dphi_leadB[1]->Fill(acos(cos(jtphi1-jtphi2)), weight);
	  if(acos(cos(jtphi1-jtphi2))>2.0944) X_leadB[1]->Fill(jtpt2/jtpt1, weight);
	}
	else if(fabs(refparton_flavorForB2)==4){
	  dphi_leadB[2]->Fill(acos(cos(jtphi1-jtphi2)), weight);
	  if(acos(cos(jtphi1-jtphi2))>2.0944) X_leadB[2]->Fill(jtpt2/jtpt1, weight);
	}
	else{
	  dphi_leadB[3]->Fill(acos(cos(jtphi1-jtphi2)), weight);
	  if(acos(cos(jtphi1-jtphi2))>2.0944) X_leadB[3]->Fill(jtpt2/jtpt1, weight);
	}
      }
      else{
	if(fabs(refparton_flavorForB1)==4){
	  dphi_leadB[3]->Fill(acos(cos(jtphi1-jtphi2)), weight);
	  if(acos(cos(jtphi1-jtphi2))>2.0944) X_leadB[3]->Fill(jtpt2/jtpt1, weight);
	}
	else{
	  dphi_leadB[4]->Fill(acos(cos(jtphi1-jtphi2)), weight);
	  if(acos(cos(jtphi1-jtphi2))>2.0944) X_leadB[4]->Fill(jtpt2/jtpt1, weight);
	}
      }
    }

    //Logic tree for Subleading B-tagged Jet
    if(jtpt1>leadjet && jtpt2>subleadjet && discr_prob2 > 0.6){
      if(fabs(refparton_flavorForB1)==5){
	if(fabs(refparton_flavorForB2)==5){
	  dphi_subleadB[0]->Fill(acos(cos(jtphi1-jtphi2)), weight);
	  if(acos(cos(jtphi1-jtphi2))>2.0944) X_subleadB[0]->Fill(jtpt2/jtpt1, weight);
	}
	else{
	  dphi_subleadB[1]->Fill(acos(cos(jtphi1-jtphi2)), weight);
	  if(acos(cos(jtphi1-jtphi2))>2.0944) X_subleadB[1]->Fill(jtpt2/jtpt1, weight);
	}
      }
      else if(fabs(refparton_flavorForB1)==4){
	if(fabs(refparton_flavorForB2)==5){
	  dphi_subleadB[1]->Fill(acos(cos(jtphi1-jtphi2)), weight);
	  if(acos(cos(jtphi1-jtphi2))>2.0944) X_subleadB[1]->Fill(jtpt2/jtpt1, weight);
	}
	else if(fabs(refparton_flavorForB2)==4){
	  dphi_subleadB[2]->Fill(acos(cos(jtphi1-jtphi2)), weight);
	  if(acos(cos(jtphi1-jtphi2))>2.0944) X_subleadB[2]->Fill(jtpt2/jtpt1, weight);
	}
	else{
	  dphi_subleadB[3]->Fill(acos(cos(jtphi1-jtphi2)), weight);
	  if(acos(cos(jtphi1-jtphi2))>2.0944) X_subleadB[3]->Fill(jtpt2/jtpt1, weight);
	}
      }
      else{
	if(fabs(refparton_flavorForB1)==4){
	  dphi_subleadB[3]->Fill(acos(cos(jtphi1-jtphi2)), weight);
	  if(acos(cos(jtphi1-jtphi2))>2.0944) X_subleadB[3]->Fill(jtpt2/jtpt1, weight);
	}
	else{
	  dphi_subleadB[4]->Fill(acos(cos(jtphi1-jtphi2)), weight);
	  if(acos(cos(jtphi1-jtphi2))>2.0944) X_subleadB[4]->Fill(jtpt2/jtpt1, weight);
	}
      }
    }

    //Logic tree for Double B-tagged Jet
    if(jtpt1>leadjet && jtpt2>subleadjet && discr_prob1 > 0.6 && discr_prob2 > 0.6){
      if(fabs(refparton_flavorForB1)==5){
	if(fabs(refparton_flavorForB2)==5){
	  dphi_doubleB[0]->Fill(acos(cos(jtphi1-jtphi2)), weight);
	  if(acos(cos(jtphi1-jtphi2))>2.0944) X_doubleB[0]->Fill(jtpt2/jtpt1, weight);
	}
	else{
	  dphi_doubleB[1]->Fill(acos(cos(jtphi1-jtphi2)), weight);
	  if(acos(cos(jtphi1-jtphi2))>2.0944) X_doubleB[1]->Fill(jtpt2/jtpt1, weight);
	}
      }
      else if(fabs(refparton_flavorForB1)==4){
	if(fabs(refparton_flavorForB2)==5){
	  dphi_doubleB[1]->Fill(acos(cos(jtphi1-jtphi2)), weight);
	  if(acos(cos(jtphi1-jtphi2))>2.0944) X_doubleB[1]->Fill(jtpt2/jtpt1, weight);
	}
	else if(fabs(refparton_flavorForB2)==4){
	  dphi_doubleB[2]->Fill(acos(cos(jtphi1-jtphi2)), weight);
	  if(acos(cos(jtphi1-jtphi2))>2.0944) X_doubleB[2]->Fill(jtpt2/jtpt1, weight);
	}
	else{
	  dphi_doubleB[3]->Fill(acos(cos(jtphi1-jtphi2)), weight);
	  if(acos(cos(jtphi1-jtphi2))>2.0944) X_doubleB[3]->Fill(jtpt2/jtpt1, weight);
	}
      }
      else{
	if(fabs(refparton_flavorForB1)==4){
	  dphi_doubleB[3]->Fill(acos(cos(jtphi1-jtphi2)), weight);
	  if(acos(cos(jtphi1-jtphi2))>2.0944) X_doubleB[3]->Fill(jtpt2/jtpt1, weight);
	}
	else{
	  dphi_doubleB[4]->Fill(acos(cos(jtphi1-jtphi2)), weight);
	  if(acos(cos(jtphi1-jtphi2))>2.0944) X_doubleB[4]->Fill(jtpt2/jtpt1, weight);
	}
      }
    }
  }

  stackHistogram(X_all,5);
  stackHistogram(dphi_all,5);
  stackHistogram(X_leadB,5);
  stackHistogram(dphi_leadB,5);
  stackHistogram(X_subleadB,5);
  stackHistogram(dphi_subleadB,5);
  stackHistogram(X_doubleB,5);
  stackHistogram(dphi_doubleB,5);

  for(int i=4; i>=0; i--){
    X_all[i]->Scale(1./X_all[0]->Integral());
    dphi_all[i]->Scale(1./dphi_all[0]->Integral());
    X_leadB[i]->Scale(1./X_leadB[0]->Integral());
    dphi_leadB[i]->Scale(1./dphi_leadB[0]->Integral());
    X_subleadB[i]->Scale(1./X_subleadB[0]->Integral());
    dphi_subleadB[i]->Scale(1./dphi_subleadB[0]->Integral());
    X_doubleB[i]->Scale(1./X_doubleB[0]->Integral());
    dphi_doubleB[i]->Scale(1./dphi_doubleB[0]->Integral());
  }

  X_all[5]->Scale(1./X_all[5]->Integral());
  dphi_all[5]->Scale(1./dphi_all[5]->Integral());
  X_leadB[5]->Scale(1./X_leadB[5]->Integral());
  dphi_leadB[5]->Scale(1./dphi_leadB[5]->Integral());
  X_subleadB[5]->Scale(1./X_subleadB[5]->Integral());
  dphi_subleadB[5]->Scale(1./dphi_subleadB[5]->Integral());
  X_doubleB[5]->Scale(1./X_doubleB[5]->Integral());
  dphi_doubleB[5]->Scale(1./dphi_doubleB[5]->Integral());
  
  TLegend *leg = new TLegend(0.2, 0.5, 0.467, 0.9);
  leg->SetFillColor(0);
  leg->AddEntry(X_all[5], "Data");
  leg->AddEntry(X_all[0], "b-bbar");
  leg->AddEntry(X_all[1], "b(bbar)+X");
  leg->AddEntry(X_all[2], "c-cbar");
  leg->AddEntry(X_all[3], "c(cbar)+light");
  leg->AddEntry(X_all[4], "udsg jets");

  TCanvas *c1 = new TCanvas("c1","",1200,600);
  c1->Divide(2,1);
  c1->cd(1);
  X_all[0]->SetTitle("No B-Jet Cut");
  X_all[0]->SetXTitle("p_{T,2}/p_{T,1}");
  X_all[0]->Draw("h");
  for(int i=1; i<6; i++){
    X_all[i]->Draw("same,h");
  }
  leg->Draw();

  c1->cd(2);
  gPad->SetLogy();
  dphi_all[0]->SetXTitle("#Delta#phi_{1,2}");
  dphi_all[0]->Draw("h");
  for(int i=1; i<6; i++){
    dphi_all[i]->Draw("same,h");
  }

  TCanvas *c2 = new TCanvas("c2","",1200,600);
  c2->Divide(2,1);
  c2->cd(1);
  X_leadB[0]->SetTitle("Leading B-Jet");
  X_leadB[0]->SetXTitle("p_{T,2}/p_{T,1}");
  X_leadB[0]->Draw("h");
  for(int i=1; i<6; i++){
    X_leadB[i]->Draw("same,h");
  }
  leg->Draw();

  c2->cd(2);
  gPad->SetLogy();
  dphi_leadB[0]->SetXTitle("#Delta#phi_{1,2}");
  dphi_leadB[0]->Draw("h");
  for(int i=1; i<6; i++){
    dphi_leadB[i]->Draw("same,h");
  }
  
  TCanvas *c3 = new TCanvas("c3","",1200,600);
  c3->Divide(2,1);
  c3->cd(1);
  X_subleadB[0]->SetTitle("Subleading B-Jet");
  X_subleadB[0]->SetXTitle("p_{T,2}/p_{T,1}");
  X_subleadB[0]->Draw("h");
  for(int i=1; i<6; i++){
    X_subleadB[i]->Draw("same,h");
  }
  leg->Draw();

  c3->cd(2);
  gPad->SetLogy();
  dphi_subleadB[0]->SetXTitle("#Delta#phi_{1,2}");
  dphi_subleadB[0]->Draw("h");
  for(int i=1; i<6; i++){
    dphi_subleadB[i]->Draw("same,h");
  }

  TCanvas *c4 = new TCanvas("c4","",1200,600);
  c4->Divide(2,1);
  c4->cd(1);
  X_doubleB[0]->SetTitle("Double B-Jet");
  X_doubleB[0]->SetXTitle("p_{T,2}/p_{T,1}");
  X_doubleB[0]->Draw("h");
  for(int i=1; i<6; i++){
    X_doubleB[i]->Draw("same,h");
  }
  leg->Draw();

  c4->cd(2);
  gPad->SetLogy();
  dphi_doubleB[0]->SetXTitle("#Delta#phi_{1,2}");
  dphi_doubleB[0]->Draw("h");
  for(int i=1; i<6; i++){
    dphi_doubleB[i]->Draw("same,h");
  }

}
