#include <iostream>
#include "TTree.h"
#include "TH1D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"

const double PI=3.14159265;

void setAllColor(TH1D *input, int color){
  input->SetFillColor(color);
  input->SetLineColor(color);
  input->SetMarkerColor(color);
}

void stackHistogram(TH1D **input){

  double intg = 0;
  for(int i=0; i<5; i++){
    for(int j=i+1; j<5; j++){
      input[i]->Add(input[j]);
    }
    if(i==0) intg = input[0]->Integral();
    input[i]->Scale(1./intg);
  }
  setAllColor(input[0], 6);
  setAllColor(input[1], 2);
  setAllColor(input[2], 3);
  setAllColor(input[3], 8);
  setAllColor(input[4], 4);
}


void plotBDijets(){

  TFile *fdata = new TFile("histos/ppdata_ppReco_Dijet_jetTrig.root");
  TFile *fQCD = new TFile("histos/ppMC_ppReco_Dijet_QCDjetTrig.root");
  TFile *fB = new TFile("histos/ppMC_ppReco_Dijet_BjetTrig.root");
  
  // 1st D = ptRatio, dphi;  2nd D = B+B, B+C, B+L, C+L, L+L
  TH1D *QCDjetData[2];
  TH1D *BjetData[2];
  TH1D *QCDjetMC[2][5];
  TH1D *BjetMC[2][5];
  QCDjetData[0] = new TH1D("QCDjetData_0","",10,0,1); QCDjetData[0]->Sumw2();
  QCDjetData[1] = new TH1D("QCDjetData_1","",30,0,PI); QCDjetData[1]->Sumw2();
  BjetData[0] = new TH1D("BjetData_0","",10,0,1); BjetData[0]->Sumw2();
  BjetData[1] = new TH1D("BjetData_1","",30,0,PI); BjetData[1]->Sumw2();
  char* histoname = new char[100];
  for(int i=0; i<5; i++){
    sprintf(histoname,"%s%d","QCDjetMC_0_",i);
    QCDjetMC[0][i] = new TH1D(histoname,"",10,0,1);
    QCDjetMC[0][i]->Sumw2();
    sprintf(histoname,"%s%d","QCDjetMC_1_",i);
    QCDjetMC[1][i] = new TH1D(histoname,"",30,0,PI);
    QCDjetMC[1][i]->Sumw2();
    sprintf(histoname,"%s%d","BjetMC_0_",i);
    BjetMC[0][i] = new TH1D(histoname,"",10,0,1);
    BjetMC[0][i]->Sumw2();
    sprintf(histoname,"%s%d","BjetMC_1_",i);
    BjetMC[1][i] = new TH1D(histoname,"",30,0,PI);
    BjetMC[1][i]->Sumw2();
  }

  TTree *dnt = (TTree*)fdata->Get("nt");
  TTree *QCDnt = (TTree*)fQCD->Get("nt");
  TTree *Bnt = (TTree*)fB->Get("nt");

  //**********************************************************//

  dnt->Draw("(jtpt1-jtpt2)/(jtpt1+jtpt2)>>QCDjetData_0","(jtpt1>80 && jtpt2>40 && pVertexFilterCutGplusUpsPP)");
  dnt->Draw("acos(cos(jtphi1-jtphi2))>>QCDjetData_1","(jtpt1>80 && jtpt2>40 && pVertexFilterCutGplusUpsPP)");
  dnt->Draw("(jtpt1-jtpt2)/(jtpt1+jtpt2)>>BjetData_0","(jtpt1>80 && jtpt2>40 && pVertexFilterCutGplusUpsPP && (discr_ssvHighEff1>2 && discr_ssvHighEff2>2))");
  dnt->Draw("acos(cos(jtphi1-jtphi2))>>BjetData_1","(jtpt1>80 && jtpt2>40 && pVertexFilterCutGplusUpsPP && (discr_ssvHighEff1>2 && discr_ssvHighEff2>2))");

  std::cout << "Finished Data" << std::endl;

  //**********************************************************//
  double jtpt1, jtpt2, jtphi1, jtphi2, refpt1, refpt2, refparton_flavorForB1, refparton_flavorForB2, weight;
  int pVertexFilterCutGplusUpsPP;
  QCDnt->SetBranchAddress("jtpt1",&jtpt1);
  QCDnt->SetBranchAddress("jtpt2",&jtpt2);
  QCDnt->SetBranchAddress("refpt1",&refpt1);
  QCDnt->SetBranchAddress("refpt2",&refpt2);
  QCDnt->SetBranchAddress("jtphi1",&jtphi1);
  QCDnt->SetBranchAddress("jtphi2",&jtphi2);
  QCDnt->SetBranchAddress("pVertexFilterCutGplusUpsPP",&pVertexFilterCutGplusUpsPP);
  QCDnt->SetBranchAddress("refparton_flavorForB1",&refparton_flavorForB1);
  QCDnt->SetBranchAddress("refparton_flavorForB2",&refparton_flavorForB2);
  QCDnt->SetBranchAddress("weight",&weight);

  for(int ie=0; ie<QCDnt->GetEntries(); ie++){
    QCDnt->GetEntry(ie);
    
    if(jtpt1>80 && jtpt2>40 && refpt1>30 && refpt2>25 && pVertexFilterCutGplusUpsPP){
      if(fabs(refparton_flavorForB1)==5 && fabs(refparton_flavorForB2)==5){
	QCDjetMC[0][0]->Fill((jtpt1-jtpt2)/(jtpt1+jtpt2), weight);
	QCDjetMC[1][0]->Fill(acos(cos(jtphi1-jtphi2)), weight);
      }
      else if(fabs(refparton_flavorForB1)==5 || fabs(refparton_flavorForB2)==5){
	QCDjetMC[0][1]->Fill((jtpt1-jtpt2)/(jtpt1+jtpt2), weight);
	QCDjetMC[1][1]->Fill(acos(cos(jtphi1-jtphi2)), weight);
      }
      else if(fabs(refparton_flavorForB1)==4 && fabs(refparton_flavorForB2)==4){
	QCDjetMC[0][2]->Fill((jtpt1-jtpt2)/(jtpt1+jtpt2), weight);
	QCDjetMC[1][2]->Fill(acos(cos(jtphi1-jtphi2)), weight);
      }
      else if(fabs(refparton_flavorForB1)==4 || fabs(refparton_flavorForB2)==4){
	QCDjetMC[0][3]->Fill((jtpt1-jtpt2)/(jtpt1+jtpt2), weight);
	QCDjetMC[1][3]->Fill(acos(cos(jtphi1-jtphi2)), weight);
      }
      else{
	QCDjetMC[0][4]->Fill((jtpt1-jtpt2)/(jtpt1+jtpt2), weight);
	QCDjetMC[1][4]->Fill(acos(cos(jtphi1-jtphi2)), weight);
      }
    }
  }
  /*QCDnt->Draw("(jtpt1-jtpt2)/(jtpt1+jtpt2)>>QCDjetMC_0_0","weight*(jtpt1>80 && jtpt2>40 && refpt1>30 && refpt2>25 && pVertexFilterCutGplusUpsPP && fabs(refparton_flavorForB1)==5 && fabs(refparton_flavorForB2)==5)");
  QCDnt->Draw("(jtpt1-jtpt2)/(jtpt1+jtpt2)>>QCDjetMC_0_1","weight*(jtpt1>80 && jtpt2>40 && refpt1>30 && refpt2>25 && pVertexFilterCutGplusUpsPP && (fabs(refparton_flavorForB1)==5 || fabs(refparton_flavorForB2)==5) && !(fabs(refparton_flavorForB1)==5 && fabs(refparton_flavorForB2)==5))");
  QCDnt->Draw("(jtpt1-jtpt2)/(jtpt1+jtpt2)>>QCDjetMC_0_2","weight*(jtpt1>80 && jtpt2>40 && refpt1>30 && refpt2>25 && pVertexFilterCutGplusUpsPP && fabs(refparton_flavorForB1)==4 && fabs(refparton_flavorForB2)==4)");
  QCDnt->Draw("(jtpt1-jtpt2)/(jtpt1+jtpt2)>>QCDjetMC_0_3","weight*(jtpt1>80 && jtpt2>40 && refpt1>30 && refpt2>25 && pVertexFilterCutGplusUpsPP && (fabs(refparton_flavorForB1)==4 || fabs(refparton_flavorForB2)==4) && !(fabs(refparton_flavorForB1)==5 || fabs(refparton_flavorForB2)==5) && !(fabs(refparton_flavorForB1)==4 && fabs(refparton_flavorForB2)==4))");
  QCDnt->Draw("(jtpt1-jtpt2)/(jtpt1+jtpt2)>>QCDjetMC_0_4","weight*(jtpt1>80 && jtpt2>40 && refpt1>30 && refpt2>25 && pVertexFilterCutGplusUpsPP && !(fabs(refparton_flavorForB1)==4 || fabs(refparton_flavorForB2)==4) && !(fabs(refparton_flavorForB1)==5 || fabs(refparton_flavorForB2)==5))");

  QCDnt->Draw("acos(cos(jtphi1-jtphi2))>>QCDjetMC_1_0","weight*(jtpt1>80 && jtpt2>40 && refpt1>30 && refpt2>25 && pVertexFilterCutGplusUpsPP && fabs(refparton_flavorForB1)==5 && fabs(refparton_flavorForB2)==5)");
  QCDnt->Draw("acos(cos(jtphi1-jtphi2))>>QCDjetMC_1_1","weight*(jtpt1>80 && jtpt2>40 && refpt1>30 && refpt2>25 && pVertexFilterCutGplusUpsPP && (fabs(refparton_flavorForB1)==5 || fabs(refparton_flavorForB2)==5) && !(fabs(refparton_flavorForB1)==5 && fabs(refparton_flavorForB2)==5))");
  QCDnt->Draw("acos(cos(jtphi1-jtphi2))>>QCDjetMC_1_2","weight*(jtpt1>80 && jtpt2>40 && refpt1>30 && refpt2>25 && pVertexFilterCutGplusUpsPP && fabs(refparton_flavorForB1)==4 && fabs(refparton_flavorForB2)==4)");
  QCDnt->Draw("acos(cos(jtphi1-jtphi2))>>QCDjetMC_1_3","weight*(jtpt1>80 && jtpt2>40 && refpt1>30 && refpt2>25 && pVertexFilterCutGplusUpsPP && (fabs(refparton_flavorForB1)==4 || fabs(refparton_flavorForB2)==4) && !(fabs(refparton_flavorForB1)==5 || fabs(refparton_flavorForB2)==5) &&  && !(fabs(refparton_flavorForB1)==4 && fabs(refparton_flavorForB2)==4))");
  QCDnt->Draw("acos(cos(jtphi1-jtphi2))>>QCDjetMC_1_4","weight*(jtpt1>80 && jtpt2>40 && refpt1>30 && refpt2>25 && pVertexFilterCutGplusUpsPP && !(fabs(refparton_flavorForB1)==4 || fabs(refparton_flavorForB2)==4) && !(fabs(refparton_flavorForB1)==5 || fabs(refparton_flavorForB2)==5))");
  */
  std::cout << "Finished QCD MC" << std::endl;

  //**********************************************************//
  
  double refparton_flavorForB3, discr_ssvHighEff1, discr_ssvHighEff2;
  Bnt->SetBranchAddress("jtpt1",&jtpt1);
  Bnt->SetBranchAddress("jtpt2",&jtpt2);
  Bnt->SetBranchAddress("refpt1",&refpt1);
  Bnt->SetBranchAddress("refpt2",&refpt2);
  Bnt->SetBranchAddress("jtphi1",&jtphi1);
  Bnt->SetBranchAddress("jtphi2",&jtphi2);
  Bnt->SetBranchAddress("pVertexFilterCutGplusUpsPP",&pVertexFilterCutGplusUpsPP);
  Bnt->SetBranchAddress("refparton_flavorForB1",&refparton_flavorForB1);
  Bnt->SetBranchAddress("refparton_flavorForB2",&refparton_flavorForB2);
  Bnt->SetBranchAddress("refparton_flavorForB3",&refparton_flavorForB3);
  Bnt->SetBranchAddress("weight",&weight);
  Bnt->SetBranchAddress("discr_ssvHighEff1",&discr_ssvHighEff1);
  Bnt->SetBranchAddress("discr_ssvHighEff2",&discr_ssvHighEff2);

  for(int ie=0; ie<Bnt->GetEntries(); ie++){
    Bnt->GetEntry(ie);

    if(jtpt1>80 && jtpt2>40 && refpt1>30 && refpt2>25 && pVertexFilterCutGplusUpsPP && discr_ssvHighEff1>2 && discr_ssvHighEff2>2){
      if(fabs(refparton_flavorForB1)==5 && fabs(refparton_flavorForB2)==5){
	BjetMC[0][0]->Fill((jtpt1-jtpt2)/(jtpt1+jtpt2), weight);
	BjetMC[1][0]->Fill(acos(cos(jtphi1-jtphi2)), weight);
      }
      else if(fabs(refparton_flavorForB1)==5 || fabs(refparton_flavorForB2)==5){
	BjetMC[0][1]->Fill((jtpt1-jtpt2)/(jtpt1+jtpt2), weight);
	BjetMC[1][1]->Fill(acos(cos(jtphi1-jtphi2)), weight);
      }
      else if(fabs(refparton_flavorForB1)==4 && fabs(refparton_flavorForB2)==4){
	BjetMC[0][2]->Fill((jtpt1-jtpt2)/(jtpt1+jtpt2), weight);
	BjetMC[1][2]->Fill(acos(cos(jtphi1-jtphi2)), weight);
      }
      else if(fabs(refparton_flavorForB1)==4 || fabs(refparton_flavorForB2)==4){
	BjetMC[0][3]->Fill((jtpt1-jtpt2)/(jtpt1+jtpt2), weight);
	BjetMC[1][3]->Fill(acos(cos(jtphi1-jtphi2)), weight);
      }
      else{
	BjetMC[0][4]->Fill((jtpt1-jtpt2)/(jtpt1+jtpt2), weight);
	BjetMC[1][4]->Fill(acos(cos(jtphi1-jtphi2)), weight);
      }
    }
  }
  /*
  Bnt->Draw("(jtpt1-jtpt2)/(jtpt1+jtpt2)>>BjetMC_0_0","weight*(jtpt1>80 && jtpt2>40 && refpt1>30 && refpt2>25 && pVertexFilterCutGplusUpsPP && fabs(refparton_flavorForB1)==5 && fabs(refparton_flavorForB2)==5)");
  Bnt->Draw("(jtpt1-jtpt2)/(jtpt1+jtpt2)>>BjetMC_0_1","weight*(jtpt1>80 && jtpt2>40 && refpt1>30 && refpt2>25 && pVertexFilterCutGplusUpsPP && (fabs(refparton_flavorForB1)==5 || fabs(refparton_flavorForB2)==5) && !(fabs(refparton_flavorForB1)==5 && fabs(refparton_flavorForB2)==5))");
  Bnt->Draw("(jtpt1-jtpt2)/(jtpt1+jtpt2)>>BjetMC_0_2","weight*(jtpt1>80 && jtpt2>40 && refpt1>30 && refpt2>25 && pVertexFilterCutGplusUpsPP && fabs(refparton_flavorForB1)==4 && fabs(refparton_flavorForB2)==4)");
  Bnt->Draw("(jtpt1-jtpt2)/(jtpt1+jtpt2)>>BjetMC_0_3","weight*(jtpt1>80 && jtpt2>40 && refpt1>30 && refpt2>25 && pVertexFilterCutGplusUpsPP && (fabs(refparton_flavorForB1)==4 || fabs(refparton_flavorForB2)==4) && !(fabs(refparton_flavorForB1)==5 || fabs(refparton_flavorForB2)==5) && !(fabs(refparton_flavorForB1)==4 && fabs(refparton_flavorForB2)==4))");
  Bnt->Draw("(jtpt1-jtpt2)/(jtpt1+jtpt2)>>BjetMC_0_4","weight*(jtpt1>80 && jtpt2>40 && refpt1>30 && refpt2>25 && pVertexFilterCutGplusUpsPP && !(fabs(refparton_flavorForB1)==4 || fabs(refparton_flavorForB2)==4) && !(fabs(refparton_flavorForB1)==5 || fabs(refparton_flavorForB2)==5))");

  Bnt->Draw("acos(cos(jtphi1-jtphi2))>>BjetMC_1_0","weight*(jtpt1>80 && jtpt2>40 && refpt1>30 && refpt2>25 && pVertexFilterCutGplusUpsPP && fabs(refparton_flavorForB1)==5 && fabs(refparton_flavorForB2)==5)");
  Bnt->Draw("acos(cos(jtphi1-jtphi2))>>BjetMC_1_1","weight*(jtpt1>80 && jtpt2>40 && refpt1>30 && refpt2>25 && pVertexFilterCutGplusUpsPP && (fabs(refparton_flavorForB1)==5 || fabs(refparton_flavorForB2)==5) && !(fabs(refparton_flavorForB1)==5 && fabs(refparton_flavorForB2)==5))");
  Bnt->Draw("acos(cos(jtphi1-jtphi2))>>BjetMC_1_2","weight*(jtpt1>80 && jtpt2>40 && refpt1>30 && refpt2>25 && pVertexFilterCutGplusUpsPP && fabs(refparton_flavorForB1)==4 && fabs(refparton_flavorForB2)==4)");
  Bnt->Draw("acos(cos(jtphi1-jtphi2))>>BjetMC_1_3","weight*(jtpt1>80 && jtpt2>40 && refpt1>30 && refpt2>25 && pVertexFilterCutGplusUpsPP && (fabs(refparton_flavorForB1)==4 || fabs(refparton_flavorForB2)==4) && !(fabs(refparton_flavorForB1)==5 || fabs(refparton_flavorForB2)==5) && !(fabs(refparton_flavorForB1)==4 && fabs(refparton_flavorForB2)==4))");
  Bnt->Draw("acos(cos(jtphi1-jtphi2))>>BjetMC_1_4","weight*(jtpt1>80 && jtpt2>40 && refpt1>30 && refpt2>25 && pVertexFilterCutGplusUpsPP && !(fabs(refparton_flavorForB1)==4 || fabs(refparton_flavorForB2)==4) && !(fabs(refparton_flavorForB1)==5 || fabs(refparton_flavorForB2)==5))");
  */
  std::cout << "Finished BJet MC" << std::endl;

  stackHistogram(BjetMC[0]);
  stackHistogram(BjetMC[1]);
  stackHistogram(QCDjetMC[0]);
  stackHistogram(QCDjetMC[1]);

  TCanvas *c10 = new TCanvas("c10","",1000,500);
  c10->SetTitle("Inclusive Jets");
  c10->Divide(2,1);
  c10->cd(1);
  QCDjetMC[0][0]->SetXTitle("A_{J}");
  QCDjetMC[0][0]->SetYTitle("evt. norm");
  for(int i=0; i<5; i++){
    if(i==0) QCDjetMC[0][0]->Draw("h");
    else QCDjetMC[0][i]->Draw("same,h");
  }
  QCDjetData[0]->Scale(1./QCDjetData[0]->Integral());
  QCDjetData[0]->Draw("same");
  TLegend *l1 = new TLegend(0.6,0.6,0.91,0.95);
  l1->AddEntry(QCDjetData[0],"pp Data");
  l1->AddEntry(QCDjetMC[0][0], "B+B");
  l1->AddEntry(QCDjetMC[0][1], "B+(C or L)");
  l1->AddEntry(QCDjetMC[0][2], "C+C");
  l1->AddEntry(QCDjetMC[0][3], "C+L");
  l1->AddEntry(QCDjetMC[0][4], "L+L");
  l1->Draw();

  c10->cd(2);
  QCDjetMC[1][0]->SetXTitle("#Delta#phi");
  QCDjetMC[1][0]->SetYTitle("evt. norm");
  for(int i=0; i<5; i++){
    if(i==0) QCDjetMC[1][0]->Draw("h");
    else QCDjetMC[1][i]->Draw("same,h");
  }
  QCDjetData[1]->Scale(1./QCDjetData[1]->Integral());
  QCDjetData[1]->Draw("same");
  
  TCanvas *c11 = new TCanvas("c11","",1000,500);
  c11->SetTitle("B-Tagged Jets");
  c11->Divide(2,1);
  c11->cd(1);
  BjetMC[0][0]->SetXTitle("A_{J}");
  BjetMC[0][0]->SetYTitle("evt. norm");
  for(int i=0; i<5; i++){
    if(i==0) BjetMC[0][0]->Draw("h");
    else BjetMC[0][i]->Draw("same,h");
  }
  BjetData[0]->Scale(1./BjetData[0]->Integral());
  BjetData[0]->Draw("same");
  TLegend *l2 = new TLegend(0.6,0.6,0.91,0.95);
  l2->AddEntry(BjetData[0],"pp Data");
  l2->AddEntry(BjetMC[0][0], "B+B");
  l2->AddEntry(BjetMC[0][1], "B+(C or L)");
  l2->AddEntry(BjetMC[0][2], "C+C");
  l2->AddEntry(BjetMC[0][3], "C+L");
  l2->AddEntry(BjetMC[0][4], "L+L");
  l2->Draw();

  c11->cd(2);
  BjetMC[1][0]->SetXTitle("#Delta#phi");
  BjetMC[1][0]->SetYTitle("evt. norm");
  for(int i=0; i<5; i++){
    if(i==0) BjetMC[1][0]->Draw("h");
    else BjetMC[1][i]->Draw("same,h");
  }
  BjetData[1]->Scale(1./BjetData[1]->Integral());
  BjetData[1]->Draw("same");

}
