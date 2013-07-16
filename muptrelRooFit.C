#include <iostream>
#include "TMath.h"
#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TROOT.h"

#ifndef __CINT__
#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooDataSet.h>
#include <RooHistPdf.h>
#include <RooPlot.h>
#include <RooAddPdf.h>
#include <RooGlobalFunc.h>
#endif

using namespace RooFit;

RooRealVar muptrelRooFit(bool fixCL=1, char *discr="discr_prob", double minXdiscr=-1, double maxXdiscr=3, char *comment="inclusive sample", char *var="muptrel", double minXvar=0, double maxXvar=5, double ptMin=60, double ptMax=500){

  TFile *fMC = new TFile("histos/ppMC.root");
  TFile *fdata = new TFile("histos/ppdata.root");
  TTree *tMC = (TTree*) fMC->Get("ntmuptrel");
  TTree *tdata = (TTree*) fdata->Get("ntmuptrel");

  TCanvas *c = new TCanvas("c","",600,600);
  
  TH1D *hB = new TH1D("hB","",50,minXvar,maxXvar);
  hB->Sumw2();
  tMC->Draw(Form("%s>>hB",var),Form("weight*(abs(refparton_flavorForB)==5&&jtpt>%f&&jtpt<%f&&%s>%f&&%s<%f)",ptMin,ptMax,discr,minXdiscr,discr,maxXdiscr));
  fixEmpty(hB);
   
  TH1D *hC = new TH1D("hC","",50,minXvar,maxXvar);
  hC->Sumw2();
  tMC->Draw(Form("%s>>hC",var),Form("weight*(abs(refparton_flavorForB)==4&&jtpt>%f&&jtpt<%f&&%s>%f&&%s<%f)",ptMin,ptMax,discr,minXdiscr,discr,maxXdiscr));
  fixEmpty(hC);

  TH1D *hL = new TH1D("hL","",50,minXvar,maxXvar);
  hL->Sumw2();
  tMC->Draw(Form("%s>>hL",var),Form("weight*(abs(refparton_flavorForB)!=5&&abs(refparton_flavorForB)!=4&&abs(refparton_flavorForB)<99&&jtpt>%f&&jtpt<%f&&%s>%f&&%s<%f)",ptMin,ptMax,discr,minXdiscr,discr,maxXdiscr));
  fixEmpty(hL);

  TH1D *hCL = new TH1D("hCL","",50,minXvar,maxXvar);
  hCL->Sumw2();
  tMC->Draw(Form("%s>>hCL",var),Form("weight*(abs(refparton_flavorForB)!=5&&abs(refparton_flavorForB)<99&&jtpt>%f&&jtpt<%f&&%s>%f&&%s<%f)",ptMin,ptMax,discr,minXdiscr,discr,maxXdiscr));
  fixEmpty(hCL);

  delete c;   

  // --- Observable ---
  RooRealVar s(var,var,0,minXvar,maxXvar);
  RooRealVar jtpt("jtpt","jtpt",0,ptMin,ptMax);
  RooRealVar discriminator(discr,discr,0,minXdiscr,maxXdiscr);
 
  // --- Build Histogram PDF ---
  RooDataHist xB("xB","xB",s,hB);
  RooHistPdf bottom("bottom","bottom PDF",s,xB);
  RooDataHist xC("xC","xC",s,hC);
  RooHistPdf charm("charm","charm PDF",s,xC);
  RooDataHist xL("xL","xL",s,hL);
  RooHistPdf light("light","light PDF",s,xL);
  RooDataHist xCL("xCL","xCL",s,hCL);
  RooHistPdf charmlight("charmlight","charmlight PDF",s,xCL);

  /*
  cout<<"hB "<<hB->Integral()<<endl;
  cout<<"hC "<<hC->Integral()<<endl;
  cout<<"hL "<<hL->Integral()<<endl;
  cout<<"hCL "<<hCL->Integral()<<endl;
  //*/

  // --- Construct signal+background PDF ---
  Double_t bInitFrac = hB->Integral()/(hB->Integral()+hCL->Integral());
  Double_t cInitFrac = hC->Integral()/(hB->Integral()+hCL->Integral());
  RooRealVar Bfraction("Bfraction","#light events",bInitFrac,0.,1);
  RooRealVar Cfraction("Cfraction","#background events",cInitFrac,0.,1); 
  if(fixCL) RooAddPdf model("model","",bottom,charmlight,Bfraction);
  else RooAddPdf model("model","",RooArgList(bottom,charm,light),RooArgList(Bfraction,Cfraction));  

  // --- Data sample ---
  RooDataSet *data = new RooDataSet("data","data",tdata,RooArgSet(s,jtpt,discriminator),Form("jtpt>%f&&jtpt<%f&&%s>%f&&%s<%f",ptMin,ptMax,discr,minXdiscr,discr,maxXdiscr));
    
  // --- Prepare canvas ---
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetLabelFont(43,"xyz");
  gStyle->SetLabelSize(14,"xyz");
  gStyle->SetTitleFont(43,"xyz");
  gStyle->SetTitleSize(16,"xyz");
  gStyle->SetTitleOffset(3.5,"x"); 
  gROOT->ForceStyle(1);
  TCanvas *cROOFIT = new TCanvas("cROOFIT",Form("Template fit of %s",var),200,10,1200,600);
  TPad* pad1 = new TPad("pad1","pad1",0,0,0.5,0.9);
  TPad*  pad2 = new TPad("pad2","pad1",0.5,0,1,0.9);
  pad1->Draw();
  pad2->Draw();
  TPaveText *header = new TPaveText(0.05,0.9,0.35,0.99);
  header->AddText(Form("Muon transverse momentum w.r.t. jet axis - %s",comment));
  header->AddText("deltaR < 0.5 ; muon pT > 5 ; jet pT > 60");
  header->SetTextSize(0.027);
  header->SetTextAlign(12);
  header->SetBorderSize(0);
  header->SetFillStyle(0);
  header->Draw();
  TPaveText *header = new TPaveText(0.4,0.9,0.7,0.99);
  header->AddText("ROOFIT ML unbinned fit");
  header->AddText(fixCL?"charm and light normalizations fixed":"3 free normalizations");
  header->SetTextSize(0.027);
  header->SetTextAlign(12);
  header->SetBorderSize(0);
  header->SetFillStyle(0);
  header->Draw();

  // --- Plot before fitting ---
  pad1->cd();
  pad1->SetLogy();
  RooPlot* sframe = s.frame();
  TH2D *htemp = new TH2D("htemp","",100,minXvar,maxXvar,100,0.5,2.5e2) ;
  htemp->SetXTitle(Form("%s %.0f < p_{T} < %.0f GeV/c",var,ptMin,ptMax));
  htemp->SetYTitle("Entries");
  htemp->Draw();
  data->plotOn(sframe,Binning(50));
  if(fixCL) {
    model.plotOn(sframe,Components(charmlight),LineStyle(kDashed),LineColor(30),LineWidth(2));
  } else {
    model.plotOn(sframe,Components(light),LineStyle(kDashed),LineColor(kBlue),LineWidth(2));
    model.plotOn(sframe,Components(charm),LineStyle(kDashed),LineColor(kGreen),LineWidth(2));
    model.plotOn(sframe,Components(RooArgSet(charm,light)),LineStyle(kDashed),LineColor(30),LineWidth(2));
  }
  model.plotOn(sframe,Components(bottom),LineStyle(kDashed),LineColor(kRed),LineWidth(2),FillColor(kRed),FillStyle(1));   
  model.plotOn(sframe,LineWidth(2),LineColor(13));
  data->plotOn(sframe,Binning(50));
  model.paramOn(sframe,Layout(0.4,0.9,0.9),Format("NEU",FixedPrecision(3)));
  sframe->Draw("same");
  TLegend *leg = new TLegend(0.61,fixCL?0.60:0.50,0.98,fixCL?0.78:0.75,comment);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry("h_data","pp @ 2.76 TeV","p");
  leg->AddEntry("model_Norm[muptrel]_Comp[bottom]","b","l");
  if(fixCL) {
    leg->AddEntry("model_Norm[muptrel]_Comp[charmlight]","c + udsg","l");
  } else {
    leg->AddEntry("model_Norm[muptrel]_Comp[charm]","c","l");
    leg->AddEntry("model_Norm[muptrel]_Comp[light]","udsg","l");
    leg->AddEntry("model_Norm[muptrel]_Comp[charm,light]","c + udsg","l");    
  }
  leg->AddEntry("model_Norm[muptrel]","b + c + udsg","l");
  leg->Draw("same");
   
  // --- Perform extended ML fit of composite PDF to data ---
  RooFitResult *fitresult = model.fitTo(*data,Save());
  
  // --- Plot after fitting ---
  pad2->cd();
  pad2->SetLogy();
  sframe = s.frame();
  htemp->Draw();
  data->plotOn(sframe,Binning(50));
  if(fixCL) {
    model.plotOn(sframe,Components(charmlight),LineStyle(kDashed),LineColor(30),LineWidth(2));
  } else {
    model.plotOn(sframe,Components(light),LineStyle(kDashed),LineColor(kBlue),LineWidth(2));
    model.plotOn(sframe,Components(charm),LineStyle(kDashed),LineColor(kGreen),LineWidth(2));
    model.plotOn(sframe,Components(RooArgSet(charm,light)),LineStyle(kDashed),LineColor(30),LineWidth(2));
  }
  model.plotOn(sframe,Components(bottom),LineStyle(kDashed),LineColor(kRed),LineWidth(2),FillColor(kRed),FillStyle(1));   
  model.plotOn(sframe,LineWidth(2),VisualizeError(*fitresult),FillColor(17));
  model.plotOn(sframe,LineWidth(2),LineColor(13));
  data->plotOn(sframe,Binning(50));
  model.paramOn(sframe,Layout(0.4,0.9,0.9),Format("NEU",FixedPrecision(3)));
  sframe->Draw("same");
  leg->Draw("same");

  // --- Save results ---
  cout <<"b jet fraction = "<<Bfraction.getVal()<<endl;
  if(!fixCL) cout <<"c jet fraction = "<<Cfraction.getVal()<<endl;
  if(fixCL) cROOFIT->SaveAs(Form("fit/%s_jtpt(%.0fto%.0f)_%s(%.2fto%.2f)_CLfixed.gif",var,ptMin,ptMax,discr,minXdiscr,maxXdiscr));
  else cROOFIT->SaveAs(Form("fit/%s_jtpt(%.0fto%.0f)_%s(%.2fto%.2f)_CLfree.gif",var,ptMin,ptMax,discr,minXdiscr,maxXdiscr));
  return Bfraction;
}

void fixEmpty(TH1 *h)
{
   for (int i=1;i<=h->GetNbinsX();i++)
   {
      if (h->GetBinContent(i)==0) h->SetBinContent(i,1e-20);
   }
}

/*
void ptDependence()
{
   const int nBins = 8;
   double ptBin[nBins+1] = {60,70,80,90,100,120,140,160,200};
   TH1D *hProb = new TH1D("hProb","",nBins,ptBin);
   TH1D *hCSV = new TH1D("hCSV","",nBins,ptBin);

   for (int n=0; n<nBins;n++)
   {
      RooRealVar f1 = bfractionFit("discr_prob",0,3.5,ptBin[n],ptBin[n+1]);
      RooRealVar f2 = bfractionFit("discr_csvSimple",0,1,ptBin[n],ptBin[n+1]);
      RooRealVar f3 = bfractionFit("discr_csvSimple",0,1,ptBin[n],ptBin[n+1]);
      hProb->SetBinContent(n+1,f1.getVal());    
      hProb->SetBinError(n+1,f1.getError());    
      hCSV->SetBinContent(n+1,f2.getVal());    
      hCSV->SetBinError(n+1,f2.getError());    
   }
   TCanvas *c2 = new TCanvas("c2","",600,600);
   hProb->SetAxisRange(0,0.05,"Y");
   hProb->SetXTitle("Jet p_{T} (GeV/c)");
   hProb->SetYTitle("b-jet fraction");
   hProb->SetTitleOffset(1.5,"Y");
   hProb->Draw();
   hCSV->SetLineColor(2);
   hCSV->SetMarkerColor(2);
   hCSV->SetMarkerStyle(24);
   hCSV->Draw("same");
   
   TLegend *leg = new TLegend(0.2,0.7,0.5,0.9);
   leg->SetBorderSize(0);
   leg->SetFillStyle(0);
   leg->SetFillColor(0);
   leg->AddEntry(hProb,"Jet Probability","pl");
   leg->AddEntry(hCSV,"CSV","pl");
   leg->Draw();
}
*/

