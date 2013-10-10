#include <iostream>
#include "TMath.h"
#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TPad.h"
#include "TRandom.h"
#include "THStack.h"
#include "TH2.h"
#include "TLatex.h"
#include "TTree.h"
#include "TPaveText.h"
#include "TFrame.h"

#include <RooGlobalFunc.h>
#include <RooFit.h>
#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooDataSet.h>
#include <RooHistPdf.h>
#include <RooPlot.h>
#include <RooAddPdf.h>
#include <RooCmdArg.h>
#include <RooFitResult.h>

using namespace RooFit;

// Numbers
struct Enumerations {

  Double_t nTaggedJetsMC;
  Double_t nUntaggedJetsMC;
  Double_t nJetsMC;
  Double_t nTaggedBjetsMC;
  Double_t nUntaggedBjetsMC;
  Double_t nBjetsMC;
  Double_t nNonBjetsMC;
  Double_t nTaggedNonBjetsMC;

  Double_t nTaggedJetsData;
  Double_t nUntaggedJetsData;

  Double_t cbForJP;
  Double_t cbForCSV;

  Double_t nTaggedJetsMCError;
  Double_t nUntaggedJetsMCError;
  Double_t nJetsMCError;
  Double_t nTaggedBjetsMCError;
  Double_t nUntaggedBjetsMCError;
  Double_t nBjetsMCError;
  Double_t nNonBjetsMCError;
  Double_t nTaggedNonBjetsMCError;

  Double_t nTaggedJetsDataError;
  Double_t nUntaggedJetsDataError;

  Double_t nJetsData;
  Double_t nJetsDataError;

};

int counter(0);
TCanvas* can1[24];
TH1D* hData[24];
TH1D* hMCC[24];
TH1D* hMCB[24];
TH1D* hMCL[24];
TH1D* hMCLC[24];
TH1D* MCTotal[24];
TH1D* MCSumGreen[24];
TH1D* hMCB_copy[24];
THStack* hs[24];
THStack* ghosths[24];
TH1D* fixCL_MCC[24];
TH1D* fixCL_MCL[24];
TH1D *fluctuateHist(TH1D* h)
{
   TH1D *hToy = (TH1D*)h->Clone(Form("hToy_%s",h->GetName()));
   for (int j=1;j<=hToy->GetNbinsX();j++) {
   	   double value = gRandom->Poisson(h->GetBinContent(j));
   	   hToy->SetBinContent(j,value);
   }
   return hToy;
}


Double_t addError(Double_t aErr, Double_t bErr) {
  // if a,b are independent
  return sqrt(aErr*aErr + bErr*bErr);
}
Double_t substractError(Double_t aErr, Double_t bErr) {
  // if b is a subset of a
  return sqrt(aErr*aErr - bErr*bErr);
}
Double_t prodError(Double_t a, Double_t b, Double_t aErr, Double_t bErr) {
  // if a,b are independent
  return sqrt(a*a*bErr*bErr + b*b*aErr*aErr);
}
Double_t fracError(Double_t a, Double_t b, Double_t aErr, Double_t bErr) {
  // error on a/(a+b) 
  // if a,b are independent
  // and aErr,bErr can be given by sqrt(sum of squares of weights)
  Double_t c=a+b;
  return sqrt( aErr*aErr * b*b/(c*c*c*c)
	      +bErr*bErr * a*a/(c*c*c*c) ) ;
}

Double_t binomError(Double_t a, Double_t b, Double_t aErr, Double_t bErr) {
  // error on a/b, where a+b=const 
  Double_t w = a/b;
  return sqrt(fabs( ( (1.-2.*w)*aErr*aErr + w*w*bErr*bErr )/(b*b) ));
}

void drawText(const std::string text, float xp, float yp);
void format1D(TH1& h1, TCanvas& c1);

class parameters{
   public:
   int fixCL;
   int cbinlo;
   int cbinhi;
   float etalo;
   float etahi;
   float ptMax;
   float ptMin;
};


TFile *fQCDMC, *fBMC, *fCMC, *fdata;
TTree *tQCDMC, *tBMC, *tCMC, *tdata;

// declare counter histos
TH1D *hTaggedLJetsMC,*hTaggedCJetsMC, *hTaggedBJetsMC, *hTaggedJetsMC, *hUntaggedLJetsMC, *hUntaggedCJetsMC, *hUntaggedBJetsMC, *hUntaggedgJetsMC, *hBjetsWithJPinfoMC, *hBjetsWithCSVinfoMC, *hTaggedJetsData, *hUntaggedJetsData, *hUntaggedJetsMC;

RooRealVar *bfractionFit(parameters p, std::string var, std::string discr, double minXdiscr, double maxXdiscr, std::string comment, double maxYaxis, bool toyMC, bool useDataBGtemplate, bool verbose);

void fillCounterHistos(std::string discr, double workingPoint, int cbinlo, int cbinhi, float etalo, float etahi);

Enumerations count(double ptMin, double ptMax);

const int nBins = 7;
double ptBin[nBins+1] = {60,70,80,90,110,130,170,250};
//const int nBins = 4;
//double ptBin[nBins+1] = {80,100,120,150,200};
//const int nBins = 1;
//double ptBin[nBins+1] = {60,70};


///XXXX
void bfractionVsJetPtPPgPlus(std::string tagger="discr_ssvHighEff", double workingPoint=2.0, int fixCL=0, std::string taggerName="SSVHE", int cbinlo=0, int cbinhi=40, float etalo=0., float etahi=2., bool useDataBGtemplate = false, bool verbose = true) {

  if(useDataBGtemplate) fixCL=1;

  parameters p;
  p.fixCL = fixCL;
  p.cbinlo = cbinlo;
  p.cbinhi = cbinhi;
  p.etalo = etalo;
  p.etahi = etahi;
  
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  gStyle->SetTextFont(42); 
  gStyle->SetLabelFont(43,"xyz");
  gStyle->SetLabelSize(28,"xyz");
  gStyle->SetTitleFont(43,"xyz");
  gStyle->SetTitleSize(28,"xyz");
  gStyle->SetTitleOffset(1.0,"x"); 
  gStyle->SetNdivisions(510,"xy");

  gROOT->ForceStyle(1);

  /*
  gStyle->SetLabelFont(43,"xyz");
  gStyle->SetLabelSize(14,"xyz");
  gStyle->SetTitleFont(43,"xyz");
  gStyle->SetTitleSize(16,"xyz");
  gStyle->SetTitleOffset(1.5,"xy"); 
  */

  fQCDMC = new TFile("histos/ppMC_ppReco_ak3PF_QCDjetTrig_noIPupperCut_9-15.root"); 
  fBMC = new TFile("histos/ppMC_ppReco_ak3PF_BjetTrig_noIPupperCut_9-15.root"); 
  fCMC = new TFile("histos/ppMC_ppReco_ak3PF_CjetTrig_noIPupperCut_9-15.root"); 
  fdata = new TFile("histos/ppdata_ppReco_ak3PF_DinkoXchk_jetTrig_noIPupperCut.root");

  tQCDMC = (TTree*) fQCDMC->Get("nt");
  tBMC = (TTree*) fBMC->Get("nt");
  tCMC = (TTree*) fCMC->Get("nt");
  tdata = (TTree*) fdata->Get("nt");

  fillCounterHistos(tagger.c_str(), workingPoint, cbinlo, cbinhi, etalo,  etahi);
  
  int doLTJP=1;
  int doLTCSV=0;
  
  Double_t bPurMC, bPurData, bEffMC, bEffDataLTJP, bEffDataLTCSV=0, taggedFracData, bFracMC, bFracData, bFracDataLTJP, bFracDataLTCSV=0, bFracJPdirect;
  Double_t bPurMCError, bPurDataError, bEffMCError, bEffDataLTJPError, bEffDataLTCSVError, taggedFracDataError, bFracMCError, bFracDataError, bFracDataLTJPError, bFracDataLTCSVError, bFracJPdirectError;
  Enumerations numbers;
  
  TH1D *hBPurityData = new TH1D("hBPurityData","hBPurityData;Jet p_{T} (GeV/c);b-Tagging purity",nBins,ptBin);
  TH1D *hBPurityMC = new TH1D("hBPurityMC","hBPurityMC;Jet p_{T} (GeV/c);b-Tagging purity",nBins,ptBin);
  TH1D *hRawBData = new TH1D("hRawBData","hRawBData;Jet p_{T} (GeV/c);raw b-jets",nBins,ptBin);
  TH1D *hRawBMC = new TH1D("hRawBMC","hRawBMC;Jet p_{T} (GeV/c);raw b-jets",nBins,ptBin);

  TH1D *hBEfficiencyMC = new TH1D("hBEfficiencyMC","hBEfficiencyMC;Jet p_{T} (GeV/c);b-Tagging efficiency",nBins,ptBin);
  TH1D *hBEfficiencyDataLTJP = new TH1D("hBEfficiencyDataLTJP","hBEfficiencyDataLTJP;Jet p_{T} (GeV/c);b-Tagging efficiency",nBins,ptBin);
  TH1D *hBEfficiencyDataLTCSV = new TH1D("hBEfficiencyDataLTCSV","hBEfficiencyDataLTCSV;Jet p_{T} (GeV/c);b-Tagging efficiency",nBins,ptBin);
  
  TH1D *hBFractionMC = new TH1D("hBFractionMC","hBFractionMC;Jet p_{T} (GeV/c);b-jet fraction",nBins,ptBin);
  TH1D *hBFractionData = new TH1D("hBFractionData","hBFractionData;Jet p_{T} (GeV/c);b-jet fraction",nBins,ptBin);
  TH1D *hBFractionDataLTJP = new TH1D("hBFractionDataLTJP","hBFractionDataLTJP;Jet p_{T} (GeV/c);b-jet fraction",nBins,ptBin);
  TH1D *hBFractionDataLTCSV = new TH1D("hBFractionDataLTCSV","hBFractionDataLTCSV;Jet p_{T} (GeV/c);b-jet fraction",nBins,ptBin);
  TH1D *hBFractionJPdirect = new TH1D("hBFractionJPdirect","hBFractionJPdirect;Jet p_{T} (GeV/c);b-jet fraction",nBins,ptBin);

  TH1D *hIncJetsMC = new TH1D("hIncJetsMC","hIncJetsMC;Jet p_{T} (GeV/c); inclusive jets in MC",nBins,ptBin);
  TH1D *hIncJetsData = new TH1D("hIncJetsData","hIncJetsData;Jet p_{T} (GeV/c); inclusive jets in Data",nBins,ptBin);
 
  int ncol=nBins;
  int nrow=1;

  if(nBins==3||nBins==2){
    ncol=nBins;
  }

  if(nBins==4){
    ncol=nBins/2;
    nrow=nBins/2;
  }
  if(nBins==6){
    ncol=3;
    nrow=2;
  }

  TCanvas *c1=new TCanvas("c1","c1",1200,600);
  c1->Divide(ncol,nrow);
  
  //TCanvas *c2=new TCanvas("c2","c2",1200,600);
  //c2->Divide(ncol,nrow);
  
  TCanvas *c3=new TCanvas("c3","c3",1200,600);
  c3->Divide(ncol,nrow);

  TCanvas *c4=new TCanvas("c4","c4",1200,600);
  c4->Divide(ncol,nrow);

  TCanvas *cCount = new TCanvas("cCount","cCount",600,600);
  bool doToyCalc = false;
  for (int n=0;n<nBins;n++) {

    cout<<"Processing jet pT bin ["<<ptBin[n]<<","<<ptBin[n+1]<<"] ..."<<endl;
    cCount->cd();
    numbers = count(ptBin[n],ptBin[n+1]);
    c1->cd(n+1);
    p.ptMin = ptBin[n];
    p.ptMax = ptBin[n+1];
    RooRealVar *fitSvtxmTag = bfractionFit(p,"svtxm",tagger,workingPoint,6,"b-tagged sample (SSVHE > 2)",9e3,doToyCalc, useDataBGtemplate, verbose);

    //c2->cd(n+1);
    //c2->GetPad(n+1)->SetLogy();
    //RooRealVar *fitJpDirect = bfractionFit(p,"discr_prob","discr_prob",0.,3.,"inclusive sample",4e5,doToyCalc);

    RooRealVar *fitJpTag=NULL;
    RooRealVar *fitJpBeforetag=NULL;
    RooRealVar *fitCsvTag=NULL;
    RooRealVar *fitCsvBeforetag=NULL;
    if (doLTJP) {
      c3->cd(n+1);
      //c3->GetPad(n+1)->SetLogy();
      fitJpBeforetag = bfractionFit(p,"discr_prob","discr_prob",0,3.,"jets with JP info",4e5,doToyCalc,false,verbose);
      c4->cd(n+1);
      //c4->GetPad(n+1)->SetLogy();
      fitJpTag = bfractionFit(p,"discr_prob",tagger,workingPoint,6,"b-tagged sample (SSVHE > 2)",4e5,doToyCalc,false,verbose);
    } 
    if (doLTCSV) {
      fitCsvBeforetag = bfractionFit(p,"discr_csvSimple",tagger,-2,10,"jets with CSV info",4e5,doToyCalc,false,verbose);
      fitCsvTag = bfractionFit(p,"discr_csvSimple",tagger,workingPoint,10,Form("b-tagged sample (%s > %.1f)",taggerName.c_str(),workingPoint),4e5,doToyCalc,false,verbose);
    } 

    taggedFracData = numbers.nTaggedJetsData / numbers.nJetsData;
    taggedFracDataError = binomError(numbers.nTaggedJetsData,numbers.nJetsData,numbers.nTaggedJetsDataError,numbers.nJetsDataError);
    
    //*  --- b-tagging purity --- 

    bPurMC = numbers.nTaggedBjetsMC / numbers.nTaggedJetsMC;
    cout<<" bPurMC "<<bPurMC<<" numbers.nTaggedBjetsMC "<<numbers.nTaggedBjetsMC<<" numbers.nTaggedJetsMC "<<numbers.nTaggedJetsMC<<endl;
    bPurMCError = binomError(numbers.nTaggedBjetsMC,numbers.nTaggedJetsMC,numbers.nTaggedBjetsMCError,numbers.nTaggedJetsMCError);
    bPurData = fitSvtxmTag->getVal();
    bPurDataError = fitSvtxmTag->getError();

    hBPurityMC->SetBinContent(n+1,bPurMC); 
    hBPurityMC->SetBinError(n+1,bPurMCError); 
    hBPurityData->SetBinContent(n+1,bPurData);    
    hBPurityData->SetBinError(n+1,bPurDataError); 

    hRawBData->SetBinContent(n+1,bPurData*numbers.nTaggedJetsData);
    hRawBData->SetBinError(n+1,bPurDataError*numbers.nTaggedJetsData);
    hRawBMC->SetBinContent(n+1,bPurMC*numbers.nTaggedJetsData);
    hRawBMC->SetBinError(n+1,bPurMCError*numbers.nTaggedJetsData);
    hIncJetsData->SetBinContent(n+1,numbers.nJetsData);
    hIncJetsData->SetBinError(n+1,numbers.nJetsDataError);
    hIncJetsMC->SetBinContent(n+1,numbers.nJetsMC);
    hIncJetsMC->SetBinError(n+1,numbers.nJetsMCError);
     

    //*/
    
    //*  --- b-tagging efficiency --- 

    bEffMC = numbers.nTaggedBjetsMC / numbers.nBjetsMC;
    bEffMCError = binomError(numbers.nTaggedBjetsMC,numbers.nBjetsMC,numbers.nTaggedBjetsMCError,numbers.nBjetsMCError);
    hBEfficiencyMC->SetBinContent(n+1,bEffMC); 
    hBEfficiencyMC->SetBinError(n+1,bEffMCError);

    if (doLTJP) {
      bEffDataLTJP = taggedFracData * numbers.cbForJP * fitJpTag->getVal() / fitJpBeforetag->getVal();
      bEffDataLTJPError = prodError(taggedFracData,fitJpTag->getVal(),taggedFracDataError,fitJpTag->getError()) * numbers.cbForJP / fitJpBeforetag->getVal(); 
      hBEfficiencyDataLTJP->SetBinContent(n+1,bEffDataLTJP);    
      hBEfficiencyDataLTJP->SetBinError(n+1,bEffDataLTJPError);
    } 

    if (doLTCSV) {
      bEffDataLTCSV = taggedFracData * numbers.cbForCSV * fitCsvTag->getVal() / fitCsvBeforetag->getVal();
      bEffDataLTCSVError = prodError(taggedFracData,fitCsvTag->getVal(),taggedFracDataError,fitCsvTag->getError()) * numbers.cbForCSV / fitCsvBeforetag->getVal(); 
      hBEfficiencyDataLTCSV->SetBinContent(n+1,bEffDataLTCSV);    
      hBEfficiencyDataLTCSV->SetBinError(n+1,bEffDataLTCSVError); 
    } 
    
    //*  --- b fraction --- 

    bFracMC = numbers.nBjetsMC / numbers.nJetsMC;
    bFracMC = numbers.nTaggedJetsMC * bPurMC / (bEffMC * numbers.nJetsMC); // for check : same as previous
    bFracMCError = binomError(numbers.nBjetsMC,numbers.nJetsMC,numbers.nBjetsMCError,numbers.nJetsMCError); 
    hBFractionMC->SetBinContent(n+1,bFracMC); 
    hBFractionMC->SetBinError(n+1,bFracMCError); 


    bFracData = taggedFracData * bPurData / bEffMC; // efficiency from MC
    bFracDataError = prodError(taggedFracData,bPurData,taggedFracDataError,bPurDataError) / bEffMC; // stat.error from purity and tagged-fraction (assumed independent)
    bFracDataError = bFracData * bPurDataError / bPurData; // stat.error only from purity
    hBFractionData->SetBinContent(n+1,bFracData);    
    hBFractionData->SetBinError(n+1,bFracDataError);

    if (doLTJP) {
      bFracDataLTJP = taggedFracData * bPurData / bEffDataLTJP ; // efficiency from LTJP method
      bFracDataLTJPError = prodError(taggedFracData,bPurData,taggedFracDataError,bPurDataError) / bEffDataLTJP; // stat.error from purity and tagged-fraction (assumed independent)
      bFracDataLTJPError = bFracDataLTJP * bPurDataError / bPurData; // stat.error only from purity
      hBFractionDataLTJP->SetBinContent(n+1,bFracDataLTJP);    
      hBFractionDataLTJP->SetBinError(n+1,bFracDataLTJPError);
    } 

    if (doLTCSV) {
      bFracDataLTCSV = taggedFracData * bPurData / bEffDataLTCSV; // efficiency from LTCSV method
      bFracDataLTCSVError = prodError(taggedFracData,bPurData,taggedFracDataError,bPurDataError) / bEffDataLTCSV; // stat.error from purity and tagged-fraction (assumed independent)
      bFracDataLTCSVError = bFracDataLTCSV * bPurDataError / bPurData; // stat.error only from purity
      hBFractionDataLTCSV->SetBinContent(n+1,bFracDataLTCSV);    
      hBFractionDataLTCSV->SetBinError(n+1,bFracDataLTCSVError);
    } 

    bFracJPdirect = fitJpBeforetag->getVal();
    bFracJPdirectError = fitJpBeforetag->getError();
    hBFractionJPdirect->SetBinContent(n+1,bFracJPdirect);   
    hBFractionJPdirect->SetBinError(n+1,bFracJPdirectError);
    //*/

    //*
    if(verbose){
      cout<<"nTaggedJetsMC "<<numbers.nTaggedJetsMC<<endl;
      cout<<"nUntaggedJetsMC "<<numbers.nUntaggedJetsMC<<endl;
      cout<<"nJetsMC "<<numbers.nJetsMC<<endl;
      cout<<"nTaggedBjetsMC "<<numbers.nTaggedBjetsMC<<endl;
      cout<<"nUntaggedBjetsMC "<<numbers.nUntaggedBjetsMC<<endl;
      cout<<"nBjetsMC "<<numbers.nBjetsMC<<endl;
      cout<<"nNonBjetsMC "<<numbers.nNonBjetsMC<<endl;
      cout<<"nTaggedNonBjetsMC "<<numbers.nTaggedNonBjetsMC<<endl;
      cout<<"nTaggedJetsData "<<numbers.nTaggedJetsData<<endl;
      cout<<"nUntaggedJetsData "<<numbers.nUntaggedJetsData<<endl;
      cout<<"bPurMC "<<bPurMC<<endl;
      cout<<"bPurData "<<bPurData<<endl;
      cout<<"bEffMC "<<bEffMC<<endl;
      cout<<"CbForJP "<<numbers.cbForJP<<endl;
      cout<<"bEffDataLTJP "<<bEffDataLTJP<<endl;
      cout<<"CbForCSV "<<numbers.cbForCSV<<endl;
      cout<<"bEffDataLTCSV "<<bEffDataLTCSV<<endl;
      cout<<"bFracMC "<<bFracMC<<endl;
      cout<<"bFracData "<<bFracData<<endl;
      cout<<"bFracDataLTJP "<<bFracDataLTJP<<endl;
      cout<<"bFracDataLTCSV "<<bFracDataLTCSV<<endl;
      cout<<"bFracJPdirect "<<bFracJPdirect<<endl;
      cout<<endl;
    }
    //*/
  }
  
  TLegend *legPur = new TLegend(0.2,0.6,0.9,0.90,Form("Purity of b-tagged sample (%s > %.1f)",taggerName.c_str(),workingPoint));
  legPur->SetBorderSize(0);
  //legPur->SetFillColor(kGray);
  legPur->SetFillStyle(0);
  legPur->AddEntry(hBPurityMC,"MC Input","pl");
  legPur->AddEntry(hBPurityData,"Data","pl");
  //legPur->SetTextSize(0.030);
  TCanvas *cBPurity = new TCanvas("cBPurity","b purity",700,600);
  cBPurity->cd();
  hBPurityMC->SetAxisRange(0,1,"Y");
  hBPurityMC->SetTitleOffset(1.3,"Y");
  hBPurityMC->SetLineColor(2);
  hBPurityMC->SetMarkerColor(2);
  hBPurityMC->SetMarkerStyle(21);
  hBPurityMC->SetMarkerSize(1.41);
  format1D(*hBPurityMC,*cBPurity);
  hBPurityMC->Draw();
  drawText("CMS Preliminary",0.15,0.965);
  drawText("#sqrt{s_{NN}} = 2.76 TeV",0.55,0.965);
  hBPurityData->SetLineColor(1);
  hBPurityData->SetMarkerColor(1);
  hBPurityData->SetMarkerStyle(20);
  hBPurityData->SetMarkerSize(1.41);
  hBPurityData->Draw("same");   
  legPur->Draw();
  //cBPurity->SaveAs("ssvhePurpp.pdf");

  TLegend *legEff = new TLegend(0.2,0.6,0.9,0.90,Form("Efficiency for tagging b-jets (%s > %.1f)",taggerName.c_str(),workingPoint));
  legEff->SetBorderSize(0);
  //legEff->SetFillColor(kGray);
  legEff->SetFillStyle(0);
  legEff->AddEntry(hBEfficiencyMC,"Simulation","pl");
  //legEff->SetTextSize(0.030);
  if (doLTJP) legEff->AddEntry(hBEfficiencyDataLTJP,"Reference Tagger","pl");
  if (doLTCSV) legEff->AddEntry(hBEfficiencyDataLTCSV,"LT method (CSV)","pl");
  TCanvas *cBEfficiency = new TCanvas("cBEfficiency","b-Tagging efficiency",700,600);
  cBEfficiency->cd();
  hBEfficiencyMC->SetAxisRange(0,1,"Y");
  hBEfficiencyMC->SetTitleOffset(1.3,"Y");
  hBEfficiencyMC->SetLineColor(2);
  hBEfficiencyMC->SetMarkerColor(2);
  hBEfficiencyMC->SetMarkerStyle(21);
  hBEfficiencyMC->SetMarkerSize(1.41);
  format1D(*hBEfficiencyMC,*cBEfficiency);
  hBEfficiencyMC->Draw();
  drawText("CMS Preliminary",0.15,0.965);
  drawText("#sqrt{s_{NN}} = 2.76 TeV",0.55,0.965);
  if (doLTJP) {
    hBEfficiencyDataLTJP->SetLineColor(kGreen+2);
    hBEfficiencyDataLTJP->SetMarkerColor(kGreen+2);
    hBEfficiencyDataLTJP->SetMarkerStyle(20);
    hBEfficiencyDataLTJP->SetMarkerSize(1.41);
    hBEfficiencyDataLTJP->Draw("same");
  }
  if (doLTCSV) {
    hBEfficiencyDataLTCSV->SetLineColor(7);
    hBEfficiencyDataLTCSV->SetMarkerColor(7);
    hBEfficiencyDataLTCSV->SetMarkerStyle(20);
    hBEfficiencyDataLTCSV->SetMarkerSize(1.41);
    hBEfficiencyDataLTCSV->Draw("same");
  }
  legEff->Draw();
  //cBEfficiency->SaveAs("ssvheEffpp.pdf");


  TLegend *legFrac = new TLegend(0.15,0.15,0.85,0.35);
  legFrac->SetBorderSize(0);
  //legFrac->SetFillColor(kGray);
  legFrac->SetFillStyle(0);
  legFrac->AddEntry(hBFractionMC,"MC Input","pl");
  //legFrac->SetTextSize(0.022);
  legFrac->AddEntry(hBFractionData,Form("%s > %.1f + pur. from SV mass + eff. from MC",taggerName.c_str(),workingPoint),"pl");
  if (doLTJP) legFrac->AddEntry(hBFractionDataLTJP,Form("%s > %.1f + pur. from SV mass + eff. from LT (JP)",taggerName.c_str(),workingPoint),"pl");
  if (doLTCSV) legFrac->AddEntry(hBFractionDataLTCSV,Form("%s > %.1f + pur. from SV mass + eff. from LT (CSV)",taggerName.c_str(),workingPoint),"pl");
  legFrac->AddEntry(hBFractionJPdirect,"Direct fit to JP","pl");
  TCanvas *cBFraction = new TCanvas("cBFraction","b-jet fraction",700,600);
  cBFraction->cd();
  hBFractionMC->SetAxisRange(0,0.03,"Y");
  hBFractionMC->SetTitleOffset(1.45,"Y");
  hBFractionMC->SetLineColor(2);
  hBFractionMC->SetMarkerColor(2);
  hBFractionMC->SetMarkerStyle(21);
  hBFractionMC->SetMarkerSize(1.41);
  format1D(*hBFractionMC,*cBFraction);
  hBFractionMC->Draw();
  drawText("CMS Preliminary",0.15,0.965);
  drawText("#sqrt{s_{NN}} = 2.76 TeV",0.55,0.965); 
  hBFractionData->SetLineColor(1);
  hBFractionData->SetMarkerColor(1);
  hBFractionData->SetMarkerStyle(20);
  hBFractionData->SetMarkerSize(1.41);
  hBFractionData->Draw("same");   
  if (doLTJP) {
    hBFractionDataLTJP->SetLineColor(kGreen+2);
    hBFractionDataLTJP->SetMarkerColor(kGreen+2);
    hBFractionDataLTJP->SetMarkerStyle(20);
    hBFractionDataLTJP->SetMarkerSize(1.41);
    hBFractionDataLTJP->Draw("same");
  }
  if (doLTCSV) {
    hBFractionDataLTCSV->SetLineColor(7);
    hBFractionDataLTCSV->SetMarkerColor(7);
    hBFractionDataLTCSV->SetMarkerStyle(20);
    hBFractionDataLTCSV->SetMarkerSize(1.41);
    hBFractionDataLTCSV->Draw("same");
  }
  hBFractionJPdirect->SetLineColor(4);
  hBFractionJPdirect->SetMarkerColor(4);
  hBFractionJPdirect->SetMarkerStyle(20);
  hBFractionJPdirect->SetMarkerSize(1.41);
  hBFractionJPdirect->Draw("same");
  legFrac->Draw();
  //cBFraction->SaveAs("ssvheFracpp.pdf");

  string mcDataString = "MC";
  if(useDataBGtemplate) mcDataString = "Data";

  TFile *fout = new TFile(Form("output/NewFormatV5_bFraction%sTemplate_pppp1_jetcut30_dinkoXchk_gPlusFilter_%sat%.1fFixCL%d_bin_%d_%d_eta_%d_%d.root",mcDataString.c_str(),taggerName.c_str(),workingPoint,fixCL,cbinlo,cbinhi,(int)etalo,(int)etahi),"recreate");

  hBFractionMC->Write();
  hBFractionData->Write();
  if (doLTJP) hBFractionDataLTJP->Write();
  if (doLTCSV) hBFractionDataLTCSV->Write();
  hBFractionJPdirect->Write();
  hBPurityMC->Write();
  hBPurityData->Write();
  hRawBMC->Write();
  hRawBData->Write();
  hBEfficiencyMC->Write();
  hBEfficiencyDataLTJP->Write();
  hIncJetsMC->Write();
  hIncJetsData->Write();
  fout->Close();

  c1->SaveAs(Form("gifs/svtxMassFit_%s.gif",fixCL?"CLfixed":"CLfree"));
  //c2->SaveAs(Form("gifs/jpDirectFit_%s.gif",fixCL?"CLfixed":"CLfree"));
  c3->SaveAs(Form("gifs/jpBeforeTag_%s.gif",fixCL?"CLfixed":"CLfree"));
  c4->SaveAs(Form("gifs/jpAfterTag_%s.gif",fixCL?"CLfixed":"CLfree"));

}





void fixEmpty(TH1 *h){
   for (int i=1;i<=h->GetNbinsX();i++){
      if (h->GetBinContent(i)==0) h->SetBinContent(i,1e-20);
   }
}

RooRealVar *bfractionFit(parameters p, std::string var, std::string discr, double minXdiscr, double maxXdiscr, std::string comment, double maxYaxis, bool toyMC, bool useDataBGtemplate, bool verbose)
{
  bool fixCL = p.fixCL;
  int cbinlo = p.cbinlo;
  int cbinhi = p.cbinhi;
  float etalo = p.etalo;
  float etahi = p.etahi;
  float ptMax = p.ptMax;
  float ptMin = p.ptMin;
  double minXvar=0.;
  double maxXvar=6.;

  if (var == "discr_prob") {
     maxXvar = 3;
  } else   if (var == "discr_probb") {
     maxXvar = 1;
  }

  // discr_prob : from (0) 0 to 3, operating point : 0.6 (1%), 0.7 
  // discr_ssvHighEff : from (-1) 1 to 6, operating point : 2 ?
  // discr_ssvHighPur : from (-1) 1 to 6, operating point : 2 ?  
  // discr_csvSimple : from (-10,-1) 0 to 1, operating point  : 0.9  
  // svtxm : from (0) 0 to 7 
  // muptrel : from (0) 0 to 5
  
  int nhistBins=30;
  if(var=="svtxm") nhistBins=24;

  double ptHatMin = 0.;
  if(ptMin>=80.&&ptMin<120.) ptHatMin = 50.;
  else if(ptMin>=120.&&ptMin<150.) ptHatMin = 65.;
  else if(ptMin>=150.&&ptMin<200.) ptHatMin = 80.;

  TH1D *hB = new TH1D("hB","hB",nhistBins,minXvar,maxXvar);
  hB->Sumw2();
  tBMC->Draw(Form("%s>>hB",var.c_str()),Form("weight*(abs(refparton_flavorForB)==5&&pVertexFilterCutGplusUpsPP&&jtpt>=%f&&jtpt<%f&&%s>=%f&&%s<%f&&bin>=%d&&bin<%d&&abs(jteta)>%f&&abs(jteta)<%f&&pthat>%f)",ptMin,ptMax,discr.c_str(),minXdiscr,discr.c_str(),maxXdiscr,cbinlo,cbinhi,etalo,etahi,ptHatMin));
  fixEmpty(hB);
  
  TH1D *hC = new TH1D("hC","hC",nhistBins,minXvar,maxXvar);
  hC->Sumw2();
  tCMC->Draw(Form("%s>>hC",var.c_str()),Form("weight*(abs(refparton_flavorForB)==4&&pVertexFilterCutGplusUpsPP&&jtpt>=%f&&jtpt<%f&&%s>=%f&&%s<%f&&bin>=%d&&bin<%d&&abs(jteta)>%f&&abs(jteta)<%f&&pthat>%f)",ptMin,ptMax,discr.c_str(),minXdiscr,discr.c_str(),maxXdiscr,cbinlo,cbinhi,etalo,etahi,ptHatMin));
  fixEmpty(hC);

  TH1D *hL = new TH1D("hL","hL",nhistBins,minXvar,maxXvar);
  hL->Sumw2();
  tQCDMC->Draw(Form("%s>>hL",var.c_str()),Form("weight*(abs(refparton_flavorForB)!=5&&pVertexFilterCutGplusUpsPP&&abs(refparton_flavorForB)!=4&&abs(refparton_flavorForB)<99&&jtpt>=%f&&jtpt<%f&&%s>=%f&&%s<%f&&bin>=%d&&bin<%d&&abs(jteta)>%f&&abs(jteta)<%f&&pthat>%f)",ptMin,ptMax,discr.c_str(),minXdiscr,discr.c_str(),maxXdiscr,cbinlo,cbinhi,etalo,etahi,ptHatMin));
  fixEmpty(hL);

  TH1D *hCL = NULL;
  if(useDataBGtemplate){
    hCL = new TH1D("hCL","hCL",nhistBins,minXvar,maxXvar);
    hCL->Sumw2();
    tdata->Draw(Form("%s>>hCL",var.c_str()),Form("weight*(discr_prob<0.5&&pVertexFilterCutGplusUpsPP&&jtpt>=%f&&jtpt<%f&&%s>=%f&&%s<%f&&bin>=%d&&bin<%d&&abs(jteta)>=%f&&abs(jteta)<%f)",ptMin,ptMax,discr.c_str(),minXdiscr,discr.c_str(),maxXdiscr,cbinlo,cbinhi,etalo,etahi));
  }
  else{
    hCL = (TH1D*) hL->Clone();
    hCL->Add(hC);
  }  
  fixEmpty(hCL);

  //*/


  // --- Observable ---
  RooRealVar s(var.c_str(),var.c_str(),0,minXvar,maxXvar);
  RooRealVar jtpt("jtpt","jtpt",0,ptMin,ptMax);
  RooRealVar discriminator(discr.c_str(),discr.c_str(),0,minXdiscr,maxXdiscr);
  RooRealVar bin("bin","bin",0,0,40); 
  RooRealVar jteta("jteta","jteta",0,-2,2); 
  RooRealVar weight("weight","weight",0,0,1e6);
  RooRealVar pVertexFilterCutGplusUpsPP("pVertexFilterCutGplusUpsPP","pVertexFilterCutGplusUpsPP",0,0,1);

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

    // fucked
  // --- Construct signal+background PDF ---
  //Double_t bInitFrac = hB->Integral()/(hB->Integral()+hCL->Integral());
  //Double_t cInitFrac = hC->Integral()/(hB->Integral()+hCL->Integral());

  double fracGuess=0.3;
  if(var != "svmtx") fracGuess =0.01;
  /*
  double bGuess = 0.03;
  double cGuess = 0.1;
  if(var != "svmtx"){
    bGuess =0.01;
    cGuess =0.05;
  }
  */
  RooRealVar *Bfraction = new RooRealVar("Bfraction","#signal events",fracGuess,0.,1);
  RooRealVar *Cfraction = new RooRealVar("Cfraction","#background events",fracGuess,0.,1); 
  RooAddPdf *model;
  if(fixCL) model = new RooAddPdf("model","",bottom,charmlight,*Bfraction);
  else model=new RooAddPdf("model","",RooArgList(bottom,charm,light),RooArgList(*Bfraction,*Cfraction));  
  //RooAddPdf model=RooAddPdf("model","",RooArgList(bottom,charm,light),RooArgList(*Bfraction,*Cfraction));  

  // --- Data sample ---
  RooDataSet *data = NULL;
  if(var == discr) data = new  RooDataSet("data","data",tdata,RooArgSet(s,jtpt,jteta,bin,weight,pVertexFilterCutGplusUpsPP),Form("pVertexFilterCutGplusUpsPP&&jtpt>=%f&&jtpt<%f&&%s>=%f&&%s<%f&&abs(jteta)>%f&&abs(jteta)<%f&&bin>=%d&&bin<%d",ptMin,ptMax,discr.c_str(),minXdiscr,discr.c_str(),maxXdiscr,etalo,etahi,cbinlo,cbinhi),"weight");
  else data = new  RooDataSet("data","data",tdata,RooArgSet(s,jtpt,jteta,bin,discriminator,weight,pVertexFilterCutGplusUpsPP),Form("pVertexFilterCutGplusUpsPP&&jtpt>=%f&&jtpt<%f&&%s>=%f&&%s<%f&&abs(jteta)>%f&&abs(jteta)<%f&&bin>=%d&&bin<%d",ptMin,ptMax,discr.c_str(),minXdiscr,discr.c_str(),maxXdiscr,etalo,etahi,cbinlo,cbinhi),"weight");
  //RooDataSet *data = new  RooDataSet("data","data",tdata,RooArgSet(s,jtpt,jteta,bin,discriminator),Form("jtpt>=%f&&jtpt<%f&&%s>=%f&&%s<%f&&abs(jteta)>%f&&abs(jteta)<%f&&bin>=%d&&bin<%d",ptMin,ptMax,discr,minXdiscr,discr,maxXdiscr,etalo,etahi,cbinlo,cbinhi));

  
  /*
     // --- Construct signal+background PDF ---
  //Double_t bInitFrac = hB->Integral()/(hB->Integral()+hCL->Integral());
  //Double_t cInitFrac = hC->Integral()/(hB->Integral()+hCL->Integral());
  RooRealVar Bfraction("Bfraction","#light events",0.3,0.,1);
  RooRealVar Cfraction("Cfraction","#background events",0.3,0.,1); hIncJetsMC
  if(fixCL) RooAddPdf model("model","",bottom,charmlight,Bfraction);
  else RooAddPdf model("model","",RooArgList(bottom,charm,light),RooArgList(Bfraction,Cfraction));  

  // --- Data sample ---
  //RooDataSet *data = new RooDataSet("data","data",tdata,RooArgSet(s,jtpt,discriminator),Form("jtpt>=%f&&jtpt<%f&&%s>=%f&&%s<%f",ptMin,ptMax,discr,minXdiscr,discr,maxXdiscr));
  //RooDataSet *data = new RooDataSet("data","data",tdata,RooArgSet(s,jtpt,discriminator),Form("jtpt>=%f&&jtpt<%f&&%s>=%f&&%s<%f&&bin>=%d&&bin<%d&&fabs(jteta)>%f&&fabs(jteta)<%f",ptMin,ptMax,discr,minXdiscr,discr,maxXdiscr,cbinlo,cbinhi,etalo,etahi));
  RooDataSet *data = new RooDataSet("data","data",tdata,RooArgSet(s,jtpt,jteta,bin,discriminator),Form("jtpt>=%f&&jtpt<%f&&%s>=%f&&%s<%f&&abs(jteta)>%f&&abs(jteta)<%f&&bin>=%d&&bin<%d",ptMin,ptMax,discr,minXdiscr,discr,maxXdiscr,etalo,etahi,cbinlo,cbinhi));
  */
  // unfucked?


  TPaveText *header = new TPaveText(0.05,0.9,0.95,0.99);
  header->AddText(Form("%s  -  ROOFIT ML unbinned fit of %s",var.c_str(),fixCL?"2 components : bottom and (charm + light)":"3 components : bottom, charm and light"));
  header->AddText(Form("Pb-Pb data - %s",comment.c_str()));
  header->AddText(Form("%s%.0f <= jet pT < %.0f",(var=="muptrel")?"deltaR < 0.5 ; muon pT > 5 ; ":"",ptMin,ptMax));
  //header->SetTextSize(0.027);
  header->SetTextSize(20);
  header->SetTextAlign(12);
  header->SetBorderSize(0);
  header->SetFillStyle(0);
  //header->Draw();


  //RooPlot* sframe = s.frame();
  TH2D *htemp = new TH2D(Form("%s%.0f%.0f",var.c_str(),ptMin,ptMax),Form("%s%.0f%.0f",var.c_str(),ptMin,ptMax),100,minXvar,maxXvar,100,0.5,maxYaxis) ;
  //htemp->SetXTitle(Form("%s %.0f < p_{T} < %.0f GeV/c",var,ptMin,ptMax));
  if(var=="svtxm")htemp->SetXTitle("SV mass (GeV/c^{2})");
  else htemp->SetXTitle("JP Disc.");
  htemp->SetYTitle("Entries");


  // --- Perform extended ML fit of composite PDF to data ---
  RooFitResult *fitresult = model->fitTo(*data,SumW2Error(kTRUE),Save(),PrintLevel(-1));
  //RooFitResult *fitresult = model->fitTo(*data,Save());
  

  RooPlot* sframe = s.frame();
  //  sframe = s.frame();

  htemp->Draw();
  if(var=="svtxm")data->plotOn(sframe,Binning(24));
  else data->plotOn(sframe,Binning(30));

  if(fixCL) {
    model->plotOn(sframe,Components(charmlight),LineStyle(kDashed),LineColor(30),LineWidth(2));
  } else {
    model->plotOn(sframe,Components(light),LineStyle(kDashed),LineColor(kBlue),LineWidth(2));
    model->plotOn(sframe,Components(charm),LineStyle(kDashed),LineColor(kGreen),LineWidth(2));
    model->plotOn(sframe,Components(RooArgSet(charm,light)),LineStyle(kDashed),LineColor(30),LineWidth(2));
  }
  model->plotOn(sframe,Components(bottom),LineStyle(kDashed),LineColor(kRed),LineWidth(2),FillColor(kRed),FillStyle(1));   
  model->plotOn(sframe,LineWidth(2),VisualizeError(*fitresult),FillColor(17));
  model->plotOn(sframe,LineWidth(2),LineColor(13));
  if(var=="svtxm")data->plotOn(sframe,Binning(24));
  else data->plotOn(sframe,Binning(30));

  model->paramOn(sframe,Layout(0.4,0.9,0.9),Format("NEU",FixedPrecision(3)));
  sframe->Draw("same");
//   TLegend *leg = new TLegend(0.61,fixCL?0.60:0.50,0.98,fixCL?0.78:0.75);
//   leg->SetBorderSize(0);
//   leg->SetFillStyle(0);
//   leg->AddEntry("h_data","PbPb data","p");
//   leg->AddEntry(Form("model_Norm[%s]_Comp[bottom]",var),"b","l");
//   if(fixCL) {
//     leg->AddEntry(Form("model_Norm[%s]_Comp[charmlight]",var),"c + udsg","l");
//   } else {
//     leg->AddEntry(Form("model_Norm[%s]_Comp[charm]",var),"c","l");
//     leg->AddEntry(Form("model_Norm[%s]_Comp[light]",var),"udsg","l");
//     leg->AddEntry(Form("model_Norm[%s]_Comp[charm,light]",var),"c + udsg","l");    
//   }
//   leg->AddEntry(Form("model_Norm[%s]",var),"b + c + udsg","l");
//   leg->Draw("same");

  //////////////////////////////////////////////////////////
  //Plot Stacked histos
  //////////////////////////////////////////////////////////
  int nXbins = 12;
  bool doLog = true;

  hMCB[counter] = new TH1D(Form("hMCB_%d",counter),Form("hMCB_%d",counter),nXbins,minXvar,maxXvar);
  hMCC[counter] = new TH1D(Form("hMCC_%d",counter),Form("hMCC_%d",counter),nXbins,minXvar,maxXvar);
  hMCL[counter] = new TH1D(Form("hMCL_%d",counter),Form("hMCL_%d",counter),nXbins,minXvar,maxXvar);
  hMCLC[counter] = new TH1D(Form("hMCLC_%d",counter),Form("hMCLC_%d",counter),nXbins,minXvar,maxXvar);
  hData[counter] = new TH1D(Form("hData_%d",counter),Form("hData_%d",counter),nXbins,minXvar,maxXvar);
  hData[counter]->Sumw2();
  hMCL[counter]->Sumw2();  hMCB[counter]->Sumw2();  hMCC[counter]->Sumw2(); hMCLC[counter]->Sumw2();

  tBMC->Draw(Form("%s>>hMCB_%d",var.c_str(),counter),Form("weight*(abs(refparton_flavorForB)==5&&pVertexFilterCutGplusUpsPP&&jtpt>=%f&&jtpt<%f&&%s>=%f&&%s<%f&&pthat>%f)",ptMin,ptMax,discr.c_str(),minXdiscr,discr.c_str(),maxXdiscr,ptHatMin),"goff");
  tCMC->Draw(Form("%s>>hMCC_%d",var.c_str(),counter),Form("weight*(abs(refparton_flavorForB)==4&&pVertexFilterCutGplusUpsPP&&jtpt>=%f&&jtpt<%f&&%s>=%f&&%s<%f&&pthat>%f)",ptMin,ptMax,discr.c_str(),minXdiscr,discr.c_str(),maxXdiscr,ptHatMin),"goff");
  tQCDMC->Draw(Form("%s>>hMCL_%d",var.c_str(),counter),Form("weight*(abs(refparton_flavorForB)!=5&&pVertexFilterCutGplusUpsPP&&abs(refparton_flavorForB)!=4&&abs(refparton_flavorForB)<99&&jtpt>=%f&&jtpt<%f&&%s>=%f&&%s<%f&&pthat>%f)",ptMin,ptMax,discr.c_str(),minXdiscr,discr.c_str(),maxXdiscr,ptHatMin),"goff");

  /*
  tBMC->Draw(Form("%s>>hMCB_%d",var,counter),Form("weight*(abs(refparton_flavorForB)==5&&jtpt>=%f&&jtpt<%f&&%s>=%f&&%s<%f)",ptMin,ptMax,discr,minXdiscr,discr,maxXdiscr),"goff");
  tCMC->Draw(Form("%s>>hMCC_%d",var,counter),Form("weight*(abs(refparton_flavorForB)==4&&jtpt>=%f&&jtpt<%f&&%s>=%f&&%s<%f)",ptMin,ptMax,discr,minXdiscr,discr,maxXdiscr),"goff");
  tQCDMC->Draw(Form("%s>>hMCL_%d",var,counter),Form("weight*(abs(refparton_flavorForB)!=5&&abs(refparton_flavorForB)!=4&&abs(refparton_flavorForB)<99&&jtpt>=%f&&jtpt<%f&&%s>=%f&&%s<%f)",ptMin,ptMax,discr,minXdiscr,discr,maxXdiscr),"goff");
  */
  tdata->Draw(Form("%s>>hData_%d",var.c_str(),counter),Form("weight*(pVertexFilterCutGplusUpsPP&&jtpt>=%f&&jtpt<%f&&%s>=%f&&%s<%f)",ptMin,ptMax,discr.c_str(),minXdiscr,discr.c_str(),maxXdiscr),"goff");
  //tdata->Draw(Form("%s>>hData_%d",var,counter),Form("(jtpt>=%f&&jtpt<%f&&%s>=%f&&%s<%f)",ptMin,ptMax,discr,minXdiscr,discr,maxXdiscr),"goff");
  hMCLC[counter]->Add( hMCL[counter]);
  hMCLC[counter]->Add( hMCC[counter]);
  fixEmpty(hMCB[counter]); fixEmpty(hMCC[counter]); fixEmpty(hMCL[counter]); fixEmpty(hMCLC[counter]); fixEmpty(hData[counter]);
  
  
  can1[counter] = new TCanvas(Form("can1_%d",counter),Form("can1_%d",counter),700,600);
  hs[counter] = new THStack(Form("hs_%d",counter),"le stack of MC histos");
  ghosths[counter] = new THStack(Form("ghosths_%d",counter),"le fake stack of MC histos");

  can1[counter]->cd();
  if (doLog) can1[counter]->cd()->SetLogy();


  if (doLog){
    hData[counter]->SetMaximum(hData[counter]->GetMaximum()*50);
    hData[counter]->SetMinimum(2);
  }
  if (!doLog){
    hData[counter]->SetMaximum(hData[counter]->GetMaximum()*1.5);
    hData[counter]->SetMinimum(0.0);
  }
  const char* yTitle = NULL;
  const char* xTitle = NULL;
  if (var=="svtxm")  xTitle = Form("Secondary vertex mass (GeV/c^{2})");
  if (var=="discr_prob") xTitle = Form("Jet probability");
  //if (var=="svtxm") yTitle = (Form("Number of Jets / %3.2f GeV",(maxXvar-minXvar)/nXbins));
  //if (var=="discr_prob") yTitle = (Form("Numbers of Jets / %3.2f",(maxXvar-minXvar)/nXbins));
  if (var=="svtxm") yTitle = ("Number of jets");
  if (var=="discr_prob") yTitle = ("Number of jets");
  hData[counter]->GetXaxis()->CenterTitle();
  hData[counter]->GetYaxis()->CenterTitle();
  hData[counter]->GetYaxis()->SetTitle(yTitle);
  hData[counter]->GetXaxis()->SetTitle(xTitle);
  hData[counter]->SetMarkerStyle(8);
  hData[counter]->SetMarkerSize(1.5);
  //hData[counter]->SetMinismum(2);
  hData[counter]->Draw();


  double Bnorm=0, Cnorm=0, Lnorm=0, LCnorm=0;
  Double_t Bfrac =Bfraction->getVal();
  Double_t Cfrac =Cfraction->getVal();
  
  //Normalize Histograms
  if(!fixCL) Lnorm = (hData[counter]->Integral(1,nXbins)/ hMCL[counter]->Integral(1,nXbins))*(1.-Cfrac-Bfrac);
  if(!fixCL) Cnorm = (hData[counter]->Integral(1,nXbins)/ hMCC[counter]->Integral(1,nXbins))*(Cfrac);
  if(!fixCL) Bnorm = (hData[counter]->Integral(1,nXbins)/ hMCB[counter]->Integral(1,nXbins))*(Bfrac);
  if(fixCL) LCnorm = (hData[counter]->Integral(1,nXbins)/hMCLC[counter]->Integral(1,nXbins))*(1.-Bfrac);
  if(fixCL)  Bnorm = (hData[counter]->Integral(1,nXbins)/ hMCB[counter]->Integral(1,nXbins))*(Bfrac);

  hMCB[counter]->SetFillColor(kRed-3);
  hMCL[counter]->SetFillColor(kBlue-3);
  hMCC[counter]->SetFillColor(kGreen-3);
  hMCLC[counter]->SetFillColor(kBlue-3);
//   hMCB[counter]->SetLineColor(kRed-4);
//   hMCL[counter]->SetLineColor(kBlue-3);
//   hMCC[counter]->SetLineColor(kGreen-3);
//   hMCLC[counter]->SetLineColor(kBlue-3);
  hMCB[counter]->SetMarkerSize(0);
  hMCC[counter]->SetMarkerSize(0);
  hMCL[counter]->SetMarkerSize(0);
  hMCLC[counter]->SetMarkerSize(0);
  //hMCB[counter]->Scale(Bnorm);
  //hMCC[counter]->Scale(Cnorm);
  //hMCL[counter]->Scale(Lnorm);
  //hMCLC[counter]->Scale(LCnorm);
  if (!fixCL){
    hMCB[counter]->Scale(Bnorm);
    hMCC[counter]->Scale(Cnorm);
    hMCL[counter]->Scale(Lnorm);
    hs[counter]->Add(hMCB[counter]);
    hs[counter]->Add(hMCC[counter]);
    hs[counter]->Add(hMCL[counter]);
  }
  if (fixCL){
    hMCB[counter]->Scale(Bnorm);
    hMCLC[counter]->Scale(LCnorm);
    hs[counter]->Add(hMCB[counter]);
    hs[counter]->Add(hMCLC[counter]);
  }

  if(!fixCL)hs[counter]->Draw("same h e");
  hData[counter]->Draw("same");


  //This is a fake THstack to also plot the constribution from the charm, 
  //even though it is supposed to be merged with the Light contribution 
  //so first we will need to keep the b constribution, AS-IS
  //then, find the relative ratio of the c to Light jets with the already-found
  //fraction (1-Bfraction) 
  if(fixCL){
    ghosths[counter]->Add(hMCB[counter]);
    double relC2LCfrac =(hMCC[counter]->Integral(1,nXbins))/(hMCC[counter]->Integral(1,nXbins) + hMCL[counter]->Integral(1,nXbins));
    double relL2LCfrac =(hMCL[counter]->Integral(1,nXbins))/(hMCC[counter]->Integral(1,nXbins) + hMCL[counter]->Integral(1,nXbins));
    fixCL_MCL[counter] = (TH1D*)hMCLC[counter]->Clone(Form("fixCL_MCL_%d",counter));
    fixCL_MCL[counter]->Scale(relL2LCfrac);
    fixCL_MCC[counter] = (TH1D*)hMCLC[counter]->Clone(Form("fixCL_MCC_%d",counter));
    fixCL_MCC[counter]->Scale(relC2LCfrac);
    fixCL_MCC[counter]->SetFillColor(kGreen-3);
    fixCL_MCL[counter]->SetFillColor(kBlue-3);  
    fixCL_MCC[counter]->SetMarkerSize(0);
    fixCL_MCL[counter]->SetMarkerSize(0);    
    ghosths[counter]->Add(fixCL_MCC[counter]);
    ghosths[counter]->Add(fixCL_MCL[counter]);
    ghosths[counter]->Draw("same h e"); 
  }
  //To obtain the overall sum histogram of the MC counts
  MCTotal[counter] = new TH1D(Form("MCTotal_%d",counter),Form("MCTotal_%d",counter),nXbins,minXvar,maxXvar);
  MCTotal[counter]->Sumw2();
  if (!fixCL){
    MCTotal[counter]->Add(hMCB[counter]);
    MCTotal[counter]->Add(hMCC[counter]);
    MCTotal[counter]->Add(hMCL[counter]);
  }
  if (fixCL){
    MCTotal[counter]->Add(hMCB[counter]);
    MCTotal[counter]->Add(hMCLC[counter]);
  }
  MCTotal[counter]->SetLineWidth(3);
  MCTotal[counter]->SetMarkerSize(0);
  MCTotal[counter]->SetMarkerColor(kGray+2);    
  MCTotal[counter]->SetLineColor(kBlue-9); 

  MCSumGreen[counter] = new TH1D(Form("MCSumGreen_%d",counter),Form("MCSumGreen_%d",counter),nXbins,minXvar,maxXvar);
  MCSumGreen[counter]->Sumw2();
  if (fixCL){
  MCSumGreen[counter]->Add(hMCB[counter]);
  MCSumGreen[counter]->Add(fixCL_MCC[counter]);
  }
  MCSumGreen[counter]->SetLineWidth(3);
  MCSumGreen[counter]->SetMarkerSize(0);
  MCSumGreen[counter]->SetMarkerColor(kGreen-7);    
  MCSumGreen[counter]->SetLineColor(kGreen-7); 

 

  //http://root.cern.ch/root/htmldoc/TH1.html#TH1:Chi2Test
  //Double_t chi2 = hData[counter]->Chi2Test(MCTotal[counter],"UW CHI2 P");
  //Double_t chi2NDF = hData[counter]->Chi2Test(MCTotal[counter],"UW CHI2/NDF P");
  Double_t chi2 = hData[counter]->Chi2Test(MCTotal[counter],"WW CHI2 P");
  Double_t chi2NDF = hData[counter]->Chi2Test(MCTotal[counter],"WW CHI2/NDF P");
  MCTotal[counter]->Draw("same,e");
  MCSumGreen[counter]->Draw("same,e");


  //  TH1F *hTaggedJetsMC = (TH1F*)hTaggedLJetsMC->Clone("hTaggedJetsMC");
  //Redraw some partial histograms to show error bars in between
  if(fixCL){
    hMCB_copy[counter] = new TH1D(Form("hMCB_copy_%d",counter),Form("hMCB_copy_%d",counter),nXbins,minXvar,maxXvar);
    hMCB_copy[counter]->Sumw2();
    hMCB_copy[counter] = (TH1D*)hMCB[counter]->Clone(Form("hMCB_copy_%d",counter));
    hMCB_copy[counter]->SetMarkerSize(0);
    hMCB_copy[counter]->SetLineWidth(3);
    hMCB_copy[counter]->SetLineColor(kRed-7); 
    hMCB_copy[counter]->SetMarkerColor(kRed-7);
    hMCB_copy[counter]->Draw("same,e");
  }
  hData[counter]->Draw("same");
  can1[counter]->GetFrame()->SetLineWidth(4);
  can1[counter]->RedrawAxis();
  TLegend *hleg = new TLegend(0.63,fixCL?0.70:0.75,0.883,fixCL?0.92:0.92);
  hleg->SetBorderSize(0);
  hleg->SetFillStyle(0);
  //hleg->SetTextSize(30);
  hleg->AddEntry(hData[counter],"pp data","lp");
  hleg->AddEntry(hMCB[counter],"b","f");
  if(!fixCL){
    hleg->AddEntry(hMCC[counter],"c","f");
    hleg->AddEntry(hMCL[counter],"udsg","f");
  }
  if(fixCL){
    hleg->AddEntry(fixCL_MCC[counter],"c","f");
    hleg->AddEntry(fixCL_MCL[counter],"usdg","f");
  }
  hleg->Draw("same");

  drawText(Form("%2.0f < p_{T} < %2.0f GeV/c",ptMin,ptMax),0.63,0.684);
  //if(fixCL)drawText(Form("%2.0f < p_{T} < %2.0f GeV/c",ptMin,ptMax),0.50,0.66);
  //drawText(Form("#chi^{2}/NDF = %3.1f / %d",chi2NDF,nXbins-1),0.53,0.55);
  drawText(Form("#chi^{2}/NDF = %3.1f / %d",chi2,(int)(chi2/chi2NDF)),0.632,0.615);
  drawText("CMS Preliminary",0.15,0.965);
  drawText("#sqrt{s_{NN}} = 2.76 TeV",0.60,0.969);
  drawText("|#eta| < 2.0",0.18,0.88);
  if(comment!="b-tagged sample (SSVHE > 2)") drawText(comment,0.18,0.80);
  if(comment=="b-tagged sample (SSVHE > 2)"){
    drawText("b-tagged sample",0.18,0.80);
    //drawText("(SSVHE > 2)",0.18,0.75);
  }

  
  bool printEach=false;
  // --- Print results ---
  //cout <<"b jet fraction in MC = "<<bInitFrac<<endl;
  if(verbose){
    cout<<"ZZZZZ TEXTcounter: "<<counter<<"   Comment: "<<comment<<endl;
    if(!fixCL) cout<<"ZZZZZ Scale factors: B= "<<Bnorm<<" C= "<<Cnorm<<" Light= "<<Lnorm<<endl;
    if(fixCL)  cout<<"ZZZZZ Scale factors: B= "<<Bnorm<<" LC= "<<LCnorm<<endl;
    cout<<"ZZZZZ Chi2: "<<chi2<<" Chi2/NDF: "<<chi2NDF<<endl;
    cout<<"ZZZZZ var: "<<var<<"  discr: "<<discr<<"    min(discr): "<<minXvar<<"    max(discr): "<<maxXvar<<endl;
    cout<<"ZZZZZ in pT [ "<<ptMin<<" , "<<ptMax<<" ]"<<endl;
    cout<<"ZZZZZ b jet fraction in data = "<<Bfraction->getVal()<<endl;
    if(!fixCL) cout <<"ZZZZZ c jet fraction = "<<Cfraction->getVal()<<endl;
  }
  std::string fixCLlabel;
  if (fixCL) fixCLlabel= "";
  if (!fixCL) fixCLlabel= "_floatC";
  std::string varLabel;
  if (var =="svtxm") varLabel = "SVmass";
  if (var =="discr_prob") varLabel = "JP";
  if (var =="discr_prob" && minXdiscr==2 ) varLabel ="JPtagged";
  std::string ptLabel;

  if (ptMin==60) ptLabel="60";
  if (ptMin==70) ptLabel="70";
  if (ptMin==80) ptLabel="80";
  if (ptMin==90) ptLabel="90";
  if (ptMin==105) ptLabel="110";
  if (ptMin==130) ptLabel="130";
  if (ptMin==170) ptLabel="170";

  //cout<<"file going to :"<<Form("PDFS/fitpp_%s_%s_%s.pdf",varLabel,fixCLlabel,ptLabel)<<endl;
  
  if (var =="svtxm"){
    if (printEach) can1[counter]->Print(Form("PDFS/fitpp_fixcl_%s_%f%s_%s.pdf",varLabel.c_str(),minXdiscr,fixCLlabel.c_str(),ptLabel.c_str()),"pdf");
    if (printEach) can1[counter]->Print(Form("MACROS/fitpp_%s_%f%s_%s.C",varLabel.c_str(),minXdiscr,fixCLlabel.c_str(),ptLabel.c_str()),"cxx");  
  }
  
  if(!fixCL){
    if (!printEach){
      if (counter<nBins*3-1) can1[counter]->Print("PDFS/bTagStackedHistos_nofixCL_NewFormatV5.pdf(","pdf");
      if (counter==nBins*3-1) can1[counter]->Print("PDFS/bTagStackedHistos_nofixCL_NewFormatV5.pdf)","pdf");
    }

  }
  if(fixCL){
    if (!printEach){
      if (counter<nBins*3-1) can1[counter]->Print("PDFS/bTagStackedHistos_fixCL_NewFormatV5.pdf(","pdf");
      if (counter==nBins*3-1) can1[counter]->Print("PDFS/bTagStackedHistos_fixCL_NewFormatV5.pdf)","pdf");
    }
  }
  counter++;
  

  
  // --- Save canvas ---
  TString path = Form("gifs/%s_jtpt%.0fto%.0f_%s%.2fto%.2f_%s.gif",var.c_str(),ptMin,ptMax,discr.c_str(),minXdiscr,maxXdiscr,fixCL?"CLfixed":"CLfree");
  //cROOFIT->SaveAs(path);

  // ========== Toy check on MC template ==========
  if (toyMC) {
     TH1F *hToyResult = new TH1F("hToyResult","",200,0,Bfraction->getVal()*2);
     char *pathToy = Form("toyMC/jtpt%.0fto%.0f_%s%.2fto%.2f_%s",ptMin,ptMax,discr.c_str(),minXdiscr,maxXdiscr,fixCL?"CLfixed":"CLfree");
     TCanvas *cToy = new TCanvas("cToy",pathToy,700,600);
     int nExp = 100;
     for (int iExp=0;iExp<nExp;iExp++) {
        TH1D* hBToy  = fluctuateHist(hB);
        TH1D* hCToy  = fluctuateHist(hC);
        TH1D* hLToy  = fluctuateHist(hL);
        TH1D* hCLToy = fluctuateHist(hCL);
	
        RooDataHist xBToy("xBToy","xBToy",s,hBToy);
        RooHistPdf bottomToy("bottomToy","bottom Toy PDF",s,xB);
        RooDataHist xCToy("xCToy","xCToy",s,hCToy);
        RooHistPdf charmToy("charmToy","charm Toy PDF",s,xCToy);
        RooDataHist xLToy("xLToy","xLToy",s,hLToy);
        RooHistPdf lightToy("lightToy","light Toy PDF",s,xLToy);
        RooDataHist xCLToy("xCLToy","xCLToy",s,hCLToy);
        RooHistPdf charmlightToy("charmlightToy","charmlight Toy PDF",s,xCLToy);

        RooRealVar *BfractionToy = new RooRealVar("BfractionToy","#light events",0.3,0.,1);
        RooRealVar CfractionToy("CfractionToy","#background events",0.3,0.,1); 
	RooAddPdf modelToy=RooAddPdf("modelToy","",bottomToy,charmlightToy,*BfractionToy);
        if(!fixCL) modelToy= RooAddPdf("modelToy","",RooArgList(bottomToy,charmToy,lightToy),RooArgList(*BfractionToy,CfractionToy));  
  
        //RooFitResult *fitresultToy = modelToy.fitTo(*data,Save(),PrintLevel(-1));
        hToyResult->Fill(BfractionToy->getVal());
	delete hBToy;
	delete hCToy;
	delete hLToy;
	delete hCLToy;
     }	
     hToyResult->Draw();
     cToy->SaveAs(Form("toyMC/jtpt_cent_%.0f-%.0f_%.0fto%.0f_%s%.2fto%.2f_%s.gif",cbinlo*2.5,cbinhi*2.5,ptMin,ptMax,discr.c_str(),minXdiscr,maxXdiscr,fixCL?"CLfixed":"CLfree"));
     cToy->SaveAs(Form("toyMC/jtpt_cent_%.0f-%.0f_%.0fto%.0f_%s%.2fto%.2f_%s.C",cbinlo*2.5,cbinhi*2.5,ptMin,ptMax,discr.c_str(),minXdiscr,maxXdiscr,fixCL?"CLfixed":"CLfree"));
     delete hToyResult;
     delete cToy;
  }  
  /*
  delete fdata;
  delete fQCDMC;
  delete fBMC;
  delete fCMC;
  */

  delete hB;
  delete hC;
  delete hL;
  delete hCL;

  return Bfraction;
  
}

void drawText(const std::string text, float xp, float yp){
  TLatex *tex = new TLatex(xp,yp,text.c_str());
  tex->SetTextFont(43);
  tex->SetTextSize(23);
  tex->SetTextColor(kBlack);
  tex->SetLineWidth(1);
  tex->SetNDC();
  tex->Draw();
}
void format1D(TH1& h1, TCanvas& c1 ){
  h1.GetXaxis()->CenterTitle();
  h1.GetYaxis()->CenterTitle();
  c1.cd();
  c1.RedrawAxis();
  return;
}

void fillCounterHistos(std::string discr, double workingPoint, int cbinlo, int cbinhi, float etalo, float etahi) {

  hTaggedLJetsMC = new TH1D("hTaggedLJetsMC","hTaggedLJetsMC",nBins,ptBin);
  hTaggedLJetsMC->Sumw2();
  tQCDMC->Draw("jtpt>>hTaggedLJetsMC",Form("weight*(pVertexFilterCutGplusUpsPP&&%s>=%f&&bin>=%d&&bin<%d&&abs(jteta)>%f&&abs(jteta)<%f&&abs(refparton_flavorForB)!=4&&abs(refparton_flavorForB)!=5)",discr.c_str(),workingPoint,cbinlo,cbinhi,etalo,etahi));

  hTaggedCJetsMC = new TH1D("hTaggedCJetsMC","hTaggedCJetsMC",nBins,ptBin);
  hTaggedCJetsMC->Sumw2();
  tCMC->Draw("jtpt>>hTaggedCJetsMC",Form("weight*(pVertexFilterCutGplusUpsPP&&%s>=%f&&bin>=%d&&bin<%d&&abs(jteta)>%f&&abs(jteta)<%f&&abs(refparton_flavorForB)==4)",discr.c_str(),workingPoint,cbinlo,cbinhi,etalo,etahi));
  
  hTaggedBJetsMC = new TH1D("hTaggedBJetsMC","hTaggedBJetsMC",nBins,ptBin);
  hTaggedBJetsMC->Sumw2();
  tBMC->Draw("jtpt>>hTaggedBJetsMC",Form("weight*(pVertexFilterCutGplusUpsPP&&%s>=%f&&bin>=%d&&bin<%d&&abs(jteta)>%f&&abs(jteta)<%f&&abs(refparton_flavorForB)==5)",discr.c_str(),workingPoint,cbinlo,cbinhi,etalo,etahi));
  
  hTaggedJetsMC = (TH1D*)hTaggedLJetsMC->Clone("hTaggedJetsMC");
  
  hTaggedJetsMC->Add(hTaggedCJetsMC);
  hTaggedJetsMC->Add(hTaggedBJetsMC);


  hUntaggedLJetsMC = new TH1D("hUntaggedLJetsMC","hUntaggedLJetsMC",nBins,ptBin);
  hUntaggedLJetsMC->Sumw2();
  tQCDMC->Draw("jtpt>>hUntaggedLJetsMC",Form("weight*(pVertexFilterCutGplusUpsPP&&%s<%f&&bin>=%d&&bin<%d&&abs(jteta)>%f&&abs(jteta)<%f&&abs(refparton_flavorForB)!=4&&abs(refparton_flavorForB)!=5&&refpt>0)",discr.c_str(),workingPoint,cbinlo,cbinhi,etalo,etahi));

  hUntaggedCJetsMC = new TH1D("hUntaggedCJetsMC","hUntaggedCJetsMC",nBins,ptBin);
  hUntaggedCJetsMC->Sumw2();
  tCMC->Draw("jtpt>>hUntaggedCJetsMC",Form("weight*(pVertexFilterCutGplusUpsPP&&%s<%f&&bin>=%d&&bin<%d&&abs(jteta)>%f&&abs(jteta)<%f&&abs(refparton_flavorForB)==4&&refpt>0)",discr.c_str(),workingPoint,cbinlo,cbinhi,etalo,etahi));
  
  hUntaggedBJetsMC = new TH1D("hUntaggedBJetsMC","hUntaggedBJetsMC",nBins,ptBin);
  hUntaggedBJetsMC->Sumw2();
  tBMC->Draw("jtpt>>hUntaggedBJetsMC",Form("weight*(pVertexFilterCutGplusUpsPP&&%s<%f&&bin>=%d&&bin<%d&&abs(jteta)>%f&&abs(jteta)<%f&&abs(refparton_flavorForB)==5&&refpt>0)",discr.c_str(),workingPoint,cbinlo,cbinhi,etalo,etahi));
  
  hUntaggedJetsMC = (TH1D*)hUntaggedLJetsMC->Clone("hUntaggedJetsMC");
  
  hUntaggedJetsMC->Add(hUntaggedCJetsMC);
  hUntaggedJetsMC->Add(hUntaggedBJetsMC);

  hBjetsWithJPinfoMC = new TH1D("hBjetsWithJPinfoMC","",nBins,ptBin);
  hBjetsWithJPinfoMC->Sumw2();
  tBMC->Draw("jtpt>>hBjetsWithJPinfoMC",Form("weight*(pVertexFilterCutGplusUpsPP&&abs(refparton_flavorForB)==5&&discr_prob>0&&bin>=%d&&bin<%d&&abs(jteta)>%f&&abs(jteta)<%f&&refpt>0)",cbinlo,cbinhi,etalo,etahi));

  hBjetsWithCSVinfoMC = new TH1D("hBjetsWithCSVinfoMC","",nBins,ptBin);
  hBjetsWithCSVinfoMC->Sumw2();
  tBMC->Draw("jtpt>>hBjetsWithCSVinfoMC",Form("weight*(pVertexFilterCutGplusUpsPP&&abs(refparton_flavorForB)==5&&discr_prob>0&&bin>=%d&&bin<%d&&abs(jteta)>%f&&abs(jteta)<%f&&refpt>0)",cbinlo,cbinhi,etalo,etahi));
 
  hTaggedJetsData = new TH1D("hTaggedJetsData","",nBins,ptBin);
  hTaggedJetsData->Sumw2();
  tdata->Draw("jtpt>>hTaggedJetsData",Form("weight*(pVertexFilterCutGplusUpsPP&&%s>=%f&&bin>=%d&&bin<%d&&abs(jteta)>%f&&abs(jteta)<%f)",discr.c_str(),workingPoint,cbinlo,cbinhi,etalo,etahi));

  hUntaggedJetsData = new TH1D("hUntaggedJetsData","",nBins,ptBin);
  hUntaggedJetsData->Sumw2();
  tdata->Draw("jtpt>>hUntaggedJetsData",Form("weight*(pVertexFilterCutGplusUpsPP&&%s<%f&&bin>=%d&&bin<%d&&abs(jteta)>%f&&abs(jteta)<%f)",discr.c_str(),workingPoint,cbinlo,cbinhi,etalo,etahi));
}

Enumerations count(double ptMin, double ptMax) {

  Enumerations res;
  int currentBin = hTaggedJetsMC->FindBin(ptMin+1.);
  if(hTaggedJetsMC->FindBin(ptMax-1.) != currentBin) cout<<" BINNING PROBLEM "<<endl;

  res.nTaggedJetsMC = hTaggedJetsMC->GetBinContent(currentBin);
  res.nUntaggedJetsMC = hUntaggedJetsMC->GetBinContent(currentBin);
  res.nJetsMC = res.nTaggedJetsMC + res.nUntaggedJetsMC;
  res.nTaggedBjetsMC = hTaggedBJetsMC->GetBinContent(currentBin);
  res.nUntaggedBjetsMC = hUntaggedBJetsMC->GetBinContent(currentBin);
  res.nBjetsMC = res.nTaggedBjetsMC + res.nUntaggedBjetsMC;
  res.nNonBjetsMC = res.nJetsMC - res.nBjetsMC;
  res.nTaggedNonBjetsMC = res.nTaggedJetsMC - res.nTaggedBjetsMC;

  res.nTaggedJetsData = hTaggedJetsData->GetBinContent(currentBin);
  res.nUntaggedJetsData = hUntaggedJetsData->GetBinContent(currentBin);

  res.cbForJP = hBjetsWithJPinfoMC->GetBinContent(currentBin) / res.nBjetsMC;
  res.cbForCSV = hBjetsWithCSVinfoMC->GetBinContent(currentBin) / res.nBjetsMC;

  res.nTaggedJetsMCError = hTaggedJetsMC->GetBinError(currentBin);
  res.nUntaggedJetsMCError = hUntaggedJetsMC->GetBinError(currentBin);
  res.nJetsMCError = addError(res.nTaggedJetsMCError,res.nUntaggedJetsMCError);
  res.nTaggedBjetsMCError = hTaggedBJetsMC->GetBinError(currentBin);
  res.nUntaggedBjetsMCError = hUntaggedBJetsMC->GetBinError(currentBin);
  res.nBjetsMCError = addError(res.nTaggedBjetsMCError,res.nUntaggedBjetsMCError);
  res.nNonBjetsMCError = substractError(res.nJetsMCError,res.nBjetsMCError);
  res.nTaggedNonBjetsMCError = substractError(res.nTaggedJetsMCError,res.nTaggedBjetsMCError);

  res.nTaggedJetsDataError = hTaggedJetsData->GetBinError(currentBin);
  res.nUntaggedJetsDataError = hUntaggedJetsData->GetBinError(currentBin);

  res.nJetsData = res.nTaggedJetsData + res.nUntaggedJetsData;
  res.nJetsDataError = addError(res.nTaggedJetsDataError,res.nUntaggedJetsDataError);
  
  return res;
}

