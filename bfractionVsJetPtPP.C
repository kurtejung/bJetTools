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
#include <RooFit.h>
#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooDataSet.h>
#include <RooHistPdf.h>
#include <RooPlot.h>
#include <RooAddPdf.h>
#include <RooGlobalFunc.h>
#endif

using namespace RooFit;

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

};
int counter(0);
TCanvas* can1[20];
TH1D* hData[20];
TH1D* hMCC[20];
TH1D* hMCB[20];
TH1D* hMCL[20];
TH1D* hMCLC[20];
TH1D* MCTotal[20];
THStack* hs[20];
THStack* fakehs[20];
void drawText(const char *text, float xp, float yp);




//void bfractionVsJetPt(char *tagger="discr_ssvHighPur", Double_t workingPoint=2, char *taggerName="SSVHP") {
void bfractionVsJetPtPP(char *tagger="discr_ssvHighEff", Double_t workingPoint=2, char *taggerName="SSVHE") {
  cout<<" OUTPUT TURNED OFF "<<endl;

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  gStyle->SetLabelFont(43,"xyz");
  gStyle->SetLabelSize(20,"xyz");
  gStyle->SetTitleFont(43,"xyz");
  gStyle->SetTitleSize(26,"xyz");
  gStyle->SetTitleOffset(1.0,"x"); 
  gStyle->SetNdivisions(510,"xy"); 

  gROOT->ForceStyle(1);


  Int_t isRecopp=0;
  Bool_t fixCL=1;
  Int_t doMuptrelFit=0;
  Int_t doLTJP=1;
  Int_t doLTCSV=0;
  Double_t CLshift=0.; // 0., -20. or 20.

  const Int_t nBins = 3;
  Double_t ptBin[nBins+1] = {80,100,120,200};

  //const Int_t nBins = 4;
//Double_t ptBin[nBins+1] = {60,80,100,120,200};
  //const Int_t nBins = 5;
  //Double_t ptBin[nBins+1] = {60,70,80,95,120,200};
  
  Double_t bPurMC, bPurData, bEffMCMuTrig, bEffDataMuTrig, SFb, bEffMC, bEffDataPtrel, bEffDataLTJP, bEffDataLTCSV, taggedFracData, bFracMC, bFracData, bFracDataLTJP, bFracDataLTCSV, bFracJPdirect;
  Double_t bPurMCError, bPurDataError, bEffMCError, bEffDataPtrelError, bEffDataLTJPError, bEffDataLTCSVError, taggedFracDataError, bFracMCError, bFracDataError, bFracDataLTJPError, bFracDataLTCSVError, bFracJPdirectError;
  Enumerations numbers, numbersMuTrig;
  
  TH1D *hBPurityData = new TH1D("hBPurityData","hBPurityData;Jet p_{T} (GeV/c);b-Tagging purity",nBins,ptBin);
  TH1D *hBPurityMC = new TH1D("hBPurityMC","hBPurityMC;Jet p_{T} (GeV/c);b-Tagging purity",nBins,ptBin);
  
  TH1D *hBEfficiencyMC = new TH1D("hBEfficiencyMC","hBEfficiencyMC;Jet p_{T} (GeV/c);b-Tagging efficiency",nBins,ptBin);
  TH1D *hBEfficiencyDataPtrel = new TH1D("hBEfficiencyDataPtrel","hBEfficiencyDataPtrel;Jet p_{T} (GeV/c);b-Tagging efficiency",nBins,ptBin);
  TH1D *hBEfficiencyDataLTJP = new TH1D("hBEfficiencyDataLTJP","hBEfficiencyDataLTJP;Jet p_{T} (GeV/c);b-Tagging efficiency",nBins,ptBin);
  TH1D *hBEfficiencyDataLTCSV = new TH1D("hBEfficiencyDataLTCSV","hBEfficiencyDataLTCSV;Jet p_{T} (GeV/c);b-Tagging efficiency",nBins,ptBin);
  
  TH1D *hBFractionMC = new TH1D("hBFractionMC","hBFractionMC;Jet p_{T} (GeV/c);b-jet fraction",nBins,ptBin);
  TH1D *hBFractionData = new TH1D(Form("hBFractionData_%sat%.1f_CLshift%.0f",taggerName,workingPoint,CLshift),"hBFractionData;Jet p_{T} (GeV/c);b-jet fraction",nBins,ptBin);
  TH1D *hBFractionDataLTJP = new TH1D(Form("hBFractionDataLTJP_%sat%.1f_CLshift%.0f",taggerName,workingPoint,CLshift),"hBFractionDataLTJP;Jet p_{T} (GeV/c);b-jet fraction",nBins,ptBin);
  TH1D *hBFractionDataLTCSV = new TH1D(Form("hBFractionDataLTCSV_%sat%.1f_CLshift%.0f",taggerName,workingPoint,CLshift),"hBFractionDataLTCSV;Jet p_{T} (GeV/c);b-jet fraction",nBins,ptBin);
  TH1D *hBFractionJPdirect = new TH1D(Form("hBFractionJPdirect_%sat%.1f_CLshift%.0f",taggerName,workingPoint,CLshift),"hBFractionJPdirect;Jet p_{T} (GeV/c);b-jet fraction",nBins,ptBin);
  

  for (Int_t n=0;n<nBins;n++) {
    //for (Int_t n=0;n<1;n++) {

    cout<<"Processing jet pT bin ["<<ptBin[n]<<","<<ptBin[n+1]<<"] ..."<<endl;

    numbers = count(isRecopp,0,ptBin[n],ptBin[n+1],tagger,workingPoint);
    RooRealVar fitSvtxmTag = bfractionFit(isRecopp,0,fixCL,CLshift,"svtxm",0,6,ptBin[n],ptBin[n+1],tagger,workingPoint,6,"b-tagged sample (SSVHE > 2)",5e4);
    RooRealVar fitJpDirect = bfractionFit(isRecopp,0,fixCL,CLshift,"discr_prob",0,3.,ptBin[n],ptBin[n+1],tagger,-999.,999.,"inclusive sample",5e4);
    
    if (doLTJP) {
      RooRealVar fitJpTag =       bfractionFit(isRecopp,0,fixCL,CLshift,"discr_prob",0.,3.,ptBin[n],ptBin[n+1],tagger,workingPoint,6,"b-tagged sample (SSVHE > 2)",5e4);
      RooRealVar fitJpBeforetag = bfractionFit(isRecopp,0,fixCL,CLshift,"discr_prob",0.,3.,ptBin[n],ptBin[n+1],"discr_prob",0,999,"jets with JP info",5e4);

    } 

    if (doLTCSV) {
      RooRealVar fitCsvBeforetag = bfractionFit(isRecopp,0,fixCL,CLshift,"discr_csvSimple",0,1,ptBin[n],ptBin[n+1],tagger,-999,999,"jets with CSV info",5e4);
      RooRealVar fitCsvTag = bfractionFit(isRecopp,0,fixCL,CLshift,"discr_csvSimple",0,1,ptBin[n],ptBin[n+1],tagger,workingPoint,999,"b-tagged sample (SSVHE > 2)",5e4);
    } 

    taggedFracData = numbers.nTaggedJetsData / (numbers.nTaggedJetsData+numbers.nUntaggedJetsData);
    taggedFracDataError = fracError(numbers.nTaggedJetsData,numbers.nUntaggedJetsData,numbers.nTaggedJetsDataError,numbers.nUntaggedJetsDataError);
    

    //*  --- b-tagging purity --- 

    bPurMC = numbers.nTaggedBjetsMC / numbers.nTaggedJetsMC;
    bPurMCError = fracError(numbers.nTaggedBjetsMC,numbers.nTaggedNonBjetsMC,numbers.nTaggedBjetsMCError,numbers.nTaggedNonBjetsMCError);
    bPurData = fitSvtxmTag.getVal();
    bPurDataError = fitSvtxmTag.getError();

    hBPurityMC->SetBinContent(n+1,bPurMC); 
    hBPurityMC->SetBinError(n+1,bPurMCError); 
    hBPurityData->SetBinContent(n+1,bPurData);    
    hBPurityData->SetBinError(n+1,bPurDataError); 
    //*/
    

    //*  --- b-tagging efficiency --- 

    bEffMC = numbers.nTaggedBjetsMC / numbers.nBjetsMC;
    bEffMCError = fracError(numbers.nTaggedBjetsMC,numbers.nUntaggedBjetsMC,numbers.nTaggedBjetsMCError,numbers.nUntaggedBjetsMCError);
    hBEfficiencyMC->SetBinContent(n+1,bEffMC); 
    hBEfficiencyMC->SetBinError(n+1,bEffMCError);

    if (doLTJP) {
      bEffDataLTJP = taggedFracData * numbers.cbForJP * fitJpTag.getVal() / fitJpBeforetag.getVal();
      bEffDataLTJPError = prodError(taggedFracData,fitJpTag.getVal(),taggedFracDataError,fitJpTag.getError()) * numbers.cbForJP / fitJpBeforetag.getVal(); 
      hBEfficiencyDataLTJP->SetBinContent(n+1,bEffDataLTJP);    
      hBEfficiencyDataLTJP->SetBinError(n+1,bEffDataLTJPError);
    } 

    if (doLTCSV) {
      bEffDataLTCSV = taggedFracData * numbers.cbForCSV * fitCsvTag.getVal() / fitCsvBeforetag.getVal();
      bEffDataLTCSVError = prodError(taggedFracData,fitCsvTag.getVal(),taggedFracDataError,fitCsvTag.getError()) * numbers.cbForCSV / fitCsvBeforetag.getVal(); 
      hBEfficiencyDataLTCSV->SetBinContent(n+1,bEffDataLTCSV);    
      hBEfficiencyDataLTCSV->SetBinError(n+1,bEffDataLTCSVError); 
    } 

    if (doMuptrelFit) { // muon triggered muptrel study 

      numbersMuTrig = count(isRecopp,1,ptBin[n],ptBin[n+1],tagger,workingPoint);
      bEffMCMuTrig = numbersMuTrig.nTaggedBjetsMC / numbersMuTrig.nBjetsMC;
      RooRealVar fitFbTag = bfractionFit(isRecopp,1,fixCL,CLshift,"muptrel",-0.5,5.5,ptBin[n],ptBin[n+1],tagger,workingPoint,999,Form("b-tagged sample (%s at %.1f)",taggerName,workingPoint)); 
      RooRealVar fitFbUntag = bfractionFit(isRecopp,1,fixCL,CLshift,"muptrel",-0.5,5.5,ptBin[n],ptBin[n+1],tagger,-999,workingPoint,Form("b-untagged sample (%s at %.1f)",taggerName,workingPoint)); 
      bEffDataMuTrig = (fitFbTag.getVal()*numbersMuTrig.nTaggedJetsData) / (fitFbTag.getVal()*numbersMuTrig.nTaggedJetsData + fitFbUntag.getVal()*numbersMuTrig.nUntaggedJetsData);
      SFb = bEffDataMuTrig / bEffMCMuTrig;
      bEffDataPtrel = bEffMC * SFb;
      bEffDataPtrelError = 0; 

      hBEfficiencyDataPtrel->SetBinContent(n+1,bEffDataPtrel);    
      hBEfficiencyDataPtrel->SetBinError(n+1,bEffDataPtrelError);
    }


    //*  --- b fraction --- 

    bFracMC = numbers.nBjetsMC / numbers.nJetsMC;
    //bFracMC = numbers.nTaggedJetsMC * bPurMC / (bEffMC * numbers.nJetsMC); // for check : same as previous
    bFracMCError = fracError(numbers.nBjetsMC,numbers.nNonBjetsMC,numbers.nBjetsMCError,numbers.nNonBjetsMCError); 
    hBFractionMC->SetBinContent(n+1,bFracMC); 
    hBFractionMC->SetBinError(n+1,bFracMCError);  

    bFracData = taggedFracData * bPurData / bEffMC; // efficiency from MC
    //bFracData = taggedFracData * bPurData / bEffDataPtrel; // efficiency from muptrel study
    bFracDataError = prodError(taggedFracData,bPurData,taggedFracDataError,bPurDataError) / bEffMC; // stat.error from purity and tagged-fraction (assumed independent)
    //bFracDataError = bFracData * bPurDataError / bPurData; // stat.error only from purity
    hBFractionData->SetBinContent(n+1,bFracData);    
    hBFractionData->SetBinError(n+1,bFracDataError);

    if (doLTJP) {
      bFracDataLTJP = taggedFracData * bPurData / bEffDataLTJP ; // efficiency from LTJP method
      bFracDataLTJPError = prodError(taggedFracData,bPurData,taggedFracDataError,bPurDataError) / bEffDataLTJP; // stat.error from purity and tagged-fraction (assumed independent)
      //bFracDataLTJPError = bFracDataLTJP * bPurDataError / bPurData; // stat.error only from purity
      hBFractionDataLTJP->SetBinContent(n+1,bFracDataLTJP);    
      hBFractionDataLTJP->SetBinError(n+1,bFracDataLTJPError);
    } 

    if (doLTCSV) {
      bFracDataLTCSV = taggedFracData * bPurData / bEffDataLTCSV; // efficiency from LTCSV method
      bFracDataLTCSVError = prodError(taggedFracData,bPurData,taggedFracDataError,bPurDataError) / bEffDataLTCSV; // stat.error from purity and tagged-fraction (assumed independent)
      //bFracDataLTCSVError = bFracDataLTCSV * bPurDataError / bPurData; // stat.error only from purity
      hBFractionDataLTCSV->SetBinContent(n+1,bFracDataLTCSV);    
      hBFractionDataLTCSV->SetBinError(n+1,bFracDataLTCSVError);
    } 

    bFracJPdirect = fitJpDirect.getVal();
    bFracJPdirectError = fitJpDirect.getError();
    hBFractionJPdirect->SetBinContent(n+1,bFracJPdirect);   
    hBFractionJPdirect->SetBinError(n+1,bFracJPdirectError);
    //*/


    //*
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
    cout<<"bEffMCMuTrig "<<bEffMCMuTrig<<endl;
    cout<<"bEffDataMuTrig "<<bEffDataMuTrig<<endl;
    cout<<"SFb "<<SFb<<endl;
    cout<<"bEffMC "<<bEffMC<<endl;
    cout<<"bEffDataPtrel "<<bEffDataPtrel<<endl;
    cout<<"CbForJP "<<numbers.cbForJP<<endl;
    cout<<"bEffDataLTJP "<<bEffDataLTJP<<endl;
    cout<<"CbForCSV "<<numbers.cbForCSV<<endl;
    cout<<"bEffDataLTCSV "<<bEffDataLTCSV<<endl;
    //cout<<"fbTag "<<fitFbTag.getVal()<<endl;
    //cout<<"fbUntag "<<fitFbUntag.getVal()<<endl;
    cout<<"bFracMC "<<bFracMC<<endl;
    cout<<"bFracMC "<<bFracMC<<endl;
    cout<<"bFracData "<<bFracData<<endl;
    cout<<"bFracDataLTJP "<<bFracDataLTJP<<endl;
    cout<<"bFracDataLTCSV "<<bFracDataLTCSV<<endl;
    cout<<"bFracJPdirect "<<bFracJPdirect<<endl;
    cout<<endl;
    //*/

  }

  TLegend *legPur = new TLegend(0.25,0.7,0.9,0.9,Form("Purity of b-tagged sample (%s > %.1f)",taggerName,workingPoint));
  legPur->SetBorderSize(0);
  legPur->SetFillStyle(0);
  legPur->AddEntry(hBPurityMC,"MC Input","p");
  legPur->AddEntry(hBPurityData,"Data","p");
  TCanvas *cBPurity = new TCanvas("cBPurity","b purity",600,600);
  hBPurityMC->SetAxisRange(0,1,"Y");
  hBPurityMC->SetTitleOffset(1.3,"Y");
  hBPurityMC->SetLineColor(2);
  hBPurityMC->SetMarkerColor(2);
  hBPurityMC->SetMarkerStyle(21);
  hBPurityMC->Draw("");
  hBPurityData->SetLineColor(1);
  hBPurityData->SetMarkerColor(1);
  hBPurityData->SetMarkerStyle(20);
  hBPurityData->Draw("same");   
  legPur->Draw();
  
  TLegend *legEff = new TLegend(0.25,0.7,0.9,0.9,Form("Efficiency for tagging b-jets (%s > %.1f)",taggerName,workingPoint));
  legEff->SetBorderSize(0);
  legEff->SetFillStyle(0);
  legEff->AddEntry(hBEfficiencyMC,"Simulation","p");
  if (doMuptrelFit) legEff->AddEntry(hBEfficiencyDataPtrel,"muon ptrel method","p");
  if (doLTJP) legEff->AddEntry(hBEfficiencyDataLTJP,"Reference Tagger","p");
  if (doLTCSV) legEff->AddEntry(hBEfficiencyDataLTCSV,"LT method (CSV)","p");
  TCanvas *cBEfficiency = new TCanvas("cBEfficiency","b-Tagging efficiency",600,600);
  hBEfficiencyMC->SetAxisRange(0,1,"Y");
  hBEfficiencyMC->SetTitleOffset(1.3,"Y");
  hBEfficiencyMC->SetLineColor(2);
  hBEfficiencyMC->SetMarkerColor(2);
  hBEfficiencyMC->SetMarkerStyle(21);
  hBEfficiencyMC->Draw();
  if (doMuptrelFit) {
    hBEfficiencyDataPtrel->SetLineColor(1);
    hBEfficiencyDataPtrel->SetMarkerColor(1);
    hBEfficiencyDataPtrel->SetMarkerStyle(20);
    hBEfficiencyDataPtrel->Draw("same");  
  } 
  if (doLTJP) {
    hBEfficiencyDataLTJP->SetLineColor(8);
    hBEfficiencyDataLTJP->SetMarkerColor(8);
    hBEfficiencyDataLTJP->SetMarkerStyle(20);
    hBEfficiencyDataLTJP->Draw("same");
  }
  if (doLTCSV) {
    hBEfficiencyDataLTCSV->SetLineColor(7);
    hBEfficiencyDataLTCSV->SetMarkerColor(7);
    hBEfficiencyDataLTCSV->SetMarkerStyle(20);
    hBEfficiencyDataLTCSV->Draw("same");
  }
  legEff->Draw();
  

  TLegend *legFrac = new TLegend(0.35,0.15,0.9,0.35);
  legFrac->SetBorderSize(1);
  legFrac->SetFillColor(18);
  legFrac->AddEntry(hBFractionMC,"MC","pl");
  legFrac->AddEntry(hBFractionJPdirect,"direct fit of JP","pl");
  legFrac->AddEntry(hBFractionData,Form("%s at %.1f + pur. from SV mass + eff. from MC",taggerName,workingPoint),"pl");
  if (doLTJP) legFrac->AddEntry(hBFractionDataLTJP,Form("%s at %.1f + pur. from SV mass + eff. from LT (JP)",taggerName,workingPoint),"pl");
  if (doLTCSV) legFrac->AddEntry(hBFractionDataLTCSV,Form("%s at %.1f + pur. from SV mass + eff. from LT (CSV)",taggerName,workingPoint),"pl");
  TCanvas *cBFraction = new TCanvas("cBFraction","b-jet fraction",600,600);
  hBFractionMC->SetAxisRange(0,0.03,"Y");
  hBFractionMC->SetTitleOffset(1.8,"Y");
  hBFractionMC->SetLineColor(2);
  hBFractionMC->SetMarkerColor(2);
  hBFractionMC->SetMarkerStyle(21);
  hBFractionMC->Draw("e1"); 
  hBFractionData->SetLineColor(1);
  hBFractionData->SetMarkerColor(1);
  hBFractionData->SetMarkerStyle(20);
  hBFractionData->Draw("e1same");   
  if (doLTJP) {
    hBFractionDataLTJP->SetLineColor(8);
    hBFractionDataLTJP->SetMarkerColor(8);
    hBFractionDataLTJP->SetMarkerStyle(20);
    hBFractionDataLTJP->Draw("e1same");
  }
  if (doLTCSV) {
    hBFractionDataLTCSV->SetLineColor(7);
    hBFractionDataLTCSV->SetMarkerColor(7);
    hBFractionDataLTCSV->SetMarkerStyle(20);
    hBFractionDataLTCSV->Draw("e1same");
  }
  hBFractionJPdirect->SetLineColor(4);
  hBFractionJPdirect->SetMarkerColor(4);
  hBFractionJPdirect->SetMarkerStyle(20);
  hBFractionJPdirect->Draw("e1same");
  legFrac->Draw();

  cout<<" NOT WRITING OUTPUT "<<endl;

  /*  TFile *fout = new TFile(Form("histos/bFraction_%sat%.1f.root",taggerName,workingPoint),"update");
  hBFractionMC->Write();
  hBFractionData->Write();
  if (doLTJP) hBFractionDataLTJP->Write();
  if (doLTCSV) hBFractionDataLTCSV->Write();
  hBFractionJPdirect->Write();
  fout->Close();
*/
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




void fixEmpty(TH1 *h){
   for (Int_t i=1;i<=h->GetNbinsX();i++){
      if (h->GetBinContent(i)==0) h->SetBinContent(i,1e-20);
   }
}


RooRealVar bfractionFit(Int_t isRecopp=0, Int_t isMuTrig=0, Bool_t fixCL=1, Double_t CLshift=0, char *var="muptrel", Double_t minXvar=0, Double_t maxXvar=5, Double_t ptMin=60, Double_t ptMax=500, 
// by default, no b-tagging :
char *discr="discr_prob", Double_t minXdiscr=-1, Double_t maxXdiscr=3, char *comment="inclusive sample", 
Double_t maxYaxis=5e2)
{


  // discr_prob : from (0) 0 to 3, operating point : 0.5, 0.7 
  // discr_ssvHighEff : from (-1) 1 to 6, operating point : 2 ?
  // discr_ssvHighPur : from (-1) 1 to 6, operating point : 2 ?  
  // discr_csvSimple : from (-10,-1) 0 to 1, operating point  : 0.9  
  // svtxm : from (0) 0 to 7 
  // muptrel : from (0) 0 to 5

  TFile *fMC, *fdata;

  if        ( isRecopp&& isMuTrig) { // pp reco, muon triggered
    fMC = new TFile("bFractionpp290512/histos/ppMC_ppReco_muTrig.root");
    fdata = new TFile("bFractionpp290512/histos/ppdata_ppReco_muTrig.root");
  } else if ( isRecopp&&!isMuTrig) { // pp reco, jet triggered
    fMC = new TFile("bFractionpp290512/histos/ppMC_ppReco_jetTrig.root");
    fdata = new TFile("bFractionpp290512/histos/ppdata_ppReco_jetTrig.root");
  } else if (!isRecopp&& isMuTrig) { // hi reco, muon triggered
    fMC = new TFile("bFractionpp290512/histos/ppMC_hiReco_muTrig.root");
    fdata = new TFile("bFractionpp290512/histos/ppdata_hiReco_muTrig.root");
  } else if (!isRecopp&&!isMuTrig) { // hi reco, jet triggered
    fMC = new TFile("bFractionpp290512/histos/ppMC_hiReco_jetTrig.root");
    fdata = new TFile("bFractionpp290512/histos/ppdata_hiReco_jetTrig.root");

  }

  TTree *tMC, *tdata;
  if (var=="muptrel") { 
    //* muon requirement
    tMC = (TTree*) fMC->Get("ntMuReq");
    tdata = (TTree*) fdata->Get("ntMuReq");
    //*/
    /* muon + away jet requirement
    tMC = (TTree*) fMC->Get("ntAwayjetReq");
    tdata = (TTree*) fdata->Get("ntAwayjetReq");
    //*/
  } else {
    tMC = (TTree*) fMC->Get("nt");
    tdata = (TTree*) fdata->Get("nt");
  } 

  //TCanvas *c = new TCanvas("c","",600,600);
  
  TH1D *hB = new TH1D("hB","",50,minXvar,maxXvar);
  hB->Sumw2();
  tMC->Draw(Form("%s>>hB",var),Form("weight*(abs(refparton_flavorForB)==5&&jtpt>=%f&&jtpt<%f&&%s>%f&&%s<%f)",ptMin,ptMax,discr,minXdiscr,discr,maxXdiscr),"goff");
  fixEmpty(hB);
  
  TH1D *hC = new TH1D("hC","",50,minXvar,maxXvar);
  hC->Sumw2();
  tMC->Draw(Form("%s>>hC",var),Form("weight*(abs(refparton_flavorForB)==4&&jtpt>=%f&&jtpt<%f&&%s>%f&&%s<%f)",ptMin,ptMax,discr,minXdiscr,discr,maxXdiscr),"goff");
  fixEmpty(hC);

  TH1D *hL = new TH1D("hL","",50,minXvar,maxXvar);
  hL->Sumw2();
  tMC->Draw(Form("%s>>hL",var),Form("weight*(abs(refparton_flavorForB)!=5&&abs(refparton_flavorForB)!=4&&abs(refparton_flavorForB)<99&&jtpt>=%f&&jtpt<%f&&%s>%f&&%s<%f)",ptMin,ptMax,discr,minXdiscr,discr,maxXdiscr),"goff");
  fixEmpty(hL);

  /*
  TH1D *hCL = new TH1D("hCL","",50,minXvar,maxXvar);
  hCL->Sumw2();
  tMC->Draw(Form("%s>>hCL",var),Form("weight*(abs(refparton_flavorForB)!=5&&abs(refparton_flavorForB)<99&&jtpt>=%f&&jtpt<%f&&%s>%f&&%s<%f)",ptMin,ptMax,discr,minXdiscr,discr,maxXdiscr));
  fixEmpty(hCL);
  //*/

  TH1D *hCL = hL->Clone();
  hCL->Add(hC,1+CLshift/100);

  //delete c;   

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
  RooDataSet *data = new RooDataSet("data","data",tdata,RooArgSet(s,jtpt,discriminator),Form("jtpt>=%f&&jtpt<%f&&%s>%f&&%s<%f",ptMin,ptMax,discr,minXdiscr,discr,maxXdiscr));
    

  TCanvas *cROOFIT = new TCanvas("cROOFIT",Form("Template fit of %s",var),200,10,1200,600);
  TPad* pad1 = new TPad("pad1","pad1",0,0,0.5,0.9);
  TPad* pad2 = new TPad("pad2","pad1",0.5,0,1,0.9);
  pad1->Draw();
  pad2->Draw();
  TPaveText *header = new TPaveText(0.05,0.9,0.95,0.99);
  header->AddText(Form("%s  -  ROOFIT ML unbinned fit of %s",var,fixCL?"2 templates : bottom and (charm + light)":"3 templates : bottom, charm and light"));
  header->AddText(Form("%s - %s - %s",isRecopp?"pp reco":"HI reco",isMuTrig?"Mu3 triggered":"Jet40 triggered",comment));
  header->AddText(Form("%s%.0f <= jet pT < %.0f",(var=="muptrel")?"deltaR < 0.5 ; muon pT > 5 ; ":"",ptMin,ptMax));
  header->SetTextSize(0.027);
  header->SetTextAlign(12);
  header->SetBorderSize(0);
  header->SetFillStyle(0);
  header->Draw();

  // --- Plot before fitting ---
  pad1->cd();
  pad1->SetLogy();
  RooPlot* sframe = s.frame();
  TH2D *htemp = new TH2D("htemp","",100,minXvar,maxXvar,100,0.5,maxYaxis) ;
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
//   TLegend *leg = new TLegend(0.61,fixCL?0.60:0.50,0.98,fixCL?0.78:0.75/*,comment*/);
//   leg->SetBorderSize(0);
//   leg->SetFillStyle(0);
//   leg->AddEntry("h_data","pp @ 2.76 TeV","p");
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
   
  // --- Perform extended ML fit of composite PDF to data ---
  RooFitResult *fitresult = model.fitTo(*data,Save(),PrintLevel(-1));
  // use this to get detailed results of fits
  //RooFitResult *fitresult = model.fitTo(*data,Save());
  
  // --- Plot after fitting ---
  pad2->cd();
  pad2->SetLogy();
  sframe = s.frame();
  htemp->Draw();
  data->plotOn(sframe,Name("data"),Binning(50));
  if(fixCL) {
    model.plotOn(sframe,Components(charmlight),LineStyle(kDashed),LineColor(30),LineWidth(2));
  } else {
    model.plotOn(sframe,Components(light),LineStyle(kDashed),LineColor(kBlue),LineWidth(2));
    model.plotOn(sframe,Components(charm),LineStyle(kDashed),LineColor(kGreen),LineWidth(2));
    model.plotOn(sframe,Components(RooArgSet(charm,light)),LineStyle(kDashed),LineColor(30),LineWidth(2));
  }
  model.plotOn(sframe,Components(bottom),LineStyle(kDashed),LineColor(kRed),LineWidth(2),FillColor(kRed),FillStyle(1));   
  model.plotOn(sframe,Name("model"),LineWidth(2),VisualizeError(*fitresult),FillColor(17));
  model.plotOn(sframe,LineWidth(2),LineColor(13));
  data->plotOn(sframe,Binning(50));
  model.paramOn(sframe,Layout(0.4,0.9,0.9),Format("NEU",FixedPrecision(3)));
  sframe->Draw("same");
  //leg->Draw("same");
  Double_t minNll = fitresult->minNll();
  cout<<" MIN NLL "<<minNll<<endl;

  // --- Get chi2/dof ---
//   Double_t chi2 = sframe->chiSquare(2);
//   cout<<" chi2/dof = "<<chi2<<endl;
//   cROOFIT->cd();
//   TPaveText *PaveChi2 = new TPaveText(0.75,0.90,0.85,0.94);
//   PaveChi2->AddText(Form("chi2/dof = %.3f",chi2));
//   PaveChi2->SetBorderSize(0);
//   PaveChi2->SetFillStyle(0);
//   PaveChi2->Draw();



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
  tMC->Draw(Form("%s>>hMCB_%d",var,counter),Form("weight*(abs(refparton_flavorForB)==5&&jtpt>=%f&&jtpt<%f&&%s>=%f&&%s<%f)",ptMin,ptMax,discr,minXdiscr,discr,maxXdiscr),"goff");
  tMC->Draw(Form("%s>>hMCC_%d",var,counter),Form("weight*(abs(refparton_flavorForB)==4&&jtpt>=%f&&jtpt<%f&&%s>=%f&&%s<%f)",ptMin,ptMax,discr,minXdiscr,discr,maxXdiscr),"goff");
  tMC->Draw(Form("%s>>hMCL_%d",var,counter),Form("weight*(abs(refparton_flavorForB)!=5&&abs(refparton_flavorForB)!=4&&abs(refparton_flavorForB)<99&&jtpt>=%f&&jtpt<%f&&%s>=%f&&%s<%f)",ptMin,ptMax,discr,minXdiscr,discr,maxXdiscr),"goff");
  tdata->Draw(Form("%s>>hData_%d",var,counter),Form("jtpt>=%f&&jtpt<%f&&%s>=%f&&%s<%f",ptMin,ptMax,discr,minXdiscr,discr,maxXdiscr),"goff");
  hMCLC[counter]->Add( hMCL[counter]);
  hMCLC[counter]->Add( hMCC[counter]);
  fixEmpty(hMCB[counter]); fixEmpty(hMCC[counter]); fixEmpty(hMCL[counter]); fixEmpty(hMCLC[counter]); fixEmpty(hData[counter]);
  
  
  can1[counter] = new TCanvas(Form("can1_%d",counter),Form("can1_%d",counter),700,600);
  hs[counter] = new THStack(Form("hs_%d",counter),"le stack of MC histos");
  fakehs[counter] = new THStack(Form("fakehs_%d",counter),"le fake stack of MC histos");

  can1[counter]->cd();
  if (doLog) can1[counter]->cd()->SetLogy();


  if (doLog){
    hData[counter]->SetMaximum(hData[counter]->GetMaximum()*50);
    hData[counter]->SetMinimum(0.5);
  }
  if (!doLog){
    hData[counter]->SetMaximum(hData[counter]->GetMaximum()*1.5);
    hData[counter]->SetMinimum(0.0);
  }
  const char* yTitle;
  const char* xTitle;
  if (var=="svtxm")  xTitle = Form("Secondary Vertex Mass [GeV/c^{2}]");
  if (var=="discr_prob") xTitle = Form("Jet Probability");
  //if (var=="svtxm") yTitle = (Form("Number of Jets / %3.2f GeV",(maxXvar-minXvar)/nXbins));
  //if (var=="discr_prob") yTitle = (Form("Numbers of Jets / %3.2f",(maxXvar-minXvar)/nXbins));
  if (var=="svtxm") yTitle = ("Number of Jets");
  if (var=="discr_prob") yTitle = ("Number of Jets");
  hData[counter]->GetXaxis()->CenterTitle();
  hData[counter]->GetYaxis()->CenterTitle();
  hData[counter]->GetYaxis()->SetTitle(yTitle);
  hData[counter]->GetXaxis()->SetTitle(xTitle);
  hData[counter]->Draw();
  double Bnorm, Cnorm, Lnorm, LCnorm;
  Double_t Bfrac =Bfraction.getVal();
  Double_t Cfrac =Cfraction.getVal();
  
  //Normalize Histograms
  //if(!fixCL) Lnorm = (hData[counter]->Integral(1,nXbins)/hMCL[counter]->Integral(1,nXbins))*(1/(Cfrac+Bfrac+1));
  //if(!fixCL) Cnorm = (hData[counter]->Integral(1,nXbins)/hMCC[counter]->Integral(1,nXbins))*(Cfrac/(Cfrac+Bfrac+1));
  //if(!fixCL) Bnorm = (hData[counter]->Integral(1,nXbins)/hMCB[counter]->Integral(1,nXbins))*(Bfrac/(Cfrac+Bfrac+1));
  //if(fixCL) LCnorm = (hData[counter]->Integral(1,nXbins)/hMCLC[counter]->Integral(1,nXbins))*(1/(Bfrac+1));
  //if(fixCL)  Bnorm = (hData[counter]->Integral(1,nXbins)/hMCB[counter]->Integral(1,nXbins))*(Bfrac/(Bfrac+1));
  if(!fixCL) Lnorm = (hData[counter]->Integral(1,nXbins)/hMCL[counter]->Integral(1,nXbins))*(1.-Cfrac-Bfrac);
  if(!fixCL) Cnorm = (hData[counter]->Integral(1,nXbins)/hMCC[counter]->Integral(1,nXbins))*(Cfrac);
  if(!fixCL) Bnorm = (hData[counter]->Integral(1,nXbins)/hMCB[counter]->Integral(1,nXbins))*(Bfrac);
  if(fixCL) LCnorm = (hData[counter]->Integral(1,nXbins)/hMCLC[counter]->Integral(1,nXbins))*(1.-Bfrac);
  if(fixCL)  Bnorm = (hData[counter]->Integral(1,nXbins)/hMCB[counter]->Integral(1,nXbins))*(Bfrac);
  
  hMCB[counter]->SetFillColor(kRed+2);
  //hMCB[counter]->SetLineWidth(4);
  hMCL[counter]->SetFillColor(kBlue+1);
  hMCC[counter]->SetFillColor(kGreen+2);
  hMCLC[counter]->SetFillColor(kBlue+2);
  hMCB[counter]->SetMarkerSize(0);
  hMCC[counter]->SetMarkerSize(0);
  hMCL[counter]->SetMarkerSize(0);
  hMCLC[counter]->SetMarkerSize(0);
  hMCB[counter]->Scale(Bnorm);
  hMCC[counter]->Scale(Cnorm);
  hMCL[counter]->Scale(Lnorm);
  hMCLC[counter]->Scale(LCnorm);
  if (!fixCL){
    hs[counter]->Add(hMCB[counter]);
    hs[counter]->Add(hMCC[counter]);
    hs[counter]->Add(hMCL[counter]);
  }
  if (fixCL){
    hs[counter]->Add(hMCB[counter]);
    hs[counter]->Add(hMCLC[counter]);
  }

  //hs[counter]->Draw("same hE2");
  hs[counter]->Draw("same h e");
  hData[counter]->Draw("same");
  can1[counter]->GetFrame()->SetLineWidth(4);
  can1[counter]->RedrawAxis();
  //checkBins(hMCB[counter]); checkBins(hMCC[counter]); checkBins(hMCL[counter]); checkBins(hData[counter]);

  //To obtain the ovaerll sum histogram of the MC counts
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
  MCTotal[counter]->SetLineWidth(3.0);
  MCTotal[counter]->SetMarkerSize(0);
  MCTotal[counter]->SetMarkerColor(kGray+2);
  //MCTotal[counter]->SetLineColor(kAzure-3);    
  MCTotal[counter]->SetLineColor(kBlue-9);    
  //MCTotal[counter]->SetLineColor(kGray+2);
  MCTotal[counter]->Draw("same e");
  hData[counter]->Draw("same");
  //http://root.cern.ch/root/htmldoc/TH1.html#TH1:Chi2Test
  Double_t chi2 = hData[counter]->Chi2Test(MCTotal[counter],"UW CHI2 P");
  Double_t chi2NDF = hData[counter]->Chi2Test(MCTotal[counter],"UW CHI2/NDF P");

  //This is a fake THstack to also plot the constribution from the charm, 
  //even though it is supposed to be merged with the Light constribution 
  if(fixCL){
    //fakehs[counter]->Add(hMCB[counter]);
    //hMCC[counter]->Scale();
    
  }

  //Redraw some partial histograms to show error bars in between
  if(fixCL){
    //hMCB[counter]->SetLineColor(kBlack);
    //hMCB[counter]->SetLineColor(kRed);
    hMCB[counter]->SetMarkerSize(0);
    hMCB[counter]->Draw("same e");
  }
  /*  // Jorge, I don't know what you're doing here, but there's a bug.  The error bars are drawn in the wrong place
  if(!fixCL){  
   hMCB[counter]->SetLineColor(kBlack);
   hMCB[counter]->SetMarkerSize(0);
   hMCB[counter]->Draw("same e");
   hMCB[counter]->Add(hMCC[counter]);
   //hMCB[counter]->SetLineColor(kRed);
   //hMCB[counter]->SetLineColor(kBlack);
   hMCB[counter]->SetMarkerSize(0);
   hMCB[counter]->Draw("same e");
  }
  */
  TLegend *hleg = new TLegend(0.5,fixCL?0.70:0.67,0.90,fixCL?0.92:0.92);
  hleg->SetBorderSize(0);
  hleg->SetFillStyle(0);
  hleg->AddEntry(hData[counter],"pp data","lp");
  hleg->AddEntry(hMCB[counter],"b","f");
  if(!fixCL){
    hleg->AddEntry(hMCC[counter],"c","f");
    hleg->AddEntry(hMCL[counter],"udsg","f");
  }
  if(fixCL){
    hleg->AddEntry(hMCLC[counter],"c + usdg","f");
  }
  hleg->Draw("same");
  if(!fixCL)drawText(Form("%2.0f < p_{T} < %2.0f GeV/c",ptMin,ptMax),0.51,0.62);
  if(fixCL)drawText(Form("%2.0f < p_{T} < %2.0f GeV/c",ptMin,ptMax),0.51,0.66);
  //drawText(Form("#chi^{2}/NDF = %3.1f",chi2NDF),0.51,0.55);
  //drawText(Form("#chi^{2}/NDF = %3.1f / %d",chi2NDF,nXbins-1),0.51,0.55);
  drawText(Form("#chi^{2}/NDF = %3.1f / %d",chi2,(int)(chi2/chi2NDF)),0.51,0.55);
  //drawText(Form("min Nll = %3.1f",minNll),0.51,0.45);
  drawText("CMS Preliminary",0.15,0.965);
  drawText("#sqrt{s_{NN}} = 2.76 TeV",0.60,0.965);
  drawText("|#eta| < 2.0",0.18,0.88);
  if(comment!="b-tagged sample (SSVHE > 2)") drawText(comment,0.18,0.80);
  if(comment=="b-tagged sample (SSVHE > 2)"){
    drawText("b-tagged sample",0.18,0.80);
    drawText("(SSVHE > 2)",0.18,0.75);
  }
  
  TString ptMinLabel;
  if(ptMin==80) ptMinLabel = Form("0%2.0f",ptMin);
  else ptMinLabel = Form("%3.0f",ptMin);
  
  bool printEach=true;
  // --- Print results ---
  //cout <<"b jet fraction in MC = "<<bInitFrac<<endl;
  cout<<"ZZZZZ TEXTcounter: "<<counter<<"   Comment: "<<comment<<endl;
  if(!fixCL) cout<<"ZZZZZ Scale factors: B= "<<Bnorm<<" C= "<<Cnorm<<" Light= "<<Lnorm<<endl;
  if(fixCL)  cout<<"ZZZZZ Scale factors: B= "<<Bnorm<<" LC= "<<LCnorm<<endl;
  cout<<"ZZZZZ Chi2: "<<chi2<<" Chi2/NDF: "<<chi2NDF<<endl;
  cout<<"ZZZZZ var: "<<var<<"  discr: "<<discr<<"    min(discr): "<<minXvar<<"    max(discr): "<<maxXvar<<endl;
  cout<<"ZZZZZ in pT [ "<<ptMin<<" , "<<ptMax<<" ]"<<endl;
  cout<<"ZZZZZ b jet fraction in data = "<<Bfraction.getVal()<<endl;
  if(!fixCL) cout <<"ZZZZZ c jet fraction = "<<Cfraction.getVal()<<endl;
  if(!fixCL){
    if (!printEach){
      if (counter<11) can1[counter]->Print("bTagStackedHistos_nofixCL.pdf(","pdf");
      if (counter==11) can1[counter]->Print("bTagStackedHistos_nofixCL.pdf)","pdf");
    }
    if (printEach) can1[counter]->Print(Form("PDFS/bStack_%sPt%s_%3.0f_nofixCL.pdf",var,ptMinLabel,ptMax),"pdf");
  }
  if(fixCL){
    if (!printEach){
      if (counter<11) can1[counter]->Print("bTagStackedHistos_fixCL.pdf(","pdf");
      if (counter==11) can1[counter]->Print("bTagStackedHistos_fixCL.pdf)","pdf");
    }
    if (printEach) can1[counter]->Print(Form("PDFS/bStack_%sPt%s_%3.0f_fixCL.pdf",var,ptMinLabel,ptMax),"pdf");
  }
  counter++;
  




  // --- Print results ---
  //cout <<"b jet fraction in MC = "<<bInitFrac<<endl;
  cout <<"b jet fraction in data = "<<Bfraction.getVal()<<endl;
  if(!fixCL) cout <<"c jet fraction = "<<Cfraction.getVal()<<endl;

  // --- Save canvas ---
  TString path = "";
  path.Append("gifs_").Append(isRecopp?"ppReco":"hiReco").Append("_").Append(isMuTrig?"muTrig":"jetTrig").Append(Form("/%s_jtpt%.0fto%.0f_%s%.2fto%.2f_",var,ptMin,ptMax,discr,minXdiscr,maxXdiscr)).Append(fixCL?"CLfixed":"CLfree").Append(".gif");
  cROOFIT->SaveAs(path);

  return Bfraction;
}
void drawText(const char *text, float xp, float yp){
  TLatex *tex = new TLatex(xp,yp,text);
  tex->SetTextFont(63);
  tex->SetTextSize(25);
  //tex->SetTextSize(0.05);
  tex->SetTextColor(kBlack);
  tex->SetLineWidth(1);
  tex->SetNDC();
  tex->Draw();
}




Enumerations count(Int_t isRecopp, Int_t isMuTrig, Double_t ptMin, Double_t ptMax, char *discr="discr_prob", Double_t workingPoint) {

  TFile *fMC, *fdata;

  if        ( isRecopp&& isMuTrig) { // pp reco, muon triggered
    fMC = new TFile("bFractionpp290512/histos/ppMC_ppReco_muTrig.root");
    fdata = new TFile("bFractionpp290512/histos/ppdata_ppReco_muTrig.root");
  } else if ( isRecopp&&!isMuTrig) { // pp reco, jet triggered
    fMC = new TFile("bFractionpp290512/histos/ppMC_ppReco_jetTrig.root");
    fdata = new TFile("bFractionpp290512/histos/ppdata_ppReco_jetTrig.root");
  } else if (!isRecopp&& isMuTrig) { // hi reco, muon triggered
    fMC = new TFile("bFractionpp290512/histos/ppMC_hiReco_muTrig.root");
    fdata = new TFile("bFractionpp290512/histos/ppdata_hiReco_muTrig.root");
  } else if (!isRecopp&&!isMuTrig) { // hi reco, jet triggered
    fMC = new TFile("bFractionpp290512/histos/ppMC_hiReco_jetTrig.root");
    fdata = new TFile("bFractionpp290512/histos/ppdata_hiReco_jetTrig.root");
    //fMC = new TFile("/home/llr/cms/mnguyen/bTagging442p5/CMSSW_4_4_2_patch5/src/bTaggingMacros/histos/ppMC_hiReco_jetTrig.root");
    //fdata = new TFile("/home/llr/cms/mnguyen/bTagging442p5/CMSSW_4_4_2_patch5/src/bTaggingMacros/histos/ppdata_hiReco_jetTrig.root");
  }

  TTree *tMC, *tdata;
  if (isMuTrig) { 
    //* muon requirement
    tMC = (TTree*) fMC->Get("ntMuReq");
    tdata = (TTree*) fdata->Get("ntMuReq");
    //*/
    /* muon + away jet requirement
    tMC = (TTree*) fMC->Get("ntAwayjetReq");
    tdata = (TTree*) fdata->Get("ntAwayjetReq");
    //*/
  } else {
    tMC = (TTree*) fMC->Get("nt");
    tdata = (TTree*) fdata->Get("nt");
  } 

  //TCanvas *c = new TCanvas("c","",600,600);
  //c->cd();

  TH1D *hTaggedJetsMC = new TH1D("hTaggedJetsMC","hTaggedJetsMC",1,ptMin,ptMax);
  hTaggedJetsMC->Sumw2();
  tMC->Draw("jtpt>>hTaggedJetsMC",Form("weight*(%s>=%f)",discr,workingPoint),"goff");
  cout<<"entries of hTaggedJetsMC"<<hTaggedJetsMC->GetEntries()<<endl;

  TH1D *hUntaggedJetsMC = new TH1D("hUntaggedJetsMC","",1,ptMin,ptMax);
  hUntaggedJetsMC->Sumw2();
  tMC->Draw("jtpt>>hUntaggedJetsMC",Form("weight*(%s<%f)",discr,workingPoint),"goff");
  
  TH1D *hTaggedBjetsMC = new TH1D("hTaggedBjetsMC","",1,ptMin,ptMax);
  hTaggedBjetsMC->Sumw2();
  tMC->Draw("jtpt>>hTaggedBjetsMC",Form("weight*(abs(refparton_flavorForB)==5&&%s>=%f)",discr,workingPoint),"goff");  

  TH1D *hUntaggedBjetsMC = new TH1D("hUntaggedBjetsMC","",1,ptMin,ptMax);
  hUntaggedBjetsMC->Sumw2();
  tMC->Draw("jtpt>>hUntaggedBjetsMC",Form("weight*(abs(refparton_flavorForB)==5&&%s<%f)",discr,workingPoint),"goff");

  TH1D *hBjetsWithJPinfoMC = new TH1D("hBjetsWithJPinfoMC","",1,ptMin,ptMax);
  hBjetsWithJPinfoMC->Sumw2();
  tMC->Draw("jtpt>>hBjetsWithJPinfoMC","weight*(abs(refparton_flavorForB)==5&&discr_prob>0)","goff");

  TH1D *hBjetsWithCSVinfoMC = new TH1D("hBjetsWithCSVinfoMC","",1,ptMin,ptMax);
  hBjetsWithCSVinfoMC->Sumw2();
  tMC->Draw("jtpt>>hBjetsWithCSVinfoMC","weight*(abs(refparton_flavorForB)==5&&discr_prob>0)","goff");
 
  TH1D *hTaggedJetsData = new TH1D("hTaggedJetsData","",1,ptMin,ptMax);
  hTaggedJetsData->Sumw2();
  tdata->Draw("jtpt>>hTaggedJetsData",Form("weight*(%s>=%f)",discr,workingPoint),"goff");

  TH1D *hUntaggedJetsData = new TH1D("hUntaggedJetsData","",1,ptMin,ptMax);
  hUntaggedJetsData->Sumw2();
  tdata->Draw("jtpt>>hUntaggedJetsData",Form("weight*(%s<%f)",discr,workingPoint),"goff");

  Enumerations res;

  res.nTaggedJetsMC = hTaggedJetsMC->GetBinContent(1);
  res.nUntaggedJetsMC = hUntaggedJetsMC->GetBinContent(1);
  res.nJetsMC = res.nTaggedJetsMC + res.nUntaggedJetsMC;
  res.nTaggedBjetsMC = hTaggedBjetsMC->GetBinContent(1);
  res.nUntaggedBjetsMC = hUntaggedBjetsMC->GetBinContent(1);
  res.nBjetsMC = res.nTaggedBjetsMC + res.nUntaggedBjetsMC;
  res.nNonBjetsMC = res.nJetsMC - res.nBjetsMC;
  res.nTaggedNonBjetsMC = res.nTaggedJetsMC - res.nTaggedBjetsMC;

  res.nTaggedJetsData = hTaggedJetsData->GetBinContent(1);
  res.nUntaggedJetsData = hUntaggedJetsData->GetBinContent(1);

  res.cbForJP = hBjetsWithJPinfoMC->GetBinContent(1) / (res.nTaggedBjetsMC+res.nUntaggedBjetsMC);
  res.cbForCSV = hBjetsWithCSVinfoMC->GetBinContent(1) / (res.nTaggedBjetsMC+res.nUntaggedBjetsMC);

  res.nTaggedJetsMCError = hTaggedJetsMC->GetBinError(1);
  res.nUntaggedJetsMCError = hUntaggedJetsMC->GetBinError(1);
  res.nJetsMCError = addError(res.nTaggedJetsMCError,res.nUntaggedJetsMCError);
  res.nTaggedBjetsMCError = hTaggedBjetsMC->GetBinError(1);
  res.nUntaggedBjetsMCError = hUntaggedBjetsMC->GetBinError(1);
  res.nBjetsMCError = addError(res.nTaggedBjetsMCError,res.nUntaggedBjetsMCError);
  res.nNonBjetsMCError = substractError(res.nJetsMCError,res.nBjetsMCError);
  res.nTaggedNonBjetsMCError = substractError(res.nTaggedJetsMCError,res.nTaggedBjetsMCError);

  res.nTaggedJetsDataError = hTaggedJetsData->GetBinError(1);
  res.nUntaggedJetsDataError = hUntaggedJetsData->GetBinError(1);

  return res;
 
}




