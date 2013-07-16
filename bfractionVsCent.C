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




void bfractionVsCent(char *tagger="discr_ssvHighEff", double workingPoint=2., int fixCL=0, char *taggerName="ssvHighEff", float ptlo=80, float pthi=100, float etalo=0., float etahi=2.) {

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetLabelFont(43,"xyz");
  gStyle->SetLabelSize(20,"xyz");
  gStyle->SetTitleFont(43,"xyz");
  gStyle->SetTitleSize(26,"xyz");
  gStyle->SetTitleOffset(1.0,"xy"); 
  gROOT->ForceStyle(1);
  
  int doLTJP=1;
  int doLTCSV=0;

  const int nBins = 2;
  double centBin[nBins+1] = {0,30,100};
  //const int nBins = 3;
  //double centBin[nBins+1] = {0,20,50,100};
 
  
  Double_t bPurMC, bPurData, bEffMC, bEffDataLTJP, bEffDataLTCSV, taggedFracData, bFracMC, bFracData, bFracDataLTJP, bFracDataLTCSV, bFracJPdirect;
  Double_t bPurMCError, bPurDataError, bEffMCError, bEffDataLTJPError, bEffDataLTCSVError, taggedFracDataError, bFracMCError, bFracDataError, bFracDataLTJPError, bFracDataLTCSVError, bFracJPdirectError;
  Enumerations numbers;
  
  TH1D *hBPurityData = new TH1D("hBPurityData","hBPurityData;Centrality;B-Tagging purity",nBins,centBin);
  TH1D *hBPurityMC = new TH1D("hBPurityMC","hBPurityMC;Centrality;B-Tagging purity",nBins,centBin);
  
  TH1D *hBEfficiencyMC = new TH1D("hBEfficiencyMC","hBEfficiencyMC;Centrality;B-Tagging efficiency",nBins,centBin);
  TH1D *hBEfficiencyDataLTJP = new TH1D("hBEfficiencyDataLTJP","hBEfficiencyDataLTJP;Centrality;B-Tagging efficiency",nBins,centBin);
  TH1D *hBEfficiencyDataLTCSV = new TH1D("hBEfficiencyDataLTCSV","hBEfficiencyDataLTCSV;Centrality;B-Tagging efficiency",nBins,centBin);
  
  TH1D *hBFractionMC = new TH1D("hBFractionMC","hBFractionMC;Centrality;B-jet fraction",nBins,centBin);
  TH1D *hBFractionData = new TH1D("hBFractionData","hBFractionData;Centrality;B-jet fraction",nBins,centBin);
  TH1D *hBFractionDataLTJP = new TH1D("hBFractionDataLTJP","hBFractionDataLTJP;Centrality;B-jet fraction",nBins,centBin);
  TH1D *hBFractionDataLTCSV = new TH1D("hBFractionDataLTCSV","hBFractionDataLTCSV;Centrality;B-jet fraction",nBins,centBin);
  TH1D *hBFractionJPdirect = new TH1D("hBFractionJPdirect","hBFractionJPdirect;Centrality;B-jet fraction",nBins,centBin);
  
  int ncol=1;
  int nrow=1;

  if(nBins==3||nBins==2){
    ncol=nBins;
  }
  if(nBins==4){
    ncol=nBins/2;
    nrow=nBins/2;
  }

  TCanvas *c1=new TCanvas("c1","c1",1200,600);
  //c1->Divide(ncol,nrow,0,0);
  c1->Divide(ncol,nrow);
  TCanvas *c2=new TCanvas("c2","c2",1200,600);
  //c2->Divide(ncol,nrow,0,0);
  c2->Divide(ncol,nrow);
  TCanvas *c3=new TCanvas("c3","c3",1200,600);
  //c3->Divide(ncol,nrow,0,0);
  c3->Divide(ncol,nrow);
  TCanvas *c4=new TCanvas("c4","c4",1200,600);
  //c4->Divide(ncol,nrow,0,0);
  c4->Divide(ncol,nrow);

  TCanvas *cCount = new TCanvas("cCount","cCount",600,600);

  for (int n=0;n<nBins;n++) {

    cout<<"Processing jet centrality bin ["<<centBin[n]<<","<<centBin[n+1]<<"] ..."<<endl;
    cCount->cd();
    cout<<"centBin[n]: "<<centBin[n]<<" centBin[n+1]: "<<centBin[n+1]<<" tagger: "<<tagger<<" workingPoint: "<<workingPoint<<" ptlo: "<<ptlo<<" pthi: "<<pthi<<" etalo: "<<etalo<<" etahi: "<<etahi<<endl;
    numbers = count(centBin[n],centBin[n+1],tagger,workingPoint,ptlo,pthi,etalo,etahi);
    
    c1->cd(n+1);
    c1->GetPad(n+1)->SetLogy();
    RooRealVar fitSvtxmTag = bfractionFit(fixCL,"svtxm",0,6,centBin[n],centBin[n+1],ptlo,pthi,etalo,etahi,tagger,workingPoint,6,"b-tagged sample (SSVHE > 2)",9e3);
    //RooRealVar fitSvtxmTag = bfractionFit(fixCL,"svtxm",0,6,centBin[n],centBin[n+1],ptlo,pthi,etalo,etahi,tagger,workingPoint,10,Form("b-tagged sample (%s at %.1f)",taggerName,workingPoint));
    //RooRealVar fitJpDirect = bfractionFit(fixCL,"discr_prob",0,3,centBin[n],centBin[n+1],ptlo,pthi,etalo,etahi,tagger,-2,10,"inclusive sample",5e4);
    c2->cd(n+1);
    c2->GetPad(n+1)->SetLogy();
    RooRealVar fitJpDirect = bfractionFit(fixCL,"discr_prob",0.0,3.,centBin[n],centBin[n+1],ptlo,pthi,etalo,etahi,"discr_prob",0.,3.,"inclusive sample",4e5);

    if (doLTJP) {
      c3->cd(n+1);
      c3->GetPad(n+1)->SetLogy();
      RooRealVar fitJpBeforetag = bfractionFit(fixCL,"discr_prob",0.0,3.,centBin[n],centBin[n+1],ptlo,pthi,etalo,etahi,"discr_prob",0,3.,"jets with JP info",4e5);
      c4->cd(n+1);
      c4->GetPad(n+1)->SetLogy();
      RooRealVar fitJpTag = bfractionFit(fixCL,"discr_prob",0.0,3.,centBin[n],centBin[n+1],ptlo,pthi,etalo,etahi,tagger,workingPoint,6,"b-tagged sample (SSVHE > 2)",4e5);
    } 
    if (doLTCSV) {
      RooRealVar fitCsvBeforetag = bfractionFit(fixCL,"discr_csvSimple",0,1,centBin[n],centBin[n+1],ptlo,pthi,etalo,etahi,tagger,-2,10,"jets with CSV info",4e5);
      RooRealVar fitCsvTag = bfractionFit(fixCL,"discr_csvSimple",0,1,centBin[n],centBin[n+1],ptlo,pthi,etalo,etahi,tagger,workingPoint,10,"b-tagged sample (SSVHE > 2)",4e5);
    } 

    taggedFracData = numbers.nTaggedJetsData / (numbers.nTaggedJetsData+numbers.nUntaggedJetsData);
    taggedFracDataError = fracError(numbers.nTaggedJetsData,numbers.nUntaggedJetsData,numbers.nTaggedJetsDataError,numbers.nUntaggedJetsDataError);
    
    //*  --- b-tagging purity --- 

    bPurMC = numbers.nTaggedBjetsMC / numbers.nTaggedJetsMC;
    cout<<" bPurMC "<<bPurMC<<" numbers.nTaggedBjetsMC "<<numbers.nTaggedBjetsMC<<" numbers.nTaggedJetsMC "<<numbers.nTaggedJetsMC<<endl;
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
    
    //*  --- b fraction --- 

    bFracMC = numbers.nBjetsMC / numbers.nJetsMC;
    //bFracMC = numbers.nTaggedJetsMC * bPurMC / (bEffMC * numbers.nJetsMC); // for check : same as previous
    bFracMCError = fracError(numbers.nBjetsMC,numbers.nNonBjetsMC,numbers.nBjetsMCError,numbers.nNonBjetsMCError); 
    hBFractionMC->SetBinContent(n+1,bFracMC); 
    hBFractionMC->SetBinError(n+1,bFracMCError); 


    bFracData = taggedFracData * bPurData / bEffMC; // efficiency from MC
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
    //*/
  }
  
  TLegend *legPur = new TLegend(0.4,0.15,0.85,0.3,Form("Purity of b-tagged sample (%s at %.1f)",taggerName,workingPoint));
  legPur->SetBorderSize(0);
  //legPur->SetFillColor(kGray);
  legPur->SetFillStyle(0);
  legPur->AddEntry(hBPurityMC,"MC Input","pl");
  legPur->AddEntry(hBPurityData,"Data","pl");
  TCanvas *cBPurity = new TCanvas("cBPurity","b purity",600,600);
  hBPurityMC->SetAxisRange(0,1,"Y");
  hBPurityMC->SetTitleOffset(1.3,"Y");
  hBPurityMC->SetLineColor(2);
  hBPurityMC->SetMarkerColor(2);
  hBPurityMC->SetMarkerStyle(21);
  hBPurityMC->Draw();
  hBPurityData->SetLineColor(1);
  hBPurityData->SetMarkerColor(1);
  hBPurityData->SetMarkerStyle(20);
  hBPurityData->Draw("same");   
  legPur->Draw();
  //cBPurity->SaveAs("purity.gif");

  TLegend *legEff = new TLegend(0.4,0.65,0.85,0.8,Form("Efficiency for tagging b-jets (%s at %.1f)",taggerName,workingPoint));
  legEff->SetBorderSize(0);
  //legEff->SetFillColor(kGray);
  legEff->SetFillStyle(0);
  legEff->AddEntry(hBEfficiencyMC,"MC Efficiency","pl");
  if (doLTJP) legEff->AddEntry(hBEfficiencyDataLTJP,"LT method (JP)","pl");
  if (doLTCSV) legEff->AddEntry(hBEfficiencyDataLTCSV,"LT method (CSV)","pl");
  TCanvas *cBEfficiency = new TCanvas("cBEfficiency","B-Tagging efficiency",600,600);
  hBEfficiencyMC->SetAxisRange(0,1,"Y");
  hBEfficiencyMC->SetTitleOffset(1.3,"Y");
  hBEfficiencyMC->SetLineColor(2);
  hBEfficiencyMC->SetMarkerColor(2);
  hBEfficiencyMC->SetMarkerStyle(21);
  hBEfficiencyMC->Draw();
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
  //cBEfficiency->SaveAs("efficiency.gif");


  TLegend *legFrac = new TLegend(0.25,0.15,0.85,0.3);
  legFrac->SetBorderSize(0);
  //legFrac->SetFillColor(kGray);
  legFrac->SetFillStyle(0);
  legFrac->AddEntry(hBFractionMC,"MC Input","pl");
  legFrac->AddEntry(hBFractionData,Form("%s at %.1f + pur. from SV mass + eff. from MC",taggerName,workingPoint),"pl");
  if (doLTJP) legFrac->AddEntry(hBFractionDataLTJP,Form("%s at %.1f + pur. from SV mass + eff. from LT (JP)",taggerName,workingPoint),"pl");
  if (doLTCSV) legFrac->AddEntry(hBFractionDataLTCSV,Form("%s at %.1f + pur. from SV mass + eff. from LT (CSV)",taggerName,workingPoint),"pl");
  legFrac->AddEntry(hBFractionJPdirect,"Direct fit to JP","pl");
  TCanvas *cBFraction = new TCanvas("cBFraction","B-jet fraction",600,600);
  hBFractionMC->SetAxisRange(0,0.03,"Y");
  hBFractionMC->SetTitleOffset(1.8,"Y");
  hBFractionMC->SetLineColor(2);
  hBFractionMC->SetMarkerColor(2);
  hBFractionMC->SetMarkerStyle(21);
  hBFractionMC->Draw(); 
  hBFractionData->SetLineColor(1);
  hBFractionData->SetMarkerColor(1);
  hBFractionData->SetMarkerStyle(20);
  hBFractionData->Draw("same");   
  if (doLTJP) {
    hBFractionDataLTJP->SetLineColor(8);
    hBFractionDataLTJP->SetMarkerColor(8);
    hBFractionDataLTJP->SetMarkerStyle(20);
    hBFractionDataLTJP->Draw("same");
  }
  if (doLTCSV) {
    hBFractionDataLTCSV->SetLineColor(7);
    hBFractionDataLTCSV->SetMarkerColor(7);
    hBFractionDataLTCSV->SetMarkerStyle(20);
    hBFractionDataLTCSV->Draw("same");
  }
  hBFractionJPdirect->SetLineColor(4);
  hBFractionJPdirect->SetMarkerColor(4);
  hBFractionJPdirect->SetMarkerStyle(20);
  hBFractionJPdirect->Draw("same");
  legFrac->Draw();
  //cBFraction->SaveAs("bfraction.gif");


  TFile *fout = new TFile(Form("output/bFractionMerged_%sat%.1fFixCL%d_centDep.root",taggerName,workingPoint,fixCL),"recreate");
  hBFractionMC->Write();
  hBFractionData->Write();
  if (doLTJP) hBFractionDataLTJP->Write();
  if (doLTCSV) hBFractionDataLTCSV->Write();
  hBFractionJPdirect->Write();
  fout->Close();

  //c1->SaveAs(Form("gifs/svtxMassFit_%s.gif",fixCL?"CLfixed":"CLfree"));
  //c2->SaveAs(Form("gifs/jpDirectFit_%s.gif",fixCL?"CLfixed":"CLfree"));
  //c3->SaveAs(Form("gifs/jpBeforeTag_%s.gif",fixCL?"CLfixed":"CLfree"));
  //c4->SaveAs(Form("gifs/jpAfterTag_%s.gif",fixCL?"CLfixed":"CLfree"));


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
   for (int i=1;i<=h->GetNbinsX();i++){
      if (h->GetBinContent(i)==0) h->SetBinContent(i,1e-20);
   }
}
//RooRealVar fitSvtxmTag = bfractionFit(fixCL,"svtxm",0,6,centBin[n],centBin[n+1],ptlo,pthi,etalo,etahi,tagger,workingPoint,6,"b-tagged sample (ssvHighEff at 2)",9e3);
RooRealVar bfractionFit(bool fixCL=1, char *var="discr_prob", double minXvar=0, double maxXvar=3, int centlo, int centhi, float ptMin, float ptMax, float etalo, float etahi,
// by default, no b-tagging :
char *discr="discr_prob", double minXdiscr=-999, double maxXdiscr=999, char *comment="inclusive sample", 
double maxYaxis=1e3)
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetLabelFont(43,"xyz");
  gStyle->SetLabelSize(20,"xyz");
  gStyle->SetTitleFont(43,"xyz");
  gStyle->SetTitleSize(26,"xyz");
  gStyle->SetTitleOffset(1.0,"x"); 
  gROOT->ForceStyle(1);


  int cbinlo = (int)centlo/2.5;
  int cbinhi = (int)centhi/2.5;

  // discr_prob : from (0) 0 to 3, operating point : 0.6 (1%), 0.7 
  // discr_ssvHighEff : from (-1) 1 to 6, operating point : 2 ?
  // discr_ssvHighPur : from (-1) 1 to 6, operating point : 2 ?  
  // discr_csvSimple : from (-10,-1) 0 to 1, operating point  : 0.9  
  // svtxm : from (0) 0 to 7 
  // muptrel : from (0) 0 to 5

  //*
  TFile *fQCDMC = new TFile("histos/PbPbQCDMC.root"); 
  TFile *fBMC = new TFile("histos/PbPbBMC.root"); 
  TFile *fCMC = new TFile("histos/PbPbCMC.root"); 
  TFile *fdata = new TFile("histos/PbPbdata.root");
  //*/

  /*
  TFile *fQCDMC = new TFile("histos/ppMC_hiReco_jetTrig.root"); 
  TFile *fBMC = new TFile("histos/ppMC_hiReco_jetTrig.root"); 
  TFile *fCMC = new TFile("histos/ppMC_hiReco_jetTrig.root"); 
  TFile *fdata = new TFile("histos/ppdata_hiReco_jetTrig.root");
  */

  TTree *tQCDMC = (TTree*) fQCDMC->Get("nt");
  TTree *tBMC = (TTree*) fBMC->Get("nt");
  TTree *tCMC = (TTree*) fCMC->Get("nt");
  TTree *tdata = (TTree*) fdata->Get("nt");

  
  int nhistBins=30;
  if(var=="svtxm") nhistBins=24;

  TH1D *hB = new TH1D("hB","",nhistBins,minXvar,maxXvar);
  hB->Sumw2();
  tBMC->Draw(Form("%s>>hB",var),Form("weight*(abs(refparton_flavorForB)==5&&jtpt>=%f&&jtpt<%f&&%s>=%f&&%s<%f&&bin>=%d&&bin<%d&&abs(jteta)>%f&&abs(jteta)<%f)",ptMin,ptMax,discr,minXdiscr,discr,maxXdiscr,cbinlo,cbinhi,etalo,etahi));
  fixEmpty(hB);
  
  TH1D *hC = new TH1D("hC","",nhistBins,minXvar,maxXvar);
  hC->Sumw2();
  tCMC->Draw(Form("%s>>hC",var),Form("weight*(abs(refparton_flavorForB)==4&&jtpt>=%f&&jtpt<%f&&%s>=%f&&%s<%f&&bin>=%d&&bin<%d&&abs(jteta)>%f&&abs(jteta)<%f)",ptMin,ptMax,discr,minXdiscr,discr,maxXdiscr,cbinlo,cbinhi,etalo,etahi));
  fixEmpty(hC);

  TH1D *hL = new TH1D("hL","",nhistBins,minXvar,maxXvar);
  hL->Sumw2();
  tQCDMC->Draw(Form("%s>>hL",var),Form("weight*(abs(refparton_flavorForB)!=5&&abs(refparton_flavorForB)!=4&&abs(refparton_flavorForB)<99&&jtpt>=%f&&jtpt<%f&&%s>=%f&&%s<%f&&bin>=%d&&bin<%d&&abs(jteta)>%f&&abs(jteta)<%f)",ptMin,ptMax,discr,minXdiscr,discr,maxXdiscr,cbinlo,cbinhi,etalo,etahi));
  fixEmpty(hL);
  /*
  TH1D *hCaux = new TH1D("hCaux","",nhistBins,minXvar,maxXvar);
  hCaux->Sumw2();
  tQCDMC->Draw(Form("%s>>hCaux",var),Form("weight*(abs(refparton_flavorForB)==4&&jtpt>=%f&&jtpt<%f&&%s>=%f&&%s<%f&&bin>=%d&&bin<%d&&abs(jteta)>%f&&abs(jteta)<%f)",ptMin,ptMax,discr,minXdiscr,discr,maxXdiscr,cbinlo,cbinhi,etalo,etahi));
  fixEmpty(hCaux);
  */
  /*
  TH1D *hCL = new TH1D("hCL","",nhistBins,minXvar,maxXvar);
  hCL->Sumw2();
  tQCDMC->Draw(Form("%s>>hCL",var),Form("weight*(abs(refparton_flavorForB)!=5&&abs(refparton_flavorForB)<99&&jtpt>=%f&&jtpt<%f&&%s>=%f&&%s<%f&&bin>=%d&&bin<%d&&abs(jteta)>%f&&abs(jteta)<%f)",ptMin,ptMax,discr,minXdiscr,discr,maxXdiscr,cbinlo,cbinhi,etalo,etahi));
  fixEmpty(hCL);
  //*/
  //*

  TH1D *hCL = hL->Clone();
  //Double_t cCoef = hCaux->Integral()/hC->Integral();
  hCL->Add(hC);

  //*/


  // --- Observable ---
  RooRealVar s(var,var,0,minXvar,maxXvar);
  RooRealVar jtpt("jtpt","jtpt",0,ptMin,ptMax);
  RooRealVar discriminator(discr,discr,0,minXdiscr,maxXdiscr);
  RooRealVar bin("bin","bin",0,0,40); 
  RooRealVar jteta("jteta","jteta",0,-2,2); 

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
  //Double_t bInitFrac = hB->Integral()/(hB->Integral()+hCL->Integral());
  //Double_t cInitFrac = hC->Integral()/(hB->Integral()+hCL->Integral());
  RooRealVar Bfraction("Bfraction","#light events",0.3,0.,1);
  RooRealVar Cfraction("Cfraction","#background events",0.3,0.,1); 
  if(fixCL) RooAddPdf model("model","",bottom,charmlight,Bfraction);
  else RooAddPdf model("model","",RooArgList(bottom,charm,light),RooArgList(Bfraction,Cfraction));  

  // --- Data sample ---
  //RooDataSet *data = new RooDataSet("data","data",tdata,RooArgSet(s,jtpt,discriminator),Form("jtpt>=%f&&jtpt<%f&&%s>=%f&&%s<%f",ptMin,ptMax,discr,minXdiscr,discr,maxXdiscr));
  //RooDataSet *data = new RooDataSet("data","data",tdata,RooArgSet(s,jtpt,discriminator),Form("jtpt>=%f&&jtpt<%f&&%s>=%f&&%s<%f&&bin>=%d&&bin<%d&&fabs(jteta)>%f&&fabs(jteta)<%f",ptMin,ptMax,discr,minXdiscr,discr,maxXdiscr,cbinlo,cbinhi,etalo,etahi));
  RooDataSet *data = new RooDataSet("data","data",tdata,RooArgSet(s,jtpt,jteta,bin,discriminator),Form("jtpt>=%f&&jtpt<%f&&%s>=%f&&%s<%f&&abs(jteta)>%f&&abs(jteta)<%f&&bin>=%d&&bin<%d",ptMin,ptMax,discr,minXdiscr,discr,maxXdiscr,etalo,etahi,cbinlo,cbinhi));

    


  TPaveText *header = new TPaveText(0.05,0.9,0.95,0.99);
  header->AddText(Form("%s  -  ROOFIT ML unbinned fit of %s",var,fixCL?"2 components : bottom and (charm + light)":"3 components : bottom, charm and light"));
  header->AddText(Form("Pb-Pb data - %s",comment));
  header->AddText(Form("%s%.0f <= jet pT < %.0f",(var=="muptrel")?"deltaR < 0.5 ; muon pT > 5 ; ":"",ptMin,ptMax));
  header->SetTextSize(0.027);
  header->SetTextAlign(12);
  header->SetBorderSize(0);
  header->SetFillStyle(0);
  //header->Draw();


  //RooPlot* sframe = s.frame();
  TH2D *htemp = new TH2D(Form("%s%.0f%.0f",var,ptMin,ptMax),Form("%s%.0f%.0f",var,ptMin,ptMax),100,minXvar,maxXvar,100,0.5,maxYaxis) ;
  //htemp->SetXTitle(Form("%s %.0f < p_{T} < %.0f GeV/c",var,ptMin,ptMax));
  if(var=="svtxm")htemp->SetXTitle("SV mass (GeV/c^{2})");
  else htemp->SetXTitle("JP Disc.");
  htemp->SetYTitle("Entries");
  /*
  htemp->Draw();
  data->plotOn(sframe,Binning(25));
  if(fixCL) {
    model.plotOn(sframe,Components(charmlight),LineStyle(kDashed),LineColor(30),LineWidth(2));
  } else {
    model.plotOn(sframe,Components(light),LineStyle(kDashed),LineColor(kBlue),LineWidth(2));
    model.plotOn(sframe,Components(charm),LineStyle(kDashed),LineColor(kGreen),LineWidth(2));
    model.plotOn(sframe,Components(RooArgSet(charm,light)),LineStyle(kDashed),LineColor(30),LineWidth(2));
  }
  model.plotOn(sframe,Components(bottom),LineStyle(kDashed),LineColor(kRed),LineWidth(2),FillColor(kRed),FillStyle(1));   
  model.plotOn(sframe,LineWidth(2),LineColor(13));
  data->plotOn(sframe,Binning(25));
  model.paramOn(sframe,Layout(0.4,0.9,0.9),Format("NEU",FixedPrecision(3)));
  sframe->Draw("same");
  */

  // --- Perform extended ML fit of composite PDF to data ---
  //RooFitResult *fitresult = model.fitTo(*data,Save(),PrintLevel(-1));
  RooFitResult *fitresult = model.fitTo(*data,Save());
  

  RooPlot* sframe = s.frame();
  //  sframe = s.frame();

  htemp->Draw();
  if(var=="svtxm")data->plotOn(sframe,Binning(24));
  else data->plotOn(sframe,Binning(30));

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
  if(var=="svtxm")data->plotOn(sframe,Binning(24));
  else data->plotOn(sframe,Binning(30));

  model.paramOn(sframe,Layout(0.4,0.9,0.9),Format("NEU",FixedPrecision(3)));
  sframe->Draw("same");
  TLegend *leg = new TLegend(0.61,fixCL?0.60:0.50,0.98,fixCL?0.78:0.75);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry("h_data","PbPb data","p");
  leg->AddEntry(Form("model_Norm[%s]_Comp[bottom]",var),"b","l");
  if(fixCL) {
    leg->AddEntry(Form("model_Norm[%s]_Comp[charmlight]",var),"c + udsg","l");
  } else {
    leg->AddEntry(Form("model_Norm[%s]_Comp[charm]",var),"c","l");
    leg->AddEntry(Form("model_Norm[%s]_Comp[light]",var),"udsg","l");
    leg->AddEntry(Form("model_Norm[%s]_Comp[charm,light]",var),"c + udsg","l");    
  }
  leg->AddEntry(Form("model_Norm[%s]",var),"b + c + udsg","l");
  leg->Draw("same");

  cout<<Form("weight*(abs(refparton_flavorForB)==5&&jtpt>=%f&&jtpt<%f&&%s>=%f&&%s<%f&&bin>=%d&&bin<%d&&abs(jteta)>%f&&abs(jteta)<%f)",ptMin,ptMax,discr,minXdiscr,discr,maxXdiscr,cbinlo,cbinhi,etalo,etahi)<<"\n\n"<<endl;
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
  tBMC->Draw(Form("%s>>hMCB_%d",var,counter),Form("weight*(abs(refparton_flavorForB)==5&&jtpt>=%f&&jtpt<%f&&%s>=%f&&%s<%f&&bin>=%d&&bin<%d&&abs(jteta)>%f&&abs(jteta)<%f)",ptMin,ptMax,discr,minXdiscr,discr,maxXdiscr,cbinlo,cbinhi,etalo,etahi),"goff");
  tCMC->Draw(Form("%s>>hMCC_%d",var,counter),Form("weight*(abs(refparton_flavorForB)==4&&jtpt>=%f&&jtpt<%f&&%s>=%f&&%s<%f&&bin>=%d&&bin<%d&&abs(jteta)>%f&&abs(jteta)<%f)",ptMin,ptMax,discr,minXdiscr,discr,maxXdiscr,cbinlo,cbinhi,etalo,etahi),"goff");
  tQCDMC->Draw(Form("%s>>hMCL_%d",var,counter),Form("weight*(abs(refparton_flavorForB)!=5&&abs(refparton_flavorForB)!=4&&abs(refparton_flavorForB)<99&&jtpt>=%f&&jtpt<%f&&%s>=%f&&%s<%f&&bin>=%d&&bin<%d&&abs(jteta)>%f&&abs(jteta)<%f)",ptMin,ptMax,discr,minXdiscr,discr,maxXdiscr,cbinlo,cbinhi,etalo,etahi),"goff");
  tdata->Draw(Form("%s>>hData_%d",var,counter),Form("jtpt>=%f&&jtpt<%f&&%s>=%f&&%s<%f&&abs(jteta)>%f&&abs(jteta)<%f&&bin>=%d&&bin<%d",ptMin,ptMax,discr,minXdiscr,discr,maxXdiscr,etalo,etahi,cbinlo,cbinhi),"goff");

  hMCLC[counter]->Add( hMCL[counter]);
  hMCLC[counter]->Add( hMCC[counter]);
  fixEmpty(hMCB[counter]); fixEmpty(hMCC[counter]); fixEmpty(hMCL[counter]); fixEmpty(hMCLC[counter]); fixEmpty(hData[counter]);
  
  
  can1[counter] = new TCanvas(Form("can1_%d",counter),Form("can1_%d",counter),700,600);
  hs[counter] = new THStack(Form("hs_%d",counter),"le stack of MC histos");
  fakehs[counter] = new THStack(Form("fakehs_%d",counter),"le fake stack of MC histos");

  can1[counter]->cd();
  if (doLog) can1[counter]->cd()->SetLogy();


  if (doLog){
    hData[counter]->SetMaximum(hData[counter]->GetMaximum()*100);
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
  hleg->AddEntry(hData[counter],"PbPb data","lp");
  hleg->AddEntry(hMCB[counter],"b","f");
  if(!fixCL){
    hleg->AddEntry(hMCC[counter],"c","f");
    hleg->AddEntry(hMCL[counter],"udsg","f");
  }
  if(fixCL){
    hleg->AddEntry(hMCLC[counter],"c + usdg","f");
  }
//   int cbinlo = (double)centlo/2.5;
//   int cbinhi = (double)centhi/2.5
//   double lowcent = (cbinlo)*2.5;
//   double hicent  = (cbinhi)*2.5;
  hleg->Draw("same");
  if(!fixCL)drawText(Form("%d - %d %%",centlo,centhi),0.52,0.62);
  if(fixCL) drawText(Form("%d - %d %%",centlo,centhi),0.52,0.66);
  drawText(Form("#chi^{2}/NDF = %3.1f / %d",chi2,(int)(chi2/chi2NDF)),0.51,0.55);
  drawText("CMS Preliminary",0.15,0.965);
  drawText("#sqrt{s_{NN}} = 2.76 TeV",0.60,0.965);
  drawText("|#eta| < 2.0",0.18,0.88);
  if(comment!="b-tagged sample (SSVHE > 2)") drawText(comment,0.18,0.80);
  if(comment=="b-tagged sample (SSVHE > 2)"){
    drawText("b-tagged sample",0.18,0.80);
    drawText("(SSVHE > 2)",0.18,0.75);
  }

  
  bool printEach=false;
  // --- Print results ---
  //cout <<"b jet fraction in MC = "<<bInitFrac<<endl;
  cout<<"ZZZZZ TEXTcounter: "<<counter<<"   Comment: "<<comment<<endl;
  if(!fixCL) cout<<"ZZZZZ Scale factors: B= "<<Bnorm<<" C= "<<Cnorm<<" Light= "<<Lnorm<<endl;
  if(fixCL)  cout<<"ZZZZZ Scale factors: B= "<<Bnorm<<" LC= "<<LCnorm<<endl;
  cout<<"ZZZZZ Chi2: "<<chi2<<" Chi2/NDF: "<<chi2NDF<<endl;
  cout<<"ZZZZZ var: "<<var<<"  discr: "<<discr<<"    min(discr): "<<minXvar<<"    max(discr): "<<maxXvar<<endl;
  cout<<"ZZZZZ in cent [ "<<centlo<<" - "<<centhi<<" ] %"<<endl;
  cout<<"ZZZZZ b jet fraction in data = "<<Bfraction.getVal()<<endl;
  if(!fixCL) cout <<"ZZZZZ c jet fraction = "<<Cfraction.getVal()<<endl;
  if(!fixCL){
    if (!printEach){
      if (counter<11) can1[counter]->Print("bTagStackedHistosInCentBin_nofixCL.pdf(","pdf");
      if (counter==11) can1[counter]->Print("bTagStackedHistosInCentBin_nofixCL.pdf)","pdf");
    }
    if (printEach) can1[counter]->Print(Form("PDFS/bStack_%sCent%d_%d_nofixCL.pdf",var,cbinlo,cbinhi),"pdf");
  }
  if(fixCL){
    if (!printEach){
      if (counter<11) can1[counter]->Print("bTagStackedHistosInCentBin_fixCL.pdf(","pdf");
      if (counter==11) can1[counter]->Print("bTagStackedHistosInCentBin_fixCL.pdf)","pdf");
    }
    if (printEach) can1[counter]->Print(Form("PDFS/bStack_%sCent%d_%d_fixCL.pdf",var,cbinlo,cbinhi),"pdf");
  }
  counter++;




  // --- Print results ---
  //cout <<"b jet fraction in MC = "<<bInitFrac<<endl;
  cout <<"b jet fraction in data = "<<Bfraction.getVal()<<endl;
  if(!fixCL) cout <<"c jet fraction = "<<Cfraction.getVal()<<endl;

  // --- Save canvas ---
  TString path = Form("gifs/%s_jtpt%.0fto%.0f_%s%.2fto%.2f_%s.gif",var,ptMin,ptMax,discr,minXdiscr,maxXdiscr,fixCL?"CLfixed":"CLfree");
  //cROOFIT->SaveAs(path);

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


Enumerations count(double centMin, double centMax, char *discr, double workingPoint, float ptlo, float pthi, float etalo, float etahi) {

  //*
  TFile *fQCDMC = new TFile("histos/PbPbQCDMC.root"); 
  TFile *fBMC = new TFile("histos/PbPbBMC.root"); 
  TFile *fCMC = new TFile("histos/PbPbCMC.root"); 
  TFile *fdata = new TFile("histos/PbPbdata.root");
  //*/

  /*
  //TFile *fMC = new TFile("histos/ppMC_hiReco_jetTrig.root"); 
  //TFile *fdata = new TFile("histos/ppdata_hiReco_jetTrig.root");
  */

  TTree *tQCDMC = (TTree*) fQCDMC->Get("nt");
  TTree *tBMC = (TTree*) fBMC->Get("nt");
  TTree *tCMC = (TTree*) fCMC->Get("nt");
  TTree *tdata = (TTree*) fdata->Get("nt");

  //TCanvas *c = new TCanvas("c","",600,600);
 
  TH1D *hTaggedLJetsMC = new TH1D("hTaggedLJetsMC","hTaggedLJetsMC",1,ptlo,pthi);
  hTaggedLJetsMC->Sumw2();
  tQCDMC->Draw("jtpt>>hTaggedLJetsMC",Form("weight*(%s>=%f&&jtpt>%f&&jtpt<%f&&abs(jteta)>%f&&abs(jteta)<%f&&abs(refparton_flavorForB)!=4&&abs(refparton_flavorForB)!=5)",discr,workingPoint,ptlo,pthi,etalo,etahi));

  TH1D *hTaggedCJetsMC = new TH1D("hTaggedCJetsMC","hTaggedCJetsMC",1,ptlo,pthi);
  hTaggedCJetsMC->Sumw2();
  tCMC->Draw("jtpt>>hTaggedCJetsMC",Form("weight*(%s>=%f&&jtpt>%f&&jtpt<%f&&abs(jteta)>%f&&abs(jteta)<%f&&abs(refparton_flavorForB)==4)",discr,workingPoint,ptlo,pthi,etalo,etahi));
  
  TH1D *hTaggedBJetsMC = new TH1D("hTaggedBJetsMC","hTaggedBJetsMC",1,ptlo,pthi);
  hTaggedBJetsMC->Sumw2();
  tBMC->Draw("jtpt>>hTaggedBJetsMC",Form("weight*(%s>=%f&&jtpt>%f&&jtpt<%f&&abs(jteta)>%f&&abs(jteta)<%f&&abs(refparton_flavorForB)==5)",discr,workingPoint,ptlo,pthi,etalo,etahi));
  
  TH1F *hTaggedJetsMC = (TH1F*)hTaggedLJetsMC->Clone("hTaggedJetsMC");
  
  hTaggedJetsMC->Add(hTaggedCJetsMC);
  hTaggedJetsMC->Add(hTaggedBJetsMC);
 

  TH1D *hUntaggedLJetsMC = new TH1D("hUntaggedLJetsMC","hUntaggedLJetsMC",1,ptlo,pthi);
  hUntaggedLJetsMC->Sumw2();
  tQCDMC->Draw("jtpt>>hUntaggedLJetsMC",Form("weight*(%s<%f&&jtpt>%f&&jtpt<%f&&abs(jteta)>%f&&abs(jteta)<%f&&abs(refparton_flavorForB)!=4&&abs(refparton_flavorForB)!=5)",discr,workingPoint,ptlo,pthi,etalo,etahi));

  TH1D *hUntaggedCJetsMC = new TH1D("hUntaggedCJetsMC","hUntaggedCJetsMC",1,ptlo,pthi);
  hUntaggedCJetsMC->Sumw2();
  tCMC->Draw("jtpt>>hUntaggedCJetsMC",Form("weight*(%s<%f&&jtpt>%f&&jtpt<%f&&abs(jteta)>%f&&abs(jteta)<%f&&abs(refparton_flavorForB)==4)",discr,workingPoint,ptlo,pthi,etalo,etahi));
  
  TH1D *hUntaggedBJetsMC = new TH1D("hUntaggedBJetsMC","hUntaggedBJetsMC",1,ptlo,pthi);
  hUntaggedBJetsMC->Sumw2();
  tBMC->Draw("jtpt>>hUntaggedBJetsMC",Form("weight*(%s<%f&&jtpt>%f&&jtpt<%f&&abs(jteta)>%f&&abs(jteta)<%f&&abs(refparton_flavorForB)==5)",discr,workingPoint,ptlo,pthi,etalo,etahi));
  
  TH1F *hUntaggedJetsMC = (TH1F*)hUntaggedLJetsMC->Clone("hUntaggedJetsMC");
  
  hUntaggedJetsMC->Add(hUntaggedCJetsMC);
  hUntaggedJetsMC->Add(hUntaggedBJetsMC);



  TH1D *hBjetsWithJPinfoMC = new TH1D("hBjetsWithJPinfoMC","",1,ptlo,pthi);
  hBjetsWithJPinfoMC->Sumw2();
  tBMC->Draw("jtpt>>hBjetsWithJPinfoMC",Form("weight*(abs(refparton_flavorForB)==5&&discr_prob>0&&jtpt>%f&&jtpt<%f&&abs(jteta)>%f&&abs(jteta)<%f)",ptlo,pthi,etalo,etahi));

  TH1D *hBjetsWithCSVinfoMC = new TH1D("hBjetsWithCSVinfoMC","",1,ptlo,pthi);
  hBjetsWithCSVinfoMC->Sumw2();
  tBMC->Draw("jtpt>>hBjetsWithCSVinfoMC",Form("weight*(abs(refparton_flavorForB)==5&&discr_prob>0&&jtpt>%f&&jtpt<%f&&abs(jteta)>%f&&abs(jteta)<%f)",ptlo,pthi,etalo,etahi));
 
  TH1D *hTaggedJetsData = new TH1D("hTaggedJetsData","",1,ptlo,pthi);
  hTaggedJetsData->Sumw2();
  tdata->Draw("jtpt>>hTaggedJetsData",Form("weight*(%s>=%f&&jtpt>%f&&jtpt<%f&&abs(jteta)>%f&&abs(jteta)<%f)",discr,workingPoint,ptlo,pthi,etalo,etahi));

  TH1D *hUntaggedJetsData = new TH1D("hUntaggedJetsData","",1,ptlo,pthi);
  hUntaggedJetsData->Sumw2();
  tdata->Draw("jtpt>>hUntaggedJetsData",Form("weight*(%s<%f&&jtpt>%f&&jtpt<%f&&abs(jteta)>%f&&abs(jteta)<%f)",discr,workingPoint,ptlo,pthi,etalo,etahi));

  Enumerations res;
 cout<<"hTaggedJetsMC entries : "<<hTaggedJetsMC->GetEntries()<<endl;
 cout<<"hTaggedJetsMC entries bin1 : "<<hTaggedJetsMC->GetBinContent(1)<<endl;
  res.nTaggedJetsMC = hTaggedJetsMC->GetBinContent(1);
  res.nUntaggedJetsMC = hUntaggedJetsMC->GetBinContent(1);
  res.nJetsMC = res.nTaggedJetsMC + res.nUntaggedJetsMC;
  res.nTaggedBjetsMC = hTaggedBJetsMC->GetBinContent(1);
  res.nUntaggedBjetsMC = hUntaggedBJetsMC->GetBinContent(1);
  res.nBjetsMC = res.nTaggedBjetsMC + res.nUntaggedBjetsMC;
  res.nNonBjetsMC = res.nJetsMC - res.nBjetsMC;
  res.nTaggedNonBjetsMC = res.nTaggedJetsMC - res.nTaggedBjetsMC;

  res.nTaggedJetsData = hTaggedJetsData->GetBinContent(1);
  res.nUntaggedJetsData = hUntaggedJetsData->GetBinContent(1);

  res.cbForJP = hBjetsWithJPinfoMC->GetBinContent(1) / res.nBjetsMC;
  res.cbForCSV = hBjetsWithCSVinfoMC->GetBinContent(1) / res.nBjetsMC;

  res.nTaggedJetsMCError = hTaggedJetsMC->GetBinError(1);
  res.nUntaggedJetsMCError = hUntaggedJetsMC->GetBinError(1);
  res.nJetsMCError = addError(res.nTaggedJetsMCError,res.nUntaggedJetsMCError);
  res.nTaggedBjetsMCError = hTaggedBJetsMC->GetBinError(1);
  res.nUntaggedBjetsMCError = hUntaggedBJetsMC->GetBinError(1);
  res.nBjetsMCError = addError(res.nTaggedBjetsMCError,res.nUntaggedBjetsMCError);
  res.nNonBjetsMCError = substractError(res.nJetsMCError,res.nBjetsMCError);
  res.nTaggedNonBjetsMCError = substractError(res.nTaggedJetsMCError,res.nTaggedBjetsMCError);

  res.nTaggedJetsDataError = hTaggedJetsData->GetBinError(1);
  res.nUntaggedJetsDataError = hUntaggedJetsData->GetBinError(1);

  return res;
 
}

