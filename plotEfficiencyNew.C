#include "TROOT.h"
#include <iostream>
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include <string>

using namespace std;

void plotEfficiencyNew(int ppPbPb=0, int doCent=0, int doBoth=1){

  if(!ppPbPb) doCent=0;
  if(doBoth) doCent=0;


  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetNdivisions(404,"y");
  gStyle->SetGridColor(kGray+1);
  if(doCent) gStyle->SetNdivisions(5,"x");

  //gStyle->SetTextSize(0.1);

  gStyle->SetMarkerSize(0.7);
  gStyle->SetTextFont(43); 
  gStyle->SetLabelFont(43,"xy"); 
  gStyle->SetTitleFont(43,"xy"); 
  gStyle->SetTitleOffset(1.3,"xy"); 
  gStyle->SetEndErrorSize(0);
  gROOT->ForceStyle(1);
  
  int centMax=4;
  if(!doCent) centMax=1;
  
  gStyle->SetOptTitle(0);
  if(doCent){
    gStyle->SetPadBottomMargin(0.12);
    gStyle->SetPadLeftMargin(0.2);
    gStyle->SetPadRightMargin(0.12);
    gStyle->SetPadTopMargin(0.12);
  }
  else{
    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadRightMargin(0.15);
    gStyle->SetPadTopMargin(0.15);
  }
  if(doCent){
    gStyle->SetTitleSize(0.06,"xy");
    gStyle->SetLabelSize(0.06,"xy");
  }
  else{
    gStyle->SetTextFont(43); 
    gStyle->SetTitleSize(28,"xy");
    gStyle->SetLabelSize(28,"xy");
  }
  gStyle->SetTitleW(0.2);
  gStyle->SetTitleH(0.1);
  
  TFile *f = new TFile();
  TFile *fC = new TFile();
  TFile *fB = new TFile();
  TH2F *hBottom_csvSimple[3], *hBottom_prob[3], *hBottom_ssvHighEff[3], *hBottom_ssvHighPur[3];
  TH2F *hCharm_csvSimple[3], *hCharm_prob[3], *hCharm_ssvHighEff[3], *hCharm_ssvHighPur[3];
  TH2F *hLight_csvSimple[3], *hLight_prob[3], *hLight_ssvHighEff[3], *hLight_ssvHighPur[3];
  TH1F *hBottom_csvSimple_proj[4], *hBottom_prob_proj[4];//, *hBottom_probb_proj[4], *hBottom_tcHighEff_proj[4], *hBottom_tcHighPur_proj[4];
  TH1F *hCharm_csvSimple_proj[4], *hCharm_prob_proj[4];//, *hCharm_probb_proj[4], *hCharm_tcHighEff_proj[4], *hCharm_tcHighPur_proj[4];
  TH1F *hLight_csvSimple_proj[4], *hLight_prob_proj[4];//, *hLight_probb_proj[4], *hLight_tcHighEff_proj[4], *hLight_tcHighPur_proj[4];
  
  TH1F *hBottom_ssvHighEff_proj[4],*hCharm_ssvHighEff_proj[4],*hLight_ssvHighEff_proj[4];
  TH1F *hBottom_ssvHighPur_proj[4],*hCharm_ssvHighPur_proj[4],*hLight_ssvHighPur_proj[4];
  
  TF1 *fx = new TF1("fx","x",0,1);
  fx->SetLineStyle(7);
  
  TGraph *gBcsvSimple[4], *gBprob[4], /* *gBprobb[4], *gBtcHighEff[4], *gBtcHighPur[4],*/ *gBssvHighEff[4], *gBssvHighPur[4];
  TGraph *gBwerkingPoint[4];

  TGraph *gCcsvSimple[4], *gCprob[4], /**gCprobb[4], *gCtcHighEff[4], *gCtcHighPur[4],*/ *gCssvHighEff[4], *gCssvHighPur[4];
  TGraph *gCwerkingPoint[4];
  
  char name[500];
  TH1F *hFrame[4];
  TText *centText[4];

  int nLoop;
  if(doBoth) nLoop = 3;
  else nLoop = 1;
  for(int ic=0; ic<nLoop; ic++){
    if(ppPbPb || (doBoth && ic==0)){
      f=TFile::Open("hist_PbPb_Light.root");
      fC= TFile::Open("hist_PbPb_Charm.root");
      fB=TFile::Open("hist_PbPb_Bottom.root");
    }
    else if(!ppPbPb || (doBoth && ic==1)){
      f=TFile::Open("hist_pp_Light.root");
      fC=TFile::Open("hist_pp_Charm.root");
      fB=TFile::Open("hist_pp_Bottom.root");
    }
    else if(doBoth && ic==2){
      f=TFile::Open("hist_pPb_Light.root");
      fC=TFile::Open("hist_pPb_Charm.root");
      fB=TFile::Open("hist_pPb_Bottom.root");
    }
    hBottom_csvSimple[ic]=(TH2F*)fB->Get("hBottom_csvSimple")->Clone("hBottom_csvSimple");
    hBottom_prob[ic]=(TH2F*)fB->Get("hBottom_prob")->Clone("hBottom_prob");
    // hBottom_probb=(TH2F*)fB->Get("hBottom_probb");
    //hBottom_tcHighEff=(TH2F*)fB->Get("hBottom_logtcHighEff");
    // hBottom_tcHighPur=(TH2F*)fB->Get("hBottom_logtcHighPur");
    hBottom_ssvHighEff[ic]=(TH2F*)fB->Get("hBottom_ssvHighEff")->Clone("hBottom_ssvHighEff");
    hBottom_ssvHighPur[ic]=(TH2F*)fB->Get("hBottom_ssvHighPur")->Clone("hBottom_ssvHighPur");

    hCharm_csvSimple[ic]=(TH2F*)fC->Get("hCharm_csvSimple")->Clone("hCharm_csvSimple");
    hCharm_prob[ic]=(TH2F*)fC->Get("hCharm_prob")->Clone("hCharm_prob");
    //hCharm_probb=(TH2F*)fC->Get("hCharm_probb");
    // hCharm_tcHighEff=(TH2F*)fC->Get("hCharm_logtcHighEff");
    // hCharm_tcHighPur=(TH2F*)fC->Get("hCharm_logtcHighPur");
    hCharm_ssvHighEff[ic]=(TH2F*)fC->Get("hCharm_ssvHighEff")->Clone("hCharm_ssvHighEff");
    hCharm_ssvHighPur[ic]=(TH2F*)fC->Get("hCharm_ssvHighPur")->Clone("hCharm_ssvHighPur");

    hLight_csvSimple[ic]=(TH2F*)f->Get("hLight_csvSimple")->Clone("hLight_csvSimple");
    hLight_prob[ic]=(TH2F*)f->Get("hLight_prob")->Clone("hLight_prob");
    //hLight_probb=(TH2F*)f->Get("hLight_probb");
    //hLight_tcHighEff=(TH2F*)f->Get("hLight_logtcHighEff");
    //hLight_tcHighPur=(TH2F*)f->Get("hLight_logtcHighPur");
    hLight_ssvHighEff[ic]=(TH2F*)f->Get("hLight_ssvHighEff")->Clone("hLight_ssvHighEff");
    hLight_ssvHighPur[ic]=(TH2F*)f->Get("hLight_ssvHighPur")->Clone("hLight_ssvHighPur");

    int binbound[5]={0,8,16,24,40};
    if(!doCent){
      if(ic==0){
	binbound[0]=0;
	binbound[1]=40;
      }
      if(ic==1){
	binbound[1]=0;
	binbound[2]=40;
      }
      if(ic==2){
	binbound[2]=0;
	binbound[3]=40;
      }
    }

    hFrame[ic]=new TH1F(Form("hFrame_%d",ic),Form("hFrame_%d",ic),1,0,1);
    hFrame[ic]->SetMinimum(0.99e-4);
    hFrame[ic]->SetMaximum(5.);
  
    if(doCent){
      if(ic==0)hFrame[ic]->SetTitle("0-20%");
      if(ic==1)hFrame[ic]->SetTitle("20-40%");
      if(ic==2)hFrame[ic]->SetTitle("40-60%");
      if(ic==3)hFrame[ic]->SetTitle("60-100%");
      if(ic==4)hFrame[ic]->SetTitle("50-70%");
      if(ic==5)hFrame[ic]->SetTitle("70-100%");
    }
    else hFrame[ic]->SetTitle("80 < p_{T} 120 GeV/c RelVal");
    hFrame[ic]->GetXaxis()->SetTitle("b-jet efficiency");
    hFrame[ic]->GetXaxis()->CenterTitle();
    hFrame[ic]->GetYaxis()->CenterTitle();
    hFrame[ic]->GetYaxis()->SetTitle("udsg-jet efficiency");
  
  
    sprintf(name,"hBottom_csvSimple_prof_%d",ic);        
    hBottom_csvSimple_proj[ic] = (TH1F*)hBottom_csvSimple[ic]->ProjectionX(name,binbound[ic]+1,binbound[ic+1]);
    sprintf(name,"hBottom_prob_prof_%d",ic);        
    hBottom_prob_proj[ic] = (TH1F*)hBottom_prob[ic]->ProjectionX(name,binbound[ic]+1,binbound[ic+1]);
    //sprintf(name,"hBottom_probb_prof_%d",ic);        
    //hBottom_probb_proj[ic] = (TH1F*)hBottom_probb->ProjectionX(name,binbound[ic]+1,binbound[ic+1]);
    // sprintf(name,"hBottom_tcHighEff_prof_%d",ic);        
    //hBottom_tcHighEff_proj[ic] = (TH1F*)hBottom_tcHighEff->ProjectionX(name,binbound[ic]+1,binbound[ic+1]);
    //sprintf(name,"hBottom_tcHighPur_prof_%d",ic);        
    //hBottom_tcHighPur_proj[ic] = (TH1F*)hBottom_tcHighPur->ProjectionX(name,binbound[ic]+1,binbound[ic+1]);
    sprintf(name,"hBottom_ssvHighEff_prof_%d",ic);        
    hBottom_ssvHighEff_proj[ic] = (TH1F*)hBottom_ssvHighEff[ic]->ProjectionX(name,binbound[ic]+1,binbound[ic+1]);
    sprintf(name,"hBottom_ssvHighPur_prof_%d",ic);        
    hBottom_ssvHighPur_proj[ic] = (TH1F*)hBottom_ssvHighPur[ic]->ProjectionX(name,binbound[ic]+1,binbound[ic+1]);
  
  
    sprintf(name,"hCharm_csvSimple_prof_%d",ic);        
    hCharm_csvSimple_proj[ic] = (TH1F*)hCharm_csvSimple[ic]->ProjectionX(name,binbound[ic]+1,binbound[ic+1]);
    sprintf(name,"hCharm_prob_prof_%d",ic);        
    hCharm_prob_proj[ic] = (TH1F*)hCharm_prob[ic]->ProjectionX(name,binbound[ic]+1,binbound[ic+1]);
    //sprintf(name,"hCharm_probb_prof_%d",ic);        
    //hCharm_probb_proj[ic] = (TH1F*)hCharm_probb->ProjectionX(name,binbound[ic]+1,binbound[ic+1]);
    //sprintf(name,"hCharm_tcHighEff_prof_%d",ic);        
    // hCharm_tcHighEff_proj[ic] = (TH1F*)hCharm_tcHighEff->ProjectionX(name,binbound[ic]+1,binbound[ic+1]);
    // sprintf(name,"hCharm_tcHighPur_prof_%d",ic);        
    // hCharm_tcHighPur_proj[ic] = (TH1F*)hCharm_tcHighPur->ProjectionX(name,binbound[ic]+1,binbound[ic+1]);
    sprintf(name,"hCharm_ssvHighEff_prof_%d",ic);        
    hCharm_ssvHighEff_proj[ic] = (TH1F*)hCharm_ssvHighEff[ic]->ProjectionX(name,binbound[ic]+1,binbound[ic+1]);
    sprintf(name,"hCharm_ssvHighPur_prof_%d",ic);        
    hCharm_ssvHighPur_proj[ic] = (TH1F*)hCharm_ssvHighPur[ic]->ProjectionX(name,binbound[ic]+1,binbound[ic+1]);


    sprintf(name,"hLight_csvSimple_prof_%d",ic);        
    hLight_csvSimple_proj[ic] = (TH1F*)hLight_csvSimple[ic]->ProjectionX(name,binbound[ic]+1,binbound[ic+1]);
    sprintf(name,"hLight_prob_prof_%d",ic);        
    hLight_prob_proj[ic] = (TH1F*)hLight_prob[ic]->ProjectionX(name,binbound[ic]+1,binbound[ic+1]);
    //sprintf(name,"hLight_probb_prof_%d",ic);        
    // hLight_probb_proj[ic] = (TH1F*)hLight_probb->ProjectionX(name,binbound[ic]+1,binbound[ic+1]);
    // sprintf(name,"hLight_tcHighEff_prof_%d",ic);        
    //hLight_tcHighEff_proj[ic] = (TH1F*)hLight_tcHighEff->ProjectionX(name,binbound[ic]+1,binbound[ic+1]);
    //sprintf(name,"hLight_tcHighPur_prof_%d",ic);        
    //hLight_tcHighPur_proj[ic] = (TH1F*)hLight_tcHighPur->ProjectionX(name,binbound[ic]+1,binbound[ic+1]);   
    sprintf(name,"hLight_ssvHighEff_prof_%d",ic);        
    hLight_ssvHighEff_proj[ic] = (TH1F*)hLight_ssvHighEff[ic]->ProjectionX(name,binbound[ic]+1,binbound[ic+1]);
    sprintf(name,"hLight_ssvHighPur_prof_%d",ic);        
    hLight_ssvHighPur_proj[ic] = (TH1F*)hLight_ssvHighPur[ic]->ProjectionX(name,binbound[ic]+1,binbound[ic+1]);    
  

    gBcsvSimple[ic]=new TGraph(hBottom_csvSimple_proj[ic]->GetNbinsX());
    gBprob[ic]=new TGraph(hBottom_prob_proj[ic]->GetNbinsX());
    //gBprobb[ic]=new TGraph(hBottom_probb_proj[ic]->GetNbinsX());
    //gBtcHighEff[ic]=new TGraph(hBottom_tcHighEff_proj[ic]->GetNbinsX());
    //gBtcHighPur[ic]=new TGraph(hBottom_tcHighPur_proj[ic]->GetNbinsX());
    gBssvHighEff[ic]=new TGraph(hBottom_ssvHighEff_proj[ic]->GetNbinsX());
    gBssvHighPur[ic]=new TGraph(hBottom_ssvHighPur_proj[ic]->GetNbinsX());

    sprintf(name,"gBcsvSimple%d",ic);
    gBcsvSimple[ic]->SetName(name);
    sprintf(name,"gBprob%d",ic);
    gBprob[ic]->SetName(name);
    // sprintf(name,"gBprobb%d",ic);
    // gBprobb[ic]->SetName(name);
    // sprintf(name,"gBtcHighEff%d",ic);
    // gBtcHighEff[ic]->SetName(name);
    // sprintf(name,"gBtcHighPur%d",ic);
    // gBtcHighPur[ic]->SetName(name);
    sprintf(name,"gBssvHighEff%d",ic);
    gBssvHighEff[ic]->SetName(name);
    sprintf(name,"gBssvHighPur%d",ic);
    gBssvHighPur[ic]->SetName(name);

    // be sure to get overflow bins
    float totalBottom = hBottom_csvSimple_proj[ic]->Integral(0,hBottom_csvSimple_proj[ic]->GetNbinsX()+1);
    float totalLight = hLight_csvSimple_proj[ic]->Integral(0,hLight_csvSimple_proj[ic]->GetNbinsX()+1);
    for(int ib=0;ib<hBottom_csvSimple_proj[ic]->GetNbinsX();ib++){
      gBcsvSimple[ic]->SetPoint(ib,hBottom_csvSimple_proj[ic]->Integral(ib+1,hBottom_csvSimple_proj[ic]->GetNbinsX()+1)/totalBottom,hLight_csvSimple_proj[ic]->Integral(ib+1,hLight_csvSimple_proj[ic]->GetNbinsX()+1)/totalLight);
    }
    totalBottom = hBottom_prob_proj[ic]->Integral(0,hBottom_prob_proj[ic]->GetNbinsX()+1);
    totalLight = hLight_prob_proj[ic]->Integral(0,hLight_prob_proj[ic]->GetNbinsX()+1);
    for(int ib=0;ib<hBottom_prob_proj[ic]->GetNbinsX();ib++){
      gBprob[ic]->SetPoint(ib,hBottom_prob_proj[ic]->Integral(ib+1,hBottom_prob_proj[ic]->GetNbinsX()+1)/totalBottom,hLight_prob_proj[ic]->Integral(ib+1,hLight_prob_proj[ic]->GetNbinsX()+1)/totalLight);
    }
    /* float totalBottom = hBottom_probb_proj[ic]->Integral(0,hBottom_probb_proj[ic]->GetNbinsX()+1);
       float totalLight = hLight_probb_proj[ic]->Integral(0,hLight_probb_proj[ic]->GetNbinsX()+1);
       for(int ib=0;ib<hBottom_probb_proj[ic]->GetNbinsX();ib++){
       gBprobb[ic]->SetPoint(ib,hBottom_probb_proj[ic]->Integral(ib+1,hBottom_probb_proj[ic]->GetNbinsX()+1)/totalBottom,hLight_probb_proj[ic]->Integral(ib+1,hLight_probb_proj[ic]->GetNbinsX()+1)/totalLight);
       }
       float totalBottom = hBottom_tcHighEff_proj[ic]->Integral(0,hBottom_tcHighEff_proj[ic]->GetNbinsX()+1);
       float totalLight = hLight_tcHighEff_proj[ic]->Integral(0,hLight_tcHighEff_proj[ic]->GetNbinsX()+1);*/
    cout<<" totalBottom "<<totalBottom<<" total on in histo "<<hBottom_prob_proj[ic]->Integral(0,hBottom_prob_proj[ic]->GetNbinsX()+1)<<endl;
       /*for(int ib=0;ib<hBottom_tcHighEff_proj[ic]->GetNbinsX();ib++){
       gBtcHighEff[ic]->SetPoint(ib,hBottom_tcHighEff_proj[ic]->Integral(ib+1,hBottom_tcHighEff_proj[ic]->GetNbinsX()+1)/totalBottom,hLight_tcHighEff_proj[ic]->Integral(ib+1,hLight_tcHighEff_proj[ic]->GetNbinsX()+1)/totalLight);
       }
       float totalBottom = hBottom_tcHighPur_proj[ic]->Integral(0,hBottom_tcHighPur_proj[ic]->GetNbinsX()+1);
       float totalLight = hLight_tcHighPur_proj[ic]->Integral(0,hLight_tcHighPur_proj[ic]->GetNbinsX()+1);    
       for(int ib=0;ib<hBottom_tcHighPur_proj[ic]->GetNbinsX();ib++){
       gBtcHighPur[ic]->SetPoint(ib,hBottom_tcHighPur_proj[ic]->Integral(ib+1,hBottom_tcHighPur_proj[ic]->GetNbinsX()+1)/totalBottom,hLight_tcHighPur_proj[ic]->Integral(ib+1,hLight_tcHighPur_proj[ic]->GetNbinsX()+1)/totalLight);
       }*/
    totalBottom = hBottom_ssvHighEff_proj[ic]->Integral(0,hBottom_ssvHighEff_proj[ic]->GetNbinsX()+1);
    totalLight = hLight_ssvHighEff_proj[ic]->Integral(0,hLight_ssvHighEff_proj[ic]->GetNbinsX()+1);
    for(int ib=0;ib<hBottom_ssvHighEff_proj[ic]->GetNbinsX();ib++){
      gBssvHighEff[ic]->SetPoint(ib,hBottom_ssvHighEff_proj[ic]->Integral(ib+1,hBottom_ssvHighEff_proj[ic]->GetNbinsX()+1)/totalBottom,hLight_ssvHighEff_proj[ic]->Integral(ib+1,hLight_ssvHighEff_proj[ic]->GetNbinsX()+1)/totalLight);
    }
    totalBottom = hBottom_ssvHighPur_proj[ic]->Integral(0,hBottom_ssvHighPur_proj[ic]->GetNbinsX()+1);
    totalLight = hLight_ssvHighPur_proj[ic]->Integral(0,hLight_ssvHighPur_proj[ic]->GetNbinsX()+1);
    for(int ib=0;ib<hBottom_ssvHighPur_proj[ic]->GetNbinsX();ib++){
      gBssvHighPur[ic]->SetPoint(ib,hBottom_ssvHighPur_proj[ic]->Integral(ib+1,hBottom_ssvHighPur_proj[ic]->GetNbinsX()+1)/totalBottom,hLight_ssvHighPur_proj[ic]->Integral(ib+1,hLight_ssvHighPur_proj[ic]->GetNbinsX()+1)/totalLight);
    }

     gCcsvSimple[ic]=new TGraph(hCharm_csvSimple_proj[ic]->GetNbinsX());
    gCprob[ic]=new TGraph(hCharm_prob_proj[ic]->GetNbinsX());
    // gCprobb[ic]=new TGraph(hCharm_probb_proj[ic]->GetNbinsX());
    //  gCtcHighEff[ic]=new TGraph(hCharm_tcHighEff_proj[ic]->GetNbinsX());
    // gCtcHighPur[ic]=new TGraph(hCharm_tcHighPur_proj[ic]->GetNbinsX());
    gCssvHighEff[ic]=new TGraph(hCharm_ssvHighEff_proj[ic]->GetNbinsX());
    gCssvHighPur[ic]=new TGraph(hCharm_ssvHighPur_proj[ic]->GetNbinsX());

    sprintf(name,"gCcsvSimple%d",ic);
    gCcsvSimple[ic]->SetName(name);
    sprintf(name,"gCprob%d",ic);
    gCprob[ic]->SetName(name);
    // sprintf(name,"gCprobb%d",ic);
    // gCprobb[ic]->SetName(name);
    // sprintf(name,"gCtcHighEff%d",ic);
    //gCtcHighEff[ic]->SetName(name);
    // sprintf(name,"gCtcHighPur%d",ic);
    // gCtcHighPur[ic]->SetName(name);
    sprintf(name,"gCssvHighEff%d",ic);
    gCssvHighEff[ic]->SetName(name);
    sprintf(name,"gCssvHighPur%d",ic);
    gCssvHighPur[ic]->SetName(name);
    

    totalBottom = hBottom_csvSimple_proj[ic]->Integral(0,hBottom_csvSimple_proj[ic]->GetNbinsX()+1);
    float totalCharm = hCharm_csvSimple_proj[ic]->Integral(0,hCharm_csvSimple_proj[ic]->GetNbinsX()+1);
    for(int ib=0;ib<hBottom_csvSimple_proj[ic]->GetNbinsX();ib++){
      gCcsvSimple[ic]->SetPoint(ib,hBottom_csvSimple_proj[ic]->Integral(ib+1,hBottom_csvSimple_proj[ic]->GetNbinsX()+1)/totalBottom,hCharm_csvSimple_proj[ic]->Integral(ib+1,hCharm_csvSimple_proj[ic]->GetNbinsX()+1)/totalCharm);

    }
    totalBottom = hBottom_prob_proj[ic]->Integral(0,hBottom_prob_proj[ic]->GetNbinsX()+1);
    totalCharm = hCharm_prob_proj[ic]->Integral(0,hCharm_prob_proj[ic]->GetNbinsX()+1);
    for(int ib=0;ib<hBottom_prob_proj[ic]->GetNbinsX();ib++){
      gCprob[ic]->SetPoint(ib,hBottom_prob_proj[ic]->Integral(ib+1,hBottom_prob_proj[ic]->GetNbinsX()+1)/totalBottom,hCharm_prob_proj[ic]->Integral(ib+1,hCharm_prob_proj[ic]->GetNbinsX()+1)/totalCharm);
    }
    /* float totalBottom = hBottom_probb_proj[ic]->Integral(0,hBottom_probb_proj[ic]->GetNbinsX()+1);
       float totalCharm = hCharm_probb_proj[ic]->Integral(0,hCharm_probb_proj[ic]->GetNbinsX()+1);
       for(int ib=0;ib<hBottom_probb_proj[ic]->GetNbinsX();ib++){
       gCprobb[ic]->SetPoint(ib,hBottom_probb_proj[ic]->Integral(ib+1,hBottom_probb_proj[ic]->GetNbinsX()+1)/totalBottom,hCharm_probb_proj[ic]->Integral(ib+1,hCharm_probb_proj[ic]->GetNbinsX()+1)/totalCharm);
       }
       float totalBottom = hBottom_tcHighEff_proj[ic]->Integral(0,hBottom_tcHighEff_proj[ic]->GetNbinsX()+1);
       float totalCharm = hCharm_tcHighEff_proj[ic]->Integral(0,hCharm_tcHighEff_proj[ic]->GetNbinsX()+1);
       for(int ib=0;ib<hBottom_tcHighEff_proj[ic]->GetNbinsX();ib++){
       gCtcHighEff[ic]->SetPoint(ib,hBottom_tcHighEff_proj[ic]->Integral(ib+1,hBottom_tcHighEff_proj[ic]->GetNbinsX()+1)/totalBottom,hCharm_tcHighEff_proj[ic]->Integral(ib+1,hCharm_tcHighEff_proj[ic]->GetNbinsX()+1)/totalCharm);
       }
       float totalBottom = hBottom_tcHighPur_proj[ic]->Integral(0,hBottom_tcHighPur_proj[ic]->GetNbinsX()+1);
       float totalCharm = hCharm_tcHighPur_proj[ic]->Integral(0,hCharm_tcHighPur_proj[ic]->GetNbinsX()+1);
       for(int ib=0;ib<hBottom_tcHighPur_proj[ic]->GetNbinsX();ib++){
       gCtcHighPur[ic]->SetPoint(ib,hBottom_tcHighPur_proj[ic]->Integral(ib+1,hBottom_tcHighPur_proj[ic]->GetNbinsX()+1)/totalBottom,hCharm_tcHighPur_proj[ic]->Integral(ib+1,hCharm_tcHighPur_proj[ic]->GetNbinsX()+1)/totalCharm);
       }*/
    totalBottom = hBottom_ssvHighEff_proj[ic]->Integral(0,hBottom_ssvHighEff_proj[ic]->GetNbinsX()+1);
    totalCharm = hCharm_ssvHighEff_proj[ic]->Integral(0,hCharm_ssvHighEff_proj[ic]->GetNbinsX()+1);
    for(int ib=0;ib<hBottom_ssvHighEff_proj[ic]->GetNbinsX();ib++){
      // cout<<"ib: "<<ib<<" ic: "<<ic<<" hCharm_ssvHighEff_proj[ic]->Integral(ib+1,hBottom_ssvHighEff_proj[ic]->GetNbinsX()+1): "<<hCharm_ssvHighEff_proj[ic]->Integral(ib+1,hBottom_ssvHighEff_proj[ic]->GetNbinsX()+1)<<"   hCharm_ssvHighEff_proj[ic]->Integral(ib+1,hCharm_ssvHighEff_proj[ic]->GetNbinsX()+1): "<<hCharm_ssvHighEff_proj[ic]->Integral(ib+1,hCharm_ssvHighEff_proj[ic]->GetNbinsX()+1)<<endl;
      gCssvHighEff[ic]->SetPoint(ib,hBottom_ssvHighEff_proj[ic]->Integral(ib+1,hBottom_ssvHighEff_proj[ic]->GetNbinsX()+1)/totalBottom,hCharm_ssvHighEff_proj[ic]->Integral(ib+1,hCharm_ssvHighEff_proj[ic]->GetNbinsX()+1)/totalCharm);
    }
    totalBottom = hBottom_ssvHighPur_proj[ic]->Integral(0,hBottom_ssvHighPur_proj[ic]->GetNbinsX()+1);
    totalCharm = hCharm_ssvHighPur_proj[ic]->Integral(0,hCharm_ssvHighPur_proj[ic]->GetNbinsX()+1);
    for(int ib=0;ib<hBottom_ssvHighPur_proj[ic]->GetNbinsX();ib++){
      gCssvHighPur[ic]->SetPoint(ib,hBottom_ssvHighPur_proj[ic]->Integral(ib+1,hBottom_ssvHighPur_proj[ic]->GetNbinsX()+1)/totalBottom,hCharm_ssvHighPur_proj[ic]->Integral(ib+1,hCharm_ssvHighPur_proj[ic]->GetNbinsX()+1)/totalCharm);
    }
    
    /* if(doCent){
      c->cd(ic+1);
      c->GetPad(ic+1)->SetLogy();
      c->GetPad(ic+1)->SetGridy();
      c->GetPad(ic+1)->SetGridx();
    }
    else{
      c->SetLogy();
      c->SetGridy();
      c->SetGridx();
      }*/
    
    
    // hFrame[ic]->Draw();
    
    if(doCent>0){
      if(ic==0) centText[ic]=new TText(0.1,1.1,"0-20%");
      if(ic==1) centText[ic]=new TText(0.1,1.1,"20-40%");
      if(ic==2) centText[ic]=new TText(0.1,1.1,"40-60%");
      if(ic==3) centText[ic]=new TText(0.1,1.1,"60-100%");
      
      centText[ic]->Draw();
      
    }


    //f->Close();
    //fC->Close();
    //fB->Close();
  }



  //fx->Draw("same");
  TLatex *prel = new TLatex(0.01,6,"CMS Simulation");
  if(!doCent)prel->SetTextSize(23);
  if(doCent)prel->SetTextSize(17);
  prel->SetTextFont(43);
  prel->SetTextColor(kBlack);
  //prel->Draw();

  TLatex *roots;
  if(ppPbPb)roots = new TLatex(0.7,6,"#sqrt{s_{NN}} = 2.76 TeV");
  if(!ppPbPb)roots = new TLatex(0.7,6,"#sqrt{s_{NN}} = 5.02 TeV");
  if(!doCent)roots->SetTextSize(23);
  if(doCent)roots->SetTextSize(17);
  roots->SetTextFont(43);
  roots->SetTextColor(kBlack);
  //roots->Draw();
 
  TLegend *leg = new TLegend(0.2,0.55,0.6,0.75);
  //if(doCent) leg->SetTextSize(20);
  //else leg->SetTextSize(23);
  if(ppPbPb&&doCent==0)leg->SetHeader("Pythia+Hydjet (PbPb), 0-100%");
  else if(ppPbPb&&doCent==1)leg->SetHeader("Pythia+Hydjet (PbPb)");
  else leg->SetHeader("Pythia (pp)");
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  //leg->AddEntry(gBcsvSimple[0],"CSV","l");
  //leg->AddEntry(gBprob[0],"JPB","l");
  leg->AddEntry(gBprob[0],"JP","lp");
  //leg->AddEntry(gBtcHighEff[0],"TCHE","l");
  //leg->AddEntry(gBtcHighPur[0],"TCHP","l");
  leg->AddEntry(gBssvHighEff[0],"SSVHE","l");
  //leg->AddEntry(gBssvHighPur[0],"SSVHP","l");
  //leg->Draw();
  //c->cd(2);



  /*

    TCanvas *c2=new TCanvas("c2","cent dependence",1);
    gcsvSimple[0]->SetLineColor(1);
    gcsvSimple[0]->Draw("al");
    for(int i=1;i<6;i++){
    gcsvSimple[i]->SetLineColor(i+1);
    gcsvSimple[i]->Draw("l");
    
    }
  */

  TCanvas *c2;
  if(doCent)c2=new TCanvas("c2","Bottom vs Bg",900,800);
  else c2=new TCanvas("c2","Bottom vs Bg",700,600);
  if(!doCent){
    //c2->cd();
    c2->SetLogy();
    c2->SetGridy();
    c2->SetGridx();
  }

  gCssvHighEff[0]->Draw();
  gCssvHighEff[1]->Draw("same");
  c2->SaveAs("test.pdf");

  TH1F *hFrame2[4];
  hFrame2[0] = (TH1F *)hFrame[1]->Clone(Form("hFrame_%d",1));
  hFrame2[0]->SetYTitle("Misidentification Probability");
  hFrame2[0]->Draw();
  
  if(doCent)c2->Divide(2,2);

    //     if(doCent){
    //       if(ic==0) centText[0]=new TText(0.1,1.1,"0-22%");
    //       if(ic==1) centText[1]=new TText(0.1,1.1,"20-44%");
    //       if(ic==2) centText[2]=new TText(0.1,1.1,"40-60%");
    //       if(ic==3) centText[3]=new TText(0.1,1.1,"60-100%");
    //     }
    // //       if(ic==0)hFrame[ic]->SetTitle("0-20%");
    // //       if(ic==1)hFrame[ic]->SetTitle("20-40%");
    // //       if(ic==2)hFrame[ic]->SetTitle("40-60%");
    // //       if(ic==3)hFrame[ic]->SetTitle("60-100%");
    // // //       if(ic==4)hFrame[ic]->SetTitle("50-70%");
    // // //       if(ic==5)hFrame[ic]->SetTitle("70-100%");
     
    //     //centText[ic]->Draw();
    //     cout<<"centText ["<<ic<<"] : "<<centText[ic]<<endl;
	
    //       // }
    


  if(!doCent) c2->cd();
  for(int ic=0; ic<nLoop; ic++){
    if(doCent){
      c2->cd(ic+1);
      c2->GetPad(ic+1)->SetLogy();
      c2->GetPad(ic+1)->SetGridy();
      c2->GetPad(ic+1)->SetGridx();
    }
    if(doCent>0){
      if(ic==0) centText[ic]=new TText(0.1,1.1,"0-20%");
      if(ic==1) centText[ic]=new TText(0.1,1.1,"20-40%");
      if(ic==2) centText[ic]=new TText(0.1,1.1,"40-60%");
      if(ic==3) centText[ic]=new TText(0.1,1.1,"60-100%");
      
      centText[ic]->Draw();
      
    }

    /* gCprobb[ic]->SetLineWidth(3);
       gCprobb[ic]->SetLineColor(kMagenta);
       gCprobb[ic]->SetMarkerColor(kMagenta);
       gCprobb[ic]->SetMarkerStyle(8);
       gCprobb[ic]->SetMarkerSize(1.1);
       gCprobb[ic]->Draw("Cp");*/

    gCprob[ic]->SetLineWidth(3);
    if(ic==0){
      gCprob[ic]->SetLineColor(kMagenta);
      gCprob[ic]->SetMarkerColor(kMagenta);
    }
    if(ic==1){
      gCprob[ic]->SetLineColor(kOrange);
      gCprob[ic]->SetMarkerColor(kOrange);
    }
    gCprob[ic]->SetMarkerStyle(20+ic);
    gCprob[ic]->SetMarkerSize(1.1);
    //gCprob[ic]->Draw("Cp");
     
    gCcsvSimple[ic]->SetLineWidth(3);
    gCprob[ic]->SetLineWidth(3);

    /* gCtcHighEff[ic]->SetLineWidth(3);
       gCtcHighEff[ic]->SetLineColor(kOrange);
       gCtcHighEff[ic]->SetMarkerColor(kOrange);
       gCtcHighEff[ic]->SetMarkerStyle(8);
       //gCtcHighEff[ic]->Draw("lp");
       gCtcHighPur[ic]->SetLineWidth(3);
       gCtcHighPur[ic]->SetLineColor(kRed);
       gCtcHighPur[ic]->SetMarkerColor(kRed);
       gCtcHighPur[ic]->SetMarkerStyle(8);
       //gCtcHighPur[ic]->Draw("lp");*/

    gCssvHighEff[ic]->SetLineWidth(3);
    if(ic==0){
      gCssvHighEff[ic]->SetLineColor(4);
      gCssvHighEff[ic]->SetMarkerColor(4);
    }
    else if(ic==1){
      gCssvHighEff[ic]->SetLineColor(kGreen+2);
      gCssvHighEff[ic]->SetMarkerColor(kGreen+2);
    }
    else if(ic==2){
      gCssvHighEff[ic]->SetLineColor(kMagenta+2);
      gCssvHighEff[ic]->SetMarkerColor(kMagenta+2);
    }
    gCssvHighEff[ic]->SetLineStyle(7);
    gCssvHighEff[ic]->SetMarkerStyle(20+ic);
    gCssvHighEff[ic]->SetMarkerSize(1.1);
    gCssvHighEff[ic]->Draw("Cp");  //UNCOMMENT ME FOR CHARM JETS
    gCssvHighPur[ic]->SetLineWidth(3);
    gCssvHighPur[ic]->SetLineColor(kGreen);
    gCssvHighPur[ic]->SetMarkerColor(kGreen);
    //gCssvHighPur[ic]->SetMarkerStyle(27);
    //gCssvHighPur[ic]->Draw("lp");

        
    gBcsvSimple[ic]->GetXaxis()->SetRangeUser(0.,1);
    gBcsvSimple[ic]->GetYaxis()->SetRangeUser(1e-4,1);
    
    gBcsvSimple[ic]->SetLineWidth(3);
    gBcsvSimple[ic]->SetLineColor(1);
    gBcsvSimple[ic]->SetMarkerColor(1);
    gBcsvSimple[ic]->SetMarkerSize(1.1);
    gBcsvSimple[ic]->SetMarkerStyle(8);
    //gBcsvSimple[ic]->Draw("lp");

    /*gBprobb[ic]->GetXaxis()->CenterTitle();
      gBprobb[ic]->GetYaxis()->CenterTitle();
      gBprobb[ic]->SetLineWidth(3);
      gBprobb[ic]->SetLineColor(kMagenta);
      gBprobb[ic]->SetMarkerColor(kMagenta);
      gBprobb[ic]->SetMarkerSize(1.1);
      gBprobb[ic]->SetMarkerStyle(8);
      //gBprobb[ic]->Draw("Cp");*/

    gBprob[ic]->SetLineWidth(3);
    if(ic==0){
      gBprob[ic]->SetLineColor(kMagenta);
      gBprob[ic]->SetMarkerColor(kMagenta);
    }
    if(ic==1){
      gBprob[ic]->SetLineColor(kOrange);
      gBprob[ic]->SetMarkerColor(kOrange);
    }
    gBprob[ic]->SetMarkerSize(1.1);
    gBprob[ic]->SetMarkerStyle(20+ic);
    //gBprob[ic]->Draw("Cp");
     
    gBcsvSimple[ic]->SetLineWidth(3);
    gBprob[ic]->SetLineWidth(3);

    /*gBtcHighEff[ic]->SetLineWidth(3);
      gBtcHighEff[ic]->SetLineColor(kOrange);
      gBtcHighEff[ic]->SetMarkerColor(kOrange);
      gBtcHighEff[ic]->SetMarkerSize(1.1);
      gBtcHighEff[ic]->SetMarkerStyle(8);
      //gBtcHighEff[ic]->Draw("lp");
      gBtcHighPur[ic]->SetLineWidth(3);
      gBtcHighPur[ic]->SetLineColor(kRed);
      gBtcHighPur[ic]->SetMarkerColor(kRed);
      gBtcHighPur[ic]->SetMarkerSize(1.1);
      gBtcHighPur[ic]->SetMarkerStyle(8);
      //gBtcHighPur[ic]->Draw("lp");*/
    gBssvHighEff[ic]->GetXaxis()->CenterTitle();
    gBssvHighEff[ic]->GetYaxis()->CenterTitle();
    gBssvHighEff[ic]->SetLineWidth(3);
    if(ic==0){
      gBssvHighEff[ic]->SetLineColor(4);
      gBssvHighEff[ic]->SetMarkerColor(4);
    }
    else if(ic==1){
      gBssvHighEff[ic]->SetLineColor(kGreen+2);
      gBssvHighEff[ic]->SetMarkerColor(kGreen+2);
    }
    else if(ic==2){
      gBssvHighEff[ic]->SetLineColor(kMagenta+2);
      gBssvHighEff[ic]->SetMarkerColor(kMagenta+2);
    }
    gBssvHighEff[ic]->SetMarkerStyle(20+ic);
    gBssvHighEff[ic]->SetMarkerSize(1.1);
    gBssvHighEff[ic]->Draw("Cp");

    gBssvHighPur[ic]->SetLineWidth(3);
    gBssvHighPur[ic]->SetLineColor(7);
    gBssvHighPur[ic]->SetMarkerColor(7);
    //gBssvHighPur[ic]->SetMarkerStyle(27);
    gBssvHighPur[ic]->SetMarkerSize(1.1);
    //gBssvHighPur[ic]->Draw("lp");

    float bEffWerk = 0.512;
    if(ic==0) bEffWerk = 0.467;
    if(ic==1) bEffWerk = 0.668;
    if(ppPbPb){
      if(!doCent){
	bEffWerk=0.4544;
      }
      else{
	if(ic==0) bEffWerk=0.4231;
	if(ic==1) bEffWerk=0.4436;
	if(ic==2) bEffWerk=0.4803;
	if(ic==3) bEffWerk=0.4914;	
      }
      
    }

    gBwerkingPoint[ic]=new TGraph(1);
    gBwerkingPoint[ic]->SetPoint(0,bEffWerk,gBssvHighEff[ic]->Eval(bEffWerk));


    gBwerkingPoint[ic]->SetMarkerSize(2.0);
    gBwerkingPoint[ic]->SetMarkerStyle(20+ic);
    gBwerkingPoint[ic]->SetMarkerColor(kRed+2);
    gBwerkingPoint[ic]->Draw("p");

    if(ppPbPb){
      if(!doCent){
	bEffWerk=0.4544;
      }
      else{
	if(ic==0) bEffWerk=0.4231;
	if(ic==1) bEffWerk=0.4436;
	if(ic==2) bEffWerk=0.4803;
	if(ic==3) bEffWerk=0.4914;	
      }

    }

    float cEffWerk = 0.0909;
    if(ic==1) cEffWerk = 0.173;
    gCwerkingPoint[ic]=new TGraph(1);
    gCwerkingPoint[ic]->SetPoint(0,bEffWerk,gCssvHighEff[ic]->Eval(bEffWerk));


    gCwerkingPoint[ic]->SetMarkerSize(2.0);
    gCwerkingPoint[ic]->SetMarkerStyle(20+ic);
    gCwerkingPoint[ic]->SetMarkerColor(kRed+2);
    gCwerkingPoint[ic]->Draw("p");
  }

  prel->Draw();
  roots->Draw();  
    
  //c2->cd(2);

  TLatex *lat1 = new TLatex(0.64,1.64,"Pythia (+Hydjet)");
  lat1->SetTextSize(22);
  lat1->SetTextFont(43);
  lat1->SetTextColor(kBlack);
  lat1->Draw("same");
 
  std::string collType;
  std::string taggerType;
  std::string probType;
  TLegend *leg2 = new TLegend(0.165,0.628,0.549,0.829);
  leg2->SetHeader("Tagging efficiency");
  leg2->SetFillColor(0);
  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0);
  for(int ic=0; ic<nLoop; ic++){
    if(ic==0) collType="PbPb";
    if(ic==1) collType="pp";
    if(ic==2) collType="pPb";
    taggerType = "SV udsg jets, "+collType;
    probType = "JP udsg jets, "+collType;
    leg2->AddEntry(gBssvHighEff[ic], taggerType.c_str(),"lp");
    //leg2->AddEntry(gBprob[ic], probType.c_str(), "lp");
    taggerType = "SV charm jets, "+collType;
    probType = "JP charm jets, "+collType;
    leg2->AddEntry(gCssvHighEff[ic], taggerType.c_str(),"lp");
    //leg2->AddEntry(gCprob[ic], probType.c_str(),"lp");
  }
  leg2->Draw("same");
  //TLatex *prel = new TLatex(0.5,1.1,"CMS simulation");
  //prel->Draw();
  TFile *fout=new TFile("temp.root","recreate");
  fout->cd();
  for(int ic=0;ic<centMax;ic++){
    gCssvHighEff[ic]->SetName("gCssvHighEff_0");
    gCssvHighEff[ic]->Write();
  }

  std::string typeLabel;

  if(doCent==1 && ppPbPb==1) typeLabel = "CentDep";
  if(ppPbPb==1 && doCent==0) typeLabel = "PbPb";
  if(ppPbPb==0 && doCent==0) typeLabel = "pp";
  //cout << "typeLabel: " << typeLabel << endl;
  if(ppPbPb==0 && doCent==1) typeLabel = "dontTrustMe";
  if(doBoth==1) typeLabel = "Both";


  // c->SaveAs(Form("bVsL_%s.pdf",typeLabel.c_str()));
  c2->SaveAs(Form("bVsX_%s.pdf",typeLabel.c_str()));
  // c->SaveAs(Form("bVsL_%s.C",typeLabel.c_str()));
  c2->SaveAs(Form("bVsX_%s.C",typeLabel.c_str()));

}
