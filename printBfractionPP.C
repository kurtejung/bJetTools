#include <iostream.h>
#include <math.h>
#include "TROOT.h"
#include "TStyle.h"
#include "TH1.h"
#include "TF1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"



Double_t abs (Double_t a);
Double_t max (Double_t a, Double_t b);
Double_t norm (Double_t a, Double_t b);
Double_t norm (Double_t a, Double_t b, Double_t c);
Double_t norm (Double_t a, Double_t b, Double_t c, Double_t d, Double_t e);

void correct2(TH1* h);
void setMeanPt(TGraphAsymmErrors *g, TH1F *h, int isData);

//void printBfraction(char *tagger="discr_ssvHighPur", Double_t workingPoint=2, char *taggerName="SSVHP") {
void printBfractionPP(char *tagger="discr_ssvHighEff", Double_t workingPoint=2, char *taggerName="SSVHE", int saveOutput=0) {


  gROOT->ForceStyle(1);

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  gStyle->SetTextFont(43);
  gStyle->SetLabelFont(43,"XYZ");
  gStyle->SetTitleFont(43,"XYZ");

  // hack
  gStyle->SetLabelSize(28,"xy");
  gStyle->SetTitleSize(28,"xy");
  gStyle->SetTitleOffset(1.1,"x");
  gStyle->SetTitleOffset(1.5,"y");
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadTopMargin(0.07);
  gStyle->SetNdivisions(408,"y");
  gROOT->ForceStyle(1);

  TFile *fin = new TFile(Form("./bFractionpp290512/histos/bFraction_%sat%.1f.root",taggerName,workingPoint));
  //TFile *fin = new TFile(Form("./histos/bFraction_regPFforJets_%sat%.1f.root",taggerName,workingPoint));
  //TFile *fin = new TFile("output/bFractionMerged_ppPbPb0_ssvHighEffat2.0FixCL1_bin_0_40_eta_0_2.root");

  TH1F *hBFractionMC = (TH1F*) fin->Get("hBFractionMC");
  TH1F *hBFractionMCRefLevel = (TH1F*) fin->Get("hBFractionMCRefLevel");

  TH1F *hBFractionData = (TH1F*) fin->Get(Form("hBFractionData_%sat%.1f_CLshift%.0f",taggerName,workingPoint,0.));
  TH1F *hBFractionDataMoreC = (TH1F*) fin->Get(Form("hBFractionData_%sat%.1f_CLshift%.0f",taggerName,workingPoint,20.));
  TH1F *hBFractionDataLessC = (TH1F*) fin->Get(Form("hBFractionData_%sat%.1f_CLshift%.0f",taggerName,workingPoint,-20.));

  TH1F *hBFractionDataLTJP = (TH1F*) fin->Get(Form("hBFractionDataLTJP_%sat%.1f_CLshift%.0f",taggerName,workingPoint,0.));
  TH1F *hBFractionDataLTJPMoreC = (TH1F*) fin->Get(Form("hBFractionDataLTJP_%sat%.1f_CLshift%.0f",taggerName,workingPoint,20.));
  TH1F *hBFractionDataLTJPLessC = (TH1F*) fin->Get(Form("hBFractionDataLTJP_%sat%.1f_CLshift%.0f",taggerName,workingPoint,-20.));

  TH1F *hBFractionJPdirect = (TH1F*) fin->Get(Form("hBFractionJPdirect_%sat%.1f_CLshift%.0f",taggerName,workingPoint,0.));
  TH1F *hBFractionJPdirectMoreC = (TH1F*) fin->Get(Form("hBFractionJPdirect_%sat%.1f_CLshift%.0f",taggerName,workingPoint,20.));
  TH1F *hBFractionJPdirectLessC = (TH1F*) fin->Get(Form("hBFractionJPdirect_%sat%.1f_CLshift%.0f",taggerName,workingPoint,-20.));

  /*  --- correction due to Jet Energy Scale (calcul) ---
  correct(hBFractionMC);
  correct(hBFractionData);
  correct(hBFractionDataMoreC);
  correct(hBFractionDataLessC);
  correct(hBFractionDataLTJP);
  correct(hBFractionDataLTJPMoreC);
  correct(hBFractionDataLTJPLessC);
  correct(hBFractionJPdirect);
  correct(hBFractionJPdirectMoreC);
  correct(hBFractionJPdirectLessC);
  //*/

  ///*  --- correction due to Jet Energy Scale (by hand) ---
  correct2(hBFractionMC);
  correct2(hBFractionData);
  correct2(hBFractionDataMoreC);
  correct2(hBFractionDataLessC);
  correct2(hBFractionDataLTJP);
  correct2(hBFractionDataLTJPMoreC);
  correct2(hBFractionDataLTJPLessC);
  correct2(hBFractionJPdirect);
  correct2(hBFractionJPdirectMoreC);
  correct2(hBFractionJPdirectLessC);
  //*/

  //  --- plots with variation of charm ---

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TCanvas *cBFraction1 = new TCanvas("cBFraction1","b-jet fraction",700,600);
  hBFractionMC->SetLineColor(2);
  hBFractionMC->SetLineWidth(2);
  hBFractionMC->SetMarkerColor(2);
  hBFractionMC->SetAxisRange(0,0.06,"Y");
  hBFractionMC->SetTitleOffset(1,"X");
  hBFractionMC->GetYaxis()->SetTitle("b-jet fraction");;

  TH1F *hBFractionMC_dummy= (TH1F*)hBFractionMC->Clone("hBFractionMC_dummy");
  hBFractionMC_dummy->SetLineWidth(0);
  hBFractionMC_dummy->SetMarkerSize(0);
  hBFractionMC_dummy->SetMarkerColor(0);
  hBFractionMC_dummy->SetLineColor(0);
  hBFractionMC_dummy->Draw(); 
  hBFractionMC->Draw("hist,same"); 

  //  hBFractionMC->Draw("e1"); 

  hBFractionData->SetMarkerStyle(8);
  hBFractionJPdirect->SetMarkerStyle(8);
  hBFractionDataLTJP->SetMarkerStyle(8);
  hBFractionDataLTJPMoreC->SetMarkerStyle(22);
  hBFractionDataLTJPLessC->SetMarkerStyle(23);

  hBFractionDataLTJP->SetLineColor(kGreen-2); 
  hBFractionDataLTJP->SetMarkerColor(kGreen-2); 
  hBFractionJPdirect->SetMarkerColor(kBlue);
  hBFractionJPdirect->SetLineColor(kBlue);
  hBFractionDataLTJPMoreC->SetLineColor(kGreen-2); 
  hBFractionDataLTJPMoreC->SetMarkerColor(kGreen-2); 
  hBFractionDataLTJPLessC->SetLineColor(kGreen-2); 
  hBFractionDataLTJPLessC->SetMarkerColor(kGreen-2); 


  hBFractionData->Draw("e1same");   
  hBFractionDataLTJP->Draw("e1same");
  hBFractionJPdirect->Draw("e1same");

  hBFractionDataMoreC->SetLineStyle(2);
  //hBFractionDataMoreC->Draw("e1same");   
  hBFractionDataLTJPMoreC->SetLineStyle(2);
  hBFractionDataLTJPMoreC->Draw("e1same");
  hBFractionJPdirectMoreC->SetLineStyle(2);
  //hBFractionJPdirectMoreC->Draw("e1same");

  hBFractionDataLessC->SetLineStyle(3);
  //hBFractionDataLessC->Draw("e1same");   
  hBFractionDataLTJPLessC->SetLineStyle(3);
  hBFractionDataLTJPLessC->Draw("e1same");
  hBFractionJPdirectLessC->SetLineStyle(3);
  //hBFractionJPdirectLessC->Draw("e1same");
  


  TLegend *legFrac1 = new TLegend(0.15,0.65,0.87,0.95);
  legFrac1->SetBorderSize(0);
  legFrac1->SetFillStyle(0);
  legFrac1->SetTextSize(23);
  legFrac1->SetHeader("pp, #sqrt{s} = 2.76 TeV");
  legFrac1->AddEntry(hBFractionDataLTJP,"SSVHE, LT method","pl");
  legFrac1->AddEntry(hBFractionDataLTJPMoreC,"SSVHE, LT method, Charm * 1.2","pl");
  legFrac1->AddEntry(hBFractionDataLTJPLessC,"SSVHE, LT method, Charm * 0.8","pl");
  legFrac1->AddEntry(hBFractionData,"SSVHE, MC eff.","pl");
  legFrac1->AddEntry(hBFractionJPdirect,"Jet Probability","pl");
  //legFrac1->Draw();


  //  --- plots of LT method with syst. uncertainty ---

  TCanvas *cBFraction2 = new TCanvas("cBFraction2","b-jet fraction",700,600);

  //*
  TH1F *hBFractionMC2 = (TH1F *)hBFractionMC->Clone("hBFractionMC2");
  hBFractionMC2->GetXaxis()->SetRangeUser(80,200);
  hBFractionMC2->SetMarkerSize(0);
  hBFractionMC2->SetMaximum(0.06);
  hBFractionMC2->SetMinimum(0.0);
  hBFractionMC2->GetXaxis()->CenterTitle();
  hBFractionMC2->GetYaxis()->CenterTitle();
  hBFractionMC2->Draw("hist");

  TGraphAsymmErrors *gBFractionMC2 = new TGraphAsymmErrors(hBFractionMC);
  setMeanPt(gBFractionMC2,hBFractionMC,0);
  gBFractionMC2->GetXaxis()->SetRangeUser(80,200);
  
  TLatex *prel = new TLatex(81,0.0615,"CMS Preliminary");
  prel->SetTextSize(23);
  prel->Draw();

  TLatex *roots = new TLatex(167,0.0615,"#sqrt{s} = 2.76 TeV");
  roots->SetTextSize(23);
  roots->Draw();

  //TLatex *ptlabel = new TLatex(90,0.0075,"Jet p_{T} > 80 GeV/c");
  //TLatex *ptlabel = new TLatex(90,0.055,"Jet p_{T} > 80 GeV/c");
  //ptlabel->Draw();

  //*/
  /*
  TH1F *hBFractionMCRefLevel2 = (TH1F *)hBFractionMCRefLevel->Clone();
  hBFractionMCRefLevel2->SetAxisRange(0,0.05,"Y");
  hBFractionMCRefLevel2->SetTitleOffset(1,"Y");
  hBFractionMCRefLevel2->SetLineColor(2);
  hBFractionMCRefLevel2->SetMarkerColor(2);
  hBFractionMCRefLevel2->SetMarkerStyle(21);
  hBFractionMCRefLevel2->Draw();
  //*/


  TGraphErrors *gSyst = new TGraphErrors(3);
  TGraphErrors *gSystCL = new TGraphErrors(3);
  TGraphErrors *gSystMeth = new TGraphErrors(3);
  TGraphErrors *gSystJES = new TGraphErrors(3);
  TGraphErrors *gSystTemp = new TGraphErrors(3);
  TGraphErrors *gSystWP = new TGraphErrors(3);

  gSyst->SetName("gSyst");
  gSystCL->SetName("gSystCL");
  gSystMeth->SetName("gSystMeth");
  gSystJES->SetName("gSystJES");
  gSystTemp->SetName("gSystTemp");
  gSystWP->SetName("gSystWP");


  Double_t errCLratio, errMethod, totalSystErr;

  //float binCentData[3]={87.88,108.2,140.6};
  //float binBound[4]={80.,100.,120.,200.};

  // hard code purity systematics
  //https://dl.dropbox.com/u/48183843/relativeShift_HIN12003.pdf
  double errDataMCtemplate[3]= {0.36, 0.11, 0.06};
  double errWorkingPoint[3]={0.18, 0.18, 0.22};


  for(Int_t i=1;i<=hBFractionDataLTJP->GetNbinsX();i++) {

    float xVal = hBFractionDataLTJP->GetBinCenter(i);

    gSyst->SetPoint(i-1,xVal,hBFractionDataLTJP->GetBinContent(i));
    gSystCL->SetPoint(i-1,xVal,hBFractionDataLTJP->GetBinContent(i));
    gSystMeth->SetPoint(i-1,xVal,hBFractionDataLTJP->GetBinContent(i));
    gSystJES->SetPoint(i-1,xVal,hBFractionDataLTJP->GetBinContent(i));
    gSystTemp->SetPoint(i-1,xVal,hBFractionDataLTJP->GetBinContent(i));
    gSystWP->SetPoint(i-1,xVal,hBFractionDataLTJP->GetBinContent(i));


//     errCLratio = max(abs(hBFractionDataLTJP->GetBinContent(i)-hBFractionDataLTJPMoreC->GetBinContent(i)),abs(hBFractionDataLTJP->GetBinContent(i)-hBFractionDataLTJPLessC->GetBinContent(i)));
//     errMethod = max(abs(hBFractionDataLTJP->GetBinContent(i)-hBFractionData->GetBinContent(i)),abs(hBFractionDataLTJP->GetBinContent(i)-hBFractionJPdirect->GetBinContent(i)));

//     double errJES = 0.07*hBFractionDataLTJP->GetBinContent(i);
//     float val = hBFractionDataLTJP->GetBinContent(i);

//     totalSystErr = norm(errCLratio,errMethod,errJES);

//     gSyst->SetPointError(i-1,hBFractionDataLTJP->GetBinWidth(i)/2,totalSystErr);
//     gSystCL->SetPointError(i-1,hBFractionDataLTJP->GetBinWidth(i)/2,errCLratio);
//     gSystMeth->SetPointError(i-1,hBFractionDataLTJP->GetBinWidth(i)/2,errMethod);
//     gSystJES->SetPointError(i-1,hBFractionDataLTJP->GetBinWidth(i)/2,errJES);

    cout<<" central value "<<hBFractionDataLTJP->GetBinContent(i)<<endl;
    errCLratio = max(abs(hBFractionDataLTJP->GetBinContent(i)-hBFractionDataLTJPMoreC->GetBinContent(i)),abs(hBFractionDataLTJP->GetBinContent(i)-hBFractionDataLTJPLessC->GetBinContent(i)));
    errMethod  = abs(hBFractionDataLTJP->GetBinContent(i)-hBFractionData->GetBinContent(i));
    double errJES = 0.07*hBFractionDataLTJP->GetBinContent(i);
    float val     = hBFractionDataLTJP->GetBinContent(i);
    totalSystErr = norm(errCLratio,errMethod,errJES,val*errDataMCtemplate[i-1],val*errWorkingPoint[i-1]);

    float xErr = hBFractionDataLTJP->GetBinWidth(i)/2;
    gSyst->SetPointError(i-1,xErr,totalSystErr);
    gSystCL->SetPointError(i-1,xErr,errCLratio);
    gSystMeth->SetPointError(i-1,xErr,errMethod);
    gSystJES->SetPointError(i-1,xErr,errJES);
    gSystTemp->SetPointError(i-1,xErr,val*errDataMCtemplate[i-1]);
    gSystWP->SetPointError(i-1,xErr,val*errWorkingPoint[i-1]);

    cout<<" rel METH "<<errMethod/val<<endl;
    cout<<" rel CL "<<errCLratio/val<<endl;
    cout<<" rel JES "<<errJES/val<<endl;
    cout<<" rel sys error "<<totalSystErr/val<<endl;
    cout<<" absolute sys error "<<totalSystErr<<endl;

    // add in MC template uncertainties
    float origStatErr = hBFractionDataLTJP->GetBinError(i);
    int statBin=i-1;


    //gSyst->SetPointError(i-1,abs(binCentData[i-1]-binBound[i-1]),abs(binCentData[i-1]-binBound[i]),totalSystErr,totalSystErr);

  }


  gSyst->SetFillColor(5);
  gSyst->Draw("2");

  gBFractionMC2->SetMarkerSize(0);
  gBFractionMC2->Draw("Z,p,same");
  hBFractionMC2->Draw("hist,same");

  TGraphAsymmErrors *gBFractionDataLTJP2 = new TGraphAsymmErrors(hBFractionDataLTJP);

  setMeanPt(gBFractionDataLTJP2,hBFractionDataLTJP,1);
  gBFractionDataLTJP2->SetLineColor(1);
  gBFractionDataLTJP2->SetMarkerColor(1);
  gBFractionDataLTJP2->SetMarkerSize(1.5);
  gBFractionDataLTJP2->Draw("p,e1,same");


  //hBFractionMCRefLevel2->Draw("e1same");
  
  TLegend *legFrac2 = new TLegend(0.3,0.15,0.8,0.34);
  legFrac2->SetHeader("#int L dt = 231 nb^{-1}");
  legFrac2->SetBorderSize(0);
  legFrac2->SetFillStyle(0);
  legFrac2->SetTextSize(23);
  legFrac2->AddEntry(gBFractionDataLTJP2,"pp Data","p");
  legFrac2->AddEntry(hBFractionMC2,"Pythia","l");
  legFrac2->AddEntry(gSyst,"Exp. uncertainty","f");
  //legFrac2->Draw();



  cBFraction2->RedrawAxis();
  if(saveOutput){
    cBFraction2->Print("bFractionPPVsPt.pdf","pdf");
    cBFraction2->Print("bFractionPPVsPt.gif","gif");
    
    TFile *fout = new TFile("pp_bFraction.root","recreate");
    gBFractionDataLTJP2->Write();
    gSyst->Write();
    gSystCL->Write();
    gSystMeth->Write();
    gSystJES->Write();
    gSystTemp->Write();
    gSystWP->Write();
    gBFractionMC2->Write();
    fout->Close();
  }
}


Double_t abs (Double_t a) {return (a>0)?a:-a;}
Double_t max (Double_t a, Double_t b) {return (a>b)?a:b;}
Double_t norm (Double_t a, Double_t b) {return sqrt(a*a+b*b);}
Double_t norm (Double_t a, Double_t b, Double_t c) {return sqrt(a*a+b*b+c*c);}
Double_t norm (Double_t a, Double_t b, Double_t c, Double_t d, Double_t e) {return sqrt(a*a+b*b+c*c+d*d+e*e);}



Double_t NtrueOverNmeas (Double_t ptMin, Double_t ptMax) {

  //Double_t n = 6;
  //Double_t n = 5.78121e+00; // from ref level
  Double_t n = 5.95950e+00; // from reco level
  Double_t c1 = 5.23359e-01;
  Double_t c2 = 1.15226e-01;
  Double_t c3 = -7.96949e-03;

  TF1 *fTrue = new TF1("fMeas","([1]+[2]*log(x)+[3]*log(x)*log(x)+[2]+2*[3]*log(x))/exp([0]*log(x*([1]+[2]*log(x)+[3]*log(x)*log(x))))",0,500); 
  fTrue->SetParameters(n,c1,c2,c3);
  TF1 *fMeas = new TF1("fTrue","1/exp([0]*log(x))",0,500);
  fMeas->SetParameter(0,n);

  return fTrue->Integral(ptMin,ptMax)/fMeas->Integral(ptMin,ptMax);

}

void correct(TH1* h) {

  Double_t coef;

  for (Int_t i=1;i<=h->GetNbinsX();i++) {

    coef = NtrueOverNmeas(h->GetBinLowEdge(i),h->GetBinLowEdge(i)+h->GetBinWidth(i));

    cout<<"  correction for bin "<<i<<" is "<<coef<<endl;

    h->SetBinContent(i,coef*h->GetBinContent(i));
    h->SetBinError(i,coef*h->GetBinError(i));

  }

}

void correct2(TH1* h) {
  /* // Simon's original values
  h->SetBinContent(1,1.88055*h->GetBinContent(1));
  h->SetBinContent(2,1.84526*h->GetBinContent(2));
  h->SetBinContent(3,1.80878*h->GetBinContent(3));
  h->SetBinContent(4,1.75974*h->GetBinContent(4));
  h->SetBinContent(5,1.68418*h->GetBinContent(5));
  h->SetBinError(1,1.88055*h->GetBinError(1));
  h->SetBinError(2,1.84526*h->GetBinError(2));
  h->SetBinError(3,1.80878*h->GetBinError(3));
  h->SetBinError(4,1.75974*h->GetBinError(4));
  h->SetBinError(5,1.68418*h->GetBinError(5));
  */

  //updated with n = 5.95

  //corrs[0]=1.84564; corrs[1]=1.78241; corrs[2]=1.72693; corrs[3]=1.65839;


  h->SetBinContent(1,1.84564*h->GetBinContent(1));
  h->SetBinContent(2,1.78241*h->GetBinContent(2));
  h->SetBinContent(3,1.70204*h->GetBinContent(3));
  h->SetBinError(1,1.84564*h->GetBinError(1));
  h->SetBinError(2,1.78241*h->GetBinError(2));
  h->SetBinError(3,1.70204*h->GetBinError(3));

  


}

void setMeanPt(TGraphAsymmErrors *g, TH1F *h, int isData){

  float meanPtData[3]={87.88,108.2,140.6};
  float meanPtMC[3]={87.83,108.2,140.1};
  
  
  for(int i=0;i<h->GetNbinsX();i++){
    if(isData){
      g->SetPoint(i,meanPtData[i],h->GetBinContent(i+1));
      g->SetPointError(i,0,0,h->GetBinError(i+1), h->GetBinError(i+1));

    }
    else{
      g->SetPoint(i,meanPtMC[i],h->GetBinContent(i+1));
      g->SetPointError(i,meanPtMC[i]-h->GetBinLowEdge(i+1),h->GetBinLowEdge(i+1)+h->GetBinWidth(i+1)-meanPtMC[i],h->GetBinError(i+1), h->GetBinError(i+1));
    }

  }
  
}
