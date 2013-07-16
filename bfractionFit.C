#include <TFile.h>
#include <TH1.h>
#include <TTree.h>
#include <TCanvas.h>

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

void fixEmpty(TH1 *h)
{
   for (int i=1;i<=h->GetNbinsX();i++)
   {
      if (h->GetBinContent(i)==0) h->SetBinContent(i,1e-20);
   }
}

RooRealVar bfractionFit(char *var = "discr_csvSimple",double minX = 0,double maxX = 1,double ptMin = 60, double ptMax = 500, bool floatC = 0)
{
   // Prepare a canvas
   TCanvas *c = new TCanvas("c","",600,600);
   c->SetLogy();
   

   // MC shape file
   TFile *infQCD = new TFile("histos/ppMC_hiReco_jetTrigQCD.root");
   TTree *tQCD = (TTree*) infQCD->Get("nt");
   TFile *infB = new TFile("histos/ppMC_hiReco_jetTrigB.root");
   TTree *tB = (TTree*) infB->Get("nt");
   TFile *infC = new TFile("histos/ppMC_hiReco_jetTrigC.root");
   TTree *tC = (TTree*) infC->Get("nt");

   // b-jet signal shape
   TH1D *hB = new TH1D("hB","",20,minX,maxX);
   hB->Sumw2();
   tB->Draw(Form("%s>>hB",var),Form("weight*(abs(refparton_flavorForB)==5&&jtpt>%f&&jtpt<%f&&pthat>80)",ptMin,ptMax));
   fixEmpty(hB);

   TH1D *hC = new TH1D("hC","",20,minX,maxX);
   hC->Sumw2();
   tC->Draw(Form("%s>>hC",var),Form("weight*(abs(refparton_flavorForB)==4&&jtpt>%f&&jtpt<%f&&pthat>80)",ptMin,ptMax));
   fixEmpty(hC);
   
   // b-jet background shape
   TH1D *hOtherFlavor = new TH1D("hOtherFlavor","",20,minX,maxX);
   TH1D *hQCDC = new TH1D("hQCDC","",20,minX,maxX);
   hOtherFlavor->Sumw2();
   hC->Sumw2();
   tQCD->Draw(Form("%s>>hOtherFlavor",var),Form("weight*(abs(refparton_flavorForB)!=5&&abs(refparton_flavorForB)!=4&&jtpt>%f&&jtpt<%f&&pthat>80)",ptMin,ptMax));
   tQCD->Draw(Form("%s>>hQCDC",var),Form("weight*(abs(refparton_flavorForB)==4&&jtpt>%f&&jtpt<%f&&pthat>80)",ptMin,ptMax));   fixEmpty(hOtherFlavor);
   
   if (!floatC) {
      hC->Scale(1./hC->Integral(0,100)*hQCDC->Integral(0,100));
      hOtherFlavor->Add(hC);
      cout <<"fracC = " <<hC->Integral(0,100)/hOtherFlavor->Integral(0,100)<<endl;
   }
   // data sample   
   TFile *infData = new TFile("histos/ppdata_hiReco_jetTrig.root");
   TTree *tData = (TTree*) infData->Get("nt");
   
   // --- Observable ---
   RooRealVar s(var,var,0,minX,maxX);
   RooRealVar jtpt("jtpt","jtpt",0,ptMin,ptMax);
 
   // --- Parameters ---
 
   // --- Build Histogram PDF ---
   RooDataHist xB("xB","xB",s,hB);
   RooDataHist xC("xC","xC",s,hC); // not used if doesn't float C
   RooHistPdf signal("signal","signal PDF",s,xB);
   RooHistPdf signalC("signalC","signalC PDF",s,xC); // not used if doesn't float C

   RooDataHist xOtherFlavor("xOtherFlavor","xOtherFlavor",s,hOtherFlavor);
   RooHistPdf background("background","Background PDF",s,xOtherFlavor);

   // --- Construct signal+background PDF ---
   // signal = bjets
   // background = all the other flavors
   //   RooRealVar nsig("nsig","#signal events",6000,0.,1e7) ;
   //   RooRealVar nbkg("nbkg","#background events",1e5,0.,1e7) ;
   //   RooAddPdf model("model","g+a",RooArgList(signal,background),RooArgList(nsig,nbkg)) ;
   RooRealVar frac("frac","#background events",0.1,0.,1) ;
   RooRealVar fracC("fracC","#background events",0.1,0.,1) ; // not sured if doesn't float C
   RooAddPdf modelBck("modelBck","g+a",signalC,background,fracC) ;
   RooAddPdf *model;
   if (!floatC) {
     model = new RooAddPdf("model","g+a",signal,background,frac) ;
   } else {
     model = new RooAddPdf("model","g+a",signal,modelBck,frac) ;
   }

   // data sample
   //RooDataSet *data = new RooDataSet("data","data",RooArgSet(s),Import(*tData));
   RooDataSet *data = new RooDataSet("data","data",tData,RooArgSet(s,jtpt),Form("jtpt>%f&&jtpt<%f",ptMin,ptMax));
   
   // --- Perform extended ML fit of composite PDF to data ---
   model->fitTo(*data) ;
   // --- Plot data and composite PDF overlaid ---
   RooPlot* sframe = s.frame() ;
   TH2D *htemp = new TH2D("htemp","",100,minX,maxX,100,1,1e6) ;
   htemp->SetXTitle(Form("%s %.0f < p_{T} < %.0f GeV/c",var,ptMin,ptMax));
   htemp->SetYTitle("Entries");
   htemp->Draw();
   cout <<"Min "<<sframe->GetMinimum()<<endl;
   data->plotOn(sframe,Binning(20)) ;
   sframe->SetTitle("");
   model->plotOn(sframe,Components(background),LineStyle(kDashed),LineColor(kBlue)) ;   
   model->plotOn(sframe,Components(signalC),LineStyle(kDashed),LineColor(kGreen)) ;   
   model->plotOn(sframe) ;
   model->plotOn(sframe,Components(signal),LineStyle(kDashed),LineColor(kRed),FillColor(kRed),FillStyle(1)) ;   

   cout <<"b jet fraction = "<<frac.getVal()<<endl;
   sframe->Draw("same");
   c->SaveAs(Form("fit/%s-%.0f-%.0f.gif",var,ptMin,ptMax));
   delete model;
   return frac;
}

void ptDependence()
{
   TFile *inf = new TFile("histos/ppMC.root");
   TTree *t = (TTree*) inf->Get("nt");

   const int nBins = 4;
   double ptBin[nBins+1] = {100,120,140,160,200};
//   const int nBins = 1;
//   double ptBin[nBins+1] = {100,400};
   
   TH1D *hProb = new TH1D("hProb","",nBins,ptBin);
   TH1D *hCSV = new TH1D("hCSV","",nBins,ptBin);
   TH1D *hSVTXM = new TH1D("hSVTXM","",nBins,ptBin);
   TProfile *pGen = new TProfile("pGen","",nBins,ptBin);
   
   for (int n=0; n<nBins;n++)
   {
      RooRealVar f1 = bfractionFit("discr_prob",0,3.5,ptBin[n],ptBin[n+1]);
      RooRealVar f2 = bfractionFit("discr_csvSimple",0,1,ptBin[n],ptBin[n+1]);
      RooRealVar f3 = bfractionFit("svtxm",0,6,ptBin[n],ptBin[n+1]);
      hProb->SetBinContent(n+1,f1.getVal());    
      hProb->SetBinError(n+1,f1.getError());    
      hCSV->SetBinContent(n+1,f2.getVal());    
      hCSV->SetBinError(n+1,f2.getError());    
      hSVTXM->SetBinContent(n+1,f3.getVal());    
      hSVTXM->SetBinError(n+1,f3.getError());    
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
   hSVTXM->SetLineColor(4);
   hSVTXM->SetMarkerColor(4);
   hSVTXM->SetMarkerStyle(24);
//   hSVTXM->Draw("same");
   t->Draw("abs(refparton_flavorForB)==5:jtpt","","prof same");
   
   TLegend *leg = new TLegend(0.2,0.7,0.5,0.9);
   leg->SetBorderSize(0);
   leg->SetFillStyle(0);
   leg->SetFillColor(0);
   leg->AddEntry(hProb,"Jet Probability","pl");
   leg->AddEntry(hCSV,"CSV","pl");
//   leg->AddEntry(hSVTXM,"SV mass","pl");
   leg->Draw();
}
