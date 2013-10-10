#include <iostream>
#include "TFile.h"
#include "TH1.h"
#include "TNtuple.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TPad.h"
#include "TPaveStats.h"
#include "TCanvas.h"
#include "TString.h"

void plotDijetAsymm(int pt1norm=1, string disc = "discr_prob", float discCut1=0.6, float discCut2=0.6, int applyEffCorr=1, int drawStats=1, int plotPtRatio=0, int doSub=1, bool ppB=1){

  int tagLead = 0;
  int tagSub = 0;

  if(discCut1>0.0) tagLead=1;
  if(discCut2>0.0) tagSub=1;

  cout<<"cutting on "<<disc.c_str()<<" > "<<discCut1<<" / "<<discCut2<<endl;

  gROOT->Reset();
  gStyle->SetOptTitle(0);
  if(drawStats) gStyle->SetOptStat(1002210);
  else gStyle->SetOptTitle(0);
  gStyle->SetNdivisions(505,"x");
  //gROOT->ForceStyle(1);

  TFile *fB = new TFile("histos/ppMC_ppReco_Dijet_BjetTrig.root");
  TFile *fC = new TFile("histos/ppMC_ppReco_Dijet_CjetTrig.root");
  TFile *fQCD = new TFile("histos/ppMC_ppReco_Dijet_QCDjetTrig.root");
  TFile *fdata = new TFile("histos/ppdata_ppReco_Dijet_jetTrig.root");

  TNtuple *ntB = (TNtuple *) fB->Get("nt");
  TNtuple *ntC = (TNtuple *) fC->Get("nt");
  TNtuple *ntQCD = (TNtuple *) fQCD->Get("nt");
  TNtuple *ntdata = (TNtuple *) fdata->Get("nt");

  TH1F *hPt1B=new TH1F("hPt1B","hPt1B",120,65.,200.);
  TH1F *hPt1C=new TH1F("hPt1C","hPt1C",120,65.,200.);
  TH1F *hPt1QCD=new TH1F("hPt1QCD","hPt1QCD",120,65.,200.);
  TH1F *hPt1Data=new TH1F("hPt1Data","hPt1Data",120,65.,200.);



  float a, b, c;
  float d, e, f;
  if(disc=="ssvHE"){
    a = 4.67356e-01;
    b =  -9.90006e-02;
    c = 9.44735e-03;
  }
  else if(disc=="jetProb"){
    a = 4.82642e-01;
    b = -2.29807e-01;
    c =  7.32250e-03;
    // test inclusive corrections
    d = 4.42378e-01;
    e = -1.06758e-01;
    f =  7.02023e-03;
  }
  else if(ppB){
    a = 6.66338e-01;
    b = -1.37254e-01;
    c =  7.96509e-03;
    d = 6.67839e-01;
    e = -1.20443e-01;
    f =  8.42849e-03;
  }
  else cout<<" no disc "<<endl;

  char  *effCorr1 = Form("%f+%f*pow(log(jtpt1*%f),2)",a,b,c);
  char  *effCorr2 = Form("%f+%f*pow(log(jtpt2*%f),2)",a,b,c);
  char  *dataEffCorr1 = Form("%f+%f*pow(log(jtpt1*%f),2)",d,e,f);
  char  *dataEffCorr2 = Form("%f+%f*pow(log(jtpt2*%f),2)",d,e,f);
  
  if(!applyEffCorr){
    effCorr1 = Form("1.");
    effCorr2 = Form("1.");
    dataEffCorr1 = Form("1.");
    dataEffCorr2 = Form("1.");
  }
  else {
    if(!tagLead){
      effCorr1 = Form("1.");
      dataEffCorr1 = Form("1.");
    }
    if(!tagSub){
      effCorr2 = Form("1.");
      dataEffCorr2 = Form("1.");
    }
  }

  cout<<dataEffCorr1<<endl;

  float epsilon = 1e-4;
  if(applyEffCorr){
    if((disc=="discr_prob" && fabs(discCut1-0.6) < epsilon) || (disc=="ssvHE" && fabs(discCut1-2.0) < epsilon) ){
      cout<<" applying efficiency corrections "<<endl;
    }
    else if((disc=="discr_prob" && fabs(discCut2-0.6) < epsilon) || (disc=="ssvHE" && fabs(discCut2-2.0) < epsilon) ){
      cout<<" applying efficiency corrections "<<endl;
    }
    else{
      cout<<disc<<endl;
      cout<<discCut1<<endl;
      cout<<discCut2<<endl;
      cout<< "efficiency not evaluated for this choice of discriminator cuts "<<endl;
      cout<<" try w/o tagging efficiency correction "<<endl;
      return;
    }
  }
  
  char *pt1Choice = Form("jtpt1");
  char *pt2Choice = Form("jtpt2");

  if(tagLead) pt1Choice = Form("jtpt1");
  if(tagSub) pt2Choice = Form("jtpt2");

  if(pt1norm){
    ntB->Draw(Form("%s>>hPt1B",pt1Choice),Form("weight/(%s)*(pthat>50&&%s>70&&%s<200&&abs(refparton_flavorForB1)==5&&%s1>=%f&&refpt1>0)",effCorr1,pt1Choice,pt1Choice,disc.c_str(),discCut1));
    ntC->Draw(Form("%s>>hPt1C",pt1Choice),Form("weight/(%s)*(pthat>50&&%s>70&&%s<200&&abs(refparton_flavorForB1)==4&&%s1>=%f&&refpt1>0)",effCorr1,pt1Choice,pt1Choice,disc.c_str(),discCut1));
    ntQCD->Draw(Form("%s>>hPt1QCD",pt1Choice),Form("weight/(%s)*(pthat>50&&%s>70&&%s<200&&abs(refparton_flavorForB1)!=5&&abs(refparton_flavorForB1)!=4&&%s1>=%f&&refpt1>0)",effCorr1,pt1Choice,pt1Choice,disc.c_str(),discCut1));    
    //ntdata->Draw(Form("%s>>hPt1Data",pt1Choice),Form("1./(%s)*(pthat>50&&%s>70&&%s<200&&%s1>=%f)",effCorr1,pt1Choice,pt1Choice,disc.c_str(),discCut1));    
    ntdata->Draw(Form("%s>>hPt1Data",pt1Choice),Form("weight/(%s)*(%s>70&&%s<200&&%s1>=%f)",effCorr1,pt1Choice,pt1Choice,disc.c_str(),discCut1));    
  }
  
  TH1F *hPt1MC = (TH1F*)hPt1B->Clone("hPt1MC");
  hPt1MC->Add(hPt1C);
  hPt1MC->Add(hPt1QCD);



  cout<<dataEffCorr1<<endl;
  
  TH1F *hBBar=new TH1F("hBBar","hBBar",10,0,1);
  TH1F *hBX=new TH1F("hBX","hBX",10,0,1);
  TH1F *hCCbar=new TH1F("hCCbar","hCCbar",10,0,1);
  TH1F *hCX=new TH1F("hCX","hCX",10,0,1);
  TH1F *hQCD=new TH1F("hQCD","hQCD",10,0,1);
  TH1F *hData=new TH1F("hData","hData",10,0,1);
  hData->Sumw2();

  TH1F *hBgBBar=new TH1F("hBgBBar","hBgBBar",10,0,1);
  TH1F *hBgBX=new TH1F("hBgBX","hBgBX",10,0,1);
  TH1F *hBgCCbar=new TH1F("hBgCCbar","hBgCCbar",10,0,1);
  TH1F *hBgCX=new TH1F("hBgCX","hBgCX",10,0,1);
  TH1F *hBgQCD=new TH1F("hBgQCD","hBgQCD",10,0,1);
  TH1F *hBgData=new TH1F("hBgData","hBgData",10,0,1);
  hBgData->Sumw2();


  string  ptVar;
  if(plotPtRatio) ptVar = Form("%s/%s",pt2Choice,pt1Choice);
  else ptVar = Form("(%s-%s)/(%s+%s)",pt1Choice,pt2Choice,pt1Choice,pt2Choice);

  cout<<" ptVar "<<ptVar<<endl;

  char *weightString = Form("(weight/(%s)/(%s))",effCorr1,effCorr2);
  char *weightStringData = Form("(weight/(%s)/(%s))",dataEffCorr1,dataEffCorr2);
  //char *cut = Form("pthat>50&&%s>70&&%s<200&&%s>40&&%s<200&&acos(cos(jtphi1-jtphi2))>2./3.*acos(-1.)&%s1>=%f&&%s2>=%f",effCorr1,effCorr2,pt1Choice,pt1Choice,pt2Choice,pt2Choice,disc.c_str(),discCut1,disc.c_str(),discCut2));

  cout<<" weightString "<<weightString<<endl;
  cout<<" weightStringData "<<weightStringData<<endl;
  cout<<"pt1Choice "<<pt1Choice<<endl;
  cout<<"pt2Choice "<<pt2Choice<<endl;
  cout<<"disc "<<disc<<endl;
  cout<<"discCut1 "<<discCut1<<endl;
  cout<<"discCut2 "<<discCut2<<endl;
  
  char drawBBar1[5000], drawBX1[5000], drawCCbar1[5000], drawCX1[5000], drawQCD1[5000], drawData1[5000];
  char drawBBar2[5000], drawBX2[5000], drawCCbar2[5000], drawCX2[5000], drawQCD2[5000], drawData2[5000];

  sprintf(drawBBar1,"%s>>hBBar",ptVar.c_str());
  sprintf(drawBBar2,"%s*(pthat>50&&%s>70&&%s<200&&%s>40&&%s<200&&acos(cos(jtphi1-jtphi2))>2./3.*acos(-1.)&&abs(refparton_flavorForB1)==5&&abs(refparton_flavorForB2)==5&&%s1>=%f&&%s2>=%f&&refpt1>0&&refpt2>0)",weightString,pt1Choice,pt1Choice,pt2Choice,pt2Choice,disc.c_str(),discCut1,disc.c_str(),discCut2);

  sprintf(drawBX1,"%s>>hBX",ptVar.c_str());
  sprintf(drawBX2,"%s*(pthat>50&&%s>70&&%s<200&&%s>40&&%s<200&&acos(cos(jtphi1-jtphi2))>2./3.*acos(-1.)&&(!(abs(refparton_flavorForB1)==5&&abs(refparton_flavorForB2)==5))&&(abs(refparton_flavorForB1)==5||abs(refparton_flavorForB2)==5)&&%s1>=%f&&%s2>=%f&&refpt1>0&&refpt2>0)",weightString,pt1Choice,pt1Choice,pt2Choice,pt2Choice,disc.c_str(),discCut1,disc.c_str(),discCut2);

  sprintf(drawCCbar1,"%s>>hCCbar",ptVar.c_str());
  sprintf(drawCCbar2,"%s*(pthat>50&&%s>70&&%s<200&&%s>40&&%s<200&&acos(cos(jtphi1-jtphi2))>2./3.*acos(-1.)&&abs(refparton_flavorForB1)==4&&abs(refparton_flavorForB2)==4&&%s1>=%f&&%s2>=%f&&refpt1>0&&refpt2>0)",weightString,pt1Choice,pt1Choice,pt2Choice,pt2Choice,disc.c_str(),discCut1,disc.c_str(),discCut2);

  sprintf(drawCX1,"%s>>hCX",ptVar.c_str());
  sprintf(drawCX2,"%s*(pthat>50&&%s>70&&%s<200&&%s>40&&%s<200&&acos(cos(jtphi1-jtphi2))>2./3.*acos(-1.)&&abs(refparton_flavorForB1)!=5&&abs(refparton_flavorForB2)!=5&&(abs(refparton_flavorForB1)==4||abs(refparton_flavorForB2)==4)&&(!(abs(refparton_flavorForB1)==4&&abs(refparton_flavorForB2)==4))&&%s1>=%f&&%s2>=%f&&refpt1>0&&refpt2>0)",weightString,pt1Choice,pt1Choice,pt2Choice,pt2Choice,disc.c_str(),discCut1,disc.c_str(),discCut2);

  sprintf(drawQCD1,"%s>>hQCD",ptVar.c_str());
  sprintf(drawQCD2,"%s*(pthat>50&&%s>70&&%s<200&&%s>40&&%s<200&&acos(cos(jtphi1-jtphi2))>2./3.*acos(-1.)&&abs(refparton_flavorForB1)!=5&&abs(refparton_flavorForB2)!=5&&abs(refparton_flavorForB1)!=4&&abs(refparton_flavorForB2)!=4&&%s1>=%f&&%s2>=%f&&refpt1>0&&refpt2>0)",weightString,pt1Choice,pt1Choice,pt2Choice,pt2Choice,disc.c_str(),discCut1,disc.c_str(),discCut2);

  sprintf(drawData1,"%s>>hData",ptVar.c_str());
  sprintf(drawData2,"%s*(%s>70&&%s<200&&%s>40&&%s<200&&acos(cos(jtphi1-jtphi2))>2./3.*acos(-1.)&&%s1>=%f&&%s2>=%f)",weightString,pt1Choice,pt1Choice,pt2Choice,pt2Choice,disc.c_str(),discCut1,disc.c_str(),discCut2);

  ntB->Draw(drawBBar1,drawBBar2);
  ntB->Draw(drawBX1,drawBX2,"same");  
  ntC->Draw(drawCCbar1,drawCCbar2,"same");
  ntC->Draw(drawCX1,drawCX2,"same");
  ntQCD->Draw(drawQCD1,drawQCD2,"same");
  ntdata->Draw(drawData1,drawData2,"sames");

  if(doSub){

    sprintf(drawBBar1,"%s>>hBgBBar",ptVar.c_str());
    sprintf(drawBBar2,"%s*(pthat>50&&%s>70&&%s<200&&%s>40&&%s<200&&acos(cos(jtphi1-jtphi2))<1./3.*acos(-1.)&&abs(refparton_flavorForB1)==5&&abs(refparton_flavorForB2)==5&&%s1>=%f&&%s2>=%f&&refpt1>0&&refpt2>0)",weightString,pt1Choice,pt1Choice,pt2Choice,pt2Choice,disc.c_str(),discCut1,disc.c_str(),discCut2);
    
    sprintf(drawBX1,"%s>>hBgBX",ptVar.c_str());
    sprintf(drawBX2,"%s*(pthat>50&&%s>70&&%s<200&&%s>40&&%s<200&&acos(cos(jtphi1-jtphi2))<1./3.*acos(-1.)&&(!(abs(refparton_flavorForB1)==5&&abs(refparton_flavorForB2)==5))&&(abs(refparton_flavorForB1)==5||abs(refparton_flavorForB2)==5)&&%s1>=%f&&%s2>=%f&&refpt1>0&&refpt2>0)",weightString,pt1Choice,pt1Choice,pt2Choice,pt2Choice,disc.c_str(),discCut1,disc.c_str(),discCut2);
    
    sprintf(drawCCbar1,"%s>>hBgCCbar",ptVar.c_str());
    sprintf(drawCCbar2,"%s*(pthat>50&&%s>70&&%s<200&&%s>40&&%s<200&&acos(cos(jtphi1-jtphi2))<1./3.*acos(-1.)&&abs(refparton_flavorForB1)==4&&abs(refparton_flavorForB2)==4&&%s1>=%f&&%s2>=%f&&refpt1>0&&refpt2>0)",weightString,pt1Choice,pt1Choice,pt2Choice,pt2Choice,disc.c_str(),discCut1,disc.c_str(),discCut2);
    
    sprintf(drawCX1,"%s>>hBgCX",ptVar.c_str());
    sprintf(drawCX2,"%s*(pthat>50&&%s>70&&%s<200&&%s>40&&%s<200&&acos(cos(jtphi1-jtphi2))<1./3.*acos(-1.)&&abs(refparton_flavorForB1)!=5&&abs(refparton_flavorForB2)!=5&&(abs(refparton_flavorForB1)==4||abs(refparton_flavorForB2)==4)&&(!(abs(refparton_flavorForB1)==4&&abs(refparton_flavorForB2)==4))&&%s1>=%f&&%s2>=%f&&refpt1>0&&refpt2>0)",weightString,pt1Choice,pt1Choice,pt2Choice,pt2Choice,disc.c_str(),discCut1,disc.c_str(),discCut2);
    
    sprintf(drawQCD1,"%s>>hBgQCD",ptVar.c_str());
    sprintf(drawQCD2,"%s*(pthat>50&&%s>70&&%s<200&&%s>40&&%s<200&&acos(cos(jtphi1-jtphi2))<1./3.*acos(-1.)&&abs(refparton_flavorForB1)!=5&&abs(refparton_flavorForB2)!=5&&abs(refparton_flavorForB1)!=4&&abs(refparton_flavorForB2)!=4&&%s1>=%f&&%s2>=%f&&refpt1>0&&refpt2>0)",weightString,pt1Choice,pt1Choice,pt2Choice,pt2Choice,disc.c_str(),discCut1,disc.c_str(),discCut2);
    
    sprintf(drawData1,"%s>>hBgData",ptVar.c_str());
    sprintf(drawData2,"%s*(%s>70&&%s<200&&%s>40&&%s<200&&acos(cos(jtphi1-jtphi2))<1./3.*acos(-1.)&&%s1>=%f&&%s2>=%f)",weightString,pt1Choice,pt1Choice,pt2Choice,pt2Choice,disc.c_str(),discCut1,disc.c_str(),discCut2);
    
    ntB->Draw(drawBBar1,drawBBar2);
    ntB->Draw(drawBX1,drawBX2,"same");  
    ntC->Draw(drawCCbar1,drawCCbar2,"same");
    ntC->Draw(drawCX1,drawCX2,"same");
    ntQCD->Draw(drawQCD1,drawQCD2,"same");
    ntdata->Draw(drawData1,drawData2,"sames");

    hBBar->Add(hBgBBar,-1);    
    hBX->Add(hBgBX,-1);
    hCCbar->Add(hBgCCbar,-1);
    hCX->Add(hBgCX,-1);
    hQCD->Add(hBgQCD,-1);
    hData->Add(hBgData,-1);
  }
  
  hData->SetMarkerStyle(8);

  hBBar->SetLineColor(1);
  hBX->SetLineColor(1);
  hCCbar->SetLineColor(1);
  hCX->SetLineColor(1);
  hQCD->SetLineColor(1);

  hBBar->SetFillColor(kMagenta);
  hBX->SetFillColor(kRed);
  hCCbar->SetFillColor(kGreen);
  hCX->SetFillColor(kGreen-2);
  hQCD->SetFillColor(kAzure);
  
  hCX->Add(hQCD);
  hCCbar->Add(hCX);
  hBX->Add(hCCbar);
  hBBar->Add(hBX);

  float norm=0.;
 
  if(pt1norm){
    norm=hPt1MC->Integral(1,121);  
    cout<<" MC:  fraction of leading jets w/ dijet = "<<hBBar->Integral()/hPt1MC->Integral()<<endl;
  }
  else norm=hBBar->Integral();

  hBBar->Scale(1./norm);
  hBX->Scale(1./norm);
  hCCbar->Scale(1./norm);
  hCX->Scale(1./norm);
  hQCD->Scale(1./norm);


  if(pt1norm){
    cout<<" Data:  fraction of leading jets w/ dijet = "<<hData->Integral()/hPt1Data->Integral()<<endl;
    hData->Scale(1./hPt1Data->Integral(1,121));
  }
  else hData->Scale(1./hData->Integral());

  TCanvas *c1=new TCanvas("c1","c1",600,600);
  c1->cd();

  hBX->SetStats(0);
  hCCbar->SetStats(0);
  hCX->SetStats(0);
  hQCD->SetStats(0);

  hBBar->Draw();
  hBX->Draw("same");
  hCCbar->Draw("same");
  hCX->Draw("same");
  hQCD->Draw("same");
  hQCD->Draw("sameaxis");
  hData->Draw("sames");
  if(plotPtRatio){
    hBBar->SetXTitle("p_{T,2} / p_{T,1}");
    if(pt1norm)hBBar->SetYTitle("1/n_{leading} dN/d(p_{T,2}/p_{T,1})");
    else hBBar->SetYTitle("unit norm.");
  }
  else{
    hBBar->SetXTitle("A_{J}");
    if(pt1norm)hBBar->SetYTitle("1/n_{leading} dN/dA_{J}");
    else hBBar->SetYTitle("unit norm.");
  }
  //hBBar->GetYaxis()->SetTitleOffset(1.2);
  //hBBar->GetXaxis()->SetNdivisions(505);

  hBBar->SetFillStyle(1001);
  hBX->SetFillStyle(1001);
  hCCbar->SetFillStyle(1001);
  hCX->SetFillStyle(1001);
  hQCD->SetFillStyle(1001);

  if(drawStats){
    gPad->Update();
    TPaveStats *stat1 = (TPaveStats*) hData->FindObject("stats");
    stat1->SetLabel("Data");
    if(plotPtRatio){
      stat1->SetY1NDC(stat1->GetY1NDC()+0.05);
      stat1->SetY2NDC(stat1->GetY2NDC()-0.07);
      stat1->SetX1NDC(stat1->GetX1NDC()-0.58);
      stat1->SetX2NDC(stat1->GetX2NDC()-0.58);
    }
    else{
      stat1->SetY1NDC(stat1->GetY1NDC()+0.03);
      stat1->SetY2NDC(stat1->GetY2NDC()-0.07);
      stat1->SetX1NDC(stat1->GetX1NDC()-0.275);
      stat1->SetX2NDC(stat1->GetX2NDC()-0.175);
    }
    stat1->Draw();

    TPaveStats *stat2 = (TPaveStats*) hBBar->FindObject("stats");
    stat2->SetLabel("MC");
    if(plotPtRatio){
      stat2->SetY1NDC(stat1->GetY1NDC()-0.28);
      stat2->SetY2NDC(stat1->GetY2NDC()-0.28);
      stat2->SetX1NDC(stat2->GetX1NDC()-0.58);
      stat2->SetX2NDC(stat2->GetX2NDC()-0.58);
    }
    else{
      stat2->SetY1NDC(stat1->GetY1NDC()-0.3);
      stat2->SetY2NDC(stat1->GetY2NDC()-0.3);
      stat2->SetX1NDC(stat2->GetX1NDC()-0.275);
      stat2->SetX2NDC(stat2->GetX2NDC()-0.175);
    }

    stat2->Draw();
  }




  hQCD->SetLineColor(1);

  TLegend *leg = NULL;
  if(plotPtRatio){
    if(drawStats) leg = new TLegend(0.5,0.2,0.8,0.5);
    else leg = new TLegend(0.19,0.655,0.49,0.955);
  }
  else{
    if(drawStats) leg = new TLegend(0.2,0.175,0.5,0.475);
    else leg = new TLegend(0.5,0.6,0.8,0.9);
  }
  leg->SetBorderSize(0.);
  leg->SetFillStyle(0.);
  leg->AddEntry(hData,"Data","pl");
  leg->AddEntry(hBBar,"b-#bar{b}","f");
  leg->AddEntry(hBX,"b(#bar{b})+X","f");
  leg->AddEntry(hCCbar,"c-#bar{c}","f");
  leg->AddEntry(hCX,"c(#bar{c})+light","f");
  leg->AddEntry(hQCD,"light jets","f");
  leg->Draw();
  
}
