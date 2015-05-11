#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TF1.h"
#include <stdlib.h>
#include <math.h>

using namespace std;

TH1* killZeroes(TH1* h1, double min){
  int j=0, k=1;
  int binwidths[h1->GetNbinsX()];
  for(int i=1; i<=h1->GetNbinsX(); i++){
    cout << "Bin " << i << " content: "<< h1->GetBinContent(i) << endl;
    if(h1->GetBinContent(i)>min || h1->GetBinContent(i) != h1->GetBinContent(i) ){
      cout << "is nonzero! "<< endl;
      binwidths[j++] = h1->GetBinWidth(i);
    }
  }
  cout << j << " nonzero bins found..." << endl;
  TH1D *hnew = new TH1D("hnew","",j,h1->GetBinLowEdge(1),h1->GetBinLowEdge(h1->GetNbinsX())+h1->GetBinWidth(h1->GetNbinsX()));
  hnew->Sumw2();
  for(int i=1; i<=h1->GetNbinsX(); i++){
    if(h1->GetBinContent(i)>min || h1->GetBinContent(i) != h1->GetBinContent(i)){
      hnew->SetBinContent(k, h1->GetBinContent(i));
      hnew->SetBinError(k++, h1->GetBinError(i));
    }
  }
  for(int i=1; i<=hnew->GetNbinsX(); i++){
    cout << "hnew bin: "<< i << " is " << hnew->GetBinContent(i) << endl;
  }
  return hnew;
}

void drawBjetSpectrum(int ppPbPb=1, int pApp=0, int useUnfolded=1, int doInclusiveJet=0, float etalo=-2., float etahi=2., int centLo=0, int centHi=100, bool useSVD=1){

  double TpAarr[4] = {6.9, 14.3, 10.2, 3.86};
  
  //for(int ieta=-2; ieta<2; ieta++){
  //etalo=(float)ieta;
  //etahi=(float)(ieta+1);
  //ROOT::gStyle->SetOptTitle(0);

  //etalo+=0.465;
  //etahi+=0.465;

  if(!ppPbPb) {
    cout<<" find Kurt :-) "<<endl;
    return;
  }
  double ppetalo=etalo+0.465;
  double ppetahi=etahi+0.465;

  //if((pApp || doInclusiveJet) && useUnfolded){
  //  cout << "warning! Selecting pp or inclusive Jets with unfolded, which is not supported! Selecting non-unfolded jets..." << endl;
  //   useUnfolded=0;
  // }
  TFile *fin = NULL;
  TFile *fin2 = NULL;
  TFile *finRev = NULL;
  if(pApp && doInclusiveJet) fin2 = new TFile(Form("output/ppMC_akPu3PF_InclJetForUnfolding_gsp0_eta%.1fTo%.1f.root",ppetalo,ppetahi));
  else if(pApp) fin2= new TFile(Form("output/ppMC_akPu3PF_bJet_matchBin0_gsp0_discrSSVHEGT2_eta%.1fTo%.1f.root",ppetalo,ppetahi)); //pp (from MC)
  else  fin2= new TFile("output/NewFormatV5_bFractionMCTemplate_pPbpp1_jetptcut30_SSVHEat2.0FixCL0_bin_0_40_eta_0_2.root"); //pA
  

  if(!useUnfolded){
    fin = new TFile("officialMCoutput/NewFormatV31_officialMC_FixedTrgMergingFinal_akPu3PF_FirstHalf_consistentEta_ssvhe20_addlJetCuts_fixBin_bFractionMCTemplate_pPbpp1_gsp0_jetptcut30_SSVHEat2.0FixCL0_bin_0_100_eta_-2.00_2.00.root");
    finRev = new TFile("officialMCoutput/NewFormatV31_officialMC_FixedTrgMergingFinal_akPu3PF_SecondHalfDataWithpPbMC_JPRetrain_v17JEC_consistentEta_ssvhe20_addlJetCuts_fixBin_bFractionMCTemplate_pPbpp1_gsp0_jetptcut30_SSVHEat2.0FixCL0_bin_0_100_eta_-2.00_2.00.root");
  }
  else if(doInclusiveJet) fin = new TFile(Form("pPb_Unfo_akPu3PF_akPu3PF_noGplus_FirstHalf_jtpt30_Inc_clo0_chi100_v8_eta_%.1fTo%.1f_.root",etalo,etahi));
  else{
    fin = new TFile(Form("pPb_Unfo_inCM_v31_officialMC_ak3PF_akPu3PF_noGplus_FirstHalfOnly_Converged_usedParameterizedUnfold0_jtpt35_bJets_clo%d_chi%d_v8_eta_%.2fTo%.2f_.root",centLo,centHi,etalo,etahi));
    finRev = new TFile(Form("pPb_Unfo_inCM_v31_officialMC_Reverse_WithResCorr_ak3PF_akPu3PF_noGplus_SecondHalfOnly_Converged_usedParameterizedUnfold0_jtpt35_bJets_clo%d_chi%d_v8_eta_%.2fTo%.2f_.root",centLo,centHi,etalo,etahi));
  }
  TFile *finpp = new TFile(Form("pPb_Unfo_inCM_v31_officialMC_ak3PF_akPu3PF_noGplus_FirstHalfOnly_Converged_usedParameterizedUnfold0_jtpt35_bJets_clo%d_chi%d_v8_eta_%.2fTo%.2f_.root",centLo,centHi,etalo,etahi));

  if(!fin){ cout << "file not found! " << endl; exit(0); }

  TH1D *hRawBData;
  TH1D *hRawBDataRev;
  TH1D *hRawpp;
  if(!useUnfolded){
    if(doInclusiveJet){
      if(pApp) hRawBData = (TH1D*)fin2->Get("hIncJetsMC");
      else hRawBData = (TH1D*)fin->Get("hIncJetsData");
    }
    else{
      if(pApp) hRawBData = (TH1D*)fin2->Get("hRawBMC");
      else{
	hRawBData = (TH1D*)fin->Get("hRawBData");
	hRawBDataRev = (TH1D*)finRev->Get("hRawBData");
	if(!pApp) hRawpp = (TH1D*)finpp->Get("hGen_cent1");
      }
    }
  }
  else{
    if(pApp) hRawBData = (TH1D*)fin->Get("hGen_cent1");
    else if(!useSVD) hRawBData = (TH1D*)fin->Get("hReco0"); //should be hReco0 for bayesian unfolding
    else hRawBData = (TH1D*)fin->Get("hRecoSVD0"); // this one is SVD unfolding
    //if(pApp) hRawBDataRev = (TH1D*)finRev->Get("hGen_cent1");
    if(!pApp && !useSVD) hRawBDataRev = (TH1D*)finRev->Get("hReco0"); //should be hReco0 for bayesian
    else if(!pApp) hRawBDataRev = (TH1D*)finRev->Get("hRecoSVD0");
    if(!pApp) hRawpp = (TH1D*)fin->Get("hGen_cent1");
  }

  if(!pApp) hRawBData->Add(hRawBDataRev); //add two datasets together!!

  if(useUnfolded){
    //  if(pApp) hRawBData = (TH1D*)killZeroes(hRawBData,1E-12);
    //  else hRawBData = (TH1D*)killZeroes(hRawBData,1.);
  }

  // divide out the bin-width
  double val1,err1,width1;
  if(!useUnfolded ||/* pApp ||*/ (useSVD && !pApp)){
    for(int i=0;i<hRawBData->GetNbinsX();i++){
      val1 =     hRawBData->GetBinContent(i+1);
      err1 =     hRawBData->GetBinError(i+1);
      width1 =   hRawBData->GetBinWidth(1+1);
      hRawBData->SetBinContent(i+1,val1/hRawBData->GetBinWidth(i+1));
      hRawBData->SetBinError(i+1,err1/hRawBData->GetBinWidth(i+1));
    }
  }
  //eff corrections taken care of in unfolding
  /*  if((!doInclusiveJet && useUnfolded && pApp)) {//|| (pApp && !doInclusiveJet)){
     double val,err,eff,effErr;
     TH1D *hBEfficiency = (TH1D*)fin2->Get("hBEfficiencyMC");
     for(int i=8;i<hRawBData->GetNbinsX();i++){
       val =     hRawBData->GetBinContent(i+1);
       err =     hRawBData->GetBinError(i+1);
       eff =     hBEfficiency->GetBinContent(i-7);
       effErr =     hBEfficiency->GetBinError(i-7);
       cout << "data bin center: "<< hRawBData->GetBinCenter(i+1) << " eff bin center: "<< hBEfficiency->GetBinCenter(i-7) << endl;
       
       if(hBEfficiency->GetBinContent(i-7)==0){
	 hRawBData->SetBinContent(i+1,0);
	 hRawBData->SetBinError(i+1,0);
       }
       else{
	 cout << "hd: " << hRawBData->GetBinCenter(i+1) << " he: " << hBEfficiency->GetBinCenter(i-7) << endl;
	 hRawBData->SetBinContent(i+1,hRawBData->GetBinContent(i+1)/hBEfficiency->GetBinContent(i-7));
	 hRawBData->SetBinError(i+1, val/eff *sqrt(effErr*effErr/eff/eff + err*err/val/val));
       }
     }
     }*/


  TCanvas *c1=new TCanvas("c1","c1",1);
  c1->SetLogy();

  hRawBData->SetMarkerStyle(8);
  hRawBData->SetYTitle("dN/dp_{T} (GeV/c)^{-1}");
  hRawBData->SetAxisRange(40.,399.,"X");


  /* TF1 *fpow = new TF1("fpow","[0]*pow(x,[1])",45,300);
 
  TH1D *hRawBDataOrig=(TH1D*)hRawBData->Clone("hRawBDataOrig");

  for(int iter=0;iter<4;iter++){

    cout<<" iteration # "<<iter<<endl;

    hRawBData->Fit(fpow);
    
    for(int i=0;i<hRawBData->GetNbinsX();i++){
      float meanVal = fpow->Integral(hRawBData->GetBinLowEdge(i+1),hRawBData->GetBinLowEdge(i+1)+hRawBData->GetBinWidth(i+1))/hRawBData->GetBinWidth(i+1);
      float centVal = fpow->Eval(hRawBData->GetBinCenter(i+1));
      float binShiftCorr = centVal/meanVal;
      cout<<" i "<<i<<" corr "<<binShiftCorr<<endl;
      
      float val2 =     hRawBDataOrig->GetBinContent(i+1);
      float err2 =     hRawBDataOrig->GetBinError(i+1);
      cout << "hrawBdata orig: " << hRawBData->GetBinContent(i+1) << endl;
      hRawBData->SetBinContent(i+1,val2*binShiftCorr);
      hRawBData->SetBinError(i+1,err2*binShiftCorr);
      cout << "modified: "<< hRawBData->GetBinContent(i+1) << endl;
      
    }
    }*/
  TFile *fout;
  if(pApp){
    if(doInclusiveJet) fout = new TFile(Form("raa_pp_denomForpA-inc_eta%.1fTo%.1f.root",etalo,etahi),"recreate");
    else if(!useUnfolded) fout = new TFile(Form("raa_pp_denomForpA-Bjet_NoUnfolding_etaCM_v2_bin%d_%d_eta%.1fTo%.1f.root",centLo,centHi,etalo,etahi),"recreate");
    else fout = new TFile(Form("raa_pp_denomForpA-Bjet_etaCM_v2_bin%d_%d_eta%.1fTo%.1f.root",centLo,centHi,etalo,etahi),"recreate");
    fout->cd();
    //hRawBData->Scale(1./(5.3E12*65E-3)); //scale by luminosity (5.3 pb^-1) and cross-section (65 mb)
    //hRawBData->Scale(1./70.); //conversion to mb from barns since its MC
    hRawBData->Draw();
    hRawBData->Write();
  }
  else{
    if(doInclusiveJet) fout = new TFile(Form("raa_pPb_numerator-inc_eta%.1fTo%.1f.root",etalo,etahi),"recreate");
    else if(!useUnfolded) fout = new TFile(Form("raa_pPb_numerator-Bjet_NoUnfolding_FullDataset_ssvhe20_etaCM_v2_bin%d_%d_eta%.1fTo%.1f.root",centLo,centHi,etalo,etahi),"recreate");
    else fout = new TFile(Form("raa_pPb_numerator-Bjet_FullDataset_ssvhe20_etaCM_v2_bin%d_%d_eta%.1fTo%.1f.root",centLo,centHi,etalo,etahi),"recreate");
    fout->cd();
    int arrnum=-1;
    if(centLo==0 && centHi==100) arrnum = 0;
    else if(centLo==0) arrnum = 1;
    if(centLo==20) arrnum = 2;
    if(centLo==50) arrnum = 3;
    cout << "TpA correction: " << TpAarr[arrnum] << endl;
    //hRawBData->Scale(1./(20.97E9*2110E-3*0.85*TpAarr[arrnum])); //scale by lumi (32 nb^-1) and cross-section (2110 mb) and account for evt selection (0.85)
    hRawBData->Scale(1./(35.092E9*(2110E-3)*0.85*TpAarr[arrnum])); //full dataset
    hRawBData->Scale(70.); //scale by sigmapp/sigmapA
    //hRawBData->Scale(1./(14.03E9*2110E-3*0.85)); //second half
    hRawBData->Draw();
    hRawpp->SetMarkerColor(2);
    //hRawpp->Scale(1./70.);
    hRawpp->Draw("same");
    hRawBData->SetNameTitle("hReco0","hReco");
    hRawBData->Write();
  }
}
//}
