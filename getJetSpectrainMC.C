//Code to generate inclusive jet and b-jet spectra from the MC histograms for pp @ 5 TeV
//Options to match the pt-binning scheme to the data one.
//Kurt Jung, Purdue University: Jan 7, 2014

#include "TFile.h"
#include "TH1D.h"
#include "TTree.h"

#include <iostream>

using namespace std;

Double_t binomError(Double_t a, Double_t b, Double_t aErr, Double_t bErr) {
  // error on a/b, where a+b=const 
  Double_t w = a/b;
  return sqrt(fabs( ( (1.-2.*w)*aErr*aErr + w*w*bErr*bErr )/(b*b) ));
}

void getJetSpectrainMC(bool doMatchingToData=0, int gsp=0){

  int nBins = 8; //not const because these can get modified to match to data
  float *ptBin = NULL;
  int ptBinDefaultBinning[9] = {40,55,70,90,110,140,170,220,400};

  //int nBins = 29;
  //double ptBinDefaultBinning[30] = {22,27,33,39,47,55,64,74,84,97,114,133,153,174,196,220,245,272,300,330,362,395,430,468,507,548,592,638,790,967};
  
  //******* VARIABLES TO PLAY WITH ******** //
  const string discriminator = "discr_ssvHighEff";
  const double discrCut = 2.0;
  const double maxDiscrCut = 6.0;
  //const double etaCut = 2.0;
  float etalo=-2., etahi=2.;
  //*************************************** //

  //fix because the pp MC got boosted somehow
  etalo+=0.465;
  etahi+=0.465;
  //crap....


  TH1D *hBFractionMC;
  TH1D *hBPurityMC, *hRawBMC, *hBEfficiencyMC;
  TH1D *hIncJetsMC;
  TH1D *crossCheck;

  //Match the pt binning to data here, or else use default values for binning scheme
  TFile *fData;
  if(doMatchingToData){
    fData = new TFile("output/NewFormatV13_FixedTrgMergingSTAR_ak3PF_fixBin_bFractionMCTemplate_pPbpp1_jetptcut30_SSVHEat2.0FixCL0_bin_0_40_eta_0_2.root");
    TH1D *htmp = (TH1D*)fData->Get("hBFractionMC")->Clone("htmp");
    nBins = htmp->GetNbinsX();
    ptBin = new float[nBins+1];
    for(int i=1; i<=nBins; i++){
      ptBin[i-1] = htmp->GetBinLowEdge(i);
    }
    ptBin[nBins] = htmp->GetBinLowEdge(nBins) + htmp->GetBinWidth(nBins);
    std::cout << "Matching Bin Scheme to Data..." << std::endl;
  }
  else{
    ptBin = new float[nBins+1];
    for(int ii=0; ii<=nBins; ii++){
      ptBin[ii] = ptBinDefaultBinning[ii];
    }
  }
  std::cout << "Using " << nBins << " pt bins of: " << std::endl << "{ ";
  for(int i=0; i<=nBins; i++){ std::cout << ptBin[i] << " ";}
  std::cout << "}" << std::endl;

  hBPurityMC = new TH1D("hBPurityMC","hBPurityMC;Jet p_{T} (GeV/c);b-Tagging purity",nBins,ptBin); hBPurityMC->Sumw2();
  hRawBMC = new TH1D("hRawBMC","hRawBMC;Jet p_{T} (GeV/c);raw b-jets",nBins,ptBin); hRawBMC->Sumw2();
  hBEfficiencyMC = new TH1D("hBEfficiencyMC","hBEfficiencyMC;Jet p_{T} (GeV/c);b-Tagging efficiency",nBins,ptBin); hBEfficiencyMC->Sumw2();
  hBFractionMC = new TH1D("hBFractionMC","hBFractionMC;Jet p_{T} (GeV/c);b-jet fraction",nBins,ptBin); hBFractionMC->Sumw2();
  hIncJetsMC = new TH1D("hIncJetsMC","hIncJetsMC;Jet p_{T} (GeV/c); inclusive jets in MC",nBins,ptBin); hIncJetsMC->Sumw2();
  crossCheck = new TH1D("crossCheck","",nBins, ptBin); crossCheck->Sumw2();
    
  //Get pp MC and fill bins pt bin by pt bin
  TFile *fin = NULL;
  //if(gsp==0) fin = new TFile("histos/ppMC_ppReco_akPu3PF_QCDjetTrig_etashift_MCWeightFinalWithVz_noTrgSelection_noIPupperCut_8.root");
  if(gsp==0) fin = new TFile("histos/ppMC_ppReco_akPu3PF_QCDjetTrig_etashift_Fix2Sample_MCWeightFinalWithVz_noTrgSelection_Full.root");
  else if(gsp==2) fin = new TFile("histos/pPbMC_ppReco_akPu3PF_QCDjetTrig_etashift_WeightCorr_useGSP2_noIPupperCut_7.root");
  else if (gsp==3) fin = new TFile("histos/pPbMC_ppReco_akPu3PF_QCDjetTrig_etashift_WeightCorr_useGSP3_noIPupperCut_7.root");
  TTree *nt = (TTree*)fin->Get("nt");

  TH1D *hTagBJets = new TH1D("hTagBJets","",nBins,ptBin); hTagBJets->Sumw2();
  TH1D *hTagJets = new TH1D("hTagJets","",nBins,ptBin); hTagJets->Sumw2();
  TH1D *hJets = new TH1D("hJets","",nBins,ptBin); hJets->Sumw2();
  TH1D *hBJets = new TH1D("hBJets","",nBins,ptBin); hBJets->Sumw2();
  TH1D *hTaggedBJetsMC = new TH1D("hTaggedBJetsMC","",nBins,ptBin); hTaggedBJetsMC->Sumw2();
  TH1D *hUntaggedBJetsMC = new TH1D("hUntaggedBJetsMC","",nBins,ptBin); hUntaggedBJetsMC->Sumw2();

  int cbinlo=0, cbinhi=100;
    
  nt->Draw("jtpt>>hTagBJets",Form("weight*(%s>=%f && %s<%f && abs(refparton_flavorForB)==5 && jteta>%f && jteta<%f)",discriminator.c_str(),discrCut,discriminator.c_str(),maxDiscrCut,etalo,etahi));
  nt->Draw("jtpt>>hTagJets",Form("weight*(%s>=%f && %s<%f && jteta>%f && jteta<%f)",discriminator.c_str(),discrCut,discriminator.c_str(),maxDiscrCut,etalo,etahi));
  nt->Draw("jtpt>>hJets",Form("weight*(jteta>%f && jteta<%f && abs(refparton_flavorForB)<99)",etalo,etahi));
  nt->Draw("jtpt>>hBJets",Form("weight*(abs(refparton_flavorForB)==5 && jteta>%f && jteta<%f)",etalo,etahi));
  nt->Draw("jtpt>>hTaggedBJetsMC",Form("weight*(%s>=%f&&bin>=%d&&bin<%d&&jteta>%f&&jteta<%f&&abs(refparton_flavorForB)==5)",discriminator.c_str(),discrCut,cbinlo,cbinhi,etalo,etahi));
  nt->Draw("jtpt>>hUntaggedBJetsMC",Form("weight*(%s<%f&&bin>=%d&&bin<%d&&jteta>%f&&jteta<%f&&abs(refparton_flavorForB)==5)",discriminator.c_str(),discrCut,cbinlo,cbinhi,etalo,etahi));

  for(int ibin=0; ibin<nBins; ibin++){
    //cout << "Starting bin " << ptBin[ibin] << " to " << ptBin[ibin+1] << "..." << endl;
    
    hBPurityMC->SetBinContent(ibin+1,hTagBJets->GetBinContent(ibin+1)/hTagJets->GetBinContent(ibin+1));
    hBPurityMC->SetBinError(ibin+1,binomError(hTagBJets->GetBinContent(ibin+1),hTagJets->GetBinContent(ibin+1),hTagBJets->GetBinError(ibin+1),hTagJets->GetBinError(ibin+1)));

    hBEfficiencyMC->SetBinContent(ibin+1,hTagBJets->GetBinContent(ibin+1)/hBJets->GetBinContent(ibin+1));
    hBEfficiencyMC->SetBinError(ibin+1,binomError(hTagBJets->GetBinContent(ibin+1),hBJets->GetBinContent(ibin+1),hTagBJets->GetBinError(ibin+1),hBJets->GetBinError(ibin+1)));

    hRawBMC->SetBinContent(ibin+1,hBJets->GetBinContent(ibin+1));
    hRawBMC->SetBinError(ibin+1,hBJets->GetBinError(ibin+1));

    hBFractionMC->SetBinContent(ibin+1,hBJets->GetBinContent(ibin+1)/hJets->GetBinContent(ibin+1));
    hBFractionMC->SetBinError(ibin+1,binomError(hBJets->GetBinContent(ibin+1),hJets->GetBinContent(ibin+1),hBJets->GetBinError(ibin+1),hJets->GetBinError(ibin+1)));

    hIncJetsMC->SetBinContent(ibin+1,hJets->GetBinContent(ibin+1));
    hIncJetsMC->SetBinError(ibin+1,hJets->GetBinError(ibin+1));

    crossCheck->SetBinContent(ibin+1,hTaggedBJetsMC->GetBinContent(ibin+1)+hUntaggedBJetsMC->GetBinContent(ibin+1));

  }

  TFile *fout = new TFile(Form("output/ppMC_akPu3PF_bJet_TEST_matchBin%d_gsp%d_discrSSVHEGT%g_eta%.1fTo%.1f.root",doMatchingToData,gsp,discrCut,etalo,etahi),"recreate");
  fout->cd();
  
  cout << "before: " << hRawBMC->GetBinContent(3) << endl;
  //WARNING! IMPORTANT STEP TAKEN HERE!
  hRawBMC->Divide(hBEfficiencyMC); //transforms "nBjets" into "nBtaggedBJets", which is what the next step expects.
  cout << "after: "<< hRawBMC->GetBinContent(3) << endl;

  hBPurityMC->Write();
  hRawBMC->Write();
  hBEfficiencyMC->Write();
  hBFractionMC->Write();
  hIncJetsMC->Write();
  crossCheck->Write();
}
