#include <iostream>
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TF1.h"
#include "TH1.h"
#include "TMath.h"

void weightPtHatBinsBC(int LCB=1, int useGSP=1, int useQCD=1){

  // LCB =
  // 0:  Light jets
  // 1:  charm jets
  // 2:  bottom jets

    
  gROOT->Reset();
  
  Int_t start=0;
  Int_t N=8;

  cout<<" opening input files "<<endl;

  TFile *fin = NULL;
  if(LCB==0) fin = new TFile("/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/merged_bjetAnalyzers_hiReco_offPV_pt30by3_oldHydjet_restrictMixTripletA_ipHICalibCentWeight_qcd.root");
  else if(LCB==1){
    if(useQCD) fin = new TFile("/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/merged_bjetAnalyzers_hiReco_offPV_pt30by3_newHydjet_restrictMixTripletA_ipHICalibCentWeight_cJetPlusQCD.root");
    else fin = new TFile("/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/merged_bjetAnalyzers_hiReco_offPV_pt30by3_newHydjet_restrictMixTripletA_ipHICalibCentWeight_cJet.root");
  }
  else if(LCB==2){
    if(useQCD) fin = new TFile("/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/merged_bjetAnalyzers_hiReco_offPV_pt30by3_newHydjet_restrictMixTripletA_ipHICalibCentWeight_bJetPlusQCD.root");
    else fin = new TFile("/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/merged_bjetAnalyzers_hiReco_offPV_pt30by3_newHydjet_restrictMixTripletA_ipHICalibCentWeight_bJet.root");
  }
  else {
    cout<<" You suck at picking input files "<<endl;
    return;
  }

  cout<<" opened input files "<<endl;

  TTree *tr_out, *tr_out_skim;

  TTree *tr_in = (TTree*)fin->Get("akPu3PFJetAnalyzer/t");
  TTree *tr_in_skim = (TTree*)fin->Get("skimanalysis/HltTree");

  cout<<tr_in->GetMaxTreeSize()<<endl;
  cout<<tr_in->GetEntries()<<endl;
  cout<<tr_in->GetMaxVirtualSize()<<endl;

  Int_t bounds[8] = {30,50,65,80,100,120,170,200};

  Double_t xSections[8]={(0)};
  
  // from the twiki
  xSections[0] = 1.079e-02;  // 30
  xSections[1] = 1.021e-03; // 50
  xSections[2] = 1.021e-03; // 50
  xSections[3] = 9.913e-05; // 80
  xSections[4] = 9.913e-05; // 80
  xSections[5] = 1.128e-05; //120
  xSections[6] = 1.470e-06; //170
  xSections[7] = 5.310e-07; // 200


  // truncate bins
  xSections[0] *= 89849./99173.;   //30-50
  //xSections[1] *= 138481./153301.;   //50-80  
  xSections[1] *= 110902./153301.;  //50-65
  xSections[2] *= 27579./153301.;    //65-80
  //xSections[3] *= 126859./143194.;    //80-120
  xSections[3] *= 99035./143194.;    //80-100
  xSections[4] *= 27824./143194.;    //100-120
  //xSections[5] *= 34457./36104.;    //120-200
  xSections[5] *= 124953./143611.;    //120-170
  xSections[6] *= 40943./64264.;    //170-200


  //These are with vz weighting  
  double cPerEventQCD[8]={
    0.0694539,
    0.110722 ,
    0.127494 ,
    0.13941  ,
    0.150601 ,
    0.158814 ,
    0.164293 ,
    0.166625 };
  double cPerEventC[8]={
    0.6977  ,
    0.914041,
    0.97952	,
    1.02432	,
    1.05478	,
    1.07694	,
    1.10489	,
    1.13004};
  
  double bPerEventB[8] = {
    0.555707,
    0.866712,
    0.962954,
    1.0185  ,
    1.06156 ,
    1.0978  ,
    1.15189 ,
    1.15061 };
    
  double bPerEventQCD[8]={
    0.0260441,
    0.0453556,
    0.0548796,
    0.0584761,
    0.0626118,
    0.0699029,
    0.0718176,
    0.0715977};


  TF1 *fCent = new TF1("fCent","pol7",0,40);
  fCent->SetParameters(14781.9,-1641.19,127.245,-8.87318,0.41423,-0.011089,0.000154744,-8.76427e-07);

  TH1F *hDataCent=new TH1F("hDataCent","hDataCent",40,-0.5,39.5);
  for(int ibin=0;ibin<40;ibin++){
    double centIntegral = fCent->Integral(ibin,ibin+1);
    if(centIntegral>0)hDataCent->SetBinContent(ibin+1,centIntegral);
    else hDataCent->SetBinContent(ibin+1,0);
  }
  hDataCent->Scale(1./hDataCent->Integral());

  //TFile *fData = new TFile("/grid_mnt/vol__vol1__u/llr/cms/mnguyen/bTagging442p5/CMSSW_4_4_2_patch5/src/bTaggingMacros/histos/PbPbdata.root");
  //TFile *fData = new TFile("/grid_mnt/vol__vol1__u/llr/cms/mnguyen/bTagging442p5/CMSSW_4_4_2_patch5/src/bTaggingMacros/histos/PbPbdata_regPFforJets.root");
  TFile *fData = new TFile("/grid_mnt/vol__vol1__u/llr/cms/mnguyen/bTagging442p5/CMSSW_4_4_2_patch5/src/bTaggingMacros/histos/PbPbdata_pt30by3_jpHICalibRepass_withDup_PU.root");

  TH1F *hDataVz = (TH1F *)fData->Get("hvz");
  hDataVz->Rebin(4);
  hDataVz->Scale(1./hDataVz->Integral());

  // preliminary
  //TFile *fMC = new TFile("/grid_mnt/vol__vol1__u/llr/cms/mnguyen/bTagging442p5/CMSSW_4_4_2_patch5/src/bTaggingMacros/histos/PbPbQCDMC.root");
  //TFile *fMC = new TFile("/grid_mnt/vol__vol1__u/llr/cms/mnguyen/bTagging442p5/CMSSW_4_4_2_patch5/src/bTaggingMacros/histos/PbPbQCDMC_regPFforJets.root");
  TFile *fMCQCD=new TFile("/grid_mnt/vol__vol1__u/llr/cms/mnguyen/bTagging442p5/CMSSW_4_4_2_patch5/src/bTaggingMacros/histos/PbPbQCDMC_pt30by3_ipHICalibCentWeight.root");
  TFile *fMCBorC = NULL;
  if(LCB==1) fMCBorC=new TFile("/grid_mnt/vol__vol1__u/llr/cms/mnguyen/bTagging442p5/CMSSW_4_4_2_patch5/src/bTaggingMacros/histos/PbPbCMC_pt30by3_ipHICalibCentWeight.root");
  if(LCB==2) fMCBorC=new TFile("/grid_mnt/vol__vol1__u/llr/cms/mnguyen/bTagging442p5/CMSSW_4_4_2_patch5/src/bTaggingMacros/histos/PbPbBMC_pt30by3_ipHICalibCentWeight.root");
  
  TH1F *hMCVz[2], *hMCCent[2];


  hMCVz[0] = (TH1F *)fMCQCD->Get("hvz");
  hMCCent[0] =  (TH1F*)fMCQCD->Get("hbin");
  hMCVz[0]->Rebin(4);
  hMCVz[0]->Scale(1./hMCVz[0]->Integral());
  hMCCent[0]->Scale(1./hMCCent[0]->Integral());
  

  if(LCB>0){
    hMCVz[1] = (TH1F *)fMCBorC->Get("hvz");
    hMCCent[1] = (TH1F *)fMCBorC->Get("hbin");    
    hMCVz[1]->Rebin(4);
    hMCVz[1]->Scale(1./hMCVz[1]->Integral());
    hMCCent[1]->Scale(1./hMCCent[1]->Integral());
  }

  //char filename[500] = "merged_bjetAnalyzers_hiReco_offPV_pt30by3_newHydjet_restrictMixTripletA_ipHICalibCentWeight";

  Double_t weight, xSecWeight, centWeight, vzWeight;
  Double_t xSecWeights[8];

  // count individual inputs one by one
  
  double fQCDentries[8]= {(0)};
  double fBorCentries[8] = {(0)};    
  double totQCDentries=0.;
  double totBorCentries=0.;

  if(useQCD && LCB>0){
    cout<<" grabbing qcd "<<endl;
    TFile *finQCDtemp = new TFile("/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/merged_bjetAnalyzers_hiReco_offPV_pt30by3_oldHydjet_restrictMixTripletA_ipHICalibCentWeight_qcd.root");
    TTree *tr_inQCDtemp = (TTree*)finQCDtemp->Get("akPu3PFJetAnalyzer/t");    
    
    for (Int_t it=start; it<N; it++) {      
      char cutname[100];
      if (it<N-1) sprintf(cutname,"pthat>%d&&pthat<%d",bounds[it],bounds[it+1]);
      else sprintf(cutname,"pthat>%d",bounds[it]);      
      fQCDentries[it] = (Double_t)tr_inQCDtemp->GetEntries(cutname);          
      totQCDentries+=fQCDentries[it];
    }
    delete tr_inQCDtemp;
    finQCDtemp->Close();
    delete finQCDtemp;  
    
    
    cout<<" determining b or c normaliztion "<<endl;
    TFile *finBorCtemp = NULL;
    if(LCB==1)finBorCtemp = new TFile("/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/merged_bjetAnalyzers_hiReco_offPV_pt30by3_newHydjet_restrictMixTripletA_ipHICalibCentWeight_cJet.root");
    if(LCB==2)finBorCtemp = new TFile("/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/merged_bjetAnalyzers_hiReco_offPV_pt30by3_newHydjet_restrictMixTripletA_ipHICalibCentWeight_bJet.root");
    TTree *tr_inBorCtemp = (TTree*)finBorCtemp->Get("akPu3PFJetAnalyzer/t");    
    
    for (Int_t it=start; it<N; it++) {      
      cout<<"it "<<it<<endl;
      char cutname[200];
      if (it<N-1) sprintf(cutname,"pthat>%d&&pthat<%d",bounds[it],bounds[it+1]);
      else sprintf(cutname,"pthat>%d",bounds[it]);      
      cout<<" cutname "<<cutname<<endl;
      fBorCentries[it] = (Double_t)tr_inBorCtemp->GetEntries(cutname);          
      totBorCentries+=fBorCentries[it];
      cout<<fBorCentries[it]<<endl;
    }
    delete tr_inBorCtemp;
    finBorCtemp->Close();
    delete finBorCtemp;  
    
  }
  

  cout<<" print b and c cross sections "<<endl;
  // determine x-section for each bin
  for (Int_t it=start; it<N; it++) {

    char cutname[100];
    if (it<N-1) sprintf(cutname,"pthat>%d&&pthat<%d",bounds[it],bounds[it+1]);
    else sprintf(cutname,"pthat>%d",bounds[it]);

    Double_t nEffEntries= 0.;
    if(useQCD&&LCB>0) nEffEntries = (double)fBorCentries[it];
    else nEffEntries = (Double_t)tr_in->GetEntries(cutname);
    if(LCB==1) nEffEntries*=cPerEventC[it]/cPerEventQCD[it];
    if(LCB==2) nEffEntries*=bPerEventB[it]/bPerEventQCD[it];
    
    cout<<" nEffEntries from flavor-filtered "<<nEffEntries<<endl;
    if(useQCD&&LCB>0){
      nEffEntries += fQCDentries[it];
      cout<<" nEffEntries w/ QCD "<<nEffEntries<<endl;
    }
    xSecWeights[it] = xSections[it]/ nEffEntries;

    if(it<N-1) cout<<bounds[it]<<" < pT,hat < "<<bounds[it+1]<<", nEffEntries = "<<nEffEntries<<endl;
    else cout<<" pT,hat > "<<bounds[it]<<", nEffEntries = "<<nEffEntries<<endl;


  }


  



  //Declaration of leaves types
  Int_t           evt;
  Float_t         b;
  Float_t         hf;
  Int_t           bin;
  Int_t           nref;
  Float_t         vx,vy,vz;
  Float_t         rawpt[300];
  Float_t         jtpt[300];
  Float_t         jteta[300];
  Float_t         jty[300];
  Float_t         jtphi[300];
  Float_t         jtpu[300];
  Float_t         discr_ssvHighEff[300];
  Float_t         discr_ssvHighPur[300];
  //Float_t         discr_csvMva[300];
  Float_t         discr_csvSimple[300];
  //Float_t         discr_muByIp3[300];
  Float_t         discr_muByPt[300];
  Float_t         discr_prob[300];
  Float_t         discr_probb[300];
  Float_t         discr_tcHighEff[300];
  Float_t         discr_tcHighPur[300];
  Int_t           nsvtx[300];
  Int_t           svtxntrk[300];
  Float_t         svtxdl[300];
  Float_t         svtxdls[300];
  Float_t         svtxm[300];
  Float_t         svtxpt[300];
  Int_t           nIPtrk[300];
  Int_t           nselIPtrk[300];
  Int_t   nIP;
  Int_t   ipJetIndex[10000];
  Float_t ipPt[10000];
  Float_t ipProb0[10000];
  //Float_t ipProb1[10000];
  Float_t ip2d[10000];
  Float_t ip2dSig[10000];
  Float_t ip3d[10000];
  Float_t ip3dSig[10000];
  Float_t ipDist2Jet[10000];
  //Float_t ipDist2JetSig[10000];
  Float_t ipClosest2Jet[10000];
  Float_t         mue[300];
  Float_t         mupt[300];
  Float_t         mueta[300];
  Float_t         muphi[300];
  Float_t         mudr[300];
  Float_t         muptrel[300];
  Int_t           muchg[300];
  Float_t         pthat;
  Int_t         beamId1;
  Int_t         beamId2;
  Float_t         refpt[300];
  Float_t         refeta[300];
  Float_t         refy[300];
  Float_t         refphi[300];
  Float_t         refdphijt[300];
  Float_t         refdrjt[300];
  Float_t         refparton_pt[300];
  Int_t           refparton_flavor[300];
  Int_t           refparton_flavorForB[300];
  Bool_t           refparton_isGSP[300];
  /*
    Int_t           ngen;
    Int_t           genmatchindex[100];
    Float_t         genpt[100];
    Float_t         geneta[100];
    Float_t         geny[100];
    Float_t         genphi[100];
    Float_t         gendphijt[100];
    Float_t         gendrjt[100];
  */
  
  Float_t chargedMax[300];
  Float_t chargedSum[300];
  Int_t chargedN[300];
  Float_t chargedHardSum[300];
  Int_t chargedHardN[300];
  Float_t photonMax[300];
  Float_t photonSum[300];
  Int_t photonN[300];
  Float_t photonHardSum[300];
  Int_t photonHardN[300];
  Float_t neutralMax[300];
  Float_t neutralSum[300];
  Int_t neutralN[300];
  Float_t eMax[300];
  Float_t eSum[300];
  Int_t eN[300];
  Float_t muMax[300];
  Float_t muSum[300];
  Int_t muN[300];
  
  
  int nHLTBit;
  bool hltBit[12];
  
  
  
  //Declaration of leaves types
  Int_t           pvSel;
  Int_t           hbheNoiseSel;
  Int_t           spikeSel;
  Int_t           collSell;
  Int_t           hltAna;
  //Int_t           superFilterPath;
  


  // Set branch addresses.
  tr_in->SetBranchAddress("evt",&evt);
  tr_in->SetBranchAddress("b",&b);
  //tr_in->SetBranchAddress("hf",&hf);
  tr_in->SetBranchAddress("bin",&bin);
  tr_in->SetBranchAddress("vx",&vx);
  tr_in->SetBranchAddress("vy",&vy);
  tr_in->SetBranchAddress("vz",&vz);
  tr_in->SetBranchAddress("nref",&nref);
  tr_in->SetBranchAddress("rawpt",rawpt);
  tr_in->SetBranchAddress("jtpt",jtpt);
  tr_in->SetBranchAddress("jteta",jteta);
  tr_in->SetBranchAddress("jty",jty);
  tr_in->SetBranchAddress("jtphi",jtphi);
  tr_in->SetBranchAddress("jtpu",jtpu);
  tr_in->SetBranchAddress("discr_ssvHighEff",discr_ssvHighEff);
  tr_in->SetBranchAddress("discr_ssvHighPur",discr_ssvHighPur);
  //tr_in->SetBranchAddress("discr_csvMva",discr_csvMva);
  tr_in->SetBranchAddress("discr_csvSimple",discr_csvSimple);
  //tr_in->SetBranchAddress("discr_muByIp3",discr_muByIp3);
  tr_in->SetBranchAddress("discr_muByPt",discr_muByPt);
  tr_in->SetBranchAddress("discr_prob",discr_prob);
  tr_in->SetBranchAddress("discr_probb",discr_probb);
  tr_in->SetBranchAddress("discr_tcHighEff",discr_tcHighEff);
  tr_in->SetBranchAddress("discr_tcHighPur",discr_tcHighPur);
  tr_in->SetBranchAddress("nsvtx",nsvtx);
  tr_in->SetBranchAddress("svtxntrk",svtxntrk);
  tr_in->SetBranchAddress("svtxdl",svtxdl);
  tr_in->SetBranchAddress("svtxdls",svtxdls);
  tr_in->SetBranchAddress("svtxm",svtxm);
  tr_in->SetBranchAddress("svtxpt",svtxpt);
  tr_in->SetBranchAddress("nIPtrk",nIPtrk);
  tr_in->SetBranchAddress("nIP",&nIP);
  tr_in->SetBranchAddress("nselIPtrk",nselIPtrk);
  tr_in->SetBranchAddress("ipJetIndex",ipJetIndex);
  tr_in->SetBranchAddress("ipPt",ipPt);
  tr_in->SetBranchAddress("ipProb0",ipProb0);
  //tr_in->SetBranchAddress("ipProb1",ipProb1);
  tr_in->SetBranchAddress("ip2d",ip2d);
  tr_in->SetBranchAddress("ip2dSig",ip2dSig);
  tr_in->SetBranchAddress("ip3d",ip3d);
  tr_in->SetBranchAddress("ip3dSig",ip3dSig);
  tr_in->SetBranchAddress("ipDist2Jet",ipDist2Jet);
  //tr_in->SetBranchAddress("ipDist2JetSig",ipDist2JetSig);
  tr_in->SetBranchAddress("ipClosest2Jet",ipClosest2Jet);
  tr_in->SetBranchAddress("mue",mue);
  tr_in->SetBranchAddress("mupt",mupt);
  tr_in->SetBranchAddress("mueta",mueta);
  tr_in->SetBranchAddress("muphi",muphi);
  tr_in->SetBranchAddress("mudr",mudr);
  tr_in->SetBranchAddress("muptrel",muptrel);
  tr_in->SetBranchAddress("muchg",muchg);
  tr_in->SetBranchAddress("pthat",&pthat);
  tr_in->SetBranchAddress("beamId1",&beamId1);
  tr_in->SetBranchAddress("beamId2",&beamId2);
  tr_in->SetBranchAddress("refpt",refpt);
  tr_in->SetBranchAddress("refeta",refeta);
  tr_in->SetBranchAddress("refy",refy);
  tr_in->SetBranchAddress("refphi",refphi);
  tr_in->SetBranchAddress("refdphijt",refdphijt);
  tr_in->SetBranchAddress("refdrjt",refdrjt);
  tr_in->SetBranchAddress("refparton_pt",refparton_pt);
  tr_in->SetBranchAddress("refparton_flavor",refparton_flavor);
  tr_in->SetBranchAddress("refparton_flavorForB",refparton_flavorForB);
  if(useGSP)tr_in->SetBranchAddress("refparton_isGSP",refparton_isGSP);
  tr_in->SetBranchAddress("nHLTBit",&nHLTBit);
  tr_in->SetBranchAddress("hltBit",hltBit);
  
  /*
    tr_in->SetBranchAddress("ngen",&ngen);
    tr_in->SetBranchAddress("genmatchindex",genmatchindex);
    tr_in->SetBranchAddress("genpt",genpt);
    tr_in->SetBranchAddress("geneta",geneta);
    tr_in->SetBranchAddress("geny",geny);
    tr_in->SetBranchAddress("genphi",genphi);
    tr_in->SetBranchAddress("gendphijt",gendphijt);
    tr_in->SetBranchAddress("gendrjt",gendrjt);
  */
  
  
  
  tr_in->SetBranchAddress("chargedMax", chargedMax);
  tr_in->SetBranchAddress("chargedSum", chargedSum);
  tr_in->SetBranchAddress("chargedN", chargedN);
  tr_in->SetBranchAddress("chargedHardSum", chargedHardSum);
  tr_in->SetBranchAddress("chargedHardN", chargedHardN);
  
  tr_in->SetBranchAddress("photonMax", photonMax);
  tr_in->SetBranchAddress("photonSum", photonSum);
  tr_in->SetBranchAddress("photonN", photonN);
  tr_in->SetBranchAddress("photonHardSum", photonHardSum);
  tr_in->SetBranchAddress("photonHardN", photonHardN);
  
  tr_in->SetBranchAddress("neutralMax", neutralMax);
  tr_in->SetBranchAddress("neutralSum", neutralSum);
  tr_in->SetBranchAddress("neutralN", neutralN);
  
  tr_in->SetBranchAddress("eMax", eMax);
  tr_in->SetBranchAddress("eSum", eSum);
  tr_in->SetBranchAddress("eN", eN);
  
  tr_in->SetBranchAddress("muMax", muMax);
  tr_in->SetBranchAddress("muSum", muSum);
  tr_in->SetBranchAddress("muN", muN);
  
  
  // Set branch addresses.
  tr_in_skim->SetBranchAddress("pvSel",&pvSel);
  tr_in_skim->SetBranchAddress("hbheNoiseSel",&hbheNoiseSel);
  tr_in_skim->SetBranchAddress("spikeSel",&spikeSel);
  tr_in_skim->SetBranchAddress("collSell",&collSell);
  tr_in_skim->SetBranchAddress("hltAna",&hltAna);
  
  char outfile[500];
  if(LCB==0) sprintf(outfile,"/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/merged_bjetAnalyzers_hiReco_offPV_pt30by3_oldHydjet_restrictMixTripletA_ipHICalibCentWeight_weighted_qcd.root");
  else if(LCB==1){
    if(useQCD)sprintf(outfile,"/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/merged_bjetAnalyzers_hiReco_offPV_pt30by3_newHydjet_restrictMixTripletA_ipHICalibCentWeight\
_weighted_cJetPlusQCD.root");
    else sprintf(outfile,"/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/merged_bjetAnalyzers_hiReco_offPV_pt30by3_newHydjet_restrictMixTripletA_ipHICalibCentWeight\
_weighted_cJet.root");
  }
  else if(LCB==2){
    if(useQCD)sprintf(outfile,"/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/merged_bjetAnalyzers_hiReco_offPV_pt30by3_newHydjet_restrictMixTripletA_ipHICalibCentWeight\
_weighted_bJetPlusQCD.root");
    else sprintf(outfile,"/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/merged_bjetAnalyzers_hiReco_offPV_pt30by3_newHydjet_restrictMixTripletA_ipHICalibCentWeight\
_weighted_bJet.root");
  }

  TFile *fout =new TFile(outfile,"recreate");


  
  tr_out = new TTree("t","Jet Analyzer");
  tr_out_skim = new TTree("HltTree","HltTree");
  

  // Set output branch addresses.
  tr_out->Branch("evt",&evt,"evt/I");
  tr_out->Branch("b",&b,"b/F");
  tr_out->Branch("hf",&hf,"hf/F");
  tr_out->Branch("bin",&bin,"bin/I");
  tr_out->Branch("vx",&vx,"vx/F");
  tr_out->Branch("vy",&vy,"vy/F");
  tr_out->Branch("vz",&vz,"vz/F");
  tr_out->Branch("nref",&nref,"nref/I");
  tr_out->Branch("rawpt",rawpt,"rawpt[nref]/F");
  tr_out->Branch("jtpt",jtpt,"jtpt[nref]/F");
  tr_out->Branch("jteta",jteta,"jteta[nref]/F");
  tr_out->Branch("jty",jty,"jty[nref]/F");
  tr_out->Branch("jtphi",jtphi,"jtphi[nref]/F");
  tr_out->Branch("jtpu",jtpu,"jtpu[nref]/F");
  tr_out->Branch("discr_ssvHighEff",discr_ssvHighEff,"discr_ssvHighEff[nref]/F");
  tr_out->Branch("discr_ssvHighPur",discr_ssvHighPur,"discr_ssvHighPur[nref]/F");
  //tr_out->Branch("discr_csvMva",discr_csvMva,"discr_csvMva[nref]/F");
  tr_out->Branch("discr_csvSimple",discr_csvSimple,"discr_csvSimple[nref]/F");
  //tr_out->Branch("discr_muByIp3",discr_muByIp3,"discr_muByIp3[nref]/F");
  tr_out->Branch("discr_muByPt",discr_muByPt,"discr_muByPt[nref]/F");
  tr_out->Branch("discr_prob",discr_prob,"discr_prob[nref]/F");
  tr_out->Branch("discr_probb",discr_probb,"discr_probb[nref]/F");
  tr_out->Branch("discr_tcHighEff",discr_tcHighEff,"discr_tcHighEff[nref]/F");
  tr_out->Branch("discr_tcHighPur",discr_tcHighPur,"discr_tcHighPur[nref]/F");
  tr_out->Branch("nsvtx",nsvtx,"nsvtx[nref]/I");
  tr_out->Branch("svtxntrk",svtxntrk,"svtxntrk[nref]/I");
  tr_out->Branch("svtxdl",svtxdl,"svtxdl[nref]/F");
  tr_out->Branch("svtxdls",svtxdls,"svtxdls[nref]/F");
  tr_out->Branch("svtxm",svtxm,"svtxm[nref]/F");
  tr_out->Branch("svtxpt",svtxpt,"svtxpt[nref]/F");
  tr_out->Branch("nIPtrk",nIPtrk,"nIPtrk[nref]/I");
  tr_out->Branch("nselIPtrk",nselIPtrk,"nselIPtrk[nref]/I");
  tr_out->Branch("nIP",&nIP,"nIP/I");
  tr_out->Branch("ipJetIndex",ipJetIndex,"ipJetIndex[nIP]/I");
  tr_out->Branch("ipPt",ipPt,"ipPt[nIP]/F");
  tr_out->Branch("ipProb0",ipProb0,"ipProb0[nIP]/F");
  //tr_out->Branch("ipProb1",ipProb1,"ipProb1[nIP]/F");
  tr_out->Branch("ip2d",ip2d,"ip2d[nIP]/F");
  tr_out->Branch("ip2dSig",ip2dSig,"ip2dSig[nIP]/F");
  tr_out->Branch("ip3d",ip3d,"ip3d[nIP]/F");
  tr_out->Branch("ip3dSig",ip3dSig,"ip3dSig[nIP]/F");
  tr_out->Branch("ipDist2Jet",ipDist2Jet,"ipDist2Jet[nIP]/F");
  //tr_out->Branch("ipDist2JetSig",ipDist2JetSig,"ipDist2JetSig[nIP]/F");
  tr_out->Branch("ipClosest2Jet",ipClosest2Jet,"ipClosest2Jet[nIP]/F");  
  tr_out->Branch("mue",mue,"mue[nref]/F");
  tr_out->Branch("mupt",mupt,"mupt[nref]/F");
  tr_out->Branch("mueta",mueta,"mueta[nref]/F");
  tr_out->Branch("muphi",muphi,"muphi[nref]/F");
  tr_out->Branch("mudr",mudr,"mudr[nref]/F");
  tr_out->Branch("muptrel",muptrel,"muptre[nref]/F");
  tr_out->Branch("muchg",muchg,"muchg[nref]/I");
  tr_out->Branch("pthat",&pthat,"pthat/F");
  tr_out->Branch("beamId1",&beamId1,"beamId1/I");
  tr_out->Branch("beamId2",&beamId1,"beamId2/I");
  tr_out->Branch("refpt",refpt,"refpt[nref]/F");
  tr_out->Branch("refeta",refeta,"refeta[nref]/F");
  tr_out->Branch("refy",refy,"refy[nref]/F");
  tr_out->Branch("refphi",refphi,"refphi[nref]/F");
  tr_out->Branch("refdphijt",refdphijt,"refdphijt[nref]/F");
  tr_out->Branch("refdrjt",refdrjt,"refdrjt[nref]/F");
  tr_out->Branch("refparton_pt",refparton_pt,"refparton_pt[nref]/F");
  tr_out->Branch("refparton_flavor",refparton_flavor,"refparton_flavor[nref]/I");
  tr_out->Branch("refparton_flavorForB",refparton_flavorForB,"refparton_flavorForB[nref]/I");
  if(useGSP)tr_out->Branch("refparton_isGSP",refparton_isGSP,"refparton_isGSP[nref]/O");
  /*
    tr_out->Branch("ngen",&ngen,"ngen/I");
    tr_out->Branch("genmatchindex",genmatchindex,"genmatchindex[nref]/I");
    tr_out->Branch("genpt",genpt,"genpt[nref]/F");
    tr_out->Branch("geneta",geneta,"geneta[nref]/F");
    tr_out->Branch("geny",geny,"geny[nref]/F");
    tr_out->Branch("genphi",genphi,"genphi[nref]/F");
    tr_out->Branch("gendphijt",gendphijt,"gendphijt[nref]/F");
    tr_out->Branch("gendrjt",gendrjt,"gendrjt[nref]/F");
  */

  
  tr_out->Branch("nHLTBit",&nHLTBit,"nHLTBit/I");
  tr_out->Branch("hltBit",hltBit,"hltBit[nHLTBit]/O");
  
  tr_out->Branch("weight",&weight,"weight/D");
  tr_out->Branch("xSecWeight",&xSecWeight,"xSecWeight/D");
  tr_out->Branch("centWeight",&centWeight,"centWeight/D");
  tr_out->Branch("vzWeight",&vzWeight,"vzWeight/D");
  
  tr_out->Branch("chargedMax", chargedMax,"chargedMax[nref]/F");
  tr_out->Branch("chargedSum", chargedSum,"chargedSum[nref]/F");
  tr_out->Branch("chargedN", chargedN,"chargedN[nref]/I");
  tr_out->Branch("chargedHardSum", chargedHardSum,"chargedHardSum[nref]/F");
  tr_out->Branch("chargedHardN", chargedHardN,"chargedHardN[nref]/I");
  
  tr_out->Branch("photonMax", photonMax,"photonMax[nref]/F");
  tr_out->Branch("photonSum", photonSum,"photonSum[nref]/F");
  tr_out->Branch("photonN", photonN,"photonN[nref]/I");
  tr_out->Branch("photonHardSum", photonHardSum,"photonHardSum[nref]/F");
  tr_out->Branch("photonHardN", photonHardN,"photonHardN[nref]/I");
  
  tr_out->Branch("neutralMax", neutralMax,"neutralMax[nref]/F");
  tr_out->Branch("neutralSum", neutralSum,"neutralSum[nref]/F");
  tr_out->Branch("neutralN", neutralN,"neutralN[nref]/I");
  
  tr_out->Branch("eMax", eMax,"eMax[nref]/F");
  tr_out->Branch("eSum", eSum,"eSum[nref]/F");
  tr_out->Branch("eN", eN,"eN[nref]/I");
  
  tr_out->Branch("muMax", muMax,"muMax[nref]/F");
  tr_out->Branch("muSum", muSum,"muSum[nref]/F");
  tr_out->Branch("muN", muN,"muN[nref]/I");
  
  // Set output branch addresses.
  tr_out_skim->Branch("pvSel",&pvSel,"pvSel/I");
  tr_out_skim->Branch("hbheNoiseSel",&hbheNoiseSel,"hbheNoiseSel/I");
  tr_out_skim->Branch("spikeSel",&spikeSel,"spikeSel/I");
  tr_out_skim->Branch("collSell",&collSell,"collSell/I");
  tr_out_skim->Branch("hltAna",&hltAna,"hltAna/I");
  
  
  
  Long64_t nentries = tr_in->GetEntries();
  Long64_t nbytes = 0;
  Long64_t nbytes_skim = 0;


  int isFiltered = 1;
  if(!LCB) isFiltered=0;

  for (Long64_t i=0; i<nentries;i++) {
    
    if((double)i==totBorCentries&&LCB>0&&useQCD){
      cout<<" switching to QCD events "<<endl;
      isFiltered = 0;
    }

    if(i%100000==0) cout<<" i = "<<i<<" out of "<<nentries<<" ("<<(int)(100*(float)i/(float)nentries)<<"%)"<<endl;
    
    nbytes += tr_in->GetEntry(i);
    
    if(pthat<bounds[0]) continue;

    int ibin=-1;
    for(ibin=0;ibin<N-1;ibin++){
      if(pthat>bounds[ibin]&&pthat<bounds[ibin+1]) break;
    }
    if(pthat>bounds[N-1]) ibin=N-1;

    
    nbytes_skim += tr_in_skim->GetEntry(i);
    
    //centWeight = fCent->Integral(bin,bin+1)/hMCCent->GetBinContent(bin+1);
    //double centWeightData = fCent->Integral(bin,bin+1)/fCent->Integral(0,40)*40.;
    double centWeightData = hDataCent->GetBinContent(bin+1);
    
    //double centWeightMC = hMCCent->GetBinContent(bin+1)/hMCCent->Integral(1,40)*40.;
    double centWeightMC = hMCCent[isFiltered]->GetBinContent(bin+1);
    centWeight = centWeightData/centWeightMC;
    //    cout<<" centWeightData "<<centWeightData<<" centWeightMC "<<centWeightMC<<" centWeight "<<centWeight<<endl;
   if(centWeight<0) centWeight=0.;

    int vzbin = (int) TMath::Ceil(vz+15.);
    if(vzbin>0&&vzbin<=30)vzWeight = hDataVz->GetBinContent(vzbin)/hMCVz[isFiltered]->GetBinContent(vzbin);
    else vzWeight=0.;

    xSecWeight = xSecWeights[ibin];

    //fudge factor for vertex/x-sec non-factorization
    if(LCB>0&&isFiltered){
      if(LCB==1) xSecWeight/=0.93022;
      if(LCB==2) xSecWeight/=0.989157; 
    }

    /*  // for testing
   vzWeight = centWeightData;
   xSecWeight = centWeightMC;
    */
    weight=xSecWeight*centWeight*vzWeight;
    
    
    
    tr_out->Fill();
    tr_out_skim->Fill();
    
  }
  
  fin->Close();

  fout->mkdir("akPu3PFJetAnalyzer");
  fout->cd("akPu3PFJetAnalyzer");
  tr_out->Write();
  fout->mkdir("skimanalysis");
  fout->cd("skimanalysis");
  tr_out_skim->Write();
  fout->Close();



}

