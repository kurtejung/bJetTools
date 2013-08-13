#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <algorithm>
#include "TH1.h"
#include "TTree.h"
#include "TFile.h"
#include "TNtupleD.h"
#include "TNtuple.h"
#include "TROOT.h"
#include "TChain.h"

//These includes cause some complications in CMSSW_5_3_8_HI_patchX.  Commented out for pp.
//If you want to recalculate the JECs on the fly again, just uncomment everything in the updateJEC blocks

//#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
//#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
//#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h" 

using namespace std;

// ******* GLOBAL DECLARATIONS **********
const int QCDpthatBins = 11;
const int HFpthatBins = 5;
const int dataFiles = 10;
//***************************************


//**********************************************************
// Do the pthat weighting for the Heavy Flavor Jets
//**********************************************************

//(HFfile,infile,HFpthatBins,QCDpthatBins,'b',usePUsub)
double *heavyJetWeighting(std::string HFfile, std::string QCDfile, int HFnfiles, int QCDnfiles, char flavor, bool usePUsub){

  const int nDivisions = 11;
  double *HFweights = new double[nDivisions];
  const int weightBlocks[nDivisions+1] = {15,30,50,80,120,170,220,280,370,460,540,1200};

  TChain *chH = NULL;
  TChain *chQCD = NULL;
  if(usePUsub){
    chH = new TChain("akPu3PFJetAnalyzer/t");
    chQCD = new TChain("akPu3PFJetAnalyzer/t");
  }
  else{
    chH = new TChain("ak3PFJetAnalyzer/t");
    chQCD = new TChain("ak3PFJetAnalyzer/t");
  }
  std::ifstream instr(HFfile.c_str(), std::ifstream::in);
  std::string filename;
  for(int ifile=0; ifile<HFnfiles; ifile++){
    instr >> filename;
    chH->Add(filename.c_str());
  }
  std::ifstream instr2(QCDfile.c_str(), std::ifstream::in);
  for(int ifile=0; ifile<QCDnfiles; ifile++){
    instr2 >> filename;
    chQCD->Add(filename.c_str());
  }
  
  int parton_flavor=0;
  if(flavor=='b') parton_flavor=5;
  if(flavor=='c') parton_flavor=4;

  char* cutname = new char[200];
  char* cutfull = new char[200];
  TH1D *bjetEntr = new TH1D("bjetEntr","",1,0,1200);
  TH1D *QCDjetEntr = new TH1D("QCDjetEntr","",1,0,1200);

  TH1D *bjetEntrFULL = new TH1D("bjetEntrFULL","",1,0,1200);
  TH1D *QCDjetEntrFULL = new TH1D("QCDjetEntrFULL","",1,0,1200);
  if(parton_flavor==5 || parton_flavor==4){
    for(int i=0; i<nDivisions; i++){
      sprintf(cutname,"pthat>%d&&pthat<%d&&refpt>0&&abs(jteta)<2&&abs(refparton_flavorForB)==%d",weightBlocks[i],weightBlocks[i+1],parton_flavor);
      sprintf(cutfull,"pthat>%d&&pthat<%d",weightBlocks[i],weightBlocks[i+1]);
      cout << "jet cut: "<< cutname << endl;
      cout << "evt cut: "<< cutfull << endl;
      chQCD->Draw("jtpt>>QCDjetEntr",cutname);
      chH->Draw("jtpt>>bjetEntr",cutname);
      double qcd1 = (double)QCDjetEntr->Integral()/(double)chQCD->GetEntries(cutfull);
      double h1 = (double)bjetEntr->Integral()/(double)chH->GetEntries(cutfull);
      cout << "between pthat " << weightBlocks[i] << " and " << weightBlocks[i+1] << " has " << qcd1/h1 << " norm, or " << qcd1 << " qcd bjets and " << h1 << " hf jets" << endl;
      if(qcd1!=qcd1 || h1!=h1) HFweights[i]=0; //check on NaN via IEEE recommended method
      else{
	cout << "QCD bjet per event, pthat " << i << ": " << qcd1 << endl;
	cout << "HF bjet per event, pthat " << i << ": " << h1 << endl;
	HFweights[i] = h1/qcd1;
      }
      bjetEntr->Reset();
      QCDjetEntr->Reset();
    }
  }
  delete bjetEntr;
  delete QCDjetEntr;
  return HFweights;
}

//**********************************************************
// Count the MC events to appropriately weight the pthat bins
//**********************************************************

int *countMCevents(std::string infile, std::string HFfile, bool usePUsub, int isMC){

  TChain *ch = NULL;
  if(usePUsub) ch = new TChain("akPu3PFJetAnalyzer/t");
  else ch = new TChain("ak3PFJetAnalyzer/t");
  std::ifstream instr(infile.c_str(), std::ifstream::in);
  std::ifstream HFstr(HFfile.c_str(), std::ifstream::in);
  std::string filename;
  for(int ifile=0; ifile<QCDpthatBins; ifile++){
    instr >> filename;
    ch->Add(filename.c_str());
  }
  int *MCentries = new int[12];
  MCentries[0] = ch->GetEntries("pthat<15");
  MCentries[1] = ch->GetEntries("pthat>=15 && pthat<30");
  MCentries[2] = ch->GetEntries("pthat>=30 && pthat<50");
  MCentries[3] = ch->GetEntries("pthat>=50 && pthat<80");
  MCentries[4] = ch->GetEntries("pthat>=80 && pthat<120");
  MCentries[5] = ch->GetEntries("pthat>=120 && pthat<170");
  MCentries[6] = ch->GetEntries("pthat>=170 && pthat<220");
  MCentries[7] = ch->GetEntries("pthat>=220 && pthat<280");
  MCentries[8] = ch->GetEntries("pthat>=280 && pthat<370");
  MCentries[9] = ch->GetEntries("pthat>=370 && pthat<460");
  MCentries[10] = ch->GetEntries("pthat>=460 && pthat<540");
  MCentries[11] = ch->GetEntries("pthat>=540 && pthat<1200");

  for(int i=0; i<12; i++){
    cout << "QCD MCentries[" << i << "]: " << MCentries[i] << endl;
  }

  if(isMC>1){
    TChain *HFch = NULL;
    if(usePUsub) HFch = new TChain("akPu3PFJetAnalyzer/t");
    else HFch = new TChain("ak3PFJetAnalyzer/t");
    for(int ifile=0; ifile<HFpthatBins; ifile++){
      HFstr >> filename;
      HFch->Add(filename.c_str());
    }
    double *HFweight = NULL;
    if(isMC==2) HFweight = heavyJetWeighting(HFfile,infile,HFpthatBins,QCDpthatBins,'b',usePUsub);
    if(isMC==3) HFweight = heavyJetWeighting(HFfile,infile,HFpthatBins,QCDpthatBins,'c',usePUsub);
    for(int i=0; i<11; i++){
      cout << "HFweight[" << i << "]: " << HFweight[i] << endl;
    }
    
    double tempEntr[12] = {0,0,0,0,0,0,0,0,0,0,0,0};
    
    /* MCentries[0] += (double)HFch->GetEntries("pthat<15")*0;
    MCentries[1] += (double)HFch->GetEntries("pthat>=15 && pthat<30")*HFweight[0];
    MCentries[2] += (double)HFch->GetEntries("pthat>=30 && pthat<50")*HFweight[1];
    MCentries[3] += (double)HFch->GetEntries("pthat>=50 && pthat<80")*HFweight[2];
    MCentries[4] += (double)HFch->GetEntries("pthat>=80 && pthat<120")*HFweight[3];
    MCentries[5] += (double)HFch->GetEntries("pthat>=120 && pthat<170")*HFweight[4];
    MCentries[6] += (double)HFch->GetEntries("pthat>=170 && pthat<220")*HFweight[5];
    MCentries[7] += (double)HFch->GetEntries("pthat>=220 && pthat<280")*HFweight[6];
    MCentries[8] += (double)HFch->GetEntries("pthat>=280 && pthat<370")*HFweight[7];
    MCentries[9] += (double)HFch->GetEntries("pthat>=370 && pthat<460")*HFweight[8];
    MCentries[10] += (double)HFch->GetEntries("pthat>=460 && pthat<540")*HFweight[9];
    MCentries[11] += (double)HFch->GetEntries("pthat>=540 && pthat<1200")*HFweight[10];*/
    tempEntr[0] = HFch->GetEntries("pthat<15");
    tempEntr[1] = HFch->GetEntries("pthat>=15 && pthat<30");
    tempEntr[2] = HFch->GetEntries("pthat>=30 && pthat<50");
    tempEntr[3] = HFch->GetEntries("pthat>=50 && pthat<80");
    tempEntr[4] = HFch->GetEntries("pthat>=80 && pthat<120");
    tempEntr[5] = HFch->GetEntries("pthat>=120 && pthat<170");
    tempEntr[6] = HFch->GetEntries("pthat>=170 && pthat<220");
    tempEntr[7] = HFch->GetEntries("pthat>=220 && pthat<280");
    tempEntr[8] = HFch->GetEntries("pthat>=280 && pthat<370");
    tempEntr[9] = HFch->GetEntries("pthat>=370 && pthat<460");
    tempEntr[10] = HFch->GetEntries("pthat>=460 && pthat<540");
    tempEntr[11] = HFch->GetEntries("pthat>=540 && pthat<1200");

    for(int i=1; i<12; i++){
      cout << "HF entries[" << i << "]: " << tempEntr[i] << endl;
      MCentries[i] += tempEntr[i]*HFweight[i-1];
    }
    for(int i=0; i<12; i++){
      cout << "Effective MCentries[" << i << "]: " << MCentries[i] << endl;
    }
  }
  
  return MCentries;
}

//**********************************************************
// Trigger-Combine the data in order to unfold properly later
//**********************************************************

//[0] = Jet20, [1] = Jet40, [2] = Jet60, [3] = Jet80
double trigComb(bool *triggerDecision, double *pscl){
  double weight=0;
  // if(triggerDecision[0] && !triggerDecision[1] && !triggerDecision[2] && !triggerDecision[3]) weight = 1./(1./pscl[0]); //Removing finnicky Jet20 sample
  if(triggerDecision[1] && !triggerDecision[2] && !triggerDecision[3]) weight = 1./(1./pscl[1]);
  if(triggerDecision[2] && !triggerDecision[3]) weight = 1./(1./pscl[1] + 1./pscl[2] - (1./(pscl[1]*pscl[2])));
  if(triggerDecision[3]) weight = 1.;
  return weight;
}

//**********************************************************
// "get" the trigger prescales by counting trigger overlap
//**********************************************************

double* getPscls(std::string infile, int nFiles, bool usePUsub){
  
  TChain *dataCH = NULL;
  if(usePUsub){
    dataCH = new TChain("akPu3PFJetAnalyzer/t");
  }
  else dataCH = new TChain("ak3PFJetAnalyzer/t");
  TChain *dataCH2 = new TChain("hltanalysis/HltTree");
  std::ifstream instr(infile.c_str(), std::ifstream::in);
  std::string filename;
  for(int ifile=0; ifile<nFiles; ifile++){
    instr >> filename;
    dataCH->Add(filename.c_str());
    dataCH2->Add(filename.c_str());
  }
  dataCH->AddFriend(dataCH2, "hltanalysis/HltTree");
  //Set up trigger combination prescales for data
  double ov1, ov2, ov3, ov4;
  ov1 = dataCH->GetEntries("HLT_PAJet20_NoJetID_v1 && HLT_PAJet80_NoJetID_v1");
  ov2 = dataCH->GetEntries("HLT_PAJet40_NoJetID_v1 && HLT_PAJet80_NoJetID_v1");
  ov3 = dataCH->GetEntries("HLT_PAJet60_NoJetID_v1 && HLT_PAJet80_NoJetID_v1");
  ov4 = dataCH->GetEntries("HLT_PAJet80_NoJetID_v1");
  double *pscls = new double[4];
  pscls[0] = ov4/ov1;
  pscls[1] = ov4/ov2;
  pscls[2] = ov4/ov3;
  pscls[3] = 1.;
  return pscls;
}

//**********************************************************
// ~~~ MAIN PROGRAM ~~~
//**********************************************************

void analyzeTrees(int isRecopp=1, int ppPbPb=0, int isMuTrig=0, int isMC=1, int doNtuples=1, int doJets=1, int doTracks=1, int updateJEC=0, int cbin=-1,int useGSP=2, int jetTrig=0, bool ExpandedTree=false, bool usePUsub=0)
{
  // isMC=0 --> Real data, ==1 --> QCD, ==2 --> bJet, ==3 --> cJet
  Float_t minJetPt=15.;
  
  if (isMuTrig) minJetPt=30;
  Float_t maxJetEta=2;
  Float_t minMuPt=5;
  
  // cbin = -1 --> 0-100%
  // cbin = 0 --> 0-20%
  // cbin = 1 --> 20-50%
  // cbin =2 --> 50-100%
  if(!ppPbPb) cbin=-1;
  int useWeight=1;

  cout << "Analyzing Trees! Assuming " << QCDpthatBins << " QCD pthat bins, and " << HFpthatBins << " B/C pthat bins." << endl;

  int pthatbin[QCDpthatBins+1] = {15,30,50,80,120,170,220,280,370,460,540,1200};
  double w = 1.;
  double wght[QCDpthatBins+1]={0.2034, 1.075E-02, 1.025E-03, 9.865E-05, 1.129E-05, 1.465E-06, 2.837E-07, 5.323E-08, 5.934E-09, 8.125E-10, 1.467E-10};

  TFile *fin=NULL;
  std::string infile;
  std::string HFfile;
  int *MCentr = NULL;
  double *pscls = NULL;

  //PbPb File load
  if(ppPbPb){
    if(isMC==0){
      if(jetTrig==1)fin = new TFile("/data_CMS/cms/mnguyen/bTaggingOutput/PbPbData/pbpbDataJet80_hiRegitSVHighPurity_pt30by3_restrictMixTripletA_jpHICalibRepass/merged_bTagAnalyzers_all.root");
      if(jetTrig==2)fin = new TFile("/data_CMS/cms/mnguyen/bTaggingOutput/PbPbData/pbpbDataJet65_hiRegitSVHighPurity_pt30by3_restrictMixTripletA_jpHICalibRepass/merged_bTagAnalyzers_all.root");
    }
    else if(isMC==1) fin = new TFile("/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/merged_bjetAnalyzers_hiReco_offPV_pt30by3_oldHydjet_restrictMixTripletA_ipHICalibCentWeight_weighted_qcd.root");
    else if(isMC==2) fin = new TFile("/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/merged_bjetAnalyzers_hiReco_offPV_pt30by3_newHydjet_restrictMixTripletA_ipHICalibCentWeight_weighted_bJetPlusQCD.root"); 
    else if(isMC==3)fin = new TFile("/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/merged_bjetAnalyzers_hiReco_offPV_pt30by3_newHydjet_restrictMixTripletA_ipHICalibCentWeight_weighted_cJetPlusQCD.root");
  }
  
  else{ //pp File Load
    if(!isMC){ 
      infile = "ppBForestList.txt";
    }
    else if(isMC){
      infile = "pythiaMCfilelist.txt";
      if(isMC==2){
	HFfile = "pythiaBJetMClist.txt";
      }
      else if(isMC==3){
	HFfile = "pythiaCJetMClist.txt";
      }
      else if(isMC>3){ 
	cout << "I don't understand this MC number!" << endl;
	exit(0);
      }
    }
  }
  if(!isMC && !ppPbPb){
    pscls = getPscls(infile,QCDpthatBins,usePUsub);
  }
  
  int dupRuns[6] = {181912,181913,181938,181950,181985,182124};
  
  std::vector<int> usedEvents[6];
  int nDup=0;

  //Declaration of leaves types                  
  Int_t           evt;
  Int_t           run;
  Int_t           bin;
  Int_t           lumi;
  Float_t         hf;
  Float_t           vz;
  Int_t           nref;
  Float_t         rawpt[1000];
  Float_t         jtpt[1000];
  Float_t         jteta[1000];
  Float_t         jty[1000];
  Float_t         jtphi[1000];
  Float_t         jtpu[1000];
  Float_t         discr_ssvHighEff[1000];
  Float_t         discr_ssvHighPur[1000];
  //Float_t         discr_csvMva[1000];
  Float_t         discr_csvSimple[1000];
  // Float_t         discr_muByIp3[1000];
  Float_t         discr_muByPt[1000];
  Float_t         discr_prob[1000];
  Float_t         discr_probb[1000];
  Float_t         discr_tcHighEff[1000];
  Float_t         discr_tcHighPur[1000];
  Int_t           nsvtx[1000];
  Int_t           svtxntrk[1000];
  Float_t         svtxdl[1000];
  Float_t         svtxdls[1000];
  Float_t         svtxm[1000];
  Float_t         svtxpt[1000];
  Int_t           nIPtrk[1000];
  Int_t           nselIPtrk[1000];
  Int_t   nIP;
  Int_t   ipJetIndex[10000];
  Float_t ipPt[10000];
  Float_t ipProb0[10000];
  Float_t ipProb1[10000];
  Float_t ip2d[10000];
  Float_t ip2dSig[10000];
  Float_t ip3d[10000];
  Float_t ip3dSig[10000];
  Float_t ipDist2Jet[10000];
  Float_t ipDist2JetSig[10000];
  Float_t ipClosest2Jet[10000];
  //Float_t         mue[1000];
  Float_t         mupt[1000];
  Float_t         muptPF[1000];
  Float_t         mueta[1000];
  Float_t         muphi[1000];
  //Float_t         mudr[1000];
  Float_t         muptrel[1000];
  //Int_t           muchg[1000];
  Float_t         pthat;
  Int_t           beamId1;
  Int_t           beamId2;
  Float_t         refpt[1000];
  Float_t         refeta[1000];
  Float_t         refy[1000];
  Float_t         refphi[1000];
  Float_t         refdphijt[1000];
  Float_t         refdrjt[1000];
  Float_t         refparton_pt[1000];
  Int_t           refparton_flavor[1000];
  Int_t           refparton_flavorForB[1000];
  Bool_t           refparton_isGSP[1000];
  Int_t 	  HLT_PAJet20_NoJetID_v1;
  Int_t 	  HLT_PAJet40_NoJetID_v1;
  Int_t 	  HLT_PAJet60_NoJetID_v1;
  Int_t 	  HLT_PAJet80_NoJetID_v1;
  Int_t 	  HLT_PAJet100_NoJetID_v1;
  Int_t           pVertexFilterCutGplusUpsPP;
  Int_t           pPAcollisionEventSelectionPA;
  Int_t           pHBHENoiseFilter;
  Int_t           pprimaryvertexFilter;

  /*
    Int_t           ngen;
    Int_t           genmatchindex[1000];
    Float_t         genpt[1000];
    Float_t         geneta[1000];
    Float_t         geny[1000];
    Float_t         genphi[1000];
    Float_t         gendphijt[1000];
    Float_t         gendrjt[1000];
  */

  float chargedMax[1000];
  float photonMax[1000];
  float neutralMax[1000];  
  float chargedSum[1000];
  float photonSum[1000];
  float neutralSum[1000];  
  float muSum[1000];  
  float eSum[1000];  

  Double_t         weight, xSecWeight, centWeight, vzWeight;
  
  int nHLTBit;
  bool hltBit[12];

  Int_t pvSel;
  Int_t hbheNoiseSel;
  Int_t spikeSel;
  Int_t collSell;

  TFile *fout=NULL;
  if(ppPbPb){

    if(cbin==-1){
      if(jetTrig==0){
	if(isMC==0) fout = new TFile("histos/PbPbdata_pt30by3_jpHICalibRepass_withDup_PU.root","recreate");
	else if(isMC==1) fout = new TFile("histos/PbPbQCDMC_pt30by3_ipHICalibCentWeight_noTrig.root","recreate");
	else if(isMC==2) fout = new TFile("histos/PbPbBMC_pt30by3_ipHICalibCentWeight_noTrig.root","recreate");
	else if(isMC==3) fout = new TFile("histos/PbPbCMC_pt30by3_ipHICalibCentWeight_noTrig.root","recreate");
      }
      if(jetTrig==1){
	if(isMC==0) fout = new TFile("histos/PbPbdata_pt30by3_jpHICalibRepass_withDup_PU.root","recreate");
	else if(isMC==1) fout = new TFile("histos/PbPbQCDMC_pt30by3_ipHICalibCentWeight.root","recreate");
	else if(isMC==2) fout = new TFile("histos/PbPbBMC_pt30by3_ipHICalibCentWeight.root","recreate");
	else if(isMC==3) fout = new TFile("histos/PbPbCMC_pt30by3_ipHICalibCentWeight.root","recreate");
      }
      if(jetTrig==2){
	if(isMC==0) fout = new TFile("histos/PbPbdata_pt30by3_jpHICalibRepass_withDup_PU_jet65.root","recreate");
	else if(isMC==1) fout = new TFile("histos/PbPbQCDMC_pt30by3_ipHICalibCentWeight_jet65.root","recreate");
	else if(isMC==2) fout = new TFile("histos/PbPbBMC_pt30by3_ipHICalibCentWeight_jet65.root","recreate");
	else if(isMC==3) fout = new TFile("histos/PbPbCMC_pt30by3_ipHICalibCentWeight_jet65.root","recreate");
      }
    }
    else{
      if(isMC==0) fout = new TFile(Form("histos/PbPbdata_%d_pt30by3_jpHICalibRepass.root",cbin),"recreate");
      else if(isMC==1) fout = new TFile(Form("histos/PbPbQCDMC_%d_pt30by3_ipHICalibCentWeight.root",cbin),"recreate");
      else if(isMC==2) fout = new TFile(Form("histos/PbPbBMC_%d_pt30by3_ipHICalibCentWeight.root",cbin),"recreate");
      else if(isMC==3) fout = new TFile(Form("histos/PbPbCMC_%d_pt30by3_ipHICalibCentWeight.root",cbin),"recreate");

    }

  }
  else{
    if( isRecopp&& isMuTrig) { // pp reco, muon triggered
      if(isMC)fout=new TFile("histos/ppMC_ppReco_muTrig_noIPupperCut.root","recreate");
      else fout=new TFile("histos/ppdata_ppReco_muTrig_noIPupperCut.root","recreate");
    }
    else if ( isRecopp && !isMuTrig && usePUsub) { // pp reco, jet triggered
      if(isMC==1) fout=new TFile("histos/ppMC_ppReco_akPu3PF_QCDjetTrig_noIPupperCut.root","recreate");
      else if(isMC==2) fout=new TFile("histos/ppMC_ppReco_akPu3PF_BjetTrig_noIPupperCut.root","recreate");
      else if(isMC==3) fout=new TFile("histos/ppMC_ppReco_akPu3PF_CjetTrig_noIPupperCut.root","recreate");
      else fout=new TFile("histos/ppdata_ppReco_akPu3PF_jetTrig_noIPupperCut.root","recreate");
    }
    else if( isRecopp&& !isMuTrig && !usePUsub){
      if(isMC==1) fout=new TFile("histos/ppMC_ppReco_ak3PF_gsp2_QCDjetTrig_noIPupperCut.root","recreate");
      else if(isMC==2) fout=new TFile("histos/ppMC_ppReco_ak3PF_gsp2_BjetTrig_noIPupperCut.root","recreate");
      else if(isMC==3) fout=new TFile("histos/ppMC_ppReco_ak3PF_gsp2_CjetTrig_noIPupperCut.root","recreate");
      else fout=new TFile("histos/ppdata_ppReco_ak3PF_gsp2_jetTrig_noIPupperCut.root","recreate");
    }
    else if (!isRecopp&& isMuTrig) { // hi reco, muon triggered
      if(isMC)fout=new TFile("histos/ppMC_hiReco_muTrig_noIPupperCut.root","recreate");
      else fout=new TFile("histos/ppdata_hiReco_muTrig_noIPupperCut.root","recreate");
    } 
    else if (!isRecopp&&!isMuTrig) { // hi reco, jet triggered
      if(isMC)fout=new TFile("histos/ppMC_hiReco_jetTrig_addGSP_up.root","recreate");
      else fout=new TFile("histos/ppdata_hiReco_jetTrig_regPFforJets.root","recreate");
    }
  }

  TH1D *hbin = new TH1D("hbin","hbin",40,-0.5,39.5);
  TH1D *hbinw = new TH1D("hbinw","hbinw",40,-0.5,39.5);
  hbin->Sumw2(); hbinw->Sumw2(); 

  TH1D *hvz = new TH1D("hvz","hvz",120,-15.,15.);
  TH1D *hvzw = new TH1D("hvzw","hvzw",120,-15.,15.);
  hvz->Sumw2(); hvzw->Sumw2(); 

  TH1D *hjtpt = new TH1D("hjtpt","hjtpt",68,80,330);
  TH1D *hjtptB = new TH1D("hjtptB","hjtptB",68,80,330);
  TH1D *hjtptC = new TH1D("hjtptC","hjtptC",68,80,330);
  TH1D *hjtptL = new TH1D("hjtptL","hjtptL",68,80,330);
  TH1D *hjtptU = new TH1D("hjtptU","hjtptU",68,80,330);
  hjtpt->Sumw2(); hjtptB->Sumw2(); hjtptC->Sumw2(); hjtptL->Sumw2(); hjtptU->Sumw2();

  TH1D *hrawpt = new TH1D("hrawpt","hrawpt",68,80,330);
  TH1D *hrawptB = new TH1D("hrawptB","hrawptB",68,80,330);
  TH1D *hrawptC = new TH1D("hrawptC","hrawptC",68,80,330);
  TH1D *hrawptL = new TH1D("hrawptL","hrawptL",68,80,330);
  hrawpt->Sumw2(); hrawptB->Sumw2(); hrawptC->Sumw2(); hrawptL->Sumw2();

  TH1D *hjteta = new TH1D("hjteta","hjteta",40,-2,2);
  TH1D *hjtetaB = new TH1D("hjtetaB","hjtetaB",40,-2,2);
  TH1D *hjtetaC = new TH1D("hjtetaC","hjtetaC",40,-2,2);
  TH1D *hjtetaL = new TH1D("hjtetaL","hjtetaL",40,-2,2);
  hjteta->Sumw2(); hjtetaB->Sumw2(); hjtetaC->Sumw2(); hjtetaL->Sumw2(); 

  TH1D *hjtphi = new TH1D("hjtphi","hjtphi",40,-1.*acos(-1.),acos(-1.));
  TH1D *hjtphiB = new TH1D("hjtphiB","hjtphiB",40,-1.*acos(-1.),acos(-1.));
  TH1D *hjtphiC = new TH1D("hjtphiC","hjtphiC",40,-1.*acos(-1.),acos(-1.));
  TH1D *hjtphiL = new TH1D("hjtphiL","hjtphiL",40,-1.*acos(-1.),acos(-1.));
  hjtphi->Sumw2(); hjtphiB->Sumw2(); hjtphiC->Sumw2(); hjtphiL->Sumw2(); 

  TH1D *hdiscr_csvSimple = new TH1D("hdiscr_csvSimple","hdiscr_csvSimple",25,0,1);
  TH1D *hdiscr_csvSimpleB = new TH1D("hdiscr_csvSimpleB","hdiscr_csvSimpleB",25,0,1);
  TH1D *hdiscr_csvSimpleC = new TH1D("hdiscr_csvSimpleC","hdiscr_csvSimpleC",25,0,1);
  TH1D *hdiscr_csvSimpleL = new TH1D("hdiscr_csvSimpleL","hdiscr_csvSimpleL",25,0,1);
  hdiscr_csvSimple->Sumw2(); hdiscr_csvSimpleB->Sumw2(); hdiscr_csvSimpleC->Sumw2(); hdiscr_csvSimpleL->Sumw2();
  
  TH1D *hdiscr_prob = new TH1D("hdiscr_prob","hdiscr_prob",25,0,2.5);
  TH1D *hdiscr_probB = new TH1D("hdiscr_probB","hdiscr_probB",25,0,2.5);
  TH1D *hdiscr_probC = new TH1D("hdiscr_probC","hdiscr_probC",25,0,2.5);
  TH1D *hdiscr_probL = new TH1D("hdiscr_probL","hdiscr_probL",25,0,2.5);
  hdiscr_prob->Sumw2(); hdiscr_probB->Sumw2(); hdiscr_probC->Sumw2(); hdiscr_probL->Sumw2();
  
  TH1D *hdiscr_ssvHighEff = new TH1D("hdiscr_ssvHighEff","hdiscr_ssvHighEff",25,1,6);
  TH1D *hdiscr_ssvHighEffB = new TH1D("hdiscr_ssvHighEffB","hdiscr_ssvHighEffB",25,1,6);
  TH1D *hdiscr_ssvHighEffC = new TH1D("hdiscr_ssvHighEffC","hdiscr_ssvHighEffC",25,1,6);
  TH1D *hdiscr_ssvHighEffL = new TH1D("hdiscr_ssvHighEffL","hdiscr_ssvHighEffL",25,1,6);
  hdiscr_ssvHighEff->Sumw2(); hdiscr_ssvHighEffB->Sumw2(); hdiscr_ssvHighEffC->Sumw2(); hdiscr_ssvHighEffL->Sumw2();
  
  TH1D *hdiscr_ssvHighPur = new TH1D("hdiscr_ssvHighPur","hdiscr_ssvHighPur",25,1,6);
  TH1D *hdiscr_ssvHighPurB = new TH1D("hdiscr_ssvHighPurB","hdiscr_ssvHighPurB",25,1,6);
  TH1D *hdiscr_ssvHighPurC = new TH1D("hdiscr_ssvHighPurC","hdiscr_ssvHighPurC",25,1,6);
  TH1D *hdiscr_ssvHighPurL = new TH1D("hdiscr_ssvHighPurL","hdiscr_ssvHighPurL",25,1,6);
  hdiscr_ssvHighPur->Sumw2(); hdiscr_ssvHighPurB->Sumw2(); hdiscr_ssvHighPurC->Sumw2(); hdiscr_ssvHighPurL->Sumw2();

  TH1D *hdiscr_tcHighEff = new TH1D("hdiscr_tcHighEff","hdiscr_tcHighEff",25,1,6);
  TH1D *hdiscr_tcHighEffB = new TH1D("hdiscr_tcHighEffB","hdiscr_tcHighEffB",25,1,6);
  TH1D *hdiscr_tcHighEffC = new TH1D("hdiscr_tcHighEffC","hdiscr_tcHighEffC",25,1,6);
  TH1D *hdiscr_tcHighEffL = new TH1D("hdiscr_tcHighEffL","hdiscr_tcHighEffL",25,1,6);
  hdiscr_tcHighEff->Sumw2(); hdiscr_tcHighEffB->Sumw2(); hdiscr_tcHighEffC->Sumw2(); hdiscr_tcHighEffL->Sumw2();
  
  TH1D *hdiscr_tcHighPur = new TH1D("hdiscr_tcHighPur","hdiscr_tcHighPur",25,1,6);
  TH1D *hdiscr_tcHighPurB = new TH1D("hdiscr_tcHighPurB","hdiscr_tcHighPurB",25,1,6);
  TH1D *hdiscr_tcHighPurC = new TH1D("hdiscr_tcHighPurC","hdiscr_tcHighPurC",25,1,6);
  TH1D *hdiscr_tcHighPurL = new TH1D("hdiscr_tcHighPurL","hdiscr_tcHighPurL",25,1,6);
  hdiscr_tcHighPur->Sumw2(); hdiscr_tcHighPurB->Sumw2(); hdiscr_tcHighPurC->Sumw2(); hdiscr_tcHighPurL->Sumw2();
  
  TH1D *hnsvtx = new TH1D("hnsvtx","hnsvtx",6,-0.5,5.5);
  TH1D *hnsvtxB = new TH1D("hnsvtxB","hnsvtxB",6,-0.5,5.5);
  TH1D *hnsvtxC = new TH1D("hnsvtxC","hnsvtxC",6,-0.5,5.5);
  TH1D *hnsvtxL = new TH1D("hnsvtxL","hnsvtxL",6,-0.5,5.5);
  hnsvtx->Sumw2(); hnsvtxB->Sumw2(); hnsvtxC->Sumw2(); hnsvtxL->Sumw2();
  
  TH1D *hsvtxntrk = new TH1D("hsvtxntrk","hsvtxntrk",12,-0.5,11.5);
  TH1D *hsvtxntrkB = new TH1D("hsvtxntrkB","hsvtxntrkB",12,-0.5,11.5);
  TH1D *hsvtxntrkC = new TH1D("hsvtxntrkC","hsvtxntrkC",12,-0.5,11.5);
  TH1D *hsvtxntrkL = new TH1D("hsvtxntrkL","hsvtxntrkL",12,-0.5,11.5);
  hsvtxntrk->Sumw2(); hsvtxntrkB->Sumw2(); hsvtxntrkC->Sumw2(); hsvtxntrkL->Sumw2();

  TH1D *hsvtxdl = new TH1D("hsvtxdl","hsvtxdl",20,0,10);
  TH1D *hsvtxdlB = new TH1D("hsvtxdlB","hsvtxdlB",20,0,10);
  TH1D *hsvtxdlC = new TH1D("hsvtxdlC","hsvtxdlC",20,0,10);
  TH1D *hsvtxdlL = new TH1D("hsvtxdlL","hsvtxdlL",20,0,10);
  hsvtxdl->Sumw2(); hsvtxdlB->Sumw2(); hsvtxdlC->Sumw2(); hsvtxdlL->Sumw2();

  TH1D *hsvtxdls = new TH1D("hsvtxdls","hsvtxdls",40,0,80);
  TH1D *hsvtxdlsB = new TH1D("hsvtxdlsB","hsvtxdlsB",40,0,80);
  TH1D *hsvtxdlsC = new TH1D("hsvtxdlsC","hsvtxdlsC",40,0,80);
  TH1D *hsvtxdlsL = new TH1D("hsvtxdlsL","hsvtxdlsL",40,0,80);
  hsvtxdls->Sumw2(); hsvtxdlsB->Sumw2(); hsvtxdlsC->Sumw2(); hsvtxdlsL->Sumw2();
  
  TH1D *hsvtxm = new TH1D("hsvtxm","hsvtxm",32,0,8);
  TH1D *hsvtxmB = new TH1D("hsvtxmB","hsvtxmB",32,0,8);
  TH1D *hsvtxmC = new TH1D("hsvtxmC","hsvtxmC",32,0,8);
  TH1D *hsvtxmL = new TH1D("hsvtxmL","hsvtxmL",32,0,8);
  hsvtxm->Sumw2(); hsvtxmB->Sumw2(); hsvtxmC->Sumw2(); hsvtxmL->Sumw2(); 
  
  TH1D *hsvtxmSV3 = new TH1D("hsvtxmSV3","hsvtxmSV3",32,0,8);
  TH1D *hsvtxmSV3B = new TH1D("hsvtxmSV3B","hsvtxmSV3B",32,0,8);
  TH1D *hsvtxmSV3C = new TH1D("hsvtxmSV3C","hsvtxmSV3C",32,0,8);
  TH1D *hsvtxmSV3L = new TH1D("hsvtxmSV3L","hsvtxmSV3L",32,0,8);
  hsvtxmSV3->Sumw2(); hsvtxmSV3B->Sumw2(); hsvtxmSV3C->Sumw2(); hsvtxmSV3L->Sumw2(); 
  
  TH1D *hsvtxpt = new TH1D("hsvtxpt","hsvtxpt",20,0,100);
  TH1D *hsvtxptB = new TH1D("hsvtxptB","hsvtxptB",20,0,100);
  TH1D *hsvtxptC = new TH1D("hsvtxptC","hsvtxptC",20,0,100);
  TH1D *hsvtxptL = new TH1D("hsvtxptL","hsvtxptL",20,0,100);
  hsvtxpt->Sumw2(); hsvtxptB->Sumw2(); hsvtxptC->Sumw2(); hsvtxptL->Sumw2(); 
  
  TH1D *hsvtxptSV3 = new TH1D("hsvtxptSV3","hsvtxptSV3",20,0,100);
  TH1D *hsvtxptSV3B = new TH1D("hsvtxptSV3B","hsvtxptSV3B",20,0,100);
  TH1D *hsvtxptSV3C = new TH1D("hsvtxptSV3C","hsvtxptSV3C",20,0,100);
  TH1D *hsvtxptSV3L = new TH1D("hsvtxptSV3L","hsvtxptSV3L",20,0,100);
  hsvtxptSV3->Sumw2(); hsvtxptSV3B->Sumw2(); hsvtxptSV3C->Sumw2(); hsvtxptSV3L->Sumw2(); 
  
  TH1D *hnIPtrk = new TH1D("hnIPtrk","hnIPtrk",100,0,100);
  TH1D *hnIPtrkB = new TH1D("hnIPtrkB","hnIPtrkB",100,0,100);
  TH1D *hnIPtrkC = new TH1D("hnIPtrkC","hnIPtrkC",100,0,100);
  TH1D *hnIPtrkL = new TH1D("hnIPtrkL","hnIPtrkL",100,0,100);
  hnIPtrk->Sumw2(); hnIPtrkB->Sumw2(); hnIPtrkC->Sumw2(); hnIPtrkL->Sumw2(); 
  
  TH1D *hnselIPtrk = new TH1D("hnselIPtrk","hnselIPtrk",100,0,100);
  TH1D *hnselIPtrkB = new TH1D("hnselIPtrkB","hnselIPtrkB",100,0,100);
  TH1D *hnselIPtrkC = new TH1D("hnselIPtrkC","hnselIPtrkC",100,0,100);
  TH1D *hnselIPtrkL = new TH1D("hnselIPtrkL","hnselIPtrkL",100,0,100);
  hnselIPtrk->Sumw2(); hnselIPtrkB->Sumw2(); hnselIPtrkC->Sumw2(); hnselIPtrkL->Sumw2(); 
  
  TH1D *hmuptrel = new TH1D("hmuptrel","hmuptrel",40,0,4);
  TH1D *hmuptrelB = new TH1D("hmuptrelB","hmuptrelB",40,0,4);
  TH1D *hmuptrelC = new TH1D("hmuptrelC","hmuptrelC",40,0,4);
  TH1D *hmuptrelL = new TH1D("hmuptrelL","hmuptrelL",40,0,4);
  hmuptrel->Sumw2(); hmuptrelB->Sumw2(); hmuptrelC->Sumw2(); hmuptrelL->Sumw2(); 
  
  TH1D *hmuptrelSV2 = new TH1D("hmuptrelSV2","hmuptrelSV2",40,0,4);
  TH1D *hmuptrelSV2B = new TH1D("hmuptrelSV2B","hmuptrelSV2B",40,0,4);
  TH1D *hmuptrelSV2C = new TH1D("hmuptrelSV2C","hmuptrelSV2C",40,0,4);
  TH1D *hmuptrelSV2L = new TH1D("hmuptrelSV2L","hmuptrelSV2L",40,0,4);
  hmuptrelSV2->Sumw2(); hmuptrelSV2B->Sumw2(); hmuptrelSV2C->Sumw2(); hmuptrelSV2L->Sumw2(); 
  
  TH1D *hmuptrelSV3 = new TH1D("hmuptrelSV3","hmuptrelSV3",40,0,4);
  TH1D *hmuptrelSV3B = new TH1D("hmuptrelSV3B","hmuptrelSV3B",40,0,4);
  TH1D *hmuptrelSV3C = new TH1D("hmuptrelSV3C","hmuptrelSV3C",40,0,4);
  TH1D *hmuptrelSV3L = new TH1D("hmuptrelSV3L","hmuptrelSV3L",40,0,4);
  hmuptrelSV3->Sumw2(); hmuptrelSV3B->Sumw2(); hmuptrelSV3C->Sumw2(); hmuptrelSV3L->Sumw2(); 
  
  TH1D *hipPt = new TH1D("hipPt","hipPt",40,0,40);
  TH1D *hipPtB = new TH1D("hipPtB","hipPtB",40,0,40);
  TH1D *hipPtC = new TH1D("hipPtC","hipPtC",40,0,40);
  TH1D *hipPtL = new TH1D("hipPtL","hipPtL",40,0,40);
  hipPt->Sumw2(); hipPtB->Sumw2(); hipPtC->Sumw2(); hipPtL->Sumw2(); 
  
  TH1D *hipProb0 = new TH1D("hipProb0","hipProb0",40,-1,1);
  TH1D *hipProb0B = new TH1D("hipProb0B","hipProb0B",40,-1,1);
  TH1D *hipProb0C = new TH1D("hipProb0C","hipProb0C",40,-1,1);
  TH1D *hipProb0L = new TH1D("hipProb0L","hipProb0L",40,-1,1);
  hipProb0->Sumw2(); hipProb0B->Sumw2(); hipProb0C->Sumw2(); hipProb0L->Sumw2(); 
  
  TH1D *hipProb1 = new TH1D("hipProb1","hipProb1",40,-1,1);
  TH1D *hipProb1B = new TH1D("hipProb1B","hipProb1B",40,-1,1);
  TH1D *hipProb1C = new TH1D("hipProb1C","hipProb1C",40,-1,1);
  TH1D *hipProb1L = new TH1D("hipProb1L","hipProb1L",40,-1,1);
  hipProb1->Sumw2(); hipProb1B->Sumw2(); hipProb1C->Sumw2(); hipProb1L->Sumw2(); 
  
  TH1D *hip2d = new TH1D("hip2d","hip2d",40,-0.1,0.1);
  TH1D *hip2dB = new TH1D("hip2dB","hip2dB",40,-0.1,0.1);
  TH1D *hip2dC = new TH1D("hip2dC","hip2dC",40,-0.1,0.1);
  TH1D *hip2dL = new TH1D("hip2dL","hip2dL",40,-0.1,0.1);
  hip2d->Sumw2(); hip2dB->Sumw2(); hip2dC->Sumw2(); hip2dL->Sumw2(); 
  
  TH1D *hip2dSig = new TH1D("hip2dSig","hip2dSig",70,-35,35);
  TH1D *hip2dSigB = new TH1D("hip2dSigB","hip2dSigB",70,-35,35);
  TH1D *hip2dSigC = new TH1D("hip2dSigC","hip2dSigC",70,-35,35);
  TH1D *hip2dSigL = new TH1D("hip2dSigL","hip2dSigL",70,-35,35);
  hip2dSig->Sumw2(); hip2dSigB->Sumw2(); hip2dSigC->Sumw2(); hip2dSigL->Sumw2(); 

  TH1D *hip2d1 = new TH1D("hip2d1","hip2d1",40,-0.1,0.1);
  TH1D *hip2d1B = new TH1D("hip2d1B","hip2d1B",40,-0.1,0.1);
  TH1D *hip2d1C = new TH1D("hip2d1C","hip2d1C",40,-0.1,0.1);
  TH1D *hip2d1L = new TH1D("hip2d1L","hip2d1L",40,-0.1,0.1);
  hip2d1->Sumw2(); hip2d1B->Sumw2(); hip2d1C->Sumw2(); hip2d1L->Sumw2(); 

  TH1D *hip2dSig1 = new TH1D("hip2dSig1","hip2dSig1",70,-35,35);
  TH1D *hip2dSig1B = new TH1D("hip2dSig1B","hip2dSig1B",70,-35,35);
  TH1D *hip2dSig1C = new TH1D("hip2dSig1C","hip2dSig1C",70,-35,35);
  TH1D *hip2dSig1L = new TH1D("hip2dSig1L","hip2dSig1L",70,-35,35);
  hip2dSig1->Sumw2(); hip2dSig1B->Sumw2(); hip2dSig1C->Sumw2(); hip2dSig1L->Sumw2(); 

  TH1D *hip2d2 = new TH1D("hip2d2","hip2d2",40,-0.1,0.1);
  TH1D *hip2d2B = new TH1D("hip2d2B","hip2d2B",40,-0.1,0.1);
  TH1D *hip2d2C = new TH1D("hip2d2C","hip2d2C",40,-0.1,0.1);
  TH1D *hip2d2L = new TH1D("hip2d2L","hip2d2L",40,-0.1,0.1);
  hip2d2->Sumw2(); hip2d2B->Sumw2(); hip2d2C->Sumw2(); hip2d2L->Sumw2(); 

  TH1D *hip2dSig2 = new TH1D("hip2dSig2","hip2dSig2",70,-35,35);
  TH1D *hip2dSig2B = new TH1D("hip2dSig2B","hip2dSig2B",70,-35,35);
  TH1D *hip2dSig2C = new TH1D("hip2dSig2C","hip2dSig2C",70,-35,35);
  TH1D *hip2dSig2L = new TH1D("hip2dSig2L","hip2dSig2L",70,-35,35);
  hip2dSig2->Sumw2(); hip2dSig2B->Sumw2(); hip2dSig2C->Sumw2(); hip2dSig2L->Sumw2(); 

  TH1D *hip2d3 = new TH1D("hip2d3","hip2d3",40,-0.1,0.1);
  TH1D *hip2d3B = new TH1D("hip2d3B","hip2d3B",40,-0.1,0.1);
  TH1D *hip2d3C = new TH1D("hip2d3C","hip2d3C",40,-0.1,0.1);
  TH1D *hip2d3L = new TH1D("hip2d3L","hip2d3L",40,-0.1,0.1);
  hip2d3->Sumw2(); hip2d3B->Sumw2(); hip2d3C->Sumw2(); hip2d3L->Sumw2(); 

  TH1D *hip2dSig3 = new TH1D("hip2dSig3","hip2dSig3",70,-35,35);
  TH1D *hip2dSig3B = new TH1D("hip2dSig3B","hip2dSig3B",70,-35,35);
  TH1D *hip2dSig3C = new TH1D("hip2dSig3C","hip2dSig3C",70,-35,35);
  TH1D *hip2dSig3L = new TH1D("hip2dSig3L","hip2dSig3L",70,-35,35);
  hip2dSig3->Sumw2(); hip2dSig3B->Sumw2(); hip2dSig3C->Sumw2(); hip2dSig3L->Sumw2(); 
  
  TH1D *hip3d = new TH1D("hip3d","hip3d",40,-0.1,0.1);
  TH1D *hip3dB = new TH1D("hip3dB","hip3dB",40,-0.1,0.1);
  TH1D *hip3dC = new TH1D("hip3dC","hip3dC",40,-0.1,0.1);
  TH1D *hip3dL = new TH1D("hip3dL","hip3dL",40,-0.1,0.1);
  hip3d->Sumw2(); hip3dB->Sumw2(); hip3dC->Sumw2(); hip3dL->Sumw2(); 
  
  TH1D *hip3dSig = new TH1D("hip3dSig","hip3dSig",70,-35,35);
  TH1D *hip3dSigB = new TH1D("hip3dSigB","hip3dSigB",70,-35,35);
  TH1D *hip3dSigC = new TH1D("hip3dSigC","hip3dSigC",70,-35,35);
  TH1D *hip3dSigL = new TH1D("hip3dSigL","hip3dSigL",70,-35,35);
  hip3dSig->Sumw2(); hip3dSigB->Sumw2(); hip3dSigC->Sumw2(); hip3dSigL->Sumw2(); 
  
  TH1D *hip3d1 = new TH1D("hip3d1","hip3d1",40,-0.1,0.1);
  TH1D *hip3d1B = new TH1D("hip3d1B","hip3d1B",40,-0.1,0.1);
  TH1D *hip3d1C = new TH1D("hip3d1C","hip3d1C",40,-0.1,0.1);
  TH1D *hip3d1L = new TH1D("hip3d1L","hip3d1L",40,-0.1,0.1);
  hip3d1->Sumw2(); hip3d1B->Sumw2(); hip3d1C->Sumw2(); hip3d1L->Sumw2(); 

  TH1D *hip3dSig1 = new TH1D("hip3dSig1","hip3dSig1",70,-35,35);
  TH1D *hip3dSig1B = new TH1D("hip3dSig1B","hip3dSig1B",70,-35,35);
  TH1D *hip3dSig1C = new TH1D("hip3dSig1C","hip3dSig1C",70,-35,35);
  TH1D *hip3dSig1L = new TH1D("hip3dSig1L","hip3dSig1L",70,-35,35);
  hip3dSig1->Sumw2(); hip3dSig1B->Sumw2(); hip3dSig1C->Sumw2(); hip3dSig1L->Sumw2(); 

  TH1D *hip3d2 = new TH1D("hip3d2","hip3d2",40,-0.1,0.1);
  TH1D *hip3d2B = new TH1D("hip3d2B","hip3d2B",40,-0.1,0.1);
  TH1D *hip3d2C = new TH1D("hip3d2C","hip3d2C",40,-0.1,0.1);
  TH1D *hip3d2L = new TH1D("hip3d2L","hip3d2L",40,-0.1,0.1);
  hip3d2->Sumw2(); hip3d2B->Sumw2(); hip3d2C->Sumw2(); hip3d2L->Sumw2(); 

  TH1D *hip3dSig2 = new TH1D("hip3dSig2","hip3dSig2",70,-35,35);
  TH1D *hip3dSig2B = new TH1D("hip3dSig2B","hip3dSig2B",70,-35,35);
  TH1D *hip3dSig2C = new TH1D("hip3dSig2C","hip3dSig2C",70,-35,35);
  TH1D *hip3dSig2L = new TH1D("hip3dSig2L","hip3dSig2L",70,-35,35);
  hip3dSig2->Sumw2(); hip3dSig2B->Sumw2(); hip3dSig2C->Sumw2(); hip3dSig2L->Sumw2(); 

  TH1D *hip3d3 = new TH1D("hip3d3","hip3d3",40,-0.1,0.1);
  TH1D *hip3d3B = new TH1D("hip3d3B","hip3d3B",40,-0.1,0.1);
  TH1D *hip3d3C = new TH1D("hip3d3C","hip3d3C",40,-0.1,0.1);
  TH1D *hip3d3L = new TH1D("hip3d3L","hip3d3L",40,-0.1,0.1);
  hip3d3->Sumw2(); hip3d3B->Sumw2(); hip3d3C->Sumw2(); hip3d3L->Sumw2(); 

  TH1D *hip3dSig3 = new TH1D("hip3dSig3","hip3dSig3",70,-35,35);
  TH1D *hip3dSig3B = new TH1D("hip3dSig3B","hip3dSig3B",70,-35,35);
  TH1D *hip3dSig3C = new TH1D("hip3dSig3C","hip3dSig3C",70,-35,35);
  TH1D *hip3dSig3L = new TH1D("hip3dSig3L","hip3dSig3L",70,-35,35);
  hip3dSig3->Sumw2(); hip3dSig3B->Sumw2(); hip3dSig3C->Sumw2(); hip3dSig3L->Sumw2(); 

  TH1D *hipDist2Jet = new TH1D("hipDist2Jet","hipDist2Jet",40,-0.1,0);
  TH1D *hipDist2JetB = new TH1D("hipDist2JetB","hipDist2JetB",40,-0.1,0);
  TH1D *hipDist2JetC = new TH1D("hipDist2JetC","hipDist2JetC",40,-0.1,0);
  TH1D *hipDist2JetL = new TH1D("hipDist2JetL","hipDist2JetL",40,-0.1,0);
  hipDist2Jet->Sumw2(); hipDist2JetB->Sumw2(); hipDist2JetC->Sumw2(); hipDist2JetL->Sumw2(); 
  
  TH1D *hipDist2JetSig = new TH1D("hipDist2JetSig","hipDist2JetSig",40,-0.1,0.1);
  TH1D *hipDist2JetSigB = new TH1D("hipDist2JetSigB","hipDist2JetSigB",40,-0.1,0.1);
  TH1D *hipDist2JetSigC = new TH1D("hipDist2JetSigC","hipDist2JetSigC",40,-0.1,0.1);
  TH1D *hipDist2JetSigL = new TH1D("hipDist2JetSigL","hipDist2JetSigL",40,-0.1,0.1);
  hipDist2JetSig->Sumw2(); hipDist2JetSigB->Sumw2(); hipDist2JetSigC->Sumw2(); hipDist2JetSigL->Sumw2(); 
  
  TH1D *hipClosest2Jet = new TH1D("hipClosest2Jet","hipClosest2Jet",40,0,1);
  TH1D *hipClosest2JetB = new TH1D("hipClosest2JetB","hipClosest2JetB",40,0,1);
  TH1D *hipClosest2JetC = new TH1D("hipClosest2JetC","hipClosest2JetC",40,0,1);
  TH1D *hipClosest2JetL = new TH1D("hipClosest2JetL","hipClosest2JetL",40,0,1);
  hipClosest2Jet->Sumw2(); hipClosest2JetB->Sumw2(); hipClosest2JetC->Sumw2(); hipClosest2JetL->Sumw2(); 

  Double_t t_jtpt, t_jteta, t_jtphi, t_rawpt, t_refpt, t_discr_prob, t_discr_ssvHighEff, t_discr_ssvHighPur, t_discr_csvSimple, t_svtxm;
  Double_t t_pthat, t_weight;
  Int_t t_refparton_flavorForB;
  Int_t trigIndex, t_bin;

  Int_t t_nIP;
  Double_t t_ipPt[100], t_ipProb0[100];
  Int_t t_ipJetIndex[100];

  TTree *nt = new TTree("nt","");
  nt->Branch("jtpt",&t_jtpt,"jtpt/D");
  nt->Branch("jteta",&t_jteta,"jteta/D");
  nt->Branch("jtphi",&t_jtphi,"jtphi/D");
  nt->Branch("rawpt",&t_rawpt,"rawpt/D");
  nt->Branch("refpt",&t_refpt,"refpt/D");
  nt->Branch("refparton_flavorForB",&t_refparton_flavorForB,"refparton_flavorForB/I");
  nt->Branch("discr_prob",&t_discr_prob,"discr_prob/D");
  nt->Branch("discr_ssvHighEff",&t_discr_ssvHighEff,"discr_ssvHighEff/D");
  nt->Branch("discr_ssvHighPur",&t_discr_ssvHighPur,"discr_ssvHighPur/D");
  nt->Branch("discr_csvSimple",&t_discr_csvSimple,"discr_csvSimple/D");
  nt->Branch("svtxm",&t_svtxm,"svtxm/D");
  nt->Branch("bin",&t_bin,"bin/I");
  if(ppPbPb){
    nt->Branch("trigIndex",&trigIndex,"trigIndex/I");
  }
  if(ExpandedTree){
    nt->Branch("nIP",&t_nIP);
    nt->Branch("ipPt",t_ipPt,"ipPt[nIP]/D");
    nt->Branch("ipProb0",t_ipProb0,"ipProb0[nIP]/D");
    nt->Branch("ipJetIndex",t_ipJetIndex,"ipJetIndex[nIP]/I");
  }
  if(!ppPbPb){
    nt->Branch("HLT_Jet20_noJetID_v1",&HLT_PAJet20_NoJetID_v1,"HLT_Jet20_noJetID_v1/I");
    nt->Branch("HLT_Jet40_noJetID_v1",&HLT_PAJet40_NoJetID_v1,"HLT_Jet40_noJetID_v1/I");
    nt->Branch("HLT_Jet60_noJetID_v1",&HLT_PAJet60_NoJetID_v1,"HLT_Jet60_noJetID_v1/I");
    nt->Branch("HLT_Jet80_noJetID_v1",&HLT_PAJet80_NoJetID_v1,"HLT_Jet80_noJetID_v1/I");
    nt->Branch("HLT_Jet100_noJetID_v1",&HLT_PAJet100_NoJetID_v1,"HLT_Jet100_noJetID_v1/I");
    nt->Branch("pVertexFilterCutGplusUpsPP",&pVertexFilterCutGplusUpsPP,"pVertexFilterCutGplusUpsPP/I");
  }

  if(isMC) nt->Branch("pthat",&t_pthat,"pthat/D");
  nt->Branch("weight",&t_weight,"weight/D");

  TNtuple *ntMuReq;
  if(isMC) ntMuReq = new TNtuple("ntMuReq","","jtpt:jteta:rawpt:refpt:refparton_flavorForB:weight:discr_prob:discr_ssvHighEff:discr_ssvHighPur:discr_csvSimple:svtxm:muptrel");
  else ntMuReq = new TNtuple("ntMuReq","","jtpt:jteta:rawpt:refparton_flavorForB:weight:discr_prob:discr_ssvHighEff:discr_ssvHighPur:discr_csvSimple:svtxm:muptrel");

  // grab the JEC's
   
  //JetCorrectorParameters* parHI442x_l2, * parHI442x_l3;
  // vector<JetCorrectorParameters> vpar_HI442x;   
  // FactorizedJetCorrector *_JEC_HI442x=NULL;
   
  if(updateJEC){   
     
    //cout<<" updating the JECs, USING REGPF "<<endl;

    //string L2Name = "JEC/JEC_regPF_L2Relative_AK3PF.txt";
    //string L3Name = "JEC/JEC_regPF_L3Absolute_AK3PF.txt";
    string L2Name = "JEC/JEC_dijet_L2Relative_AK3PF.txt";
    string L3Name = "JEC/JEC_dijet_L3Absolute_AK3PF.txt";
     
    // parHI442x_l2 = new JetCorrectorParameters(L2Name.c_str());
    //parHI442x_l3 = new JetCorrectorParameters(L3Name.c_str());
          
    // vpar_HI442x.push_back(*parHI442x_l2);
    // vpar_HI442x.push_back(*parHI442x_l3);
    // _JEC_HI442x = new FactorizedJetCorrector(vpar_HI442x);
  }     
   
  std::ifstream instr(infile.c_str(), std::ifstream::in);
  std::ifstream HFstr(HFfile.c_str(), std::ifstream::in);
  std::string filename;
  int nFiles=0;
  if(ppPbPb) nFiles=1;
  else if(isMC){
    nFiles=QCDpthatBins;
    if(isMC>1) nFiles+=HFpthatBins;
  }
  else{
    nFiles=dataFiles;
  }
  for(int ifile=0; ifile<nFiles; ifile++){
    
    //Add b/c statistics to the HF statistics
    if(!ppPbPb){
      if((isMC && ifile<QCDpthatBins) || !isMC){
	instr >> filename;
      }
      else if(isMC && ifile>=QCDpthatBins){
	HFstr >> filename;
      }
      std::cout << "File: " << filename << std::endl;
      fin = TFile::Open(filename.c_str());
    }
    TTree *t;
    if(usePUsub) t = (TTree*) fin->Get("akPu3PFJetAnalyzer/t");
    else t = (TTree*) fin->Get("ak3PFJetAnalyzer/t");
    TTree *tSkim = (TTree*) fin->Get("skimanalysis/HltTree");
    TTree *tEvt = NULL;
    if(!ppPbPb) tEvt = (TTree*) fin->Get("hiEvtAnalyzer/HiTree");
    TTree *tHlt = NULL;
    if(!ppPbPb) tHlt = (TTree*) fin->Get("hltanalysis/HltTree");
    TTree *tmu = (TTree*) fin->Get("muonTree/HLTMuTree");
    if(!t || !tSkim || (!tEvt&&!ppPbPb) || (!tHlt&&!ppPbPb)){ cout << "Error! Can't find one of the trees!" << endl; exit(0);}
     
    if(tEvt) t->AddFriend("hiEvtAnalyzer/HiTree");    
    if(tHlt) t->AddFriend("hltanalysis/HltTree");
    if(tSkim) t->AddFriend("skimanalysis/HltTree");
     
    t->SetBranchAddress("evt",&evt);
    t->SetBranchAddress("lumi",&lumi);
    if(cbin != -1 || ppPbPb) t->SetBranchAddress("bin",&bin);
    if(!isMC) t->SetBranchAddress("run",&run);
    if(ppPbPb) t->SetBranchAddress("hf",&hf);
    t->SetBranchAddress("vz",&vz);           
    t->SetBranchAddress("nref",&nref);
    t->SetBranchAddress("rawpt",rawpt);
    t->SetBranchAddress("jtpt",jtpt);
    t->SetBranchAddress("jteta",jteta);
    t->SetBranchAddress("jtphi",jtphi);
    if(!ppPbPb){
      t->SetBranchAddress("jtpu",jtpu);
      t->SetBranchAddress("jty",jty);
    }
    t->SetBranchAddress("discr_ssvHighEff",discr_ssvHighEff);
    t->SetBranchAddress("discr_ssvHighPur",discr_ssvHighPur);
    //t->SetBranchAddress("discr_csvMva",discr_csvMva);
    t->SetBranchAddress("discr_csvSimple",discr_csvSimple);
    //t->SetBranchAddress("discr_muByIp3",discr_muByIp3);
    t->SetBranchAddress("discr_muByPt",discr_muByPt);
    t->SetBranchAddress("discr_prob",discr_prob);

    t->SetBranchAddress("discr_probb",discr_probb);
    t->SetBranchAddress("discr_tcHighEff",discr_tcHighEff);
    t->SetBranchAddress("discr_tcHighPur",discr_tcHighPur);
    t->SetBranchAddress("nsvtx",nsvtx);
    t->SetBranchAddress("svtxntrk",svtxntrk);
    t->SetBranchAddress("svtxdl",svtxdl);
    t->SetBranchAddress("svtxdls",svtxdls);
    t->SetBranchAddress("svtxpt",svtxpt);
    
    t->SetBranchAddress("svtxm",svtxm);

    t->SetBranchAddress("nIPtrk",nIPtrk);
    t->SetBranchAddress("nselIPtrk",nselIPtrk); 
    t->SetBranchAddress("nIP",&nIP);

    if(doTracks){
      t->SetBranchAddress("ipJetIndex",ipJetIndex);
      t->SetBranchAddress("ipPt",ipPt);
      t->SetBranchAddress("ipProb0",ipProb0);
      //t->SetBranchAddress("ipProb1",ipProb1);
      t->SetBranchAddress("ip2d",ip2d);
      t->SetBranchAddress("ip2dSig",ip2dSig);
      t->SetBranchAddress("ip3d",ip3d);
      t->SetBranchAddress("ip3dSig",ip3dSig);
      t->SetBranchAddress("ipDist2Jet",ipDist2Jet);
      //t->SetBranchAddress("ipDist2JetSig",ipDist2JetSig);
      t->SetBranchAddress("ipClosest2Jet",ipClosest2Jet);
    }

    t->SetBranchAddress("mupt",mupt);
    if(ppPbPb) t->SetBranchAddress("muptPF",muptPF);
    if(!ppPbPb) t->SetBranchAddress("pVertexFilterCutGplusUpsPP",&pVertexFilterCutGplusUpsPP);

    /*
      t->SetBranchAddress("mue",mue);
      t->SetBranchAddress("mueta",mueta);
      t->SetBranchAddress("muphi",muphi);
      t->SetBranchAddress("mudr",mudr);
      t->SetBranchAddress("muptrel",muptrel);
      t->SetBranchAddress("muchg",muchg);
    */
    if(isMC){
      t->SetBranchAddress("pthat",&pthat);
      t->SetBranchAddress("beamId1",&beamId1);
      t->SetBranchAddress("beamId2",&beamId2);
      t->SetBranchAddress("refpt",refpt);
      t->SetBranchAddress("refeta",refeta);
      t->SetBranchAddress("refy",refy);
      t->SetBranchAddress("refphi",refphi);
      t->SetBranchAddress("refdphijt",refdphijt);
      t->SetBranchAddress("refdrjt",refdrjt);
      t->SetBranchAddress("refparton_pt",refparton_pt);
      t->SetBranchAddress("refparton_flavor",refparton_flavor);
      t->SetBranchAddress("refparton_flavorForB",refparton_flavorForB);

      TBranch* tweight;
      if(isMC){
	tweight = t->GetBranch("weight");
	if(!tweight){
	  if(ifile==0){
	    cout << "Weight not found in Tree. Calculating..." << endl;
	    useWeight=0;
	  }
	}
	if(!ppPbPb && !useWeight && ifile==0){
	  MCentr = countMCevents(infile, HFfile, usePUsub, isMC);
	  // if(isMC>1){
	  //  for(int lm=HFpthatBins+2; lm<QCDpthatBins+1; lm++){
	  //    MCentr[HFpthatBins] += MCentr[lm]; //hack because we go to pthat bin 540 in QCD jet and only pthat bin 170 in b/c jet MC
	  //  }
	  // }
	  for(int i=0; i<10; i++){
	    cout << "MCentr["<<i<<"]: " << *(MCentr+i) << endl;
	  }
	}
      }
    }
    t->SetBranchAddress("chargedMax",chargedMax);
    t->SetBranchAddress("photonMax",photonMax);
    t->SetBranchAddress("neutralMax",neutralMax);
    t->SetBranchAddress("chargedSum",chargedSum);
    t->SetBranchAddress("photonSum",photonSum);
    t->SetBranchAddress("neutralSum",neutralSum);
    t->SetBranchAddress("muSum",muSum);
    t->SetBranchAddress("eSum",eSum);
    
    if(isMC&&useWeight){
      t->SetBranchAddress("weight",&weight);
      t->SetBranchAddress("xSecWeight",&xSecWeight);
      if(ppPbPb)t->SetBranchAddress("centWeight",&centWeight);
      t->SetBranchAddress("vzWeight",&vzWeight);
    }

    if(ppPbPb){
      t->SetBranchAddress("nHLTBit",&nHLTBit);
      t->SetBranchAddress("hltBit",hltBit);
      
      tSkim->SetBranchAddress("pvSel",&pvSel);
      tSkim->SetBranchAddress("hbheNoiseSel",&hbheNoiseSel);
      tSkim->SetBranchAddress("spikeSel",&spikeSel);
      tSkim->SetBranchAddress("collSell",&collSell);
    }
    if(!ppPbPb){
      t->SetBranchAddress("HLT_PAJet20_NoJetID_v1",&HLT_PAJet20_NoJetID_v1);
      t->SetBranchAddress("HLT_PAJet40_NoJetID_v1",&HLT_PAJet40_NoJetID_v1);
      t->SetBranchAddress("HLT_PAJet60_NoJetID_v1",&HLT_PAJet60_NoJetID_v1);
      t->SetBranchAddress("HLT_PAJet80_NoJetID_v1",&HLT_PAJet80_NoJetID_v1);
      t->SetBranchAddress("HLT_PAJet100_NoJetID_v1",&HLT_PAJet100_NoJetID_v1);
      t->SetBranchAddress("pPAcollisionEventSelectionPA",&pPAcollisionEventSelectionPA);
      t->SetBranchAddress("pHBHENoiseFilter",&pHBHENoiseFilter);
      t->SetBranchAddress("pprimaryvertexFilter",&pprimaryvertexFilter);
    }

    Long64_t nentries = t->GetEntries();

    int gspCounter=0;
    //nentries=10;
    for (Long64_t i=0; i<nentries;i++) {
       
      if (i%100000==0) cout<<" i = "<<i<<" out of "<<nentries<<" ("<<(int)(100*(float)i/(float)nentries)<<"%)"<<endl; 
      
      tSkim->GetEntry(i);
      t->GetEntry(i);
      if(ppPbPb && isMC){
	// temporarily remove cuts from MC
	if(!pvSel||!spikeSel) continue; //hbheNoise doesn't work in mixed events
      }
      else if(ppPbPb){
	//if(!pvSel||!hbheNoiseSel||!spikeSel) continue;
	// turn off spike and on coll Sel
	if(!pvSel||!hbheNoiseSel||!collSell){
	  //cout<<" selection failed, pvSel="<<pvSel<<", hbheNoiseSel="<<hbheNoiseSel<<" , collSell="<<collSell<<endl;
	  continue;
	}
      }

      //Cut to remove events that correspond to the twiki "good events" but not the golden lumi filter
      if(!isMC){
	if(((int)run==211821 && lumi>=57 && lumi<370) || ((int)run==211821 && lumi>420)) continue;
      }

      if(!ppPbPb){
        if(!isMC){
          if(!pHBHENoiseFilter || !pprimaryvertexFilter || !pPAcollisionEventSelectionPA) continue;
        }
        else{
          if(!pHBHENoiseFilter || !pPAcollisionEventSelectionPA) continue;
        }
      }

      if(!HLT_PAJet40_NoJetID_v1 && !HLT_PAJet60_NoJetID_v1 && !HLT_PAJet80_NoJetID_v1) continue;
      
      if(ppPbPb){
	if(cbin==-1){
	  // do nothing
	}
	else if(cbin==0){
	  if(bin>=8) continue;
	}
	else if(cbin==1) {
	  if(bin<8||bin>=20) continue;
	}
	else if(cbin==2){
	  if(bin<20) continue;
	}
	else {
	  cout<<" bin not defined "<<endl;
	  return;
	}
      }
      if(ppPbPb) t_bin=bin;
      else  t_bin=39;
      
      if(isMC&&!ppPbPb){
	if(beamId1==2112 || beamId2==2112)  continue;
      }
      
      if(ppPbPb){
	if(jetTrig==1&&!hltBit[10]) continue;
	if(jetTrig==2&&!hltBit[9]) continue;
      }

      if(fabs(vz)>15.) continue;
      
      // pileup rejection
      if(ppPbPb && hf>150000.){
	cout<<" rejecting pileup, "<<" hf "<<hf<<" bin "<<bin<<endl;
	for(int ij=0;ij<nref;ij++) if(jtpt[ij]>65.&&fabs(jteta[ij])<2.)cout<<" # associated tracks =  "<<nselIPtrk[ij]<<endl;
	continue;
      }
      bool isNoise=false;
      if(ppPbPb){
	for(int ij=0; ij<nref; ij++){	  
	  if(jtpt[ij]>4000&&fabs(jteta[ij])<2) cout<<" mupt "<<mupt[0]<<" muptPF "<<muptPF[0]<<endl;
	  if(jtpt[ij]>minJetPt&&fabs(jteta[ij])<2){
	    if(neutralMax[ij]/(neutralMax[ij]+chargedMax[ij]+photonMax[ij])>0.975){
	      cout<<" cleaning event with jet of  "<<jtpt[ij]<<", eta "<<jteta[ij]<<" noise = "<<neutralMax[ij]/(neutralMax[ij]+chargedMax[ij]+photonMax[ij])<<endl;
	      isNoise=true;
	    }
	    if(muptPF[ij]>10&&mupt[ij]/muptPF[ij]<0.75){
	      cout<<" cleaning event with jet of  "<<jtpt[ij]<<", eta "<<jteta[ij]<<" muptPF = "<<muptPF[ij]<<" mupt "<<mupt[ij]<<endl;
	      isNoise=true;
	    }
	    if(chargedSum[ij]+photonSum[ij]+neutralSum[ij]+muSum[ij]+eSum[ij]<0.5*rawpt[ij]){
	      cout<<" cleaning event with jet of  "<<jtpt[ij]<<", eta "<<jteta[ij]<<" sum PF pt = "<<chargedSum[ij]+photonSum[ij]+neutralSum[ij]+muSum[ij]+eSum[ij]<<endl;
	      isNoise=true;
	    }
	  }
	}
      }
      if(isNoise) continue;
      
      
      if(!isMC&&ppPbPb&&jetTrig<2){
	tmu->GetEntry(i);
	
	bool foundEvt = false;
	for(int irun=0;irun<6;irun++){       
	  if(run==dupRuns[irun]) {
	    // binary search does not give the right behavior for some reason
	    //if(binary_search(usedEvents[irun].begin(), usedEvents[irun].end(), evt)) {
	    // use the slower find instead
	    std::vector<int>::iterator it;
	    // iterator to vector element:
	    it = find (usedEvents[irun].begin(), usedEvents[irun].end(), evt);
	    if(it!=usedEvents[irun].end()){
	      
	      nDup++;
	      foundEvt = true;
	      //cout<< " duplicate event, run: "<<run<<" evt: "<<evt<<endl;
	      break;
	    }
	    usedEvents[irun].push_back(evt);
	  }
	}
	
	if(foundEvt) continue;
      }
      
      if(updateJEC){
	
	//for(int ij=0; ij<nref; ij++){	  
	//  _JEC_HI442x->setJetEta(jteta[ij]);
	//  _JEC_HI442x->setJetPt(rawpt[ij]);
	//  jtpt[ij] = rawpt[ij]*_JEC_HI442x->getCorrection(); 
	//}	
      }
      
      if(useWeight){
	if(isMC)w=weight;
	else if(ppPbPb){
	  if(jetTrig==2){
	    if(hltBit[10]) w=0.;
	    else w=1./8.93333857823361388e-01;
	  }
	}
      }
      //trigger weighting in pp data
      if(!ppPbPb && !isMC){
	bool trgDec[4] = {(bool)HLT_PAJet20_NoJetID_v1, (bool)HLT_PAJet40_NoJetID_v1, (bool)HLT_PAJet60_NoJetID_v1, (bool)HLT_PAJet80_NoJetID_v1};
	w = trigComb(trgDec, pscls);
      }

      if(ppPbPb){
	if(hltBit[10]) trigIndex=3;
	else if(hltBit[9]) trigIndex=2;
	else if(hltBit[8]) trigIndex=1;
	else trigIndex=0;
      }
      
      //Do the weighting = x-sec / Nentries, where Nentries is weighted differently for B/C jets and QCD jets
      if(isMC){
	t_pthat=pthat;
	int j=0;
	while(pthat>pthatbin[j] && j<QCDpthatBins) j++;
	//	if(isMC>1 && ifile>=QCDpthatBins){
	  // cout << "pthat: "<< pthat << endl;
	  //cout << "bin: "<< j << endl;
	  /*int k = (j<HFpthatBins+1 ? j : HFpthatBins+1); //WATCH THIS! IF YOU ADD AN HF pthat-15 sample, this must be FIXED!!
	    w = (wght[k-1]/MCentr[j]);*/
	  //cout << "weight: "<< wght[k-1] << endl;
	  //cout << "MCentr: "<< MCentr[j] << endl;
	  // w *= HFweight[k-1]; //do HF reweighting for b/c samples (changed in MCcounter - added QCD and B/C samples together)
	  //if(pthat>220) continue;
	//	}
	//	else{
	  w = (wght[j-1]/MCentr[j]); //wght[0] = pthat>15, MCentr[0] = pthat<15.  I know it's dumb - bear with me.
	  //	}
      }
      t_weight=w;	  
      
      int useEvent=0;
      
      int trackPosition =0;

      for(int ij=0;ij<nref;ij++){
      
	trackPosition+=nselIPtrk[ij];
      
	if(useGSP==2){
	  if(refparton_isGSP[ij]==1){
	    gspCounter++;
	    if(gspCounter%2==0) continue;
	  }	
	}
	if(useGSP==3){
	  if(refparton_isGSP[ij]==0){
	    gspCounter++;
	    if(gspCounter%2==0) continue;
	  }	
	}
      
	if(jtpt[ij]>minJetPt && fabs(jteta[ij])<maxJetEta){ 
	  if(doNtuples){
	  
	    t_jtpt=jtpt[ij];
	    t_jteta=jteta[ij];
	    t_jtphi=jtphi[ij];
	    t_rawpt=rawpt[ij];
	    t_refpt=refpt[ij];
	    t_refparton_flavorForB=refparton_flavorForB[ij];
	    t_discr_prob=discr_prob[ij];
	    t_discr_ssvHighEff=discr_ssvHighEff[ij];
	    t_discr_ssvHighPur=discr_ssvHighPur[ij];
	    t_discr_csvSimple=discr_csvSimple[ij];
	    t_svtxm=svtxm[ij];
	    
	    //Find jet tracks that correspond to the jet & apply proximity cuts
	    if(ExpandedTree){
	      t_nIP=nselIPtrk[ij];
	      int counter=0;
	      for(int itrk=0; itrk<nIP; itrk++){
		if(ipJetIndex[itrk] == ij){
		  t_ipProb0[counter] = ipProb0[itrk];
		  t_ipPt[counter] = ipPt[itrk];
		  t_ipJetIndex[counter] = ij;
		  counter++;
		}
	      }
	    }
	    nt->Fill();

	    if (sqrt(acos(cos(jtphi[ij]-muphi[ij]))*acos(cos(jtphi[ij]-muphi[ij]))+(jteta[ij]-mueta[ij])*(jteta[ij]-mueta[ij]))<0.5 && mupt[ij]>minMuPt) { 
	    
	      if(isMC)ntMuReq->Fill(jtpt[ij],jteta[ij],rawpt[ij],refpt[ij],refparton_flavorForB[ij],w,discr_prob[ij],discr_ssvHighEff[ij],discr_ssvHighPur[ij],discr_csvSimple[ij],svtxm[ij],muptrel[ij]); 
	      else ntMuReq->Fill(jtpt[ij],jteta[ij],rawpt[ij],refparton_flavorForB[ij],w,discr_prob[ij],discr_ssvHighEff[ij],discr_ssvHighPur[ij],discr_csvSimple[ij],svtxm[ij],muptrel[ij]); 
	    }
	  }

	
	  if(!doJets) continue;
	
	  if (isMuTrig) {
	    //muon requirement
	    if (sqrt(acos(cos(jtphi[ij]-muphi[ij]))*acos(cos(jtphi[ij]-muphi[ij]))+(jteta[ij]-mueta[ij])*(jteta[ij]-mueta[ij]))>0.5 || mupt[ij]<minMuPt) continue;
	  }
	
	  useEvent=1;
	
	  hjtpt->Fill(jtpt[ij],w);    
	  if(isMC){
	    if(abs(refparton_flavorForB[ij])==5)hjtptB->Fill(jtpt[ij],w);    
	    else if(abs(refparton_flavorForB[ij])==4)hjtptC->Fill(jtpt[ij],w);    
	    else if(abs(refparton_flavorForB[ij])<99)hjtptL->Fill(jtpt[ij],w);    
	    else hjtptU->Fill(jtpt[ij],w);    
	  }
	  hrawpt->Fill(rawpt[ij],w);    
	  if(isMC){
	    if(abs(refparton_flavorForB[ij])==5)hrawptB->Fill(rawpt[ij],w);    
	    else if(abs(refparton_flavorForB[ij])==4)hrawptC->Fill(rawpt[ij],w);    
	    else if(abs(refparton_flavorForB[ij])<99)hrawptL->Fill(rawpt[ij],w);    
	  }
	  hjteta->Fill(jteta[ij],w);    
	  if(isMC){
	    if(abs(refparton_flavorForB[ij])==5)hjtetaB->Fill(jteta[ij],w);    
	    else if(abs(refparton_flavorForB[ij])==4)hjtetaC->Fill(jteta[ij],w);    
	    else if(abs(refparton_flavorForB[ij])<99)hjtetaL->Fill(jteta[ij],w);    
	  }
	  hjtphi->Fill(jtphi[ij],w);    
	  if(isMC){
	    if(abs(refparton_flavorForB[ij])==5)hjtphiB->Fill(jtphi[ij],w);    
	    else if(abs(refparton_flavorForB[ij])==4)hjtphiC->Fill(jtphi[ij],w);    
	    else if(abs(refparton_flavorForB[ij])<99)hjtphiL->Fill(jtphi[ij],w);    
	  }
	  //*
	  hdiscr_csvSimple->Fill(discr_csvSimple[ij],w);    
	  if(isMC){
	    if(abs(refparton_flavorForB[ij])==5)hdiscr_csvSimpleB->Fill(discr_csvSimple[ij],w);    
	    else if(abs(refparton_flavorForB[ij])==4)hdiscr_csvSimpleC->Fill(discr_csvSimple[ij],w);    
	    else if(abs(refparton_flavorForB[ij])<99)hdiscr_csvSimpleL->Fill(discr_csvSimple[ij],w);    
	  }
	
	  hdiscr_prob->Fill(discr_prob[ij],w);    
	  if(isMC){
	    if(abs(refparton_flavorForB[ij])==5)hdiscr_probB->Fill(discr_prob[ij],w); 
	    else if(abs(refparton_flavorForB[ij])==4)hdiscr_probC->Fill(discr_prob[ij],w);    
	    else if(abs(refparton_flavorForB[ij])<99)hdiscr_probL->Fill(discr_prob[ij],w);    
	  }

	  hdiscr_ssvHighEff->Fill(discr_ssvHighEff[ij],w);    
	  if(isMC){
	    if(abs(refparton_flavorForB[ij])==5)hdiscr_ssvHighEffB->Fill(discr_ssvHighEff[ij],w); 
	    else if(abs(refparton_flavorForB[ij])==4)hdiscr_ssvHighEffC->Fill(discr_ssvHighEff[ij],w);    
	    else if(abs(refparton_flavorForB[ij])<99)hdiscr_ssvHighEffL->Fill(discr_ssvHighEff[ij],w);    
	  }
	
	  hdiscr_ssvHighPur->Fill(discr_ssvHighPur[ij],w);    
	  if(isMC){
	    if(abs(refparton_flavorForB[ij])==5)hdiscr_ssvHighPurB->Fill(discr_ssvHighPur[ij],w); 
	    else if(abs(refparton_flavorForB[ij])==4)hdiscr_ssvHighPurC->Fill(discr_ssvHighPur[ij],w);    
	    else if(abs(refparton_flavorForB[ij])<99)hdiscr_ssvHighPurL->Fill(discr_ssvHighPur[ij],w);    
	  }

	  hdiscr_tcHighEff->Fill(discr_tcHighEff[ij],w);    
	  if(isMC){
	    if(abs(refparton_flavorForB[ij])==5)hdiscr_tcHighEffB->Fill(discr_tcHighEff[ij],w); 
	    else if(abs(refparton_flavorForB[ij])==4)hdiscr_tcHighEffC->Fill(discr_tcHighEff[ij],w);    
	    else if(abs(refparton_flavorForB[ij])<99)hdiscr_tcHighEffL->Fill(discr_tcHighEff[ij],w);    
	  }
	
	  hdiscr_tcHighPur->Fill(discr_tcHighPur[ij],w);    
	  if(isMC){
	    if(abs(refparton_flavorForB[ij])==5)hdiscr_tcHighPurB->Fill(discr_tcHighPur[ij],w); 
	    else if(abs(refparton_flavorForB[ij])==4)hdiscr_tcHighPurC->Fill(discr_tcHighPur[ij],w);    
	    else if(abs(refparton_flavorForB[ij])<99)hdiscr_tcHighPurL->Fill(discr_tcHighPur[ij],w);    
	  }
	  //*
	  hnsvtx->Fill(nsvtx[ij],w);    
	  if(isMC){
	    if(abs(refparton_flavorForB[ij])==5)hnsvtxB->Fill(nsvtx[ij],w);    
	    else if(abs(refparton_flavorForB[ij])==4)hnsvtxC->Fill(nsvtx[ij],w);    
	    else if(abs(refparton_flavorForB[ij])<99)hnsvtxL->Fill(nsvtx[ij],w); 
	  }
	
	  if(nsvtx[ij]>0){
	  
	    hsvtxntrk->Fill(svtxntrk[ij],w);    
	    if(isMC){
	      if(abs(refparton_flavorForB[ij])==5)hsvtxntrkB->Fill(svtxntrk[ij],w);    
	      else if(abs(refparton_flavorForB[ij])==4)hsvtxntrkC->Fill(svtxntrk[ij],w);
	      else if(abs(refparton_flavorForB[ij])<99)hsvtxntrkL->Fill(svtxntrk[ij],w);
	    }
	  
	    // require at least 1 tracks as in btagging @ 7 TeV note
	    if(svtxntrk[ij]>1){	  
	    
	      hsvtxdl->Fill(svtxdl[ij],w);    
	      if(isMC){
		if(abs(refparton_flavorForB[ij])==5)hsvtxdlB->Fill(svtxdl[ij],w);    
		else if(abs(refparton_flavorForB[ij])==4)hsvtxdlC->Fill(svtxdl[ij],w);
		else if(abs(refparton_flavorForB[ij])<99)hsvtxdlL->Fill(svtxdl[ij],w);
	      }
	    
	      hsvtxdls->Fill(svtxdls[ij],w);    
	      if(isMC){
		if(abs(refparton_flavorForB[ij])==5)hsvtxdlsB->Fill(svtxdls[ij],w);    
		else if(abs(refparton_flavorForB[ij])==4)hsvtxdlsC->Fill(svtxdls[ij],w); 
		else if(abs(refparton_flavorForB[ij])<99)hsvtxdlsL->Fill(svtxdls[ij],w); 
	      }
	    
	      hsvtxm->Fill(svtxm[ij],w);    
	      if(isMC){
		if(abs(refparton_flavorForB[ij])==5)hsvtxmB->Fill(svtxm[ij],w);    
		else if(abs(refparton_flavorForB[ij])==4)hsvtxmC->Fill(svtxm[ij],w);
		else if(abs(refparton_flavorForB[ij])<99)hsvtxmL->Fill(svtxm[ij],w); 
	      }
	    
	      hsvtxpt->Fill(svtxpt[ij],w);    
	      if(isMC){
		if(abs(refparton_flavorForB[ij])==5)hsvtxptB->Fill(svtxpt[ij],w);    
		else if(abs(refparton_flavorForB[ij])==4)hsvtxptC->Fill(svtxpt[ij],w);
		else if(abs(refparton_flavorForB[ij])<99)hsvtxptL->Fill(svtxpt[ij],w);
	      }
	    
	      if(svtxntrk[ij]>=3) {
	      
		hsvtxmSV3->Fill(svtxm[ij],w);    
		if(isMC){
		  if(abs(refparton_flavorForB[ij])==5)hsvtxmSV3B->Fill(svtxm[ij],w);    
		  else if(abs(refparton_flavorForB[ij])==4)hsvtxmSV3C->Fill(svtxm[ij],w);
		  else if(abs(refparton_flavorForB[ij])<99)hsvtxmSV3L->Fill(svtxm[ij],w); 
		}
	      
		hsvtxptSV3->Fill(svtxpt[ij],w);    
		if(isMC){
		  if(abs(refparton_flavorForB[ij])==5)hsvtxptSV3B->Fill(svtxpt[ij],w);    
		  else if(abs(refparton_flavorForB[ij])==4)hsvtxptSV3C->Fill(svtxpt[ij],w);
		  else if(abs(refparton_flavorForB[ij])<99)hsvtxptSV3L->Fill(svtxpt[ij],w);
		}
	      }
	    }
	  }
	  
	  hnIPtrk->Fill(nIPtrk[ij],w);    
	  if(isMC){
	    if(abs(refparton_flavorForB[ij])==5)hnIPtrkB->Fill(nIPtrk[ij],w);    
	    else if(abs(refparton_flavorForB[ij])==4)hnIPtrkC->Fill(nIPtrk[ij],w);    
	    else if(abs(refparton_flavorForB[ij])<99)hnIPtrkL->Fill(nIPtrk[ij],w);    
	  }
	
	
	  hnselIPtrk->Fill(nselIPtrk[ij],w);    
	
	  if(isMC){
	    if(abs(refparton_flavorForB[ij])==5)hnselIPtrkB->Fill(nselIPtrk[ij],w);    
	    else if(abs(refparton_flavorForB[ij])==4)hnselIPtrkC->Fill(nselIPtrk[ij],w);
	    else if(abs(refparton_flavorForB[ij])<99)hnselIPtrkL->Fill(nselIPtrk[ij],w);
	  }
	
	  if (sqrt(acos(cos(jtphi[ij]-muphi[ij]))*acos(cos(jtphi[ij]-muphi[ij]))+(jteta[ij]-mueta[ij])*(jteta[ij]-mueta[ij]))<0.5 && mupt[ij]>minMuPt) { 
	  
	    hmuptrel->Fill(muptrel[ij],w);    
	    if(isMC){
	      if(abs(refparton_flavorForB[ij])==5)hmuptrelB->Fill(muptrel[ij],w);
	      else if(abs(refparton_flavorForB[ij])==4)hmuptrelC->Fill(muptrel[ij],w);
	      else if(abs(refparton_flavorForB[ij])<99)hmuptrelL->Fill(muptrel[ij],w);
	    }
	  
	    if(svtxntrk[ij]>=2) {  
	      hmuptrelSV2->Fill(muptrel[ij],w);    
	      if(isMC){
		if(abs(refparton_flavorForB[ij])==5)hmuptrelSV2B->Fill(muptrel[ij],w);
		else if(abs(refparton_flavorForB[ij])==4)hmuptrelSV2C->Fill(muptrel[ij],w);
		else if(abs(refparton_flavorForB[ij])<99)hmuptrelSV2L->Fill(muptrel[ij],w);
	      }
	    }
	  
	    if(svtxntrk[ij]>=3) {  
	      hmuptrelSV3->Fill(muptrel[ij],w);    
	      if(isMC){
		if(abs(refparton_flavorForB[ij])==5)hmuptrelSV3B->Fill(muptrel[ij],w);
		else if(abs(refparton_flavorForB[ij])==4)hmuptrelSV3C->Fill(muptrel[ij],w);
		else if(abs(refparton_flavorForB[ij])<99)hmuptrelSV3L->Fill(muptrel[ij],w);
	      }
	    }
	  
	  }
	  //*/
     
	  float ip2d1MostSig=-999.;
	  float ip3d1MostSig=-999.;
	  float ip2dSig1MostSig=-999.;
	  float ip3dSig1MostSig=-999.;

	  float ip2d2MostSig=-999.;
	  float ip3d2MostSig=-999.;
	  float ip2dSig2MostSig=-999.;
	  float ip3dSig2MostSig=-999.;

	  float ip2d3MostSig=-999.;
	  float ip3d3MostSig=-999.;
	  float ip2dSig3MostSig=-999.;
	  float ip3dSig3MostSig=-999.;

	  for(int it=trackPosition-nselIPtrk[ij];it<trackPosition;it++){
	  
	    if(fabs(ipDist2Jet[it])>0.07) continue;
	    if(ipClosest2Jet[it] > 5.0) continue;
	  
	    if(ip2dSig[it]>ip2dSig1MostSig){	    
	      ip2d3MostSig=ip2d2MostSig;
	      ip2dSig3MostSig=ip2dSig2MostSig;
	      ip2d2MostSig=ip2d1MostSig;
	      ip2dSig2MostSig=ip2dSig1MostSig;
	      ip2d1MostSig=ip2d[it];
	      ip2dSig1MostSig=ip2dSig[it];
	    }
	    else if(ip2dSig[it]>ip2dSig2MostSig){	    
	      ip2d3MostSig=ip2d2MostSig;
	      ip2dSig3MostSig=ip2dSig2MostSig;
	      ip2d2MostSig=ip2d[it];
	      ip2dSig2MostSig=ip2dSig[it];
	    }
	    else if(ip2dSig[it]>ip2dSig3MostSig){	    
	      ip2d3MostSig=ip2d[it];
	      ip2dSig3MostSig=ip2dSig[it];
	    }

	  
	    if(ip3dSig[it]>ip3dSig1MostSig){	    
	      ip3d3MostSig=ip3d2MostSig;
	      ip3dSig3MostSig=ip3dSig2MostSig;
	      ip3d2MostSig=ip3d1MostSig;
	      ip3dSig2MostSig=ip3dSig1MostSig;
	      ip3d1MostSig=ip3d[it];
	      ip3dSig1MostSig=ip3dSig[it];
	    }
	    else if(ip3dSig[it]>ip3dSig2MostSig){	    
	      ip3d3MostSig=ip3d2MostSig;
	      ip3dSig3MostSig=ip3dSig2MostSig;
	      ip3d2MostSig=ip3d[it];
	      ip3dSig2MostSig=ip3dSig[it];
	    }
	    else if(ip3dSig[it]>ip3dSig3MostSig){	    
	      ip3d3MostSig=ip3d[it];
	      ip3dSig3MostSig=ip3dSig[it];
	    }



	    hipPt->Fill(ipPt[it],w);    

	    if(isMC){
	      if(abs(refparton_flavorForB[ij])==5)hipPtB->Fill(ipPt[it],w);
	      else if(abs(refparton_flavorForB[ij])==4)hipPtC->Fill(ipPt[it],w); 
	      else if(abs(refparton_flavorForB[ij])<99)hipPtL->Fill(ipPt[it],w); 
	    }
	  
	    hipProb0->Fill(ipProb0[it],w);    
	    if(isMC){
	      if(abs(refparton_flavorForB[ij])==5)hipProb0B->Fill(ipProb0[it],w);
	      else if(abs(refparton_flavorForB[ij])==4)hipProb0C->Fill(ipProb0[it],w);
	      else if(abs(refparton_flavorForB[ij])<99)hipProb0L->Fill(ipProb0[it],w);
	    }
	  
	    hipProb1->Fill(ipProb1[it],w);    
	    if(isMC){
	      if(abs(refparton_flavorForB[ij])==5)hipProb1B->Fill(ipProb1[it],w);
	      else if(abs(refparton_flavorForB[ij])==4)hipProb1C->Fill(ipProb1[it],w);
	      else if(abs(refparton_flavorForB[ij])<99)hipProb1L->Fill(ipProb1[it],w);
	    }
	  
	    hip2d->Fill(ip2d[it],w);    
	    if(isMC){
	      if(abs(refparton_flavorForB[ij])==5)hip2dB->Fill(ip2d[it],w);
	      else if(abs(refparton_flavorForB[ij])==4)hip2dC->Fill(ip2d[it],w); 
	      else if(abs(refparton_flavorForB[ij])<99)hip2dL->Fill(ip2d[it],w); 
	    }
	  
	    hip2dSig->Fill(ip2dSig[it],w);    
	    if(isMC){
	      if(abs(refparton_flavorForB[ij])==5)hip2dSigB->Fill(ip2dSig[it],w);
	      else if(abs(refparton_flavorForB[ij])==4)hip2dSigC->Fill(ip2dSig[it],w);
	      else if(abs(refparton_flavorForB[ij])<99)hip2dSigL->Fill(ip2dSig[it],w);
	    }
	  
	    hip3d->Fill(ip3d[it],w);    
	    if(isMC){
	      if(abs(refparton_flavorForB[ij])==5)hip3dB->Fill(ip3d[it],w);
	      else if(abs(refparton_flavorForB[ij])==4)hip3dC->Fill(ip3d[it],w); 
	      else if(abs(refparton_flavorForB[ij])<99)hip3dL->Fill(ip3d[it],w); 
	    }
	  
	    hip3dSig->Fill(ip3dSig[it],w);    
	    if(isMC){
	      if(abs(refparton_flavorForB[ij])==5)hip3dSigB->Fill(ip3dSig[it],w);
	      else if(abs(refparton_flavorForB[ij])==4)hip3dSigC->Fill(ip3dSig[it],w);
	      else if(abs(refparton_flavorForB[ij])<99)hip3dSigL->Fill(ip3dSig[it],w);
	    }
	  
	    hipDist2Jet->Fill(ipDist2Jet[it],w);    
	    if(isMC){
	      if(abs(refparton_flavorForB[ij])==5)hipDist2JetB->Fill(ipDist2Jet[it],w);
	      else if(abs(refparton_flavorForB[ij])==4)hipDist2JetC->Fill(ipDist2Jet[it],w);
	      else if(abs(refparton_flavorForB[ij])<99)hipDist2JetL->Fill(ipDist2Jet[it],w);
	    }
	  
	    hipDist2JetSig->Fill(ipDist2JetSig[it],w);    
	    if(isMC){
	      if(abs(refparton_flavorForB[ij])==5)hipDist2JetSigB->Fill(ipDist2JetSig[it],w);
	      else if(abs(refparton_flavorForB[ij])==4)hipDist2JetSigC->Fill(ipDist2JetSig[it],w);
	      else if(abs(refparton_flavorForB[ij])<99)hipDist2JetSigL->Fill(ipDist2JetSig[it],w);
	    }
	  
	    hipClosest2Jet->Fill(ipClosest2Jet[it],w);    
	    if(isMC){
	      if(abs(refparton_flavorForB[ij])==5)hipClosest2JetB->Fill(ipClosest2Jet[it],w);
	      else if(abs(refparton_flavorForB[ij])==4)hipClosest2JetC->Fill(ipClosest2Jet[it],w);
	      else if(abs(refparton_flavorForB[ij])<99)hipClosest2JetL->Fill(ipClosest2Jet[it],w);
	    }
	  	  
	  }

	  //if(jtpt[ij]<90){
	  if(jtpt[ij]<200){
	    hip2d1->Fill(ip2d1MostSig,w);    
	    if(isMC){
	      if(abs(refparton_flavorForB[ij])==5)hip2d1B->Fill(ip2d1MostSig,w);
	      else if(abs(refparton_flavorForB[ij])==4)hip2d1C->Fill(ip2d1MostSig,w);
	      else if(abs(refparton_flavorForB[ij])<99)hip2d1L->Fill(ip2d1MostSig,w);
	    }
	    hip2dSig1->Fill(ip2dSig1MostSig,w);    
	    if(isMC){
	      if(abs(refparton_flavorForB[ij])==5)hip2dSig1B->Fill(ip2dSig1MostSig,w);
	      else if(abs(refparton_flavorForB[ij])==4)hip2dSig1C->Fill(ip2dSig1MostSig,w);
	      else if(abs(refparton_flavorForB[ij])<99)hip2dSig1L->Fill(ip2dSig1MostSig,w);
	    }
	    hip3d1->Fill(ip3d1MostSig,w);    
	    if(isMC){
	      if(abs(refparton_flavorForB[ij])==5)hip3d1B->Fill(ip3d1MostSig,w);
	      else if(abs(refparton_flavorForB[ij])==4)hip3d1C->Fill(ip3d1MostSig,w);
	      else if(abs(refparton_flavorForB[ij])<99)hip3d1L->Fill(ip3d1MostSig,w);
	    }
	    hip3dSig1->Fill(ip3dSig1MostSig,w);    
	    if(isMC){
	      if(abs(refparton_flavorForB[ij])==5)hip3dSig1B->Fill(ip3dSig1MostSig,w);
	      else if(abs(refparton_flavorForB[ij])==4)hip3dSig1C->Fill(ip3dSig1MostSig,w);
	      else if(abs(refparton_flavorForB[ij])<99)hip3dSig1L->Fill(ip3dSig1MostSig,w);
	    }

	    hip2d2->Fill(ip2d2MostSig,w);    
	    if(isMC){
	      if(abs(refparton_flavorForB[ij])==5)hip2d2B->Fill(ip2d2MostSig,w);
	      else if(abs(refparton_flavorForB[ij])==4)hip2d2C->Fill(ip2d2MostSig,w);
	      else if(abs(refparton_flavorForB[ij])<99)hip2d2L->Fill(ip2d2MostSig,w);
	    }
	    hip2dSig2->Fill(ip2dSig2MostSig,w);    
	    if(isMC){
	      if(abs(refparton_flavorForB[ij])==5)hip2dSig2B->Fill(ip2dSig2MostSig,w);
	      else if(abs(refparton_flavorForB[ij])==4)hip2dSig2C->Fill(ip2dSig2MostSig,w);
	      else if(abs(refparton_flavorForB[ij])<99)hip2dSig2L->Fill(ip2dSig2MostSig,w);
	    }
	    hip3d2->Fill(ip3d2MostSig,w);    
	    if(isMC){
	      if(abs(refparton_flavorForB[ij])==5)hip3d2B->Fill(ip3d2MostSig,w);
	      else if(abs(refparton_flavorForB[ij])==4)hip3d2C->Fill(ip3d2MostSig,w);
	      else if(abs(refparton_flavorForB[ij])<99)hip3d2L->Fill(ip3d2MostSig,w);
	    }
	    hip3dSig2->Fill(ip3dSig2MostSig,w);    
	    if(isMC){
	      if(abs(refparton_flavorForB[ij])==5)hip3dSig2B->Fill(ip3dSig2MostSig,w);
	      else if(abs(refparton_flavorForB[ij])==4)hip3dSig2C->Fill(ip3dSig2MostSig,w);
	      else if(abs(refparton_flavorForB[ij])<99)hip3dSig2L->Fill(ip3dSig2MostSig,w);
	    }

	    hip2d3->Fill(ip2d3MostSig,w);    
	    if(isMC){
	      if(abs(refparton_flavorForB[ij])==5)hip2d3B->Fill(ip2d3MostSig,w);
	      else if(abs(refparton_flavorForB[ij])==4)hip2d3C->Fill(ip2d3MostSig,w);
	      else if(abs(refparton_flavorForB[ij])<99)hip2d3L->Fill(ip2d3MostSig,w);
	    }
	    hip2dSig3->Fill(ip2dSig3MostSig,w);    
	    if(isMC){
	      if(abs(refparton_flavorForB[ij])==5)hip2dSig3B->Fill(ip2dSig3MostSig,w);
	      else if(abs(refparton_flavorForB[ij])==4)hip2dSig3C->Fill(ip2dSig3MostSig,w);
	      else if(abs(refparton_flavorForB[ij])<99)hip2dSig3L->Fill(ip2dSig3MostSig,w);
	    }
	    hip3d3->Fill(ip3d3MostSig,w);    
	    if(isMC){
	      if(abs(refparton_flavorForB[ij])==5)hip3d3B->Fill(ip3d3MostSig,w);
	      else if(abs(refparton_flavorForB[ij])==4)hip3d3C->Fill(ip3d3MostSig,w);
	      else if(abs(refparton_flavorForB[ij])<99)hip3d3L->Fill(ip3d3MostSig,w);
	    }
	    hip3dSig3->Fill(ip3dSig3MostSig,w);    
	    if(isMC){
	      if(abs(refparton_flavorForB[ij])==5)hip3dSig3B->Fill(ip3dSig3MostSig,w);
	      else if(abs(refparton_flavorForB[ij])==4)hip3dSig3C->Fill(ip3dSig3MostSig,w);
	      else if(abs(refparton_flavorForB[ij])<99)hip3dSig3L->Fill(ip3dSig3MostSig,w);
	    }	
	  }
	}
      }
      if(useEvent){
	if(isMC){
	  hbinw->Fill(bin,w);
	  hbin->Fill(bin,xSecWeight*vzWeight);
	}
	else hbin->Fill(bin);
	
	if(isMC){
	  hvzw->Fill(vz,w);
	  if(ppPbPb)hvz->Fill(vz,xSecWeight*centWeight);
	  else hvz->Fill(vz,xSecWeight);
	}
	else hvz->Fill(vz);
      }
      
    }
  }
  fout->cd();

  hbin->Write(); hbinw->Write(); hvz->Write(); hvzw->Write();
  
  hjtpt->Write();
  if(isMC) hjtptB->Write(); hjtptC->Write(); hjtptL->Write(); hjtptU->Write();
  
  hrawpt->Write();
  if(isMC) hrawptB->Write(); hrawptC->Write(); hrawptL->Write(); 
  
  hjteta->Write();
  if(isMC) hjtetaB->Write(); hjtetaC->Write(); hjtetaL->Write(); 
  
  hjtphi->Write();
  if(isMC) hjtphiB->Write(); hjtphiC->Write(); hjtphiL->Write(); 
  
  hdiscr_csvSimple->Write();
  if(isMC) hdiscr_csvSimpleB->Write(); hdiscr_csvSimpleC->Write(); hdiscr_csvSimpleL->Write(); 
  
  hdiscr_prob->Write();
  if(isMC) hdiscr_probB->Write(); hdiscr_probC->Write(); hdiscr_probL->Write(); 
  
  hdiscr_ssvHighEff->Write();
  if(isMC) hdiscr_ssvHighEffB->Write(); hdiscr_ssvHighEffC->Write(); hdiscr_ssvHighEffL->Write(); 

  hdiscr_ssvHighPur->Write();
  if(isMC) hdiscr_ssvHighPurB->Write(); hdiscr_ssvHighPurC->Write(); hdiscr_ssvHighPurL->Write(); 

  hdiscr_tcHighEff->Write();
  if(isMC) hdiscr_tcHighEffB->Write(); hdiscr_tcHighEffC->Write(); hdiscr_tcHighEffL->Write(); 

  hdiscr_tcHighPur->Write();
  if(isMC) hdiscr_tcHighPurB->Write(); hdiscr_tcHighPurC->Write(); hdiscr_tcHighPurL->Write(); 

  hnsvtx->Write();
  if(isMC) hnsvtxB->Write(); hnsvtxC->Write(); hnsvtxL->Write();

  hsvtxntrk->Write();
  if(isMC) hsvtxntrkB->Write();hsvtxntrkC->Write(); hsvtxntrkL->Write();
  
  hsvtxdl->Write();
  if(isMC) hsvtxdlB->Write(); hsvtxdlC->Write(); hsvtxdlL->Write();

  hsvtxdls->Write();
  if(isMC) hsvtxdlsB->Write(); hsvtxdlsC->Write(); hsvtxdlsL->Write();

  hsvtxm->Write();
  if(isMC) hsvtxmB->Write(); hsvtxmC->Write(); hsvtxmL->Write();

  hsvtxmSV3->Write();
  if(isMC) hsvtxmSV3B->Write(); hsvtxmSV3C->Write(); hsvtxmSV3L->Write();

  hsvtxpt->Write();
  if(isMC) hsvtxptB->Write(); hsvtxptC->Write(); hsvtxptL->Write();

  hsvtxptSV3->Write();
  if(isMC) hsvtxptSV3B->Write(); hsvtxptSV3C->Write(); hsvtxptSV3L->Write();

  hnIPtrk->Write();
  if(isMC) hnIPtrkB->Write(); hnIPtrkC->Write(); hnIPtrkL->Write();

  hnselIPtrk->Write();
  if(isMC) hnselIPtrkB->Write(); hnselIPtrkC->Write(); hnselIPtrkL->Write();

  hmuptrel->Write();
  if(isMC) hmuptrelB->Write(); hmuptrelC->Write(); hmuptrelL->Write();

  hmuptrelSV2->Write();
  if(isMC) hmuptrelSV2B->Write(); hmuptrelSV2C->Write(); hmuptrelSV2L->Write();

  hmuptrelSV3->Write();
  if(isMC) hmuptrelSV3B->Write(); hmuptrelSV3C->Write(); hmuptrelSV3L->Write();

  hipPt->Write();
  if(isMC) hipPtB->Write(); hipPtC->Write(); hipPtL->Write();

  hipProb0->Write();
  if(isMC) hipProb0B->Write(); hipProb0C->Write(); hipProb0L->Write();

  hipProb1->Write();
  if(isMC) hipProb1B->Write(); hipProb1C->Write(); hipProb1L->Write();

  hip2d->Write();
  if(isMC) hip2dB->Write(); hip2dC->Write(); hip2dL->Write();

  hip2dSig->Write();
  if(isMC) hip2dSigB->Write(); hip2dSigC->Write(); hip2dSigL->Write();

  hip2d1->Write();
  if(isMC) hip2d1B->Write(); hip2d1C->Write(); hip2d1L->Write();

  hip2dSig1->Write();
  if(isMC) hip2dSig1B->Write(); hip2dSig1C->Write(); hip2dSig1L->Write();

  hip2d2->Write();
  if(isMC) hip2d2B->Write(); hip2d2C->Write(); hip2d2L->Write();

  hip2dSig2->Write();
  if(isMC) hip2dSig2B->Write(); hip2dSig2C->Write(); hip2dSig2L->Write();

  hip2d3->Write();
  if(isMC) hip2d3B->Write(); hip2d3C->Write(); hip2d3L->Write();

  hip2dSig3->Write();
  if(isMC) hip2dSig3B->Write(); hip2dSig3C->Write(); hip2dSig3L->Write();


  hip3d->Write();
  if(isMC) hip3dB->Write(); hip3dC->Write(); hip3dL->Write();

  hip3dSig->Write();
  if(isMC) hip3dSigB->Write(); hip3dSigC->Write(); hip3dSigL->Write();

  hip3d1->Write();
  if(isMC) hip3d1B->Write(); hip3d1C->Write(); hip3d1L->Write();

  hip3dSig1->Write();
  if(isMC) hip3dSig1B->Write(); hip3dSig1C->Write(); hip3dSig1L->Write();

  hip3d2->Write();
  if(isMC) hip3d2B->Write(); hip3d2C->Write(); hip3d2L->Write();

  hip3dSig2->Write();
  if(isMC) hip3dSig2B->Write(); hip3dSig2C->Write(); hip3dSig2L->Write();

  hip3d3->Write();
  if(isMC) hip3d3B->Write(); hip3d3C->Write(); hip3d3L->Write();

  hip3dSig3->Write();
  if(isMC) hip3dSig3B->Write(); hip3dSig3C->Write(); hip3dSig3L->Write();


  hipDist2Jet->Write();
  if(isMC) hipDist2JetB->Write(); hipDist2JetC->Write(); hipDist2JetL->Write();

  hipDist2JetSig->Write();
  if(isMC) hipDist2JetSigB->Write(); hipDist2JetSigC->Write(); hipDist2JetSigL->Write();

  hipClosest2Jet->Write();
  if(isMC) hipClosest2JetB->Write(); hipClosest2JetC->Write(); hipClosest2JetL->Write();

  nt->Write();
  ntMuReq->Write();
  
  fout->Close();

}
