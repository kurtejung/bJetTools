#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TF1.h"
#include <iostream>

using namespace std;


void processPbPb(int pp_PbPb=0, int isMC=2, bool ispPb = true, bool recreate_all=true){
  //////////////////////////////////////////////////////////
  //   This file has been automatically generated 
  //     (Thu Jan 19 10:22:25 2012 by ROOT version5.30/03)
  //   from TTree t/akPu3PFpatJets Jet Analysis Tree
  //   found on file: merged_jetTrackAnalyzers_hiRegit.root
  //////////////////////////////////////////////////////////


  //Reset ROOT and connect tree file
  gROOT->Reset();
  
  //Declaration of leaves types

  bool hltBit[12];

  Int_t           hf;
  Int_t           bin;
  Double_t           vz;
  Int_t           nref;
  Double_t         rawpt;
  Double_t         jtpt;
  Double_t         jteta;
  Double_t         jtphi;
  Double_t         discr_ssvHighEff;
  Double_t         discr_ssvHighPur;
  Double_t         discr_csvSimple;
  Double_t         discr_cJetHighEff;
  Double_t         discr_cJetHighPur;
  // Double_t         discr_csvMva;
  // Double_t         discr_muByIp3;
  // Double_t         discr_muByPt;
  Double_t         discr_prob;
  // Double_t         discr_probb;
  // Double_t         discr_tcHighEff;
  //Double_t         discr_tcHighPur;
  Double_t         pthat;
  //Int_t         refparton_flavor;
  Double_t         refparton_flavorForB;
  //Double_t         refparton_flavorForB;
  Double_t         weight;

  Double_t svtxm;

  TFile *f = new TFile();
  TFile *fout;
  int nloop;
  if(recreate_all) nloop=3;
  else nloop=1;
  for(int iloop=1; iloop<=nloop; iloop++){
    if(recreate_all) isMC=iloop;
    if(pp_PbPb){
      if(isMC==1)f = TFile::Open("~/bTagTrees/histos/PbPbMC_QCD_pt30by3_ipHICalibCentWeight_noTrig.root");
      if(isMC==2)f = TFile::Open("~/bTagTrees/histos/PbPbMC_C_pt30by3_ipHICalibCentWeight_noTrig.root");
      if(isMC==3)f = TFile::Open("~/bTagTrees/histos/PbPbMC_B_pt30by3_ipHICalibCentWeight_noTrig.root");
    }
    else if(!ispPb){
      if(isMC==1)f = TFile::Open("/Users/kjung/bTagTrees/pPb/histos/ppMC_ppReco_akPu3PF_QCDjetTrig_etashift_Fix2Sample_MCWeightFinalWithVz_noTrgSelection_Full.root");
      if(isMC==2)f = TFile::Open("/Users/kjung/bTagTrees/pPb/histos/ppMC_ppReco_akPu3PF_QCDjetTrig_etashift_Fix2Sample_MCWeightFinalWithVz_noTrgSelection_Full.root");
      if(isMC==3)f = TFile::Open("/Users/kjung/bTagTrees/pPb/histos/ppMC_ppReco_akPu3PF_QCDjetTrig_etashift_Fix2Sample_MCWeightFinalWithVz_noTrgSelection_Full.root");
    }
    else{
      if(isMC==1)f = TFile::Open("/Users/kjung/charmJets/pPb/input/DMesonCJet_QCDJetOnly_pPbMC_ppReco_akPu3PF_convertToJetTree_temp.root");
      if(isMC==2)f = TFile::Open("/Users/kjung/charmJets/pPb/input/DMesonCJet_QCDJetOnly_pPbMC_ppReco_akPu3PF_convertToJetTree_temp.root");
      if(isMC==3)f = TFile::Open("/Users/kjung/charmJets/pPb/input/DMesonCJet_QCDJetOnly_pPbMC_ppReco_akPu3PF_convertToJetTree_temp.root");
    }

    cout << "opened: " << f->GetName() << endl;
     
    if(f==NULL){
      cout<<"ppPbPb not defined "<<endl;
      exit(0);
    }

    TTree *t = (TTree*)f->Get("jets");
    // Set branch addresses.
    if(pp_PbPb>0){
      t->SetBranchAddress("hltBit",hltBit);
    }
    if(pp_PbPb>0){
      t->SetBranchAddress("hf",&hf);
      t->SetBranchAddress("bin",&bin);
    }
    t->SetBranchAddress("rawpt",&rawpt);
    if(pp_PbPb) t->SetBranchAddress("jtptB",&jtpt);
    else t->SetBranchAddress("jtpt",&jtpt);
    t->SetBranchAddress("jteta",&jteta);
    t->SetBranchAddress("jtphi",&jtphi);
    t->SetBranchAddress("discr_ssvHighEff",&discr_ssvHighEff);
    t->SetBranchAddress("discr_ssvHighPur",&discr_ssvHighPur);
    t->SetBranchAddress("discr_csvSimple",&discr_csvSimple);
    t->SetBranchAddress("discr_cJetHighEff",&discr_cJetHighEff);
    t->SetBranchAddress("discr_cJetHighPur",&discr_cJetHighPur);
    t->SetBranchAddress("discr_prob",&discr_prob);
    t->SetBranchAddress("svtxm",&svtxm);
    t->SetBranchAddress("pthat",&pthat);
    t->SetBranchAddress("refparton_flavorForB",&refparton_flavorForB);
    t->SetBranchAddress("weight",&weight);

    //     This is the loop skeleton
    //       To read only selected branches, Insert statements like:
    // t->SetBranchStatus("*",0);  // disable all branches
    // TTreePlayer->SetBranchStatus("branchname",1);  // activate branchname

    Long64_t nentries = t->GetEntries();

    TH1F *hBottomPt = new TH1F("hBottomPt","hBottomPt",14,60,200);
    TH1F *hCharmPt = new TH1F("hCharmPt","hCharmPt",14,60,200);
    TH1F *hLightPt = new TH1F("hLightPt","hLightPt",14,60,200);

    TH2F *hBottom_ssvHighEff = new TH2F("hBottom_ssvHighEff","hBottom_ssvHighEff",20,1.4,5,40,-0.5,99.5);
    TH2F *hBottom_ssvHighPur = new TH2F("hBottom_ssvHighPur","hBottom_ssvHighPur",20,1.4,5,40,-0.5,99.5);
    TH2F *hBottom_cJetHighEff = new TH2F("hBottom_cJetHighEff","hBottom_cJetHighEff",20,1.4,5,40,-0.5,99.5);
    TH2F *hBottom_cJetHighPur = new TH2F("hBottom_cJetHighPur","hBottom_cJetHighPur",20,1.4,5,40,-0.5,99.5);
    //TH2F *hBottom_csvMva = new TH2F("hBottom_csvMva","hBottom_csvMva",50,0,t->GetMaximum("discr_csvMva")+0.01,40,-0.5,99.5);
    TH2F *hBottom_csvSimple = new TH2F("hBottom_csvSimple","hBottom_csvSimple",50,0,1.01,40,-0.5,99.5);
    TH2F *hBottom_svtxm = new TH2F("hBottom_svtxm","hBottom_svtxm",50,0,t->GetMaximum("svtxm")+0.01,40,-0.5,99.5);
    TH2F *hBottom_prob = new TH2F("hBottom_prob","hBottom_prob",50,0.001,t->GetMaximum("discr_prob")+0.01,40,-0.5,99.5);
    //TH2F *hBottom_probb = new TH2F("hBottom_probb","hBottom_probb",50,0.001,t->GetMaximum("discr_probb")+0.01,40,-0.5,99.5);
    //TH2F *hBottom_tcHighEff = new TH2F("hBottom_tcHighEff","hBottom_tcHighEff",50,0,40,40,-0.5,99.5);
    //TH2F *hBottom_tcHighPur = new TH2F("hBottom_tcHighPur","hBottom_tcHighPur",50,0,40,40,-0.5,99.5);
    //TH2F *hBottom_logtcHighEff = new TH2F("hBottom_logtcHighEff","hBottom_logtcHighEff",50,-1,3.5,40,-0.5,99.5);
    //TH2F *hBottom_logtcHighPur = new TH2F("hBottom_logtcHighPur","hBottom_logtcHighPur",50,-1,3.5,40,-0.5,99.5);

    TH2F *hCharm_ssvHighEff = new TH2F("hCharm_ssvHighEff","hCharm_ssvHighEff",20,1.4,5,40,-0.5,99.5);
    TH2F *hCharm_ssvHighPur = new TH2F("hCharm_ssvHighPur","hCharm_ssvHighPur",20,1.4,5,40,-0.5,99.5);
    TH2F *hCharm_cJetHighEff = new TH2F("hCharm_cJetHighEff","hCharm_cJetHighEff",20,1.4,5,40,-0.5,99.5);
    TH2F *hCharm_cJetHighPur = new TH2F("hCharm_cJetHighPur","hCharm_cJetHighPur",20,1.4,5,40,-0.5,99.5);
    //TH2F *hCharm_csvMva = new TH2F("hCharm_csvMva","hCharm_csvMva",50,0,t->GetMaximum("discr_csvMva")+0.01,40,-0.5,99.5);
    TH2F *hCharm_csvSimple = new TH2F("hCharm_csvSimple","hCharm_csvSimple",50,0,1.01,40,-0.5,99.5);
    TH2F *hCharm_svtxm = new TH2F("hCharm_svtxm","hCharm_svtxm",50,0,t->GetMaximum("svtxm")+0.01,40,-0.5,99.5);
    TH2F *hCharm_prob = new TH2F("hCharm_prob","hCharm_prob",50,0.001,t->GetMaximum("discr_prob")+0.01,40,-0.5,99.5);
    //TH2F *hCharm_probb = new TH2F("hCharm_probb","hCharm_probb",50,0.001,t->GetMaximum("discr_probb")+0.01,40,-0.5,99.5);
    //TH2F *hCharm_tcHighEff = new TH2F("hCharm_tcHighEff","hCharm_tcHighEff",50,0,40,40,-0.5,99.5);
    //TH2F *hCharm_tcHighPur = new TH2F("hCharm_tcHighPur","hCharm_tcHighPur",50,0,40,40,-0.5,99.5);
    //TH2F *hCharm_logtcHighEff = new TH2F("hCharm_logtcHighEff","hCharm_logtcHighEff",50,-1,3.5,40,-0.5,99.5);
    //TH2F *hCharm_logtcHighPur = new TH2F("hCharm_logtcHighPur","hCharm_logtcHighPur",50,-1,3.5,40,-0.5,99.5);

    TH2F *hLight_ssvHighEff = new TH2F("hLight_ssvHighEff","hLight_ssvHighEff",20,1.4,5,40,-0.5,99.5);
    TH2F *hLight_ssvHighPur = new TH2F("hLight_ssvHighPur","hLight_ssvHighPur",20,1.4,5,40,-0.5,99.5);
    TH2F *hLight_cJetHighEff = new TH2F("hLight_cJetHighEff","hLight_cJetHighEff",20,1.4,5,40,-0.5,99.5);
    TH2F *hLight_cJetHighPur = new TH2F("hLight_cJetHighPur","hLight_cJetHighPur",20,1.4,5,40,-0.5,99.5);
    //  TH2F *hLight_csvMva = new TH2F("hLight_csvMva","hLight_csvMva",50,0,t->GetMaximum("discr_csvMva")+0.01,40,-0.5,99.5);
    TH2F *hLight_csvSimple = new TH2F("hLight_csvSimple","hLight_csvSimple",50,0,1.01,40,-0.5,99.5);
    TH2F *hLight_svtxm = new TH2F("hLight_svtxm","hLight_svtxm",50,0,t->GetMaximum("svtxm")+0.01,40,-0.5,99.5);
    TH2F *hLight_prob = new TH2F("hLight_prob","hLight_prob",50,0.001,t->GetMaximum("discr_prob")+0.01,40,-0.5,99.5);
    // TH2F *hLight_probb = new TH2F("hLight_probb","hLight_probb",50,0.001,t->GetMaximum("discr_probb")+0.01,40,-0.5,99.5);
    // TH2F *hLight_tcHighEff = new TH2F("hLight_tcHighEff","hLight_tcHighEff",50,0,40,40,-0.5,99.5);
    // TH2F *hLight_tcHighPur = new TH2F("hLight_tcHighPur","hLight_tcHighPur",50,0,40,40,-0.5,99.5);
    // TH2F *hLight_logtcHighEff = new TH2F("hLight_logtcHighEff","hLight_logtcHighEff",50,-1,3.5,40,-0.5,99.5);
    // TH2F *hLight_logtcHighPur = new TH2F("hLight_logtcHighPur","hLight_logtcHighPur",50,-1,3.5,40,-0.5,99.5);


    //TF1 *weightFunc = new TF1("weightFunc","[0]+[1]*x+[2]*x*x+[3]*x*x*x",0,40);
    //weightFunc->SetParameters(1.43178e+04,-1.13159e+03,3.01112e+01,-2.69347e-01);

    Long64_t nbytes = 0;
    for (Long64_t i=0; i<nentries;i++) {
      nbytes += t->GetEntry(i);
      if(!pp_PbPb) bin=39;

      if(pthat<50) continue;

      // hard coded Jet80                                                                                                                                                     
      if(pp_PbPb){
	//	if(!hltBit[10]) continue;
      }
      //if(isMC&&ppPbPb){                                                                                                                                                          
      // if(fabs(vz)>15.) continue;

      // for(int j=0;j<nref;j++){
	
      if(jtpt<80||fabs(jteta)>2.4) continue;
      //if(discr_csvSimple>0.9) continue;
      //if(fabs(jteta[j])>2.0) continue;
      
      /*
	float weight = 1.;
	if(weightCent&&pp_PbPb) {
	weight = weightFunc->Eval(bin);
	if(weight<0) {
	cout<<" binf "<<bin<<endl;
	weight=0.;
	}
	}
      */
      //cout<<" weight "<<weight<<endl;
      
      if(abs(refparton_flavorForB)==5){
	hBottomPt->Fill(jtpt,weight);
	hBottom_ssvHighEff->Fill(discr_ssvHighEff,bin,weight);
	hBottom_ssvHighPur->Fill(discr_ssvHighPur,bin,weight);
	hBottom_cJetHighEff->Fill(discr_cJetHighEff,bin,weight);
	hBottom_cJetHighPur->Fill(discr_cJetHighPur,bin,weight);
	// hBottom_csvMva->Fill(discr_csvMva,bin,weight);
        hBottom_csvSimple->Fill(discr_csvSimple,bin,weight);
	hBottom_svtxm->Fill(svtxm,bin,weight);
	hBottom_prob->Fill(discr_prob,bin,weight);
	// hBottom_probb->Fill(discr_probb,bin,weight);
	// hBottom_tcHighEff->Fill(discr_tcHighEff,bin,weight);
	// hBottom_tcHighPur->Fill(discr_tcHighPur,bin,weight);
	// if(discr_tcHighEff>-1) hBottom_logtcHighEff->Fill(log(discr_tcHighEff+1),bin,weight);
	// else hBottom_logtcHighEff->Fill(-999.,bin,weight);
	// if(discr_tcHighPur>-1) hBottom_logtcHighPur->Fill(log(discr_tcHighPur+1),bin,weight);
	// else hBottom_logtcHighPur->Fill(-999.,bin,weight);
      }
      else if(abs(refparton_flavorForB)==4){
	hCharmPt->Fill(jtpt,weight);
	hCharm_ssvHighEff->Fill(discr_ssvHighEff,bin,weight);
	hCharm_ssvHighPur->Fill(discr_ssvHighPur,bin,weight);
	hCharm_cJetHighEff->Fill(discr_cJetHighEff,bin,weight);
	hCharm_cJetHighPur->Fill(discr_cJetHighPur,bin,weight);
	//hCharm_csvMva->Fill(discr_csvMva,bin,weight);
	hCharm_csvSimple->Fill(discr_csvSimple,bin,weight);
	hCharm_svtxm->Fill(svtxm,bin,weight);
	hCharm_prob->Fill(discr_prob,bin,weight);
	// hCharm_probb->Fill(discr_probb,bin,weight);
	// hCharm_tcHighEff->Fill(discr_tcHighEff,bin,weight);
	// hCharm_tcHighPur->Fill(discr_tcHighPur,bin,weight);
	// if(discr_tcHighEff>-1) hCharm_logtcHighEff->Fill(log(discr_tcHighEff+1),bin,weight);
	// else hCharm_logtcHighEff->Fill(-999.,bin,weight);
	//if(discr_tcHighPur>-1) hCharm_logtcHighPur->Fill(log(discr_tcHighPur+1),bin,weight);
	// else hCharm_logtcHighPur->Fill(-999.,bin,weight);
	
      }
      else if(refparton_flavorForB>=-99){
	hLightPt->Fill(jtpt,weight);
	hLight_ssvHighEff->Fill(discr_ssvHighEff,bin,weight);
	hLight_ssvHighPur->Fill(discr_ssvHighPur,bin,weight);
	hLight_cJetHighEff->Fill(discr_cJetHighEff,bin,weight);
	hLight_cJetHighPur->Fill(discr_cJetHighPur,bin,weight);
	// hLight_csvMva->Fill(discr_csvMva,bin,weight);
	hLight_csvSimple->Fill(discr_csvSimple,bin,weight);
	hLight_svtxm->Fill(svtxm,bin,weight);
	hLight_prob->Fill(discr_prob,bin,weight);
	// hLight_probb->Fill(discr_probb,bin,weight);
	// hLight_tcHighEff->Fill(discr_tcHighEff,bin,weight);
	//hLight_tcHighPur->Fill(discr_tcHighPur,bin,weight);
	// if(discr_tcHighEff>-1) hLight_logtcHighEff->Fill(log(discr_tcHighEff+1),bin,weight);
	// else hLight_logtcHighEff->Fill(-999.,bin,weight);
	// if(discr_tcHighPur>-1) hLight_logtcHighPur->Fill(log(discr_tcHighPur+1),bin,weight);
	// else hLight_logtcHighPur->Fill(-999.,bin,weight);
      }
      
      
    }
   
    if(pp_PbPb){
      if(isMC==1)fout=new TFile("hist_PbPb_Light.root","RECREATE");
      else if(isMC==2)fout=new TFile("hist_PbPb_Charm.root","RECREATE");
      else if(isMC==3)fout=new TFile("hist_PbPb_Bottom.root","RECREATE");
    }
    else if(!ispPb){
      if(isMC==1) fout=new TFile("hist_pp_Light_5Tev.root","RECREATE");
      else if(isMC==2)fout=new TFile("hist_pp_Charm_5TeV.root","RECREATE");
      else if(isMC==3)fout=new TFile("hist_pp_Bottom_5TeV.root","RECREATE");
    }
    else{
      if(isMC==1) fout=new TFile("hist_pPb_Light_cJet.root","RECREATE");
      else if(isMC==2)fout=new TFile("hist_pPb_Charm_cJet.root","RECREATE");
      else if(isMC==3)fout=new TFile("hist_pPb_Bottom_cJet.root","RECREATE");
    }
   
    hBottomPt->Write();
    hBottom_ssvHighEff->Write();
    hBottom_ssvHighPur->Write();
    hBottom_cJetHighEff->Write();
    hBottom_cJetHighPur->Write();
    //hBottom_csvMva->Write();
    hBottom_csvSimple->Write();
    hBottom_svtxm->Write();
    hBottom_prob->Write();
    // hBottom_probb->Write();
    //hBottom_tcHighEff->Write();
    // hBottom_tcHighPur->Write();
    // hBottom_logtcHighEff->Write();
    // hBottom_logtcHighPur->Write();
    hCharmPt->Write();
    hCharm_ssvHighEff->Write();
    hCharm_ssvHighPur->Write();
    hCharm_cJetHighEff->Write();
    hCharm_cJetHighPur->Write();
    // hCharm_csvMva->Write();
    hCharm_csvSimple->Write();
    hCharm_svtxm->Write();
    hCharm_prob->Write();
    //hCharm_probb->Write();
    // hCharm_tcHighEff->Write();
    //hCharm_tcHighPur->Write();
    // hCharm_logtcHighEff->Write();
    //hCharm_logtcHighPur->Write();
    hLightPt->Write();
    hLight_ssvHighEff->Write();
    hLight_ssvHighPur->Write();
    hLight_cJetHighEff->Write();
    hLight_cJetHighPur->Write();
    // hLight_csvMva->Write();
    hLight_csvSimple->Write();
    hLight_svtxm->Write();
    hLight_prob->Write();
    // hLight_probb->Write();
    // hLight_tcHighEff->Write();
    // hLight_tcHighPur->Write();
    // hLight_logtcHighEff->Write();
    // hLight_logtcHighPur->Write();

    fout->Close();
    f->Close();
  }
}
