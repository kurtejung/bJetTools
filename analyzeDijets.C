#include <iostream>
#include <iterator>
#include <vector>
#include <algorithm>
#include "TH1.h"
#include "TTree.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TROOT.h"

void analyzeDijets(int ppPbPb=1, int isMC=2, int useWeight=1, int doNtuples=1,  int jetTrig=2)
{
  // isMC=0 --> Real data, ==1 --> QCD, ==2 --> bJet, ==3 --> cJet
  Float_t minJetPt=0.;
  if(jetTrig==1) minJetPt=80.;
  if(jetTrig==2) minJetPt=65.;

  //if(!ppPbPb) minJetPt=65;
  Float_t maxJetEta=2;

  
  TFile *fin;

  if(ppPbPb){
    if(isMC==0){
      if(jetTrig==1)fin = new TFile("/data_CMS/cms/mnguyen/bTaggingOutput/PbPbData/pbpbDataJet80_hiRegitSVHighPurity_pt30by3_restrictMixTripletA_jpHICalibRepass/merged_bTagAnalyzers_all.root");
      if(jetTrig==2)fin = new TFile("/data_CMS/cms/mnguyen/bTaggingOutput/PbPbData/pbpbDataJet65_hiRegitSVHighPurity_pt30by3_restrictMixTripletA_jpHICalibRepass/merged_bTagAnalyzers_all.root");
    }
    else if(isMC==1) fin = new TFile("/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/merged_bjetAnalyzers_hiReco_offPV_pt30by3_oldHydjet_restrictMixTripletA_ipHICalibCentWeight_weighted_qcd.root");
    else if(isMC==2) fin = new TFile("/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/merged_bjetAnalyzers_hiReco_offPV_pt30by3_newHydjet_restrictMixTripletA_ipHICalibCentWeight_weighted_bJetPlusQCD.root"); 
    else if(isMC==3)fin = new TFile("/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/merged_bjetAnalyzers_hiReco_offPV_pt30by3_newHydjet_restrictMixTripletA_ipHICalibCentWeight_weighted_cJetPlusQCD.root");

  }
  else{
    cout<<" not defined "<<endl;
    return;
  }

  TTree *t = (TTree*) fin->Get("akPu3PFJetAnalyzer/t");
  //TTree *t = (TTree*) fin->Get("ak5PFJetAnalyzer/t"); //for ppReco_jetTrig
  TTree *tSkim = (TTree*) fin->Get("skimanalysis/HltTree");
  
  TTree *tmu = (TTree*) fin->Get("muonTree/HLTMuTree");
  
  int dupRuns[6] = {181912,181913,181938,181950,181985,182124};
  
  std::vector<int> usedEvents[6];
  int nDup=0;

  //Declaration of leaves types                  
  Int_t           evt;
  Int_t           run;
  Int_t           bin;
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
  Float_t         discr_csvMva[1000];
  Float_t         discr_csvSimple[1000];
  Float_t         discr_muByIp3[1000];
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
  Float_t         mue[1000];
  Float_t         mupt[1000];
  Float_t         muptPF[1000];
  Float_t         mueta[1000];
  Float_t         muphi[1000];
  Float_t         mudr[1000];
  Float_t         muptrel[1000];
  Int_t           muchg[1000];
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

  Double_t         weight, xSecWeight, centWeight, vzWeight;
  
  int nHLTBit;
  bool hltBit[12];

  Int_t pvSel;
  Int_t hbheNoiseSel;
  Int_t spikeSel;
  Int_t collSell;
  
  t->SetBranchAddress("evt",&evt);
  if(isMC==0)tmu->SetBranchAddress("Run",&run);
  t->SetBranchAddress("bin",&bin);           
  t->SetBranchAddress("hf",&hf);           
  t->SetBranchAddress("vz",&vz);           
  t->SetBranchAddress("nref",&nref);
  t->SetBranchAddress("rawpt",rawpt);
  t->SetBranchAddress("jtpt",jtpt);
  t->SetBranchAddress("jteta",jteta);
  t->SetBranchAddress("jty",jty);
  t->SetBranchAddress("jtphi",jtphi);
  t->SetBranchAddress("jtpu",jtpu);
  t->SetBranchAddress("discr_ssvHighEff",discr_ssvHighEff);
  t->SetBranchAddress("discr_ssvHighPur",discr_ssvHighPur);
  t->SetBranchAddress("discr_csvMva",discr_csvMva);
  t->SetBranchAddress("discr_csvSimple",discr_csvSimple);
  t->SetBranchAddress("discr_muByIp3",discr_muByIp3);
  t->SetBranchAddress("discr_muByPt",discr_muByPt);
  t->SetBranchAddress("discr_prob",discr_prob);
  t->SetBranchAddress("discr_probb",discr_probb);
  t->SetBranchAddress("discr_tcHighEff",discr_tcHighEff);
  t->SetBranchAddress("discr_tcHighPur",discr_tcHighPur);
  /*
  t->SetBranchAddress("nsvtx",nsvtx);
  t->SetBranchAddress("svtxntrk",svtxntrk);
  t->SetBranchAddress("svtxdl",svtxdl);
  t->SetBranchAddress("svtxdls",svtxdls);
  t->SetBranchAddress("svtxm",svtxm);
  t->SetBranchAddress("svtxpt",svtxpt);
 
  t->SetBranchAddress("nIPtrk",nIPtrk);
  t->SetBranchAddress("nselIPtrk",nselIPtrk); 
  t->SetBranchAddress("nIP",&nIP);
  */
  t->SetBranchAddress("mupt",mupt);
  t->SetBranchAddress("muptPF",muptPF);

  t->SetBranchAddress("mue",mue);
  t->SetBranchAddress("mueta",mueta);
  t->SetBranchAddress("muphi",muphi);
  t->SetBranchAddress("mudr",mudr);
  t->SetBranchAddress("muptrel",muptrel);
  t->SetBranchAddress("muchg",muchg);

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
    t->SetBranchAddress("refparton_isGSP",refparton_isGSP);
    /*
    t->SetBranchAddress("ngen",&ngen);
    t->SetBranchAddress("genmatchindex",genmatchindex);
    t->SetBranchAddress("genpt",genpt);
    t->SetBranchAddress("geneta",geneta);
    t->SetBranchAddress("geny",geny);
    t->SetBranchAddress("genphi",genphi);
    t->SetBranchAddress("gendphijt",gendphijt);
    t->SetBranchAddress("gendrjt",gendrjt);
    */    
  }


  if(isMC&&useWeight){
    t->SetBranchAddress("weight",&weight);
    t->SetBranchAddress("xSecWeight",&xSecWeight);
    if(ppPbPb)t->SetBranchAddress("centWeight",&centWeight);
    t->SetBranchAddress("vzWeight",&vzWeight);
  }
  t->SetBranchAddress("nHLTBit",&nHLTBit);
  t->SetBranchAddress("hltBit",hltBit);

  tSkim->SetBranchAddress("pvSel",&pvSel);
  tSkim->SetBranchAddress("hbheNoiseSel",&hbheNoiseSel);
  tSkim->SetBranchAddress("spikeSel",&spikeSel);
  tSkim->SetBranchAddress("collSell",&collSell);


  TFile *fout;
  if(ppPbPb){
    if(isMC==0){
      if(jetTrig==1)fout = new TFile("dijetNtuples/dijetNtuple_PbPbdata_pt30by3_jpHICalibRepass_jet80.root","recreate");
      if(jetTrig==2)fout = new TFile("dijetNtuples/dijetNtuple_PbPbdata_pt30by3_jpHICalibRepass_jet65.root","recreate");
    }
    else if(isMC==1) fout = new TFile("dijetNtuples/dijetNtuple_PbPbQCDMC_pt30by3_ipHICalibCentWeight.root","recreate");
    else if(isMC==2) fout = new TFile("dijetNtuples/dijetNtuple_PbPbBMC_pt30by3_ipHICalibCentWeight.root","recreate");
    else if(isMC==3) fout = new TFile("dijetNtuples/dijetNtuple_PbPbCMC_pt30by3_ipHICalibCentWeight.root","recreate");
  }
  else{
    cout<<" Name your file, stupid "<<endl;
  }


  TNtuple *nt;
  if(isMC) nt= new TNtuple("nt","","bin:pt1:pt2:pt3:eta1:eta2:eta3:phi1:phi2:phi3:jetProb1:jetProb2:jetProb3:ssvHE1:ssvHE2:ssvHE3:ssvHP1:ssvHP2:ssvHP3:csvSimple1:csvSimple2:csvSimple3:svtxm1:svtxm2:svtxm3:refpt1:refpt2:refpt3:flavor1:flavor2:flavor3:isGSP1:isGSP2:isGSP3:w:pthat:jet55:jet65:jet80");
  else nt= new TNtuple("nt","","bin:pt1:pt2:pt3:eta1:eta2:eta3:phi1:phi2:phi3:jetProb1:jetProb2:jetProb3:ssvHE1:ssvHE2:ssvHE3:ssvHP1:ssvHP2:ssvHP3:csvSimple1:csvSimple2:csvSimple3:svtxm1:svtxm2:svtxm3");




  

  Long64_t nentries = t->GetEntries();

  int gspCounter=0;

  for (Long64_t i=0; i<nentries;i++) {

    if (i%100000==0) cout<<" i = "<<i<<" out of "<<nentries<<" ("<<(int)(100*(float)i/(float)nentries)<<"%)"<<endl; 

    tSkim->GetEntry(i);
    if(isMC){
      // temporarily remove cuts from MC
      if(!pvSel||!spikeSel) continue; //hbheNoise doesn't work in mixed events
    }
    else{
      //if(!pvSel||!hbheNoiseSel||!spikeSel) continue;
      // turn off spike and on coll Sel
      if(!pvSel||!hbheNoiseSel||!collSell){
	//cout<<" selection failed, pvSel="<<pvSel<<", hbheNoiseSel="<<hbheNoiseSel<<" , collSell="<<collSell<<endl;
	continue;
      }
    }

    t->GetEntry(i);



    if(isMC&&!ppPbPb){
      if(beamId1==2112 || beamId2==2112)  continue;
    }

    if(ppPbPb&&isMC==0){
      if(jetTrig==1&&!hltBit[10]) continue;
      if(jetTrig==2&&!hltBit[9]) continue;
    }
    //if(isMC&&ppPbPb){

    if(fabs(vz)>15.) continue;
    
    // pileup rejection
    if(hf>150000.){
      cout<<" rejecting pileup, "<<" hf "<<hf<<" bin "<<bin<<endl;
      for(int ij=0;ij<nref;ij++) if(jtpt[ij]>65.&&fabs(jteta[ij])<2.)cout<<" # associated tracks =  "<<nselIPtrk[ij]<<endl;
      continue;
    }
    bool isNoise=false;

    for(int ij=0; ij<nref; ij++){	  
      if(jtpt[ij]>4000&&fabs(jteta[ij])<2)cout<<" hello "<<muptPF[0]<<" "<<mupt[0]<<endl; 
      if(jtpt[ij]>minJetPt&&fabs(jteta[ij])<2){
	if(neutralMax[ij]/(neutralMax[ij]+chargedMax[ij]+photonMax[ij])>0.975){
	  //cout<<" cleaning event with jet of  "<<jtpt[ij]<<", eta "<<jteta[ij]<<" noise = "<<neutralMax[ij]/(neutralMax[ij]+chargedMax[ij]+photonMax[ij])<<endl;
	  isNoise=true;
	}
	if(muptPF[ij]>10&&mupt[ij]/muptPF[ij]<0.75){
	  cout<<" cleaning event with jet of  "<<jtpt[ij]<<", eta "<<jteta[ij]<<" muptPF = "<<muptPF[ij]<<" mupt "<<mupt[ij]<<endl;
	  isNoise=true;
	}
      }
    }
    if(isNoise) continue;

    
    if(!isMC&&ppPbPb){
      tmu->GetEntry(i);
      
      bool foundEvt = false;
      for(int irun=0;irun<6;irun++){       
	if(run==dupRuns[irun]) {
	  // binary search does not give the right behavior for some reason
	  //if(binary_search(usedEvents[irun].begin(), usedEvents[irun].end(), evt)) {
	  // use the slower find instead
	  
	  // iterator to vector element:
	  std::vector<int>::iterator myIt = std::find(usedEvents[irun].begin(), usedEvents[irun].end(), evt);
	  if(myIt!=usedEvents[irun].end()){
	    
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


    
    double w=1.;
    if(isMC&&useWeight) w=weight;
    
    int useEvent=0;
    
    //find leading jet
    int index1=-1;
    float pt1=0;
    

    for(int j1=0;j1<nref;j1++){
      if(fabs(jteta[j1])<2&&jtpt[j1]>pt1){
	index1=j1;
	pt1=jtpt[j1];	
      }	
    }
    
    if(index1<0) continue;
    
    useEvent=1;
    
    //find subleading jet
    int index2=-1;
    float pt2=0;
    
    for(int j2=0;j2<nref;j2++){
      if(j2==index1||jtpt[j2]>jtpt[index1]) continue;
      if(fabs(jteta[j2])<2&&jtpt[j2]>pt2){
	index2=j2;
	pt2=jtpt[j2];
      }	
    }
    
      
      
    //find thirdleading jet
    int index3=-1;
    float pt3=0;
    
    for(int j3=0;j3<nref;j3++){
      if(j3==index2||j3==index1||jtpt[j3]>jtpt[index2]) continue;
      if(fabs(jteta[j3])<2&&jtpt[j3]>pt3){
	index3=j3;
	pt3=jtpt[j3];
      }	
    }


    

    float eta1=0; float eta2=0; float eta3=0.;
    float phi1=0; float phi2=0; float phi3=0.;
    float jetProb1=0; float jetProb2=0; float jetProb3=0.;
    float ssvHE1=0; float ssvHE2=0; float ssvHE3=0.;
    float ssvHP1=0; float ssvHP2=0; float ssvHP3=0.;
    float csvSimple1=0; float csvSimple2=0; float csvSimple3=0.;
    float svtxm1=0; float svtxm2=0; float svtxm3=0.;
    float refpt1=0; float refpt2=0; float refpt3=0.;
    float flavor1=0; float flavor2=0; float flavor3=0.;
    float isGSP1=0; float isGSP2=0; float isGSP3=0.;
    

    if(index1>=0){
      eta1=jteta[index1];
      phi1=jtphi[index1];
      jetProb1=discr_prob[index1];
      ssvHE1=discr_ssvHighEff[index1];
      ssvHP1=discr_ssvHighPur[index1];
      csvSimple1=discr_csvSimple[index1];
      svtxm1=svtxm[index1];
      if(isMC){
	flavor1=refparton_flavorForB[index1];
	isGSP1=(float)refparton_isGSP[index1];
	refpt1=refpt[index1];
      }
    }

    if(index2>=0){
      eta2=jteta[index2];
      phi2=jtphi[index2];
      jetProb2=discr_prob[index2];
      ssvHE2=discr_ssvHighEff[index2];
      ssvHP2=discr_ssvHighPur[index2];
      csvSimple2=discr_csvSimple[index2];
      svtxm2=svtxm[index2];
      if(isMC){
	flavor2=refparton_flavorForB[index2];
	isGSP2=(float)refparton_isGSP[index2];
	refpt2=refpt[index2];
      }
    }

    if(index3>=0){
      eta3=jteta[index3];
      phi3=jtphi[index3];
      jetProb3=discr_prob[index3];
      ssvHE3=discr_ssvHighEff[index3];
      ssvHP3=discr_ssvHighPur[index3];
      csvSimple3=discr_csvSimple[index3];
      svtxm3=svtxm[index3];
      if(isMC){
	flavor3=refparton_flavorForB[index3];
	isGSP3=(float)refparton_isGSP[index3];
	refpt3=refpt[index3];
      }
    }
    


    if(doNtuples){
      if(ppPbPb){
	if(isMC){
	  float dummy[39];
	  dummy[0]=bin;
	  dummy[1]=pt1; dummy[2]=pt2; dummy[3]=pt3;
	  dummy[4]=eta1; dummy[5]=eta2; dummy[6]=eta3;
	  dummy[7]=phi1; dummy[8]=phi2; dummy[9]=phi3;
	  dummy[10]=jetProb1; dummy[11]=jetProb2; dummy[12]= jetProb3;
	  dummy[13]=ssvHE1; dummy[14]=ssvHE2; dummy[15]= ssvHE3;
	  dummy[16]=ssvHP1; dummy[17]=ssvHP2; dummy[18]= ssvHP3;
	  dummy[19]=csvSimple1; dummy[20]=csvSimple2; dummy[21]= csvSimple3;
	  dummy[22]=svtxm1; dummy[23]=svtxm2; dummy[24]= svtxm3;	      
	  dummy[25]=refpt1; dummy[26]=refpt2; dummy[27]=refpt3;
	  dummy[28]=flavor1; dummy[29]=flavor2; dummy[30]= flavor3;
	  dummy[31]=isGSP1; dummy[32]=isGSP2; dummy[33]= isGSP3;
	  dummy[34]=w; dummy[35]=pthat; 	      
	  dummy[36]=hltBit[8]; dummy[37]=hltBit[9]; dummy[38]=hltBit[10];
	  if(isMC&&minJetPt>65){
	    if(pthat<50) dummy[38]=0;
	  }
	  nt->Fill(dummy);
	}
	else{
	  float dummy[25];
	  dummy[0]=bin;
	  dummy[1]=pt1; dummy[2]=pt2; dummy[3]=pt3;
	  dummy[4]=eta1; dummy[5]=eta2; dummy[6]=eta3;
	  dummy[7]=phi1; dummy[8]=phi2; dummy[9]=phi3;
	  dummy[10]=jetProb1; dummy[11]=jetProb2; dummy[12]= jetProb3;
	  dummy[13]=ssvHE1; dummy[14]=ssvHE2; dummy[15]= ssvHE3;
	  dummy[16]=ssvHP1; dummy[17]=ssvHP2; dummy[18]= ssvHP3;
	  dummy[19]=csvSimple1; dummy[20]=csvSimple2; dummy[21]= csvSimple3;
	  dummy[22]=svtxm1; dummy[23]=svtxm2; dummy[24]= svtxm3;	      
	  nt->Fill(dummy);
	}
      }	
    }
  }


  nt->Write();
  
  fout->Close();

}
