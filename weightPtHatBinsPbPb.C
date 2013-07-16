#include <iostream>
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TF1.h"
#include "TH1.h"
#include "TMath.h"

void weightPtHatBins(int LCB=0){
    
  gROOT->Reset();
  
  Int_t start=4;
  Int_t N=7;

  TFile *fin[7], *fout[7];
  TTree *tr_in[7], *tr_out[7];
  TTree *tr_in_skim[7], *tr_out_skim[7];

  Int_t bounds[7] = {15,30,50,80,120,170,200};


  Double_t xSections[7]={(0)};

  xSections[1] = 1.079e-02; //30
  xSections[2] = 1.021e-03; //50
  xSections[3] = 9.913e-05; //80         
  xSections[4] = 1.128e-05; //120
  xSections[5] = 1.470e-06; //170
  xSections[6] = 5.310e-07; //200
  
  // multiply by the fraction of events on the interval
  xSections[1] *= 89849./99173.;   //30-50
  xSections[2] *= 138481./153301.;   //50-80
  xSections[3] *= 126859./143194.;   //80-120
  //xSections[4] *= 34457./36104.;   //120-200
  xSections[4] *= 31479./36104.;    //120-170
  xSections[5] *= 31707./49715.;    //170-200

  // fit to centrality distribution
  TF1 *fCent = new TF1("fCent","pol7",0,40);
  fCent->SetParameters(14781.9,-1641.19,127.245,-8.87318,0.41423,-0.011089,0.000154744,-8.76427e-07);

  TFile *fData = new TFile("/grid_mnt/vol__vol1__u/llr/cms/mnguyen/bTagging442p5/CMSSW_4_4_2_patch5/src/bTaggingMacros/histos/PbPbdata.root");
  TH1F *hDataVz = (TH1F *)fData->Get("hvz");
  hDataVz->Rebin(4);
  hDataVz->Scale(1./hDataVz->Integral());

  TFile *fMC = new TFile("/grid_mnt/vol__vol1__u/llr/cms/mnguyen/bTagging442p5/CMSSW_4_4_2_patch5/src/bTaggingMacros/histos/PbPbQCDMC.root");
  TH1F *hMCCent =  (TH1F*)fMC->Get("hbin");
  TH1F *hMCVz = (TH1F *)fMC->Get("hvz");
  hMCVz->Rebin(4);
  hMCVz->Scale(1./hMCVz->Integral());

  char filename[500] = "merged_bjetAnalyzers_hiRecoV3_offPV_centUp_regFix";

  Double_t weight, xSecWeight, centWeight, vzWeight;
  
  for (Int_t it=start; it<N; it++) {
      
    cout<<" file # "<<it<<endl;
    TString inputPath = "/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/";
    TString outputPath = "/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/";

    if(LCB==0)inputPath.Append(Form("qcd%d/%s.root",bounds[it],filename));
    else if(LCB==1)inputPath.Append(Form("cJet%d/%s.root",bounds[it],filename));
    else if(LCB==2)inputPath.Append(Form("bJet%d/%s.root",bounds[it],filename));
    
    cout<<"   reading from "<<inputPath<<endl; 
    fin[it] = new TFile(inputPath);

    fin[it]->cd("/akPu3PFJetAnalyzer");
    tr_in[it] = (TTree*)gDirectory->Get("t");
    fin[it]->cd("/skimanalysis");
    tr_in_skim[it] = (TTree*)gDirectory->Get("HltTree");

    outputPath.Append(Form("%s_weighted_WithUpperCut_%d.root",filename,bounds[it]));
    fout[it] = new TFile(outputPath,"RECREATE");
    cout<<"   writing into "<<outputPath<<endl; 

    fout[it]->mkdir("akPu3PFJetAnalyzer");
    fout[it]->mkdir("skimanalysis");

    tr_out[it] = new TTree("t","Jet Analyzer");
    tr_out_skim[it] = new TTree("HltTree","HltTree");

    char cutname[100];
    if (it<N-1) sprintf(cutname,"pthat>%d&&pthat<%d",bounds[it],bounds[it+1]);
    else sprintf(cutname,"pthat>%d",bounds[it]);
    //sprintf(cutname,"pthat>%d",bounds[it]);
    cout<<cutname<<endl;
    Double_t fentries = (Double_t)tr_in[it]->GetEntries(cutname);
    xSecWeight = xSections[it]/(fentries);
    cout<<"  x-sec weight ("<< bounds[it] <<") : "<< xSecWeight <<" entries: "<<fentries<<endl;


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
    Float_t         discr_csvMva[300];
    Float_t         discr_csvSimple[300];
    Float_t         discr_muByIp3[300];
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
    Float_t ipProb1[10000];
    Float_t ip2d[10000];
    Float_t ip2dSig[10000];
    Float_t ip3d[10000];
    Float_t ip3dSig[10000];
    Float_t ipDist2Jet[10000];
    Float_t ipDist2JetSig[10000];
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

    int nHLTBit;
    bool hltBit[12];


    // Set branch addresses.
    tr_in[it]->SetBranchAddress("evt",&evt);
    tr_in[it]->SetBranchAddress("b",&b);
    //tr_in[it]->SetBranchAddress("hf",&hf);
    tr_in[it]->SetBranchAddress("bin",&bin);
    tr_in[it]->SetBranchAddress("vx",&vx);
    tr_in[it]->SetBranchAddress("vy",&vy);
    tr_in[it]->SetBranchAddress("vz",&vz);
    tr_in[it]->SetBranchAddress("nref",&nref);
    tr_in[it]->SetBranchAddress("rawpt",rawpt);
    tr_in[it]->SetBranchAddress("jtpt",jtpt);
    tr_in[it]->SetBranchAddress("jteta",jteta);
    tr_in[it]->SetBranchAddress("jty",jty);
    tr_in[it]->SetBranchAddress("jtphi",jtphi);
    tr_in[it]->SetBranchAddress("jtpu",jtpu);
    tr_in[it]->SetBranchAddress("discr_ssvHighEff",discr_ssvHighEff);
    tr_in[it]->SetBranchAddress("discr_ssvHighPur",discr_ssvHighPur);
    tr_in[it]->SetBranchAddress("discr_csvMva",discr_csvMva);
    tr_in[it]->SetBranchAddress("discr_csvSimple",discr_csvSimple);
    tr_in[it]->SetBranchAddress("discr_muByIp3",discr_muByIp3);
    tr_in[it]->SetBranchAddress("discr_muByPt",discr_muByPt);
    tr_in[it]->SetBranchAddress("discr_prob",discr_prob);
    tr_in[it]->SetBranchAddress("discr_probb",discr_probb);
    tr_in[it]->SetBranchAddress("discr_tcHighEff",discr_tcHighEff);
    tr_in[it]->SetBranchAddress("discr_tcHighPur",discr_tcHighPur);
    tr_in[it]->SetBranchAddress("nsvtx",nsvtx);
    tr_in[it]->SetBranchAddress("svtxntrk",svtxntrk);
    tr_in[it]->SetBranchAddress("svtxdl",svtxdl);
    tr_in[it]->SetBranchAddress("svtxdls",svtxdls);
    tr_in[it]->SetBranchAddress("svtxm",svtxm);
    tr_in[it]->SetBranchAddress("svtxpt",svtxpt);
    tr_in[it]->SetBranchAddress("nIPtrk",nIPtrk);
    tr_in[it]->SetBranchAddress("nIP",&nIP);
    tr_in[it]->SetBranchAddress("nselIPtrk",nselIPtrk);
    tr_in[it]->SetBranchAddress("ipJetIndex",ipJetIndex);
    tr_in[it]->SetBranchAddress("ipPt",ipPt);
    tr_in[it]->SetBranchAddress("ipProb0",ipProb0);
    tr_in[it]->SetBranchAddress("ipProb1",ipProb1);
    tr_in[it]->SetBranchAddress("ip2d",ip2d);
    tr_in[it]->SetBranchAddress("ip2dSig",ip2dSig);
    tr_in[it]->SetBranchAddress("ip3d",ip3d);
    tr_in[it]->SetBranchAddress("ip3dSig",ip3dSig);
    tr_in[it]->SetBranchAddress("ipDist2Jet",ipDist2Jet);
    tr_in[it]->SetBranchAddress("ipDist2JetSig",ipDist2JetSig);
    tr_in[it]->SetBranchAddress("ipClosest2Jet",ipClosest2Jet);
    //tr_in[it]->SetBranchAddress("mue",mue);
    //tr_in[it]->SetBranchAddress("mupt",mupt);
    //tr_in[it]->SetBranchAddress("mueta",mueta);
    //tr_in[it]->SetBranchAddress("muphi",muphi);
    //tr_in[it]->SetBranchAddress("mudr",mudr);
    //tr_in[it]->SetBranchAddress("muptrel",muptrel);
    //tr_in[it]->SetBranchAddress("muchg",muchg);
    tr_in[it]->SetBranchAddress("pthat",&pthat);
    tr_in[it]->SetBranchAddress("beamId1",&beamId1);
    tr_in[it]->SetBranchAddress("beamId2",&beamId2);
    tr_in[it]->SetBranchAddress("refpt",refpt);
    tr_in[it]->SetBranchAddress("refeta",refeta);
    tr_in[it]->SetBranchAddress("refy",refy);
    tr_in[it]->SetBranchAddress("refphi",refphi);
    tr_in[it]->SetBranchAddress("refdphijt",refdphijt);
    tr_in[it]->SetBranchAddress("refdrjt",refdrjt);
    tr_in[it]->SetBranchAddress("refparton_pt",refparton_pt);
    tr_in[it]->SetBranchAddress("refparton_flavor",refparton_flavor);
    tr_in[it]->SetBranchAddress("refparton_flavorForB",refparton_flavorForB);
    tr_in[it]->SetBranchAddress("nHLTBit",&nHLTBit);
    tr_in[it]->SetBranchAddress("hltBit",hltBit);

    /*
    tr_in[it]->SetBranchAddress("ngen",&ngen);
    tr_in[it]->SetBranchAddress("genmatchindex",genmatchindex);
    tr_in[it]->SetBranchAddress("genpt",genpt);
    tr_in[it]->SetBranchAddress("geneta",geneta);
    tr_in[it]->SetBranchAddress("geny",geny);
    tr_in[it]->SetBranchAddress("genphi",genphi);
    tr_in[it]->SetBranchAddress("gendphijt",gendphijt);
    tr_in[it]->SetBranchAddress("gendrjt",gendrjt);
    */
 
    // Set output branch addresses.
    tr_out[it]->Branch("evt",&evt,"evt/I");
    tr_out[it]->Branch("b",&b,"b/F");
    tr_out[it]->Branch("hf",&hf,"hf/F");
    tr_out[it]->Branch("bin",&bin,"bin/I");
    tr_out[it]->Branch("vx",&vx,"vx/F");
    tr_out[it]->Branch("vy",&vy,"vy/F");
    tr_out[it]->Branch("vz",&vz,"vz/F");
    tr_out[it]->Branch("nref",&nref,"nref/I");
    tr_out[it]->Branch("rawpt",rawpt,"rawpt[nref]/F");
    tr_out[it]->Branch("jtpt",jtpt,"jtpt[nref]/F");
    tr_out[it]->Branch("jteta",jteta,"jteta[nref]/F");
    tr_out[it]->Branch("jty",jty,"jty[nref]/F");
    tr_out[it]->Branch("jtphi",jtphi,"jtphi[nref]/F");
    tr_out[it]->Branch("jtpu",jtpu,"jtpu[nref]/F");
    tr_out[it]->Branch("discr_ssvHighEff",discr_ssvHighEff,"discr_ssvHighEff[nref]/F");
    tr_out[it]->Branch("discr_ssvHighPur",discr_ssvHighPur,"discr_ssvHighPur[nref]/F");
    tr_out[it]->Branch("discr_csvMva",discr_csvMva,"discr_csvMva[nref]/F");
    tr_out[it]->Branch("discr_csvSimple",discr_csvSimple,"discr_csvSimple[nref]/F");
    tr_out[it]->Branch("discr_muByIp3",discr_muByIp3,"discr_muByIp3[nref]/F");
    tr_out[it]->Branch("discr_muByPt",discr_muByPt,"discr_muByPt[nref]/F");
    tr_out[it]->Branch("discr_prob",discr_prob,"discr_prob[nref]/F");
    tr_out[it]->Branch("discr_probb",discr_probb,"discr_probb[nref]/F");
    tr_out[it]->Branch("discr_tcHighEff",discr_tcHighEff,"discr_tcHighEff[nref]/F");
    tr_out[it]->Branch("discr_tcHighPur",discr_tcHighPur,"discr_tcHighPur[nref]/F");
    tr_out[it]->Branch("nsvtx",nsvtx,"nsvtx[nref]/I");
    tr_out[it]->Branch("svtxntrk",svtxntrk,"svtxntrk[nref]/I");
    tr_out[it]->Branch("svtxdl",svtxdl,"svtxdl[nref]/F");
    tr_out[it]->Branch("svtxdls",svtxdls,"svtxdls[nref]/F");
    tr_out[it]->Branch("svtxm",svtxm,"svtxm[nref]/F");
    tr_out[it]->Branch("svtxpt",svtxpt,"svtxpt[nref]/F");
    tr_out[it]->Branch("nIPtrk",nIPtrk,"nIPtrk[nref]/I");
    tr_out[it]->Branch("nselIPtrk",nselIPtrk,"nselIPtrk[nref]/I");
    tr_out[it]->Branch("nIP",&nIP,"nIP/I");
    tr_out[it]->Branch("ipJetIndex",ipJetIndex,"ipJetIndex[nIP]/I");
    tr_out[it]->Branch("ipPt",ipPt,"ipPt[nIP]/F");
    tr_out[it]->Branch("ipProb0",ipProb0,"ipProb0[nIP]/F");
    tr_out[it]->Branch("ipProb1",ipProb1,"ipProb1[nIP]/F");
    tr_out[it]->Branch("ip2d",ip2d,"ip2d[nIP]/F");
    tr_out[it]->Branch("ip2dSig",ip2dSig,"ip2dSig[nIP]/F");
    tr_out[it]->Branch("ip3d",ip3d,"ip3d[nIP]/F");
    tr_out[it]->Branch("ip3dSig",ip3dSig,"ip3dSig[nIP]/F");
    tr_out[it]->Branch("ipDist2Jet",ipDist2Jet,"ipDist2Jet[nIP]/F");
    tr_out[it]->Branch("ipDist2JetSig",ipDist2JetSig,"ipDist2JetSig[nIP]/F");
    tr_out[it]->Branch("ipClosest2Jet",ipClosest2Jet,"ipClosest2Jet[nIP]/F");  
    //tr_out[it]->Branch("mue",mue,"mue[nref]/F");
    //tr_out[it]->Branch("mupt",mupt,"mupt[nref]/F");
    //tr_out[it]->Branch("mueta",mueta,"mueta[nref]/F");
    //tr_out[it]->Branch("muphi",muphi,"muphi[nref]/F");
    //tr_out[it]->Branch("mudr",mudr,"mudr[nref]/F");
    //tr_out[it]->Branch("muptrel",muptrel,"muptre[nref]/F");
    //tr_out[it]->Branch("muchg",muchg,"muchg[nref]/I");
    tr_out[it]->Branch("pthat",&pthat,"pthat/F");
    tr_out[it]->Branch("beamId1",&beamId1,"beamId1/I");
    tr_out[it]->Branch("beamId2",&beamId1,"beamId2/I");
    tr_out[it]->Branch("refpt",refpt,"refpt[nref]/F");
    tr_out[it]->Branch("refeta",refeta,"refeta[nref]/F");
    tr_out[it]->Branch("refy",refy,"refy[nref]/F");
    tr_out[it]->Branch("refphi",refphi,"refphi[nref]/F");
    tr_out[it]->Branch("refdphijt",refdphijt,"refdphijt[nref]/F");
    tr_out[it]->Branch("refdrjt",refdrjt,"refdrjt[nref]/F");
    tr_out[it]->Branch("refparton_pt",refparton_pt,"refparton_pt[nref]/F");
    tr_out[it]->Branch("refparton_flavor",refparton_flavor,"refparton_flavor[nref]/I");
    tr_out[it]->Branch("refparton_flavorForB",refparton_flavorForB,"refparton_flavorForB[nref]/I");
    /*
    tr_out[it]->Branch("ngen",&ngen,"ngen/I");
    tr_out[it]->Branch("genmatchindex",genmatchindex,"genmatchindex[nref]/I");
    tr_out[it]->Branch("genpt",genpt,"genpt[nref]/F");
    tr_out[it]->Branch("geneta",geneta,"geneta[nref]/F");
    tr_out[it]->Branch("geny",geny,"geny[nref]/F");
    tr_out[it]->Branch("genphi",genphi,"genphi[nref]/F");
    tr_out[it]->Branch("gendphijt",gendphijt,"gendphijt[nref]/F");
    tr_out[it]->Branch("gendrjt",gendrjt,"gendrjt[nref]/F");
    */
    tr_out[it]->Branch("weight",&weight,"weight/D");
    tr_out[it]->Branch("xSecWeight",&xSecWeight,"xSecWeight/D");
    tr_out[it]->Branch("centWeight",&centWeight,"centWeight/D");
    tr_out[it]->Branch("vzWeight",&vzWeight,"vzWeight/D");

    tr_out[it]->Branch("nHLTBit",&nHLTBit,"nHLTBit/I");
    tr_out[it]->Branch("hltBit",hltBit,"hltBit[nHLTBit]/O");


    //Declaration of leaves types
    Int_t           pvSel;
    Int_t           hbheNoiseSel;
    Int_t           spikeSel;
    Int_t           collSell;
    Int_t           hltAna;
    //Int_t           superFilterPath;

    // Set branch addresses.
    tr_in_skim[it]->SetBranchAddress("pvSel",&pvSel);
    tr_in_skim[it]->SetBranchAddress("hbheNoiseSel",&hbheNoiseSel);
    tr_in_skim[it]->SetBranchAddress("spikeSel",&spikeSel);
    tr_in_skim[it]->SetBranchAddress("collSell",&collSell);
    tr_in_skim[it]->SetBranchAddress("hltAna",&hltAna);

    // Set output branch addresses.
    tr_out_skim[it]->Branch("pvSel",&pvSel,"pvSel/I");
    tr_out_skim[it]->Branch("hbheNoiseSel",&hbheNoiseSel,"hbheNoiseSel/I");
    tr_out_skim[it]->Branch("spikeSel",&spikeSel,"spikeSel/I");
    tr_out_skim[it]->Branch("collSell",&collSell,"collSell/I");
    tr_out_skim[it]->Branch("hltAna",&hltAna,"hltAna/I");


    Long64_t nentries = tr_in[it]->GetEntries();
    Long64_t nbytes = 0;
    Long64_t nbytes_skim = 0;

    for (Long64_t i=0; i<nentries;i++) {

      if(i%100000==0) cout<<" i = "<<i<<" out of "<<nentries<<" ("<<(int)(100*(float)i/(float)nentries)<<"%)"<<endl;

      nbytes += tr_in[it]->GetEntry(i);
      nbytes_skim += tr_in_skim[it]->GetEntry(i);

      centWeight = fCent->Integral(bin,bin+1)/hMCCent->GetBinContent(bin+1);
      //centWeight = fCent->Integral(bin,bin+1);
      if(centWeight<0) centWeight=0.;

      int vzbin = (int) TMath::Ceil(vz+15.);
      if(vzbin>0&&vzbin<=30)vzWeight = hDataVz->GetBinContent(vzbin)/hMCVz->GetBinContent(vzbin);
      else vzWeight=0.;

      weight=xSecWeight*centWeight*vzWeight;


      if (it<N-1){
	if (pthat>bounds[it+1]) continue;
      }

      tr_out[it]->Fill();
      tr_out_skim[it]->Fill();

    }

    fout[it]->cd("akPu3PFJetAnalyzer");
    tr_out[it]->Write();
    fout[it]->cd("skimanalysis");
    tr_out_skim[it]->Write();

    fout[it]->Close();

  }

}

