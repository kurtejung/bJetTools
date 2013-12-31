

#include "TTree.h"
#include "TFile.h"
#include "TH1.h"

#include <iostream>
#include <fstream>

void TTreeCreator(){


  //for some reason if you don't make your outfile before your ttree, root hates it
  TFile *fout = new TFile("treeCreator.root","recreate");

  //Declare your temporary storage variables
  Double_t t_jtpt, t_jteta, t_jtphi, t_rawpt, t_refpt, t_refparton_flavorForB, t_discr_prob, t_discr_ssvHighEff, t_discr_ssvHighPur, t_discr_csvSimple, t_svtxm;
  Double_t t_pthat, t_bin, t_nMCentries, t_weight;
  Int_t t_HLT_Jet20, t_HLT_Jet40, t_HLT_Jet60, t_HLT_Jet80, t_HLT_Jet100;

  //Create TTree
  TTree *nt = new TTree("nt","");
  //Tell the tree what branches you want and where to look to fill those branches
  nt->Branch("jtpt",&t_jtpt);
  nt->Branch("jteta",&t_jteta);
  nt->Branch("jtphi",&t_jtphi);
  nt->Branch("rawpt",&t_rawpt);
  nt->Branch("refpt",&t_refpt);
  nt->Branch("refparton_flavorForB",&t_refparton_flavorForB);
  nt->Branch("discr_prob",&t_discr_prob);
  nt->Branch("discr_ssvHighEff",&t_discr_ssvHighEff);
  nt->Branch("discr_ssvHighPur",&t_discr_ssvHighPur);
  nt->Branch("discr_csvSimple",&t_discr_csvSimple);
  nt->Branch("svtxm",&t_svtxm);
  nt->Branch("HLT_Jet20_noJetID_v1",&t_HLT_Jet20,"HLT_Jet20_noJetID_v1/I");
  nt->Branch("HLT_Jet40_noJetID_v1",&t_HLT_Jet40,"HLT_Jet40_noJetID_v1/I");
  nt->Branch("HLT_Jet60_noJetID_v1",&t_HLT_Jet60,"HLT_Jet60_noJetID_v1/I");
  nt->Branch("HLT_Jet80_noJetID_v1",&t_HLT_Jet80,"HLT_Jet80_noJetID_v1/I");
  nt->Branch("HLT_Jet100_noJetID_v1",&t_HLT_Jet100,"HLT_Jet100_noJetID_v1/I");
  nt->Branch("pthat",&t_pthat,"pthat/D");


  //Now declare a bunch of variables that will get filled by the incoming Forest file
  //Essentially, these guys will get reset and refilled by root on each new event from the Forest.
  //The arrays are there because the branches in the Forest are arrays of jets.  1000 is just a big number so we don't overload the array.  There's never going to be 1000 jets in an event.
  Int_t           evt;
  Float_t         vz;
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
  Float_t         mue[1000];
  Float_t         mupt[1000];
  Float_t         mueta[1000];
  Float_t         muphi[1000];
  Float_t         mudr[1000];
  Float_t         muptrel[1000];
  Int_t           muchg[1000];
  Float_t         pthat;
  Float_t         refparton_pt[1000];
  Int_t           refparton_flavor[1000];
  Int_t           refparton_flavorForB[1000];
  Int_t 	  HLT_PAJet20_NoJetID_v1;
  Int_t 	  HLT_PAJet40_NoJetID_v1;
  Int_t 	  HLT_PAJet60_NoJetID_v1;
  Int_t 	  HLT_PAJet80_NoJetID_v1;
  Int_t 	  HLT_PAJet100_NoJetID_v1;

  //Now lets loop over our filelist
  TFile *fin;
  std::string infile = "ppNewJEC_BForest.txt";  //txt list of files

  std::ifstream instr(infile.c_str(), std::ifstream::in);
  std::string filename;
  for(int ifile=0; ifile<nFiles; ifile++){

    instr >> filename;
    std::cout << "File: " << filename << std::endl;
    fin = TFile::Open(filename.c_str());

    //Load in the trees from the file

    TTree *t = (TTree*) fin->Get("akPu3PFJetAnalyzer/t");
    TTree *tSkim = (TTree*) fin->Get("skimanalysis/HltTree");
    TTree *tEvt = (TTree*) fin->Get("hiEvtAnalyzer/HiTree");
    if(!t || !tSkim || !tEvt){ cout << "Warning! Can't find one of the trees!" << endl; exit(0);}

    if(tEvt) t->AddFriend("hiEvtAnalyzer/HiTree");    
    t->AddFriend("hltanalysis/HltTree"); 

    //And now tell ROOT where you want to store all the tree branches for each event (i.e. at the variables we just declared above)
    t->SetBranchAddress("evt",&evt);        
    t->SetBranchAddress("vz",&vz);           
    t->SetBranchAddress("nref",&nref);
    t->SetBranchAddress("HLT_PAJet20_NoJetID_v1",&HLT_PAJet20_NoJetID_v1);
    t->SetBranchAddress("HLT_PAJet40_NoJetID_v1",&HLT_PAJet40_NoJetID_v1);
    t->SetBranchAddress("HLT_PAJet60_NoJetID_v1",&HLT_PAJet60_NoJetID_v1);
    t->SetBranchAddress("HLT_PAJet80_NoJetID_v1",&HLT_PAJet80_NoJetID_v1);
    t->SetBranchAddress("HLT_PAJet100_NoJetID_v1",&HLT_PAJet100_NoJetID_v1);
    t->SetBranchAddress("rawpt",rawpt);
    t->SetBranchAddress("jtpt",jtpt);
    t->SetBranchAddress("jteta",jteta);
    t->SetBranchAddress("jtphi",jtphi);
    t->SetBranchAddress("jty",jty);
    t->SetBranchAddress("jtphi",jtphi);
    t->SetBranchAddress("jtpu",jtpu);
    t->SetBranchAddress("discr_ssvHighEff",discr_ssvHighEff);
    t->SetBranchAddress("discr_ssvHighPur",discr_ssvHighPur);
    t->SetBranchAddress("discr_csvSimple",discr_csvSimple);
    t->SetBranchAddress("discr_muByPt",discr_muByPt);
    t->SetBranchAddress("discr_prob",discr_prob);
    t->SetBranchAddress("discr_probb",discr_probb);
    t->SetBranchAddress("discr_tcHighEff",discr_tcHighEff);
    t->SetBranchAddress("discr_tcHighPur",discr_tcHighPur);
    t->SetBranchAddress("nsvtx",nsvtx);

    t->SetBranchAddress("mue",mue);
    t->SetBranchAddress("mupt",mupt);
    t->SetBranchAddress("mueta",mueta);
    t->SetBranchAddress("muphi",muphi);
    t->SetBranchAddress("mudr",mudr);
    t->SetBranchAddress("muptrel",muptrel);
    t->SetBranchAddress("muchg",muchg);

    if(isMC){
      t->SetBranchAddress("pthat",&pthat);
      t->SetBranchAddress("refparton_pt",refparton_pt);
      t->SetBranchAddress("refparton_flavor",refparton_flavor);
      t->SetBranchAddress("refparton_flavorForB",refparton_flavorForB);
    }

    Long64_t nentries = t->GetEntries();
    cout << "entries: "<< nentries << endl;

     for (Long64_t i=0; i<t->GetEntries(); i++) {      
      if (i%100000==0) std::cout<<" i = "<<i<<" out of "<<t->GetEntries()<<" ("<<(int)(100*(float)i/(float)t->GetEntries())<<"%)"<<std::endl; 

      //This does your heavy lifting and loads all branches into your temp variables.
      t->GetEntry(i);

      //Here's the space where you could fill a histogram, for example
      //Or, even better, set your temporary variables (declared at the beginning)

      //Now fill the information in our temporary trees to one line in the new TTree
      nt->Fill();
      //close entries loop
     }

     //close file loop
  }

  //Write your TTree to your output file
  fout->cd();
  nt->Write();
}
  
