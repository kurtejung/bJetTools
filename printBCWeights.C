void printBCWeights(int BorC=0){

  /*
    This macro determines the normalization of the b and c samples relative to the QCD samples, by comparing the number of b-jets or c-jets
    The z-vertex weighting is applied b/c the b and c samples have a slightly different vertex distribution

  */


  float N=8;
  float bounds[8]={30,50,65,80,100,120,170,200};

  //TFile *fin=new TFile("merged_bjetAnalyzers_hiReco_offPV_pt30by3_oldHydjet_restrictMixTripletA_ipHICalibCentWeight_weighted_qcd.root");    
  TFile *fin=new TFile("merged_bjetAnalyzers_hiReco_offPV_pt30by3_restrictMixTripletA_ipHICalibCentWeight_weighted_bJet.root");    
  TTree *t=(TTree*)fin->Get("akPu3PFJetAnalyzer/t");
    

  for(int it=0;it<N;it++){

    //TFile *fin=new TFile(Form("cJet%d/merged_bjetAnalyzers_hiReco_offPV_pt30by3_newHydjet_restrictMixTripletA_ipHICalibCentWeight.root",(int)bounds[it]));
    //TFile *fin=new TFile("merged_bjetAnalyzers_hiReco_offPV_pt30by3_newHydjet_restrictMixTripletA_ipHICalibCentWeight_weighted_cJet.root");    

    char cutname[100], cCut[100], bCut[100];

    if (it<N-1) sprintf(cutname,"centWeight*vzWeight*(pthat>%d&&pthat<%d)",bounds[it],bounds[it+1]);
    else sprintf(cutname,"centWeight*vzWeight*(pthat>%d)",bounds[it]);

    if (it<N-1) sprintf(cCut,"centWeight*vzWeight*(pthat>%d&&pthat<%d&&refpt>0&&abs(jteta)<2&&abs(refparton_flavorForB)==4)",bounds[it],bounds[it+1]);
    else sprintf(cCut,"centWeight*vzWeight*(pthat>%d&&refpt>0&&abs(jteta)<2&&abs(refparton_flavorForB)==4)",bounds[it]);

    if (it<N-1) sprintf(bCut,"centWeight*vzWeight*(pthat>%d&&pthat<%d&&refpt>0&&abs(jteta)<2&&abs(refparton_flavorForB)==5)",bounds[it],bounds[it+1]);
    else sprintf(bCut,"centWeight*vzWeight*(pthat>%d&&refpt>0&&abs(jteta)<2&&abs(refparton_flavorForB)==5)",bounds[it]);


    TH1D *htemp1=new TH1D("htemp1","htemp1",1,0,9999);
    TH1D *htemp2=new TH1D("htemp2","htemp2",1,0,9999);

    if(BorC){
      t->Draw("refpt>>htemp1",cCut);
      double nc = htemp1->Integral();
      //double nentc = t->GetEntries(cutname);
      t->Draw("pthat>>htemp2",cutname);
      double nentc = htemp2->Integral();
      double ratioc=0.;
      if(nentc>0)ratioc = nc/nentc;
    ///cout<<nc<<"/"<<nentc<<"="<<ratioc<<endl;
      cout<<ratioc<<endl;
    }
    else{
      t->Draw("refpt>>htemp1",bCut);
      double nb = htemp1->Integral();
      //double nentb = t->GetEntries(cutname);
      t->Draw("pthat>>htemp2",cutname);
      double nentb = htemp2->Integral();
      double ratiob=0.;
      if(nentb>0)ratiob = nb/nentb;
      ///cout<<nb<<"/"<<nentb<<"="<<ratiob<<endl;
      cout<<ratiob<<endl;
    }
    delete htemp1;
    delete htemp2;



  }



}
