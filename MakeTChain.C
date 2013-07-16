void MakeTChain(int doQCD=1, int doC = 0, int doB=0){

  if(doC&&doB){
    cout<<" you can't do that on television "<<endl;
    return;
  }

  TChain ch("t");  		   

  if(doC){
    ch.Add("/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/cJet30/merged_bjetAnalyzers_hiReco_offPV_pt30by3_newHydjet_restrictMixTripletA_ipHICalibCentWeight.root/akPu3PFJetAnalyzer/t");
    ch.Add("/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/cJet50/merged_bjetAnalyzers_hiReco_offPV_pt30by3_newHydjet_restrictMixTripletA_ipHICalibCentWeight.root/akPu3PFJetAnalyzer/t");
    ch.Add("/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/cJet65/merged_bjetAnalyzers_hiReco_offPV_pt30by3_newHydjet_restrictMixTripletA_ipHICalibCentWeight.root/akPu3PFJetAnalyzer/t");
    ch.Add("/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/cJet80/merged_bjetAnalyzers_hiReco_offPV_pt30by3_newHydjet_restrictMixTripletA_ipHICalibCentWeight.root/akPu3PFJetAnalyzer/t");
    ch.Add("/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/cJet100/merged_bjetAnalyzers_hiReco_offPV_pt30by3_newHydjet_restrictMixTripletA_ipHICalibCentWeight.root/akPu3PFJetAnalyzer/t");
    ch.Add("/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/cJet120/merged_bjetAnalyzers_hiReco_offPV_pt30by3_newHydjet_restrictMixTripletA_ipHICalibCentWeight.root/akPu3PFJetAnalyzer/t");
    ch.Add("/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/cJet200/merged_bjetAnalyzers_hiReco_offPV_pt30by3_newHydjet_restrictMixTripletA_ipHICalibCentWeight.root/akPu3PFJetAnalyzer/t");
  }
  if(doB){
    ch.Add("/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/bJet30/merged_bjetAnalyzers_hiReco_offPV_pt30by3_newHydjet_restrictMixTripletA_ipHICalibCentWeight.root/akPu3PFJetAnalyzer/t");
    ch.Add("/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/bJet50/merged_bjetAnalyzers_hiReco_offPV_pt30by3_newHydjet_restrictMixTripletA_ipHICalibCentWeight.root/akPu3PFJetAnalyzer/t");
    ch.Add("/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/bJet65/merged_bjetAnalyzers_hiReco_offPV_pt30by3_newHydjet_restrictMixTripletA_ipHICalibCentWeight.root/akPu3PFJetAnalyzer/t");
    ch.Add("/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/bJet80/merged_bjetAnalyzers_hiReco_offPV_pt30by3_newHydjet_restrictMixTripletA_ipHICalibCentWeight.root/akPu3PFJetAnalyzer/t");
    ch.Add("/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/bJet100/merged_bjetAnalyzers_hiReco_offPV_pt30by3_newHydjet_restrictMixTripletA_ipHICalibCentWeight.root/akPu3PFJetAnalyzer/t");
    ch.Add("/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/bJet120/merged_bjetAnalyzers_hiReco_offPV_pt30by3_newHydjet_restrictMixTripletA_ipHICalibCentWeight.root/akPu3PFJetAnalyzer/t");
    ch.Add("/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/bJet200/merged_bjetAnalyzers_hiReco_offPV_pt30by3_newHydjet_restrictMixTripletA_ipHICalibCentWeight.root/akPu3PFJetAnalyzer/t");
  }
  if(doQCD){
    ch.Add("/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/qcd30/merged_bjetAnalyzers_hiReco_offPV_pt30by3_restrictMixTripletA_ipHICalibCentWeight_noTrig.root/akPu3PFJetAnalyzer/t");
    ch.Add("/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/qcd50/merged_bjetAnalyzers_hiReco_offPV_pt30by3_restrictMixTripletA_ipHICalibCentWeight_noTrig.root/akPu3PFJetAnalyzer/t");
    ch.Add("/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/qcd80/merged_bjetAnalyzers_hiReco_offPV_pt30by3_restrictMixTripletA_ipHICalibCentWeight_noTrig.root/akPu3PFJetAnalyzer/t");
    ch.Add("/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/qcd120/merged_bjetAnalyzers_hiReco_offPV_pt30by3_restrictMixTripletA_ipHICalibCentWeight_noTrig.root/akPu3PFJetAnalyzer/t");
    ch.Add("/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/qcd170/merged_bjetAnalyzers_hiReco_offPV_pt30by3_restrictMixTripletA_ipHICalibCentWeight_noTrig.root/akPu3PFJetAnalyzer/t");
    ch.Add("/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/qcd200/merged_bjetAnalyzers_hiReco_offPV_pt30by3_restrictMixTripletA_ipHICalibCentWeight_noTrig.root/akPu3PFJetAnalyzer/t");
  }


  TChain ch2("HltTree");  		   
  if(doC){
    ch2.Add("/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/cJet30/merged_bjetAnalyzers_hiReco_offPV_pt30by3_newHydjet_restrictMixTripletA_ipHICalibCentWeight.root/skimanalysis/HltTree");
    ch2.Add("/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/cJet50/merged_bjetAnalyzers_hiReco_offPV_pt30by3_newHydjet_restrictMixTripletA_ipHICalibCentWeight.root/skimanalysis/HltTree");
    ch2.Add("/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/cJet65/merged_bjetAnalyzers_hiReco_offPV_pt30by3_newHydjet_restrictMixTripletA_ipHICalibCentWeight.root/skimanalysis/HltTree");
    ch2.Add("/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/cJet80/merged_bjetAnalyzers_hiReco_offPV_pt30by3_newHydjet_restrictMixTripletA_ipHICalibCentWeight.root/skimanalysis/HltTree");
    ch2.Add("/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/cJet100/merged_bjetAnalyzers_hiReco_offPV_pt30by3_newHydjet_restrictMixTripletA_ipHICalibCentWeight.root/skimanalysis/HltTree");
    ch2.Add("/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/cJet120/merged_bjetAnalyzers_hiReco_offPV_pt30by3_newHydjet_restrictMixTripletA_ipHICalibCentWeight.root/skimanalysis/HltTree");
    ch2.Add("/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/cJet200/merged_bjetAnalyzers_hiReco_offPV_pt30by3_newHydjet_restrictMixTripletA_ipHICalibCentWeight.root/skimanalysis/HltTree");
  }
  if(doB){
    ch2.Add("/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/bJet30/merged_bjetAnalyzers_hiReco_offPV_pt30by3_newHydjet_restrictMixTripletA_ipHICalibCentWeight.root/skimanalysis/HltTree");
    ch2.Add("/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/bJet50/merged_bjetAnalyzers_hiReco_offPV_pt30by3_newHydjet_restrictMixTripletA_ipHICalibCentWeight.root/skimanalysis/HltTree");
    ch2.Add("/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/bJet65/merged_bjetAnalyzers_hiReco_offPV_pt30by3_newHydjet_restrictMixTripletA_ipHICalibCentWeight.root/skimanalysis/HltTree");
    ch2.Add("/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/bJet80/merged_bjetAnalyzers_hiReco_offPV_pt30by3_newHydjet_restrictMixTripletA_ipHICalibCentWeight.root/skimanalysis/HltTree");
    ch2.Add("/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/bJet100/merged_bjetAnalyzers_hiReco_offPV_pt30by3_newHydjet_restrictMixTripletA_ipHICalibCentWeight.root/skimanalysis/HltTree");
    ch2.Add("/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/bJet120/merged_bjetAnalyzers_hiReco_offPV_pt30by3_newHydjet_restrictMixTripletA_ipHICalibCentWeight.root/skimanalysis/HltTree");
    ch2.Add("/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/bJet200/merged_bjetAnalyzers_hiReco_offPV_pt30by3_newHydjet_restrictMixTripletA_ipHICalibCentWeight.root/skimanalysis/HltTree");
  }
  if(doQCD){
    ch2.Add("/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/qcd30/merged_bjetAnalyzers_hiReco_offPV_pt30by3_restrictMixTripletA_ipHICalibCentWeight_noTrig.root/skimanalysis/HltTree");
    ch2.Add("/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/qcd50/merged_bjetAnalyzers_hiReco_offPV_pt30by3_restrictMixTripletA_ipHICalibCentWeight_noTrig.root/skimanalysis/HltTree");
    ch2.Add("/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/qcd80/merged_bjetAnalyzers_hiReco_offPV_pt30by3_restrictMixTripletA_ipHICalibCentWeight_noTrig.root/skimanalysis/HltTree");
    ch2.Add("/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/qcd120/merged_bjetAnalyzers_hiReco_offPV_pt30by3_restrictMixTripletA_ipHICalibCentWeight_noTrig.root/skimanalysis/HltTree");
    ch2.Add("/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/qcd170/merged_bjetAnalyzers_hiReco_offPV_pt30by3_restrictMixTripletA_ipHICalibCentWeight_noTrig.root/skimanalysis/HltTree");
    ch2.Add("/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/qcd200/merged_bjetAnalyzers_hiReco_offPV_pt30by3_restrictMixTripletA_ipHICalibCentWeight_noTrig.root/skimanalysis/HltTree");
  }
  TChain ch3("HltTree");  		   
  if(doC){
    ch3.Add("/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/cJet30/merged_bjetAnalyzers_hiReco_offPV_pt30by3_newHydjet_restrictMixTripletA_ipHICalibCentWeight.root/hltanalysis/HltTree");
    ch3.Add("/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/cJet50/merged_bjetAnalyzers_hiReco_offPV_pt30by3_newHydjet_restrictMixTripletA_ipHICalibCentWeight.root/hltanalysis/HltTree");
    ch3.Add("/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/cJet65/merged_bjetAnalyzers_hiReco_offPV_pt30by3_newHydjet_restrictMixTripletA_ipHICalibCentWeight.root/hltanalysis/HltTree");
    ch3.Add("/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/cJet80/merged_bjetAnalyzers_hiReco_offPV_pt30by3_newHydjet_restrictMixTripletA_ipHICalibCentWeight.root/hltanalysis/HltTree");
    ch3.Add("/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/cJet100/merged_bjetAnalyzers_hiReco_offPV_pt30by3_newHydjet_restrictMixTripletA_ipHICalibCentWeight.root/hltanalysis/HltTree");
    ch3.Add("/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/cJet120/merged_bjetAnalyzers_hiReco_offPV_pt30by3_newHydjet_restrictMixTripletA_ipHICalibCentWeight.root/hltanalysis/HltTree");
    ch3.Add("/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/cJet200/merged_bjetAnalyzers_hiReco_offPV_pt30by3_newHydjet_restrictMixTripletA_ipHICalibCentWeight.root/hltanalysis/HltTree");
  }
  if(doB){
    ch3.Add("/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/bJet30/merged_bjetAnalyzers_hiReco_offPV_pt30by3_newHydjet_restrictMixTripletA_ipHICalibCentWeight.root/hltanalysis/HltTree");
    ch3.Add("/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/bJet50/merged_bjetAnalyzers_hiReco_offPV_pt30by3_newHydjet_restrictMixTripletA_ipHICalibCentWeight.root/hltanalysis/HltTree");
    ch3.Add("/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/bJet65/merged_bjetAnalyzers_hiReco_offPV_pt30by3_newHydjet_restrictMixTripletA_ipHICalibCentWeight.root/hltanalysis/HltTree");
    ch3.Add("/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/bJet80/merged_bjetAnalyzers_hiReco_offPV_pt30by3_newHydjet_restrictMixTripletA_ipHICalibCentWeight.root/hltanalysis/HltTree");
    ch3.Add("/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/bJet100/merged_bjetAnalyzers_hiReco_offPV_pt30by3_newHydjet_restrictMixTripletA_ipHICalibCentWeight.root/hltanalysis/HltTree");
    ch3.Add("/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/bJet120/merged_bjetAnalyzers_hiReco_offPV_pt30by3_newHydjet_restrictMixTripletA_ipHICalibCentWeight.root/hltanalysis/HltTree");
    ch3.Add("/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/bJet200/merged_bjetAnalyzers_hiReco_offPV_pt30by3_newHydjet_restrictMixTripletA_ipHICalibCentWeight.root/hltanalysis/HltTree");
  }
  if(doQCD){
    ch3.Add("/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/qcd30/merged_bjetAnalyzers_hiReco_offPV_pt30by3_restrictMixTripletA_ipHICalibCentWeight_noTrig.root/hltanalysis/HltTree");
    ch3.Add("/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/qcd50/merged_bjetAnalyzers_hiReco_offPV_pt30by3_restrictMixTripletA_ipHICalibCentWeight_noTrig.root/hltanalysis/HltTree");
    ch3.Add("/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/qcd80/merged_bjetAnalyzers_hiReco_offPV_pt30by3_restrictMixTripletA_ipHICalibCentWeight_noTrig.root/hltanalysis/HltTree");
    ch3.Add("/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/qcd120/merged_bjetAnalyzers_hiReco_offPV_pt30by3_restrictMixTripletA_ipHICalibCentWeight_noTrig.root/hltanalysis/HltTree");
    ch3.Add("/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/qcd170/merged_bjetAnalyzers_hiReco_offPV_pt30by3_restrictMixTripletA_ipHICalibCentWeight_noTrig.root/hltanalysis/HltTree");
    ch3.Add("/data_CMS/cms/mnguyen/bTaggingOutput/hydjetEmbedded/qcd200/merged_bjetAnalyzers_hiReco_offPV_pt30by3_restrictMixTripletA_ipHICalibCentWeight_noTrig.root/hltanalysis/HltTree");
  }

  char filename[500];
  if(doC){
    if(doQCD) sprintf(filename,"merged_bjetAnalyzers_hiReco_offPV_pt30by3_newHydjet_restrictMixTripletA_ipHICalibCentWeight_cJetPlusQCD.root");
    else sprintf(filename,"merged_bjetAnalyzers_hiReco_offPV_pt30by3_newHydjet_restrictMixTripletA_ipHICalibCentWeight_cJet.root");
  }
  if(doB){
    if(doQCD) sprintf(filename,"merged_bjetAnalyzers_hiReco_offPV_pt30by3_newHydjet_restrictMixTripletA_ipHICalibCentWeight_bJetPlusQCD.root");
    else sprintf(filename,"merged_bjetAnalyzers_hiReco_offPV_pt30by3_newHydjet_restrictMixTripletA_ipHICalibCentWeight_bJet.root");
  }
  else if(doQCD) sprintf(filename,"merged_bjetAnalyzers_hiReco_offPV_pt30by3_oldHydjet_restrictMixTripletA_ipHICalibCentWeight_qcd.root");
  TFile *fout =new TFile(filename,"recreate");
  fout->cd();
  TDirectory *t1 = fout->mkdir("akPu3PFJetAnalyzer");
  t1->cd();
  ch.Write();

  TDirectory *t2 = fout->mkdir("skimanalysis");
  t2->cd();
  ch2.Write();

  TDirectory *t3 = fout->mkdir("hltanalysis");
  t3->cd();
  ch3.Write();

  fout->Close();


}
