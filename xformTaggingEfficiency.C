void xformTaggingEfficiency(){


  TFile *fMatrix = new TFile("output/reco2GenMatrix.root");
  TFile *fin = new TFile("output/NewFormatV5_bFractionMCTemplate_pppp1_SSVHEat2.0FixCL0_bin_0_40_eta_0_2.root");

  // reco 2 gen matrix
  TH2F *hXform = (TH2F*)fMatrix->Get("hRecoVsGenNorm");

  // grab reco binned histos  
  TH1F *hRecoEffMC = (TH1F*)fin->Get("hBEfficiencyMC");
  TH1F *hRecoEffDataLTJP = (TH1F*)fin->Get("hBEfficiencyDataLTJP");

  TH1F *hRecoPurMC = (TH1F*)fin->Get("hBPurityMC");
  TH1F *hRecoPurData = (TH1F*)fin->Get("hBPurityData");

  TH1F *hRecoSpecMC = (TH1F*)fin->Get("hRawBMC");
  TH1F *hRecoSpecData = (TH1F*)fin->Get("hRawBData");
  

  // declar gen binned histos
  
  TH1F *hGenEffMC = hRecoEffMC->Clone("hGenEffMC");
  hGenEffMC->Reset();
  TH1F *hGenEffDataLTJP = hRecoEffDataLTJP->Clone("hGenEffDataLTJP");
  hGenEffDataLTJP->Reset();

  TH1F *hGenPurMC = hRecoPurMC->Clone("hGenPurMC");
  hGenPurMC->Reset();
  TH1F *hGenPurData = hRecoPurData->Clone("hGenPurData");
  hGenPurData->Reset();

  TH1F *hGenSpecMC = hRecoSpecMC->Clone("hGenSpecMC");
  hGenSpecMC->Reset();
  TH1F *hGenSpecData = hRecoSpecData->Clone("hGenSpecData");
  hGenSpecData->Reset();

  if(hRecoEffMC->GetNbinsX() != hXform->GetNbinsX()){
    cout<<" FAIL "<<endl;
    return; 
  }


  for(int i=0;i<hXform->GetNbinsX();i++){
    
    float genEffMC = 0;
    float genEffErrMC = 0;
    float genEffDataLTJP = 0;
    float genEffErrDataLTJP = 0;

    float genPurMC = 0;
    float genPurErrMC = 0;
    float genPurData = 0;
    float genPurErrData = 0;
    
    float genSpecMC = 0;
    float genSpecErrMC = 0;
    float genSpecData = 0;
    float genSpecErrData = 0;
    
    for(int j=0;j<hXform->GetNbinsY();j++){
      
      float coeff = hXform->GetBinContent(i+1,j+1);
      
      float recoEffMC = hRecoEffMC->GetBinContent(j+1);
      float recoEffErrMC = hRecoEffMC->GetBinError(j+1);
      genEffMC += coeff * recoEffMC;
      genEffErrMC +=  coeff * recoEffErrMC * coeff * recoEffErrMC;

      float recoEffDataLTJP = hRecoEffDataLTJP->GetBinContent(j+1);
      float recoEffErrDataLTJP = hRecoEffDataLTJP->GetBinError(j+1);
      genEffDataLTJP += coeff * recoEffDataLTJP;
      genEffErrDataLTJP +=  coeff * recoEffErrDataLTJP * coeff * recoEffErrDataLTJP;
      
      float recoPurMC = hRecoPurMC->GetBinContent(j+1);
      float recoPurErrMC = hRecoPurMC->GetBinError(j+1);
      genPurMC += coeff * recoPurMC;
      genPurErrMC +=  coeff * recoPurErrMC * coeff * recoPurErrMC;

      float recoPurData = hRecoPurData->GetBinContent(j+1);
      float recoPurErrData = hRecoPurData->GetBinError(j+1);
      genPurData += coeff * recoPurData;
      genPurErrData +=  coeff * recoPurErrData * coeff * recoPurErrData;
      
      float recoSpecMC = hRecoSpecMC->GetBinContent(j+1);
      float recoSpecErrMC = hRecoSpecMC->GetBinError(j+1);
      genSpecMC += coeff * recoSpecMC;
      genSpecErrMC +=  coeff * recoSpecErrMC * coeff * recoSpecErrMC;

      float recoSpecData = hRecoSpecData->GetBinContent(j+1);
      float recoSpecErrData = hRecoSpecData->GetBinError(j+1);
      genSpecData += coeff * recoSpecData;
      genSpecErrData +=  coeff * recoSpecErrData * coeff * recoSpecErrData;
      
    }
    
    genEffErrMC = sqrt(genEffErrMC);
    hGenEffMC->SetBinContent(i+1, recoEffMC);
    hGenEffMC->SetBinError(i+1, recoEffErrMC);

    genEffErrDataLTJP = sqrt(genEffErrDataLTJP);
    hGenEffDataLTJP->SetBinContent(i+1, recoEffDataLTJP);
    hGenEffDataLTJP->SetBinError(i+1, recoEffErrDataLTJP);

    hGenEffMC->SetBinContent(i+1, genEffMC);
    hGenEffMC->SetBinError(i+1, genEffErrMC);

    hGenEffDataLTJP->SetBinContent(i+1, genEffDataLTJP);
    hGenEffDataLTJP->SetBinError(i+1, genEffErrDataLTJP);

    genPurErrMC = sqrt(genPurErrMC);
    hGenPurMC->SetBinContent(i+1, recoPurMC);
    hGenPurMC->SetBinError(i+1, recoPurErrMC);

    genPurErrData = sqrt(genPurErrData);
    hGenPurData->SetBinContent(i+1, recoPurData);
    hGenPurData->SetBinError(i+1, recoPurErrData);

    hGenPurMC->SetBinContent(i+1, genPurMC);
    hGenPurMC->SetBinError(i+1, genPurErrMC);

    hGenPurData->SetBinContent(i+1, genPurData);
    hGenPurData->SetBinError(i+1, genPurErrData);

    genSpecErrMC = sqrt(genSpecErrMC);
    hGenSpecMC->SetBinContent(i+1, recoSpecMC);
    hGenSpecMC->SetBinError(i+1, recoSpecErrMC);

    genSpecErrData = sqrt(genSpecErrData);
    hGenSpecData->SetBinContent(i+1, recoSpecData);
    hGenSpecData->SetBinError(i+1, recoSpecErrData);

    hGenSpecMC->SetBinContent(i+1, genSpecMC);
    hGenSpecMC->SetBinError(i+1, genSpecErrMC);

    hGenSpecData->SetBinContent(i+1, genSpecData);
    hGenSpecData->SetBinError(i+1, genSpecErrData);

    
  }
  
  TFile *fout=new TFile("outputTowardsFinal/genBinnedHistos.root","recreate");

  hRecoEffMC->SetXTitle("recoJet p_{T} (GeV/c)");
  hRecoEffDataLTJP->SetXTitle("recoJet p_{T} (GeV/c)");
  hGenEffMC->SetXTitle("genJet p_{T} (GeV/c)");
  hGenEffDataLTJP->SetXTitle("genJet p_{T} (GeV/c)");
  hRecoEffMC->Write();
  hRecoEffDataLTJP->Write();
  hGenEffMC->Write();
  hGenEffDataLTJP->Write();

  hRecoPurMC->SetXTitle("recoJet p_{T} (GeV/c)");
  hRecoPurData->SetXTitle("recoJet p_{T} (GeV/c)");
  hGenPurMC->SetXTitle("genJet p_{T} (GeV/c)");
  hGenPurData->SetXTitle("genJet p_{T} (GeV/c)");
  hRecoPurMC->Write();
  hRecoPurData->Write();
  hGenPurMC->Write();
  hGenPurData->Write();

  hRecoSpecMC->SetXTitle("recoJet p_{T} (GeV/c)");
  hRecoSpecData->SetXTitle("recoJet p_{T} (GeV/c)");
  hGenSpecMC->SetXTitle("genJet p_{T} (GeV/c)");
  hGenSpecData->SetXTitle("genJet p_{T} (GeV/c)");
  hRecoSpecMC->Write();
  hRecoSpecData->Write();
  hGenSpecMC->Write();
  hGenSpecData->Write();

  fout->Close();
  
}
