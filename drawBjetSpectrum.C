void drawBjetSpectrum(int ppPbPb=1){


  gStyle->SetOptTitle(0);

  if(!ppPbPb) {
    cout<<" find Kurt :-) "<<endl;
    return;
  }

  TFile *fin= new TFile("pPb/output/NewFormatV5_bFractionMCTemplate_pPbpp1_jetptcut30_SSVHEat2.0FixCL0_bin_0_40_eta_0_2.root");

  //TH1F *hRawBData = (TH1F *)fin->Get("hRawBData");
  TH1F *hRawBData = (TH1F *)fin->Get("hIncJetsData");

  // divide out the bin-width
  for(int i=0;i<hRawBData->GetNbinsX();i++){
    float val =     hRawBData->GetBinContent(i+1);
    float err =     hRawBData->GetBinError(i+1);
    float width =   hRawBData->GetBinWidth(i+1);
    hRawBData->SetBinContent(i+1,val/width);
    hRawBData->SetBinError(i+1,err/width);
  }
  
  TH1F *hBEfficiency = (TH1F *)fin->Get("hBEfficiencyMC");
  for(int i=0;i<hRawBData->GetNbinsX();i++){
    float val =     hRawBData->GetBinContent(i+1);
    float err =     hRawBData->GetBinError(i+1);
    float eff =     hBEfficiency->GetBinContent(i+1);
    float effErr =     hBEfficiency->GetBinError(i+1);   
    
    //hRawBData->SetBinContent(i+1,val/eff);
    //hRawBData->SetBinError(i+1, val/eff *sqrt(effErr*effErr/eff/eff + err*err/val/val));
  }




  TCanvas *c1=new TCanvas("c1","c1",1);
  c1->SetLogy();

  hRawBData->SetMarkerStyle(8);
  hRawBData->SetYTitle("dN/dp_{T} (GeV/c)^{-1}");
  hRawBData->Draw();


  TF1 *fpow = new TF1("fpow","[0]*pow(x,[1])",55,250);
 
  TH1F *hRawBDataOrig=hRawBData->Clone("hRawBDataOrig");

  for(int iter=0;iter<4;iter++){

    cout<<" iteration # "<<iter<<endl;

    hRawBData->Fit(fpow);
    
    for(int i=0;i<hRawBData->GetNbinsX();i++){
      float meanVal = fpow->Integral(hRawBData->GetBinLowEdge(i+1),hRawBData->GetBinLowEdge(i+1)+hRawBData->GetBinWidth(i+1))/hRawBData->GetBinWidth(i+1);
      float centVal = fpow->Eval(hRawBData->GetBinCenter(i+1));
      float binShiftCorr = centVal/meanVal;
      cout<<" i "<<i<<" corr "<<binShiftCorr<<endl;
      
      float val =     hRawBDataOrig->GetBinContent(i+1);
      float err =     hRawBDataOrig->GetBinError(i+1);
      
      hRawBData->SetBinContent(i+1,val*binShiftCorr);
      hRawBData->SetBinError(i+1,err*binShiftCorr);
      
    }
  }
  TFile *fout = new TFile("raa_inclJet_pPb_numberator.root","recreate");
  fout->cd();
  hRawBData->Scale(1./3.093E7); //scale by luminosity (5.3 pb^-1 (pp), 30.93 nb^-1 )
  hRawBData->Write();

}
