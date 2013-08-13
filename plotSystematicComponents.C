TH1F *transformReco2Gen(TH1F *h, TH2F *hXform){
  
  if(h->GetNbinsX() != hXform->GetNbinsX()){ cout << "Warning! Matrix binning and histogram binning non-compatible!" << endl; exit(0);}

  TH1F *htemp = h->Clone(h->GetName());
  for(int i=0;i<hXform->GetNbinsX();i++){
    float genMC = 0;
    float genErrMC = 0;
    for(int j=0;j<hXform->GetNbinsY();j++){
      
      float coeff = hXform->GetBinContent(i+1,j+1);
      float recoMC = h->GetBinContent(j+1);
      float recoErrMC = h->GetBinError(j+1);
      genMC += coeff * recoMC;
      genErrMC +=  coeff * recoErrMC * coeff * recoErrMC;
    }
    genErrMC = sqrt(genErrMC);
    htemp->SetBinContent(i+1, genMC);
    //htemp->SetBinError(i+1, genErrMC);
  }
  return htemp;
}

TH1F *zeroErrors(TH1F *h){
  for(int i=0;i<h->GetNbinsX();i++) h->SetBinError(i+1,0);
}

void plotSystematicComponents(bool doTransform=1, bool ppPbPb=1, bool plotSymmetrized=1)
{

  gStyle->SetErrorX(0);

  //  TFile *fMatrix = NULL;
  TH2F *xNorm=NULL;

  if(doTransform){
    if(ppPbPb){
      fMatrix= new TFile("outputTowardsFinal/reco2GenMatrix.root");
      xNorm = (TH2F*)(fMatrix->Get("hRecoVsGenNorm"))->Clone("xNorm");
    }
    else{
      fMatrix = new TFile("output/reco2GenMatrix_pp.root");
      xNorm = (TH2F*)(fMatrix->Get("hRecoVsGenNorm"))->Clone("xNorm");
    }
  }

  TFile *f0, *f1a, *f1b, *f2, *f3, *f4a, *f4b;

  if(ppPbPb){
    // Default parameters and LTJP systematic
    f0 = new TFile("outputTowardsFinal/bFractionMCTemplate_ppPbPb1_SSVHEat2.0FixCL0_bin_0_40_eta_0_2_binomErrors_jet55_wideBin_v2.root");
    // Working point variation
    f1a = new TFile("outputTowardsFinal/bFractionMCTemplate_ppPbPb1_SSVHEat2.5FixCL0_bin_0_40_eta_0_2_binomErrors_jet55_wideBin_v2.root");
    f1b = new TFile("outputTowardsFinal/bFractionMCTemplate_ppPbPb1_SSVHEat1.8FixCL0_bin_0_40_eta_0_2_binomErrors_jet55_wideBin_v2.root");
    // Data-driven template
    f2 = new TFile("outputTowardsFinal/bFractionDataTemplate_ppPbPb1_SSVHEat2.0FixCL1_bin_0_40_eta_0_2_binomErrors_jet55_wideBin_v2.root");
    // fix charm
    f3 = new TFile("outputTowardsFinal/bFractionMCTemplate_ppPbPb1_SSVHEat2.0FixCL1_bin_0_40_eta_0_2_binomErrors_jet55_wideBin_v2.root");
    f4a = new TFile("outputTowardsFinal/bFractionMCTemplate_gspUp_ppPbPb1_SSVHEat2.0FixCL0_bin_0_40_eta_0_2_binomErrors_jet55_wideBin_v2.root");
    f4b = new TFile("outputTowardsFinal/bFractionMCTemplate_gspDown_ppPbPb1_SSVHEat2.0FixCL0_bin_0_40_eta_0_2_binomErrors_jet55_wideBin_v2.root");
  }
  else{
    // Default parameters and LTJP systematic
    f0 = new TFile("output/NewFormatV5_bFractionMCTemplate_pppp1_SSVHEat2.0FixCL0_bin_0_40_eta_0_2.root");
    // Working point variation
    f1a = new TFile("output/NewFormatV5_bFractionMCTemplate_pppp1_SSVHEat2.5FixCL0_bin_0_40_eta_0_2.root");
    f1b = new TFile("output/NewFormatV5_bFractionMCTemplate_pppp1_SSVHEat1.8FixCL0_bin_0_40_eta_0_2.root");
    // Data-driven template
    f2 = new TFile("output/NewFormatV5_bFractionDataTemplate_pppp1_SSVHEat2.0FixCL1_bin_0_40_eta_0_2.root");
    // fix charm
    f3 = new TFile("output/NewFormatV5_bFractionMCTemplate_pppp1_SSVHEat2.0FixCL1_bin_0_40_eta_0_2.root");
  }
  // Divide by efficiency

  // LTJP and default
  TH1F *hRawSpec0 =  (TH1F*)f0->Get("hRawBData");
  TH1F *hEff0 =  (TH1F*)f0->Get("hBEfficiencyMC");
  TH1F *hEffLTJP =  (TH1F*)f0->Get("hBEfficiencyDataLTJP");

  TH1F *hDefault = hRawSpec0->Clone("hDefault");
  TH1F *hLTJP = hRawSpec0->Clone("hLTJP");
  
  hDefault->Divide(hEff0);
  hLTJP->Divide(hEffLTJP);

  if(doTransform){
    hDefault = transformReco2Gen(hDefault, xNorm);
    hLTJP = transformReco2Gen(hLTJP, xNorm);
  }

  // Working point variation
  TH1F *hRawSpec1a =  (TH1F*)f1a->Get("hRawBData");
  TH1F *hEff1a =  (TH1F*)f1a->Get("hBEfficiencyMC");
  TH1F *hRawSpec1b =  (TH1F*)f1b->Get("hRawBData");
  TH1F *hEff1b =  (TH1F*)f1b->Get("hBEfficiencyMC");
  
  TH1F *hWorkingPointUp = hRawSpec1a->Clone("hiWorkingPointUp");
  TH1F *hWorkingPointDown = hRawSpec1b->Clone("hWorkingPointDown");
  
  hWorkingPointUp->Divide(hEff1a);
  hWorkingPointDown->Divide(hEff1b);

 if(doTransform){
  hWorkingPointUp = transformReco2Gen(hWorkingPointUp, xNorm);
  hWorkingPointDown = transformReco2Gen(hWorkingPointDown, xNorm);
 }

  // Data-driven template
  TH1F *hRawSpec2 =  (TH1F*)f2->Get("hRawBData");
  TH1F *hEff2 =  (TH1F*)f2->Get("hBEfficiencyMC");
  
  TH1F *hDataDriven = hRawSpec2->Clone("hDataDriven");
  hDataDriven->Divide(hEff2);

 if(doTransform) hDataDriven = transformReco2Gen(hDataDriven, xNorm);

  // Fix charm
  TH1F *hRawSpec3 =  (TH1F*)f3->Get("hRawBData");
  TH1F *hEff3 =  (TH1F*)f3->Get("hBEfficiencyMC");
  
  TH1F *hCharm = hRawSpec3->Clone("hCharm");
  hCharm->Divide(hEff3);

 if(doTransform) hCharm = transformReco2Gen(hCharm, xNorm);
 
 // gluon systematics
  TH1F *hRawSpec4a =  (TH1F*)f4a->Get("hRawBData");
  TH1F *hEff4a =  (TH1F*)f4a->Get("hBEfficiencyMC");
  
  TH1F *hGlueUp = hRawSpec4a->Clone("hGlueUp");
  hGlueUp->Divide(hEff4a);

  TH1F *hRawSpec4b =  (TH1F*)f4b->Get("hRawBData");
  TH1F *hEff4b =  (TH1F*)f4b->Get("hBEfficiencyMC");
  
  TH1F *hGlueDown = hRawSpec4b->Clone("hGlueDown");
  hGlueDown->Divide(hEff4b);

  if(doTransform) hGlueUp = transformReco2Gen(hGlueUp, xNorm);
  if(doTransform) hGlueDown = transformReco2Gen(hGlueDown, xNorm);

 //plotting style stuff

 hLTJP->SetMarkerColor(2);
 hWorkingPointUp->SetMarkerColor(3);
 hWorkingPointDown->SetMarkerColor(7);
 hDataDriven->SetMarkerColor(6);
 hCharm->SetMarkerColor(4);
 hGlueUp->SetMarkerColor(kOrange+1);
 hGlueDown->SetMarkerColor(kOrange+1);
 
 hLTJP->SetLineColor(2);
 hWorkingPointUp->SetLineColor(3);
 hWorkingPointDown->SetLineColor(7);
 hDataDriven->SetLineColor(6);
 hCharm->SetLineColor(4);
 hGlueUp->SetLineColor(kOrange+1);
 hGlueDown->SetLineColor(kOrange+1);
 
 
 hDefault->Draw();
 hLTJP->Draw("same");
 
 hWorkingPointUp->Draw("same");
 hWorkingPointDown->Draw("same");
 hDataDriven->Draw("same");
 hCharm->Draw("same");
 hGlueUp->Draw("same");
 hGlueDown->Draw("same");
 

  // Now take fractional systematics

  TCanvas *c2=new TCanvas("c2","c2",600,600);

  TH1F *hFracTotal = hDefault->Clone("hFracTotal");
  TH1F *hFracLTJP = hLTJP->Clone("hFracLTJP");
  TH1F *hFracWorkingPointUp = hWorkingPointUp->Clone("hFracWorkingPointUp");
  TH1F *hFracWorkingPointDown = hWorkingPointDown->Clone("hFracWorkingPointDown");
  TH1F *hFracDataDriven = hDataDriven->Clone("hFracDataDriven");
  TH1F *hFracCharm = hCharm->Clone("hFracCharm");
  TH1F *hFracGlueUp = hGlueUp->Clone("hFracGlueUp");
  TH1F *hFracGlueDown = hGlueDown->Clone("hFracGlueDown");
  
  hFracTotal->Reset();
  hFracLTJP->Reset();
  hFracWorkingPointUp->Reset();
  hFracWorkingPointDown->Reset();
  hFracDataDriven->Reset();
  hFracCharm->Reset();
  hFracGlueUp->Reset();
  hFracGlueDown->Reset();

  TH1F *hFracWorkingPoint = hWorkingPointUp->Clone("hFracWorkingPoint");
  TH1F *hFracGlue = hGlueUp->Clone("hFracGlue");
  
  TH1F *hFracTotalInv = (TH1F*)hFracTotal->Clone("hFracTotalInv");
  TH1F *hFracLTJPInv = (TH1F*)hFracLTJP->Clone("hFracLTJPInv");
  TH1F *hFracWorkingPointInv = (TH1F*)hFracWorkingPoint->Clone("hFracWorkingPointInv");
  TH1F *hFracDataDrivenInv = (TH1F*)hFracDataDriven->Clone("hFracDataDrivenInv");
  TH1F *hFracCharmInv = (TH1F*)hFracCharm->Clone("hFracCharmInv");
  TH1F *hFracGlueInv = (TH1F*)hFracGlue->Clone("hFracGlueInv");

  for(int i =0;i < hDefault->GetNbinsX(); i++){

    float valDefault = hDefault->GetBinContent(i+1);

    float fracLTJP = hLTJP->GetBinContent(i+1)/valDefault-1.;
    float fracWorkingPointUp = hWorkingPointUp->GetBinContent(i+1)/valDefault-1.;
    float fracWorkingPointDown = hWorkingPointDown->GetBinContent(i+1)/valDefault-1.;
    float fracDataDriven = hDataDriven->GetBinContent(i+1)/valDefault-1.;
    float fracCharm = hCharm->GetBinContent(i+1)/valDefault-1.;
    float fracGlueUp = hGlueUp->GetBinContent(i+1)/valDefault-1.;
    float fracGlueDown = hGlueDown->GetBinContent(i+1)/valDefault-1.;
    
    hFracLTJP->SetBinContent(i+1,fracLTJP);
    hFracWorkingPointUp->SetBinContent(i+1,fracWorkingPointUp);
    hFracWorkingPointDown->SetBinContent(i+1,fracWorkingPointDown);
    hFracDataDriven->SetBinContent(i+1,fracDataDriven);
    hFracCharm->SetBinContent(i+1,fracCharm);
    hFracGlueUp->SetBinContent(i+1,fracGlueUp);
    hFracGlueDown->SetBinContent(i+1,fracGlueDown);

    float fracWorkingPoint = TMath::Max(fabs(fracWorkingPointUp),fabs(fracWorkingPointDown));
    float fracGlue = TMath::Max(fabs(fracGlueUp),fabs(fracGlueDown));

    hFracWorkingPoint->SetBinContent(i+1,fracWorkingPoint);
    hFracGlue->SetBinContent(i+1,fracGlue);

    hFracLTJPInv->SetBinContent(i+1,-1.*fracLTJP);
    hFracWorkingPointInv->SetBinContent(i+1,-1.*fracWorkingPoint);
    hFracDataDrivenInv->SetBinContent(i+1,-1.*fracDataDriven);
    hFracCharmInv->SetBinContent(i+1,-1.*fracCharm);
    hFracGlueInv->SetBinContent(i+1,-1.*fracGlue);
        

    float fracTotal = sqrt(fracLTJP*fracLTJP+
			   fracWorkingPoint*fracWorkingPoint+
			   fracDataDriven*fracDataDriven+
			   fracCharm*fracCharm+
			   fracGlue*fracGlue
			   );
    
    hFracTotal->SetBinContent(i+1,fracTotal);
    hFracTotalInv->SetBinContent(i+1,-1.*fracTotal);

    // Don't need errors
    /*
    float effDefault = hDefault->GetBinError(i+1);
    float relErrDefault = effDefault/valDefault;

    float relErrLTJP = hLTJP->GetBinError(i+1)/hLTJP->GetBinContent(i+1);
    float relErrWorkingPointUp = hWorkingPointUp->GetBinError(i+1)/hWorkingPointUp->GetBinContent(i+1);
    float relErrWorkingPointDown = hWorkingPointDown->GetBinError(i+1)/hWorkingPointDown->GetBinContent(i+1);
    float relErrDataDriven = hDataDriven->GetBinError(i+1)/hDataDriven->GetBinContent(i+1);
    float relErrCharm = hCharm->GetBinError(i+1)/hCharm->GetBinContent(i+1);


    float errFracLTJP = fracLTJP*sqrt(relErrDefault*relErrDefault+relErrLTJP*relErrLTJP);
    float errFracWorkingPointUp = fracWorkingPointUp*sqrt(relErrDefault*relErrDefault+relErrWorkingPointUp*relErrWorkingPointUp);
    float errFracWorkingPointDown = fracWorkingPointDown*sqrt(relErrDefault*relErrDefault+relErrWorkingPointDown*relErrWorkingPointDown);
    float errFracDataDriven = fracDataDriven*sqrt(relErrDefault*relErrDefault+relErrDataDriven*relErrDataDriven);
    float errFracCharm = fracCharm*sqrt(relErrDefault*relErrDefault+relErrCharm*relErrCharm);

    hFracLTJP->SetBinError(i+1,errFracLTJP);
    hFracWorkingPointUp->SetBinError(i+1,errFracWorkingPointUp);
    hFracWorkingPointDown->SetBinError(i+1,errFracWorkingPointDown);
    hFracDataDriven->SetBinError(i+1,errFracDataDriven);
    hFracCharm->SetBinError(i+1,errFracCharm);
    */
    


  }

  hFracLTJP->GetYaxis()->SetTitle("Relative Error");
  hFracLTJP->GetXaxis()->SetNdivisions(505);

  hFracLTJP->SetMaximum(0.5);
  hFracLTJP->SetMinimum(-0.5);
  hFracLTJP->SetFillColor(hFracLTJP->GetMarkerColor());
  hFracLTJP->SetFillStyle(3018);
  hFracLTJP->Draw("h");

  hFracLTJPInv->SetMaximum(0.5);
  hFracLTJPInv->SetMinimum(-0.5);
  hFracLTJPInv->SetFillColor(hFracLTJP->GetMarkerColor());
  hFracLTJPInv->SetFillStyle(3018);
  if(plotSymmetrized)hFracLTJPInv->Draw("h,same");


  hFracWorkingPointUp->SetFillColor(hFracWorkingPointUp->GetMarkerColor());
  hFracWorkingPointUp->SetFillStyle(3006);
  if(!plotSymmetrized)hFracWorkingPointUp->Draw("h,same");
  hFracWorkingPointDown->SetFillColor(hFracWorkingPointDown->GetMarkerColor());
  hFracWorkingPointDown->SetFillStyle(3006);
  if(!plotSymmetrized)hFracWorkingPointDown->Draw("h,same");

  hFracWorkingPoint = zeroErrors(hFracWorkingPoint);
  hFracWorkingPointInv = zeroErrors(hFracWorkingPointInv);

  hFracWorkingPoint = zeroErrors(hFracWorkingPoint);
  hFracWorkingPoint->SetFillColor(hFracWorkingPoint->GetMarkerColor());
  hFracWorkingPoint->SetFillStyle(3007);
  if(plotSymmetrized)hFracWorkingPoint->Draw("h,same");
  hFracWorkingPointInv->SetFillColor(hFracWorkingPoint->GetMarkerColor());
  hFracWorkingPointInv->SetFillStyle(3007);
  if(plotSymmetrized)hFracWorkingPointInv->Draw("h,same");



  hFracDataDriven->SetFillColor(hFracDataDriven->GetMarkerColor());
  hFracDataDriven->SetFillStyle(3004);
  hFracDataDriven->Draw("h,same");
  hFracDataDrivenInv->SetFillColor(hFracDataDriven->GetMarkerColor());
  hFracDataDrivenInv->SetFillStyle(3004);
  if(plotSymmetrized)hFracDataDrivenInv->Draw("h,same");


  hFracCharm->SetFillColor(hFracCharm->GetMarkerColor());
  hFracCharm->SetFillStyle(3005);
  hFracCharm->Draw("h,same");
  hFracCharmInv->SetFillColor(hFracCharm->GetMarkerColor());
  hFracCharmInv->SetFillStyle(3005);
  if(plotSymmetrized)hFracCharmInv->Draw("h,same");
  
  hFracGlue = zeroErrors(hFracGlue);
  hFracGlueInv = zeroErrors(hFracGlueInv);

  hFracGlue->SetFillColor(hFracGlue->GetMarkerColor());
  hFracGlue->SetFillStyle(3007);
  hFracGlue->Draw("h,same");
  hFracGlueInv->SetFillColor(hFracGlue->GetMarkerColor());
  hFracGlueInv->SetFillStyle(3007);
  if(plotSymmetrized)hFracGlueInv->Draw("h,same");

  hFracTotal->Draw("h,same");
  hFracTotalInv->Draw("h,same");

  

  TLegend *leg=new TLegend(0.5,0.7,0.85,0.95);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(hFracTotal,"Total","l");
  leg->AddEntry(hFracLTJP,"Lifetime Tagger","f");
  if(!plotSymmetrized)leg->AddEntry(hFracWorkingPointUp,"Working point up","f");
  if(plotSymmetrized)leg->AddEntry(hFracWorkingPoint,"Working point","f");
  if(!plotSymmetrized)leg->AddEntry(hFracWorkingPointDown,"Working point down","f");
  leg->AddEntry(hFracDataDriven,"Data-driven template","f");
  leg->AddEntry(hFracCharm,"Charm","f");
  leg->AddEntry(hFracGlue,"Glue","f");
  leg->Draw();

  



  TFile *fout=NULL;
  if(ppPbPb)fout = new TFile("systematics_components_PbPb.root","recreate");
  else fout = new TFile("systematics_components_pp.root","recreate");
  hFracLTJP->Write();
  hFracWorkingPointUp->Write();
  hFracWorkingPointDown->Write();
  hFracWorkingPoint->Write();
  hFracDataDriven->Write();
  hFracCharm->Write();

  hFracLTJPInv->Write();
  hFracWorkingPointInv->Write();
  hFracDataDrivenInv->Write();
  hFracCharmInv->Write();

  fout->Close();
  
}
