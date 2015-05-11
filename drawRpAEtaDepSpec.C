void drawRpAEtaDepSpec(){

  TFile *f1 = new TFile("raa_pPb_numerator-Bjet_FullDataset_ssvhe20_etaCM_v2_bin0_100_eta-2.0To-1.0.root");
  TFile *f2 = new TFile("raa_pPb_numerator-Bjet_FullDataset_ssvhe20_etaCM_v2_bin0_100_eta-1.0To0.0.root");
  TFile *f3 = new TFile("raa_pPb_numerator-Bjet_FullDataset_ssvhe20_etaCM_v2_bin0_100_eta0.0To1.0.root");
  TFile *f4 = new TFile("raa_pPb_numerator-Bjet_FullDataset_ssvhe20_etaCM_v2_bin0_100_eta1.0To2.0.root");
  TFile *f5 = new TFile("raa_pp_denomForpA-Bjet_etaCM_v2_bin0_100_eta-2.0To-1.0.root");
  TFile *f6 = new TFile("raa_pp_denomForpA-Bjet_etaCM_v2_bin0_100_eta-1.0To0.0.root");
  TFile *f7 = new TFile("raa_pp_denomForpA-Bjet_etaCM_v2_bin0_100_eta0.0To1.0.root");
  TFile *f8 = new TFile("raa_pp_denomForpA-Bjet_etaCM_v2_bin0_100_eta1.0To2.0.root");

  TFile *xchk = new TFile("histos/ppMC_ppReco_akPu3PF_QCDjetTrig_etashift_Fix2Sample_MCWeightFinalWithVz_noTrgSelection_Full.root");
  TTree *nt1 = (TTree*)xchk->Get("nt");
  double xbins2[17] = {0, 5, 10, 15, 20, 25, 30, 35, 40, 55, 70, 90, 110, 140, 170, 220, 400};
  TH1D *hf = new TH1D("hf","",16,xbins2); hf->Sumw2();
  //nt1->Project("hf","refpt","weight*(abs(refparton_flavorForB)==5 && jteta>(-2+0.465) && jteta<(-1+0.465) && jtpt>35 && rawpt>20)");
  for(int i=1; i<=hf->GetNbinsX(); i++){
    hf->SetBinError(i,hf->GetBinError(i)/hf->GetBinWidth(i));
    hf->SetBinContent(i,hf->GetBinContent(i)/hf->GetBinWidth(i));
  }
  hf->SetLineColor(kRed+2);
  hf->Scale(1E9);
  //hf->Scale(1./70.);

  TH1D *hspec[8];
  TH1D *hspec[0] = (TH1D*)f1->Get("hReco0")->Clone("hspec_0");
  TH1D *hspec[1] = (TH1D*)f2->Get("hReco0")->Clone("hspec_1");
  TH1D *hspec[2] = (TH1D*)f3->Get("hReco0")->Clone("hspec_2");
  TH1D *hspec[3] = (TH1D*)f4->Get("hReco0")->Clone("hspec_3");
  TH1D *hspec[4] = (TH1D*)f5->Get("hGen_cent1")->Clone("hspec_4");
  TH1D *hspec[5] = (TH1D*)f6->Get("hGen_cent1")->Clone("hspec_5");
  TH1D *hspec[6] = (TH1D*)f7->Get("hGen_cent1")->Clone("hspec_6");
  TH1D *hspec[7] = (TH1D*)f8->Get("hGen_cent1")->Clone("hspec_7");

  int colorArray[4] = {1,2,4,8};
  int styleArray[4] = {20,21,24,25};

  double yerrs[8] = {0.0, 0.107355, 0.173463, 0.217392, 0.245444, 0.267618,0.241417,0.198722};
  double lumiErr[8] = {0.0,0.4,0.4,0.03,0.04,0.04,0.02,0.02};
  for(int i=0; i<8; i++){ yerrs[i] = sqrt(pow(yerrs[i],2)*4+pow(lumiErr[i],2)); }
  double xbins[8] = {47.5,62.5,80,100,125,155,195,310};
  double ybins[4][8];
  double xerrs[8] = {5,5,5,10,10,10,10,40};
  TGraphErrors *specErr[4][8];
  for(int i=0; i<4; i++){
    for(int j=0; j<8; j++){
      //cout << hspec[i]->GetBinCenter(j+9) << endl;
      ybins[i][j] = hspec[i]->GetBinContent(j+9) * pow(10,(i*2)) * 1E9;
      double yerrtot = ybins[i][j]*yerrs[j];
      specErr[i][j] = new TGraphErrors(1,&xbins[j],&ybins[i][j],&xerrs[j],&yerrtot);
      specErr[i][j]->SetFillColor(kYellow);
    }
  }

  for(int i=0; i<8; i++){
    for(int ibin=1; ibin<hspec[i]->GetNbinsX(); ibin++){
      if(hspec[i]->GetBinLowEdge(ibin)<55) hspec[i]->SetBinContent(ibin,-1);
    }
  }

  /* double val1,err1,width1;
  for(int j=0; j<4; j++){
    for(int i=0;i<hspec[j]->GetNbinsX();i++){
      val1 =     hspec[j]->GetBinContent(i+1);
      err1 =     hspec[j]->GetBinError(i+1);
      width1 =   hspec[j]->GetBinWidth(1+1);
      hspec[j]->SetBinContent(i+1,hspec[j]->GetBinContent(i+1)/hspec[j]->GetBinWidth(i+1));
      hspec[j]->SetBinError(i+1,hspec[j]->GetBinError(i+1)/hspec[j]->GetBinWidth(i+1));
    }
    }*/

  gStyle->SetErrorX(0);
  TCanvas *c1 = new TCanvas("c1","",600,700);
  c1->cd();
  c1->SetLogy();
  c1->SetLeftMargin(0.20);
  c1->SetBottomMargin(0.07);
  hspec[0]->SetXTitle("b-jet p_{T} [GeV/c]");
  hspec[0]->SetYTitle("#frac{d#sigma_{b}^{pp}}{dp_{T}}, #frac{1}{T_{pA}} #frac{d^{2}#sigma_{b}^{pA}}{dp_{T} d#eta} #left(#frac{pb}{GeV/c}#right)");
  hspec[0]->GetYaxis()->SetTitleOffset(2.0);
  hspec[0]->GetXaxis()->SetRangeUser(55.1,400);
  hspec[0]->GetXaxis()->CenterTitle(0);
  hspec[0]->GetXaxis()->SetTitleSize(0.035);
  hspec[0]->GetXaxis()->SetLabelSize(0.03);
  hspec[0]->GetXaxis()->SetLabelOffset(0.003);
  hspec[0]->GetYaxis()->SetLabelSize(0.035);
  hspec[0]->GetYaxis()->SetTitleSize(0.04);
  hspec[0]->GetXaxis()->SetTitleOffset(0.9);
  hspec[0]->GetYaxis()->CenterTitle(1);
  hspec[0]->GetXaxis()->CenterTitle(1);
  hspec[0]->SetTitle("");
  //hspec[0]->Scale(1./(32E9*2110E-3));
  hspec[0]->SetMarkerColor(colorArray[0]);
  //hspec[0]->Scale(1./6.9);
  hspec[0]->Scale(1E9);
  hspec[0]->SetMaximum(1E11);
  hspec[0]->SetMinimum(5E-3);
  hspec[0]->Draw();
  for(int i=0; i<8; i++){ specErr[0][i]->Draw("2,same"); }
  hspec[0]->Draw("same");
  for(int j=1; j<4; j++){
    hspec[j]->SetMarkerColor(colorArray[j]);
    hspec[j]->SetMarkerStyle(styleArray[j]);
    hspec[j]->Scale(1E9);
    // hspec[j]->Scale(1./6.9);
    hspec[j]->Scale(pow(10,(j*2)));
    hspec[j]->Draw("same");
    for(int i=0; i<8; i++){ specErr[j][i]->Draw("2,same"); }
    hspec[j]->Draw("same");
  }
  //hspec[4]->Scale(1./4.);
  for(int j=4; j<8; j++){
    hspec[j]->SetMarkerSize(0);
    hspec[j]->Scale(1E9);
    //hspec[j]->Scale(6.9);
    hspec[j]->SetMarkerColor(colorArray[j-4]);
    hspec[j]->Scale(pow(10,((j-4)*2)));
    hspec[j]->SetLineColor(1);
    hspec[j]->Draw("h,same");
  }
  //hf->Draw("h,same");
  TLegend *l1 = new TLegend(0.44,0.73,0.90,0.93);
  l1->SetFillColor(0);
  for(int i=3; i>=0; i--){
    if(i==0) l1->AddEntry(hspec[i],Form("Data, %.1f < #eta_{CM} < %.1f",(float)(i-2),(float)(i-1)),"P");
    else l1->AddEntry(hspec[i],Form("Data, %.1f < #eta_{CM} < %.1f (x 10^{%d})",(float)(i-2),(float)(i-1),(i)*2),"P");
  }
  l1->AddEntry(hspec[4],"PYTHIA (matched #eta_{CM})","Lp");
  l1->Draw("same");
  TLatex *cmsP = new TLatex(70,1.4E11,"CMS ");
  cmsP->SetTextFont(62);
  cmsP->SetTextSize(0.05);
  cmsP->Draw("same");
  TLatex *l4 = new TLatex(252,1.45E11,"35 nb^{-1} (5.02 TeV)");//; PbPb L = 150 #mub^{-1}");
  l4->SetTextFont(43);
  l4->SetTextSize(25);
  l4->Draw("same");
  
  //TLatex *cmsP = new TLatex(255,8E11," CMS                #sqrt{s_{NN}} = 5.02 TeV             L = 35 nb^{-1}");
  //cmsP->SetTextFont(43);
  //cmsP->SetTextSize(20.5);
  // cmsP->Draw("same");

  //gStyle->SetErrorX(1);
  TCanvas *c2 = new TCanvas("c2","",600,700);
  c2->Divide(1,4,0,0.05);
  TLine *line1 = new TLine(55,1,400,1);
  TLatex *latex1[4];
  line1->SetLineStyle(2);
  TH1D *ratios[4];
  TGraphErrors *ratioErr[4][8];
  TLatex *cmsP2 = new TLatex(55,2.31,"CMS");
  TLatex *l5 = new TLatex(255,2.31,"35 nb^{-1} (5.02 TeV)");
  c2->SetLeftMargin(0.23);
  for(int i=3; i>=0; i--){
    ratios[i] = (TH1D*)hspec[i]->Clone(Form("ratios_%d",i));
    ratios[i]->Divide(hspec[i+4]);
    for(int j=0; j<8; j++){
      double ycont = ratios[i]->GetBinContent(j+9);
      double yerrtot = ratios[i]->GetBinContent(j+9)*yerrs[j];
      ratioErr[i][j] = new TGraphErrors(1,&xbins[j],&ycont,&xerrs[j],&yerrtot);
      ratioErr[i][j]->SetFillColor(kYellow);
    }
    c2->GetPad(i+1)->SetRightMargin(0.025);
    if(i!=0) c2->GetPad(i+1)->SetTopMargin(0.01);
    else c2->GetPad(i+1)->SetTopMargin(0.18);
    if(i!=3) c2->GetPad(i+1)->SetBottomMargin(0.01);
    else c2->GetPad(i+1)->SetBottomMargin(0.25);
    ratios[i]->SetMaximum(2.2);
    ratios[i]->SetMinimum(-0.2);
    ratios[i]->SetTitle("");
    ratios[i]->GetYaxis()->SetTitleSize(0.20);
    ratios[i]->GetYaxis()->SetNdivisions(8);
    if(i==0) ratios[i]->SetXTitle("b-jet p_{T} [GeV/c]");
    else ratios[i]->SetXTitle("");
    ratios[i]->GetYaxis()->SetTitleOffset(0.35);
    ratios[i]->GetXaxis()->SetTitleOffset(0.25);
    ratios[i]->GetYaxis()->SetLabelSize(0.12);
    ratios[i]->GetXaxis()->SetRangeUser(55.1,400);
    if(i==2) ratios[i]->SetYTitle("R_{pA}^{PYTHIA}");
    else ratios[i]->SetYTitle("");
    ratios[i]->GetYaxis()->CenterTitle(1);
    c2->cd(4-(i));
    ratios[i]->Draw();
    for(int j=0; j<8; j++){ ratioErr[i][j]->Draw("2,same"); }
    line1->Draw("same");
    ratios[i]->Draw("same");
    latex1[i] = new TLatex(250,1.73,Form("%.1f < #eta_{CM} < %.1f",(float)(i-2),(float)(i-1)));
    if(i!=3) latex1[i]->SetTextSize(0.125);
    else latex1[i]->SetTextSize(0.10);
    latex1[i]->Draw("same");
    if(i==3){
      cmsP2->SetTextFont(62);
      cmsP2->SetTextSize(0.18);
      cmsP2->Draw("same");
      l5->SetTextFont(43);
      l5->SetTextSize(25);
      l5->Draw("same");
    }
  }
  ratios[0]->GetXaxis()->SetLabelSize(0.1);
  ratios[0]->GetYaxis()->SetLabelSize(0.1);
  ratios[0]->GetXaxis()->SetTitleSize(0.12);
  ratios[0]->GetYaxis()->SetTitleOffset(0.45);
  ratios[0]->GetYaxis()->SetTitleSize(0.1);
  ratios[0]->GetXaxis()->SetTitleOffset(0.95);
  
}
