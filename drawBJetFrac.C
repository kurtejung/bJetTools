{

  TFile* fIncl1 = new TFile("pPb_Unfo_v20_akPu3PF_akPu3PF_noGplus_FirstHalfOnly_Converged_InclRecoPt_usedParameterizedUnfold0_jtpt35_bJets_clo0_chi100_v8_eta_-2.00To2.00_.root");
  TFile* fIncl2 = new TFile("pPb_Unfo_v20_Reverse_akPu3PF_akPu3PF_noGplus_SecondHalfOnly_Converged_InclRecoPt_usedParameterizedUnfold0_jtpt35_bJets_clo0_chi100_v8_eta_-2.00To2.00_.root");

  TH1D *InclSpec = (TH1D*)fIncl1->Get("hReco0")->Clone("InclSpec");
  TH1D *InclSpecRev = (TH1D*)fIncl2->Get("hReco0")->Clone("InclSpecRev");
  InclSpec->Add(InclSpecRev);

  cout << "nbins: "<< InclSpec->GetNbinsX();
  cout << "bin scheme:" << endl;
  for(int i=0; i<InclSpec->GetNbinsX(); i++){
    cout << InclSpec->GetBinLowEdge(i+1) << ", ";
  }
  cout << InclSpec->GetBinLowEdge(InclSpec->GetNbinsX())+InclSpec->GetBinWidth(InclSpec->GetNbinsX()) << endl;

  TFile *bj1 = new TFile("pPb_Unfo_inCM_v31_officialMC_ak3PF_akPu3PF_noGplus_FirstHalfOnly_Converged_usedParameterizedUnfold0_jtpt35_bJets_clo0_chi100_v8_eta_-2.00To2.00_.root");
  TFile *bj2 = new TFile("pPb_Unfo_inCM_v31_officialMC_Reverse_WithResCorr_ak3PF_akPu3PF_noGplus_SecondHalfOnly_Converged_usedParameterizedUnfold0_jtpt35_bJets_clo0_chi100_v8_eta_-2.00To2.00_.root");

  TH1D *bSpec = (TH1D*)bj1->Get("hRecoSVD0")->Clone("bSpec");
  TH1D *bSpecRev = (TH1D*)bj2->Get("hRecoSVD0")->Clone("bSpecRev");
  bSpec->Add(bSpecRev);

  double xbins2[17] = {0, 5, 10, 15, 20, 25, 30, 35, 40, 55, 70, 90, 110, 140, 170, 220, 400};
  bSpec = (TH1D*)bSpec->Rebin(16,bSpec->GetName(),xbins2);
  for(int ibin=0; ibin<16; ibin++){
    bSpec->SetBinContent(ibin+1,bSpec->GetBinContent(ibin+1)/bSpec->GetBinWidth(ibin+1));
    bSpec->SetBinError(ibin+1,bSpec->GetBinError(ibin+1)/bSpec->GetBinWidth(ibin+1));
  }
  //bSpec->Scale(35./20.97);

  TCanvas *c2 = new TCanvas("c2","",600,600);
  c2->cd();
  TH1D *bcln = (TH1D*)bSpec->Clone("bSpec");
  TH1D *inclCln = (TH1D*)InclSpec->Clone("inclCln");
  bcln->Draw();
  inclCln->Draw("same");

  TFile *pythia = new TFile("histos/pPbMC_ppReco_akPu3PF_QCDjetTrig_officialMC_etashift_noTrgSelection_JetCleaned_MCWeightFinal_WithCentVzWeight_4256.root");
  
  double *tmp = new double[bSpec->GetNbinsX()+1];
  for(int i=1; i<=bSpec->GetNbinsX(); i++){
    tmp[i-1] = bSpec->GetBinLowEdge(i);
  }
  tmp[bSpec->GetNbinsX()] = bSpec->GetBinLowEdge(bSpec->GetNbinsX()) + bSpec->GetBinWidth(bSpec->GetNbinsX());

  TH1D *pBFrac = new TH1D("pBFrac","",bSpec->GetNbinsX(),tmp); pBFrac->Sumw2();
  TH1D *inclForDiv = new TH1D("inclForDiv","",bSpec->GetNbinsX(),tmp); inclForDiv->Sumw2();
  TTree *nt = (TTree*)pythia->Get("nt");
  nt->Draw("refpt>>inclForDiv","weight*(abs(jteta)<2)");
  nt->Draw("refpt>>pBFrac","weight*(abs(jteta)<2 && abs(refparton_flavorForB)==5)");

  pBFrac->Divide(inclForDiv);
  bSpec->Divide(InclSpec);

  TCanvas *c3= new TCanvas("c3","",600,600);
  c3->cd();
  c3->SetTopMargin(0.06993);

  pBFrac->SetLineStyle(7);
  pBFrac->SetLineWidth(2);
  pBFrac->SetLineColor(4);
  pBFrac->SetMinimum(0);
  pBFrac->SetMaximum(0.09);
  pBFrac->SetXTitle("b-jet p_{T} [GeV/c]");
  pBFrac->SetYTitle("b-jet Fraction");
  pBFrac->GetXaxis()->SetRangeUser(55,400);
  pBFrac->Draw("hist");
  for(int i=1; i<=pBFrac->GetNbinsX(); i++){
    cout << "bin " << i << " : " << pBFrac->GetBinContent(i) << " " << bSpec->GetBinContent(i) << endl;
    if(bSpec->GetBinContent(i)>0) cout << "err: "<< 1-pBFrac->GetBinContent(i)/bSpec->GetBinContent(i) << endl;
  }

  const int nBins = 16;
  double xbins[nBins], ybins[nBins], xerr[nBins], yerr[nBins];
  double yerrTot[nBins] = {0,0,0,0,0,0,0,0,0.13729, 0.1306, 0.1235, 0.146166, 0.16266, 0.1437,0.1816,0.23813};
  for(int i=1; i<=bSpec->GetNbinsX(); i++){
    xbins[i-1] = bSpec->GetBinLowEdge(i)+bSpec->GetBinWidth(i)/2.;
    ybins[i-1] = bSpec->GetBinContent(i);
    xerr[i-1] = bSpec->GetBinWidth(i)*0.495;
    yerr[i-1] = yerrTot[i-1]*ybins[i-1];
  }
  TGraphErrors *systErr = new TGraphErrors(nBins,xbins,ybins,xerr,yerr);
  systErr->SetFillColor(kRed-7);
  systErr->Draw("2,same");
  bSpec->SetMarkerStyle(20);
  bSpec->SetMarkerColor(1);
  bSpec->Draw("same");
  pBFrac->Draw("hist,same");

  TLegend *l1 = new TLegend(0.42,0.7377,0.889,0.8898);
  l1->SetFillColor(0);
  l1->AddEntry(systErr,"pPb Data, -2 < #eta_{CM} < 2","fp");
  l1->AddEntry(pBFrac,"PYTHIA Z2 + HIJING","l");
  l1->Draw("same");

  //50.85,0.0926233
  TLatex *cmsP = new TLatex(70.9,0.0813,"CMS ");
  cmsP->SetTextFont(62);
  cmsP->SetTextSize(0.0558);
  cmsP->Draw("same");
  TLatex *l4 = new TLatex(258.5,0.0915,"35 nb^{-1} (5.02 TeV)");//; PbPb L = 150 #mub^{-1}");
  l4->SetTextFont(43);
  l4->SetTextSize(25);
  l4->Draw("same");
  
}
