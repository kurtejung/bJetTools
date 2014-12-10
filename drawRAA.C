
void formatCanvas(TCanvas *c){
  c->Divide(1,2,0.01,0.01);
  c->cd(1);
  //c->GetPad(1)->SetLogy();
  c->GetPad(1)->SetPad(0.,0.225,1.,1.);
  c->GetPad(2)->SetPad(0.,0.0,1.,0.3);
  c->GetPad(2)->SetBottomMargin(0.3);
  //c->GetPad(2)->SetGridy(1);
}


void drawRAA(bool doInclusiveJet=0, bool useUnfolded=1, bool drawTheory=1){

  float etalo=-2;
  float etahi=2;

  TH1D *dEta80To110 = new TH1D("dEta80To110","",4,-2.465,1.535);
  TH1D *dEta110To150 = new TH1D("dEta110To150","",4,-2.365,1.635);
  TH1D *dEta70To80 = new TH1D("dEta70To80","",4,-2.565,1.435);
  TH1D *dEta150To190 = new TH1D("dEta150To190","",4,-2.265,1.735);

  int nbins = 3;
  double xbins[4] = {0,20,50,100};
  double xbins2[4] = {2,22,52,102};
  double xbins3[4] = {-2,18,48,98};
  double xbins4[4] = {4,24,54,104};
  TH1D *dCent70To80 = new TH1D("dCent70To80","",nbins,xbins3);
  TH1D *dCent80To110 = new TH1D("dCent80To110","",nbins,xbins);
  TH1D *dCent110To150 = new TH1D("dCent110To150","",nbins,xbins2);
  TH1D *dCent150To190 = new TH1D("dCent150To190","",nbins,xbins4);
  
  //for(int ieta=-2; ieta<2; ieta++){
  //  etalo=(float)(ieta+0.5);
  //  etahi=(float)(ieta+1.5);

  int binLo, binHi;
  //for(int icent=0; icent<3; icent++){
  //  if(icent==0){ binLo = 0; binHi=20; }
  //  if(icent==1){ binLo = 20; binHi=50; }
  //  if(icent==2){ binLo = 50; binHi=100; }
  
    TH1D *h = new TH1D();
    TH1D *h2 = new TH1D();
    TH1D *h_den = new TH1D();
    h->Sumw2();
    h2->Sumw2();
    h_den->Sumw2();

    //grab unfolded histograms from the unfolding output (inclusive)
    TFile *f0, *f2;
    TFile *fIncl = new TFile("~/bTagTrees/pPb/raa_pPb_numerator-inc_eta-1.0To1.0.root");
    TFile *fInclB = new TFile("~/bTagTrees/pPb/raa_pp_denomForpA-inc_eta-1.0To1.0.root");
    TH1D *hinclTop = (TH1D*)fIncl->Get("hReco0")->Clone("hinclTop");
    TH1D *hinclBot = (TH1D*)fInclB->Get("hGen_cent1")->Clone("hinclBot");
    for(int i=1; i<=hinclTop->GetNbinsX(); i++){
      hinclTop->SetBinError(i,hinclTop->GetBinError(i)/hinclTop->GetBinWidth(i));
      hinclTop->SetBinContent(i,hinclTop->GetBinContent(i)/hinclTop->GetBinWidth(i));
      hinclBot->SetBinError(i,hinclBot->GetBinError(i)/hinclBot->GetBinWidth(i));
      hinclBot->SetBinContent(i,hinclBot->GetBinContent(i)/hinclBot->GetBinWidth(i));
    }
    hinclTop->Scale(1./6.9);
    hinclTop->Divide(hinclBot);
    hinclTop->SetMarkerColor(4);
    //hinclTop->Scale(0.689);

    //grab unfolded histograms from the unfolding output (bjet)
    if(doInclusiveJet) f0 = new TFile(Form("~/bTagTrees/pPb/raa_pPb_numerator-inc_eta%.1fTo%.1f.root",etalo,etahi));
    //else f0 = new TFile(Form("~/bTagTrees/pPb/raa_pPb_numerator-Incljet_FullDataset_ssvhe20_etaCM_v2_bin0_100_eta%.1fTo%.1f.root",etalo,etahi)); //incl jet unfolding
    else f0 = new TFile(Form("~/bTagTrees/pPb/raa_pPb_numerator-Bjet_FullDataset_ssvhe20_etaCM_v2_bin0_100_eta%.1fTo%.1f.root",etalo,etahi)); //bjet unfolding
    

    if(useUnfolded) h = (TH1D*)(f0->Get("hReco0"));
    else if(doInclusiveJet) h = (TH1D*)(f0->Get("hIncJetsData"));
    else h = (TH1D*)(f0->Get("hRawBData"));

    /* for(int i=1; i<=h->GetNbinsX(); i++){
       h->SetBinError(i, h->GetBinError(i)/h->GetBinWidth(i));
       h->SetBinContent(i, h->GetBinContent(i)/h->GetBinWidth(i));
       }*/
    //h->Scale(1./6.9); //pPb total N coll (L*pp x-sec * A, sigma_pp = 70 mb @ CERN)  
  
    if(doInclusiveJet) f2 = new TFile(Form("~/bTagTrees/pPb/raa_pp_denomForpA-inc_eta%.1fTo%.1f.root",etalo,etahi),"OLD");
    else f2 = new TFile(Form("~/bTagTrees/pPb/raa_pp_denomForpA-Bjet_etaCM_v2_bin0_100_eta%.1fTo%.1f.root",etalo,etahi),"OLD"); //incl jet unfolding
    //else f2 = new TFile(Form("raa_pp_denomForpA-Bjet_etaLab_eta%.1fTo%.1f.root",etalo,etahi),"OLD"); //bjet unfolding
    f2->cd();
    //h_den = (TH1D*)(f2->Get("hRawBMC"));
    if(useUnfolded) h_den = (TH1D*)f2->Get("hGen_cent1");
    else if(doInclusiveJet) h_den = (TH1D*)(f2->Get("hIncJetsMC"));
    else h_den = (TH1D*)(f2->Get("hRawBMC"));

    //Fix raw B MC histogram to have 14 bins to be compatible with the unfolded nonsense
    /*const int fixNbin = h->GetNbinsX();
    Double_t binBoundaries[15];
    for(int ibin=0; ibin<fixNbin; ibin++){
      binBoundaries[ibin] = h->GetBinLowEdge(ibin+1);
    }
    binBoundaries[fixNbin] = h->GetBinLowEdge(fixNbin)+h->GetBinWidth(fixNbin);
    TH1D *h_den_fix = new TH1D("h_den_fix","",fixNbin,binBoundaries); h_den_fix->Sumw2();
    for(int ibin=0; ibin<fixNbin-h_den->GetNbinsX(); ibin++){
      h_den_fix->SetBinContent(ibin+1,0);
      h_den_fix->SetBinError(ibin+1,0);
    }
    for(int ibin=fixNbin-h_den->GetNbinsX(); ibin<fixNbin; ibin++){
      h_den_fix->SetBinContent(ibin+1,h_den->GetBinContent(ibin-(fixNbin-h_den->GetNbinsX())+1));
      h_den_fix->SetBinError(ibin+1,h_den->GetBinError(ibin-(fixNbin-h_den->GetNbinsX())+1));
      }*/
    //end fix

    //h->Draw();
    h->SetTitle("");
    //h_den->SetMarkerColor(2);
    //h_den->Draw("same");
    //h->Divide(h_den);
    h->SetYTitle("#frac{1}{N} #frac{dN}{dp_{T}} (1/GeV/c)");
    h->SetMaximum(1E-6);
    h->SetMinimum(1E-11);
    h->Draw();
    h_den->SetMarkerColor(4);
    h_den->Draw("same");
    //c1->SetLogy();

    TLegend *t1 = new TLegend(0.5,0.697,0.879,0.906);
    t1->AddEntry(h,"pPb data b-jet spec.","p");
    t1->AddEntry(h2,"pPb MC b-jet spec.","p");
    t1->AddEntry(h_den,"pp MC b-jet spec.","p");
    t1->SetFillColor(0);
    t1->Draw("same");

    TCanvas *c2 = new TCanvas("c2","",600,600); //800,800
    //formatCanvas(c2);
    TH1D *hcln = (TH1D*)h->Clone("h_cln");
    TH1D *hdencln = (TH1D*)h_den->Clone("hdencln");
    cout << "WARNING! Scaling denominator for neutrinoless definition change!" << endl;
    hdencln->Scale(1./1.243);
    hcln->Divide(hdencln);
    //hcln->Divide(hcln);
  
    //Apply systematic errors to the inclusive eta plot
    TGraphErrors *systErr[7];
    TGraphErrors *unfoldErr[7];
    cout << "nbin: "<< hcln->GetNbinsX() << endl;
    //set systematics
    const int nBins = 16;
    double xp[nBins], yp[nBins], xerr[nBins], yerr[nBins];
    //double yerrTot[nBins] = {0,0,0,0,0,0,0,0,0.13729, 0.1306, 0.1235, 0.146166, 0.16266, 0.1437,0.1816,0.23813}; //before recalibration
    double yerrTot[nBins] = {0,0,0,0,0,0,0,0,0,0.15092,0.1084, 0.11274, 0.1450, 0.1429, 0.17421, 0.2114}; //after recalibration
    double unfoldErrTot[nBins] = {0,0,0,0,0,0,0,0,0,0.10,0.10,0.06,0.06,0.05,0.05,0.06};
    double unfoldErrTotEta1[nBins] = {0,0,0,0,0,0,0,0,0,0.15,0.08,0.06,0.06,0.08,0.04,0.04};
    //double pythiaSysErr[nBins+1] = {0,0,0,0,0,0,0,0.136,0.138,0.140,0.142,0.148,0.152,0.157,0.164,0.174,0.182}; //option B
    // double optionBExtra[nBins+1] = {0,0,0,0,0,0,0,0,0,0.10,0.111,0.078,0.119,0.108,0.053,0.122,0.168};
    double pythiaSysErr[1] = {0.22}; //option A
    double pythiaYpoints[1] = {1};
    double pythiaXpoints[1] = {10};
    double pythiaXerr[1] = {10};
    for(int i=0; i<nBins; i++){ yerrTot[i] = sqrt(pow(yerrTot[i],2)+pow(unfoldErrTot[i],2)); }
    int j=0;
    if(!doInclusiveJet){
      for(int i=1; i<=hcln->GetNbinsX(); i++){
	xp[i-1] = hcln->GetBinLowEdge(i)+hcln->GetBinWidth(i)/2.;
	yp[i-1] = hcln->GetBinContent(i);
	xerr[i-1] = hcln->GetBinWidth(i)*0.495;
	yerr[i-1] = yerrTot[i-1];
	//pythiaXpoints[i-1] = xp[i-1];
	//pythiaXerr[i-1] = xerr[i-1];

	if(yerr[i-1]>0.001 && yp[i-1]>-1){ 
	  systErr[j] = new TGraphErrors(1,&xp[i-1],&yp[i-1],&xerr[i-1],&yerrTot[i-1]);
	  unfoldErr[j] = new TGraphErrors(1,&xp[i-1],&yp[i-1],&xerr[i-1],&unfoldErrTot[i-1]);
	  systErr[j]->SetName(Form("RpA_SystErr_bin%d",j));
	  j++;
	}
      }
      //pythiaXerr[16] = 0.;
      //pythiaXpoints[16] = 400.;
      /*for(int jj=0; jj<nBins+1; jj++){
	pythiaSysErr[jj]*=(1+optionBExtra[jj]);
	}*/
      TGraphErrors *pythiaErr = new TGraphErrors(1,pythiaXpoints,pythiaYpoints,pythiaXerr,pythiaSysErr);
    }
    //TGraphErrors *systErr = new TGraphErrors(nBins,xp,yp,xerr,yerr);
    double jesx[1]={2.365};
    double jesy[1]={1};
    double jesxerr[1]={.1};
    double jesyerr[1]={0.035};
    TGraphErrors *JESNorm = new TGraphErrors(1,jesx,jesy,jesxerr,jesyerr);
    JESNorm->SetFillColor(kGreen+3);

    hcln->SetMaximum(2.5);
    hcln->SetMinimum(0);
    if(doInclusiveJet) hcln->SetYTitle("Inclusive Jet R_{pA}");
    else hcln->SetYTitle("Nuclear Modification Factor");
    //else hcln->SetYTitle("b-jet R_{pA}^{PYTHIA}");
    hcln->SetXTitle("b-jet p_{T} [GeV/c]");
    hcln->GetXaxis()->SetRangeUser(0,400);
    hcln->GetXaxis()->CenterTitle(0);
    hcln->GetYaxis()->CenterTitle(0);
    hcln->SetMarkerStyle(25);
    hcln->SetMarkerSize(1.2);
    for(int i=1; i<hcln->GetNbinsX(); i++){
      if(hcln->GetBinLowEdge(i) < 55){
	hcln->SetBinContent(i,-1);
	hcln->SetBinError(i,0.0001);
      }
    }
    hcln->Draw();
    TLine *at2 = new TLine(0,1,400,1);
    at2->SetLineStyle(7);
    at2->SetLineWidth(1);
    pythiaErr->SetFillColor(kRed-7);
    //pythiaErr->SetFillStyle(1); //3844
    if(!doInclusiveJet){
      for(int i=0; i<7; i++){
	unfoldErr[i]->SetFillColor(kCyan-7);
	systErr[i]->SetFillColor(kYellow-7);
	systErr[i]->SetMarkerStyle(25);
	systErr[i]->Draw("2,same");
	//unfoldErr[i]->Draw("2,same");
      }
    }
    pythiaErr->Draw("2,same");
    double lumix[1] = {30}; double lumiy[1] = {1}; double lumierrx[1] = {10}; double lumierry[1] = {0.035};
    TGraphErrors *lumiErr = new TGraphErrors(1,lumix,lumiy,lumierrx,lumierry);
    lumiErr->SetFillColor(kGreen+3);
    lumiErr->Draw("2,same");
    at2->Draw("same");
    //JESNorm->Draw("2,same");
    hcln->Draw("same");
    
    //JESNorm->Draw("2,same");
  TLine *phigh1 = new TLine(380,1.22,400,1.22);
  phigh1->SetLineColor(2);
  phigh1->SetLineStyle(4);
  phigh1->SetLineWidth(1);
  //phigh1->Draw("same");
  TLine *phigh2 = new TLine(380,0.78,400,0.78);
  phigh2->SetLineColor(2);
  phigh2->SetLineStyle(4);
  phigh2->SetLineWidth(1);
  //phigh2->Draw("same");
  double zeroTemp[1] = {0};
  TGraphErrors *inclJetSys = new TGraphErrors(1,zeroTemp,zeroTemp,zeroTemp,zeroTemp);
  inclJetSys->SetMarkerColor(2);
  inclJetSys->SetMarkerStyle(21);
  inclJetSys->SetFillColor(8);
  //gROOT->ProcessLine(".x ~/Downloads/UnfoldedakPu3PFIncJetRpAvsPtInterpolatedRefEtaBin7.C");
  //gROOT->ProcessLine(".x ~/Downloads/RAA_drawFinal.C");
  if(drawTheory){
    gROOT->ProcessLine(".x ~/Downloads/drawpPbTheory.C");
  }
  TGraph *raaErr = new TGraph();
  raaErr->SetFillColor(kMagenta+1);
  raaErr->SetMarkerStyle(20);
  raaErr->SetMarkerSize(1.2);

  TLegend *leg1 = new TLegend(0.5,0.729,0.88,0.886);
  leg1->SetFillColor(0);
  //leg1->AddEntry(inclJetSys,"Inclusive-Jet R_{pA}, |#eta_{CM}|<0.5","fp");
  //leg1->AddEntry(raaErr,"b-jet R_{AA}, (0-100%), |#eta|<2","fp");
  leg1->AddEntry(systErr[0],"b-jet R_{pA}^{PYTHIA}, -2.4<#eta_{CM}<1.6","fp");
  leg1->AddEntry(grepPbtheory,"Vitev, et.al. arXiv:1306.0909","lp");
  //leg1->AddEntry(pythiaErr,"Reference Unc.","f");
  //leg1->AddEntry(lumiErr,"Luminosity Unc.","f");
  //leg1->AddEntry(systErr[0],"Total Systematic Error","f");
  //leg1->AddEntry(unfoldErr[0],"Unfolding Uncertainty","f");
  leg1->Draw("same");

  TLegend *leg3 = new TLegend(0.2,0.78,0.472,0.878);
  leg3->SetFillColor(0);
  leg3->AddEntry(lumiErr,"pPb Luminosity Unc.","f");
  leg3->AddEntry(pythiaErr,"pPb Reference Unc.","f");
  leg3->Draw("same");

  grepPbtheory->Draw("C,same");

  c2->SetTopMargin(0.0734);
  c2->RedrawAxis();

  TLatex *cmsP = new TLatex(0.5,2.56,"CMS Preliminary");
  cmsP->SetTextFont(43);
  cmsP->SetTextSize(25);
  cmsP->Draw("same");
  TLatex *l1 = new TLatex(-0.5,0.18,Form("0-100%%",etalo,etahi));
  //l1->Draw("same");
  TLatex *l4 = new TLatex(160,2.56,"pPb L = 35 nb^{-1}");//; PbPb L = 150 #mub^{-1}");
  l4->SetTextFont(43);
  l4->SetTextSize(25);
  l4->Draw("same");
  TLatex *l2 = new TLatex(20.77,0.175,"pPb #sqrt{s_{NN}} = 5.02 TeV #int L dt = 35 nb^{-1}"); //35, 20.7, 14
  l2->SetTextSize(0.035);
  //l2->Draw("same");
  TLatex *l3 = new TLatex(-2,0.22,"PbPb #sqrt{s_{NN}} = 2.76 TeV #int L dt = 150 #mub^{-1}");
  l3->SetTextSize(0.035);
  //l3->Draw("same");
  //TLatex *ltag = new TLatex(65,2.2,"CSV > 0.679");
  TLatex *ltag = new TLatex(65,2.2,"SSVHE > 2.0");
  //ltag->Draw("same");
  TFile *out = new TFile("RpA_BJet_Output.root","RECREATE");
  out->cd();
  hcln->SetName("RpA");//,"b jet RpA");
  hcln->Write();
  h->Write();
  h_den->Write();
  for(int i=0; i<7; i++){
    systErr[i]->Write();
  }
  out->Close();

  // hcln->Fit("pol0","","",45,250);
  //hcln->Fit("pol1","","",45,250);

  cout << "values: { ";
  for(int i=1; i<=hcln->GetNbinsX(); i++){
    cout << hcln->GetBinContent(i) << " ";
  }
  cout << " }" << endl;

  cout << "stat err: { ";
  for(int i=1; i<=hcln->GetNbinsX(); i++){
    cout << hcln->GetBinError(i) << " ";
  }
  cout << " }" << endl;

  cout << "syst err: { ";
  for(int i=1; i<=hcln->GetNbinsX(); i++){
    cout << yerrTot[i-1]+(yerrTot[i-1]*hcln->GetBinContent(i)) << " ";
  }
  cout << " }" << endl;
}
