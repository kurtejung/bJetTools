
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

  bool drawPowheg = 0;
  
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
    TH1D *he1 = new TH1D();
    TH1D *he2 = new TH1D();
    TH1D *he3 = new TH1D();
    TH1D *he4 = new TH1D();
    he1->Sumw2(); he2->Sumw2(); he3->Sumw2(); he4->Sumw2();

    //grab unfolded histograms from the unfolding output (inclusive)
    TFile *f0, *f2, *f5,*f6,*f7,*f8;
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
    else if(useUnfolded){
      f0 = new TFile(Form("~/bTagTrees/pPb/raa_pPb_numerator-Bjet_FullDataset_ssvhe20_etaCM_v2_bin0_100_eta%.1fTo%.1f.root",etalo,etahi)); //bjet unfolding
      f5 = new TFile("~/bTagTrees/pPb/raa_pPb_numerator-Bjet_FullDataset_ssvhe20_etaCM_eta-2.0To-1.0.root");
      f6 = new TFile("~/bTagTrees/pPb/raa_pPb_numerator-Bjet_FullDataset_ssvhe20_etaCM_eta-1.0To0.0.root");
      f7 = new TFile("~/bTagTrees/pPb/raa_pPb_numerator-Bjet_FullDataset_ssvhe20_etaCM_eta0.0To1.0.root");
      f8 = new TFile("~/bTagTrees/pPb/raa_pPb_numerator-Bjet_FullDataset_ssvhe20_etaCM_eta1.0To2.0.root");
    }
    else{
      f0 = new TFile(Form("~/bTagTrees/pPb/raa_pPb_numerator-Bjet_NoUnfolding_FullDataset_ssvhe20_etaCM_v2_bin0_100_eta%.1fTo%.1f.root",etalo,etahi));
    }

    if(useUnfolded){
      h = (TH1D*)(f0->Get("hReco0"));
      he1 = (TH1D*)(f5->Get("hReco0"));
      he2 = (TH1D*)(f6->Get("hReco0"));
      he3 = (TH1D*)(f7->Get("hReco0"));
      he4 = (TH1D*)(f8->Get("hReco0"));
    }
    else if(doInclusiveJet) h = (TH1D*)(f0->Get("hIncJetsData"));
    else h = (TH1D*)(f0->Get("hReco0"));

    //h = (TH1D*)he1->Clone("h");
    //he1->Add(he2);
    //he1->Add(he3);
    //he1->Add(he4);
    //h->Scale(1./6.9);

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
    else h_den = (TH1D*)(f2->Get("hGen_cent1"));

    TFile *xchk = new TFile("histos/ppMC_ppReco_akPu3PF_QCDjetTrig_etashift_Fix2Sample_MCWeightFinalWithVz_noTrgSelection_Full.root");
    TTree *nt1 = (TTree*)xchk->Get("nt");
    double xbinsXchk[19] = {0, 5, 10, 15, 20, 25, 30, 35, 40, 55, 70, 90, 110, 140, 170, 220, 400, 700, 1000};
    TH1D *hf = new TH1D("hf","",18,xbinsXchk); hf->Sumw2();
    nt1->Project("hf","refpt","weight*(abs(refparton_flavorForB)==5 && jteta<2 && jteta>-2 && jtpt>35 && rawpt>20)");
    for(int i=1; i<=hf->GetNbinsX(); i++){
      hf->SetBinError(i,hf->GetBinError(i)/hf->GetBinWidth(i));
      hf->SetBinContent(i,hf->GetBinContent(i)/hf->GetBinWidth(i));
    }
    hf->SetLineColor(kRed+2);
    //hf->Scale(1E9);
    //hf->Scale(1./70.);
    //hf->Scale(10.);
    //h_den = hf;

    TFile *ff = new TFile("pPb_Unfo_inCM_v31_officialMC_ak3PF_akPu3PF_noGplus_FirstHalfOnly_Converged_usedParameterizedUnfold0_jtpt35_bJets_clo0_chi100_v8_eta_-2.00To2.00_.root");
    TH1D *hfCheck2 = (TH1D*)ff->Get("hGen_cent1");
    TFile *ff2 = new TFile("pPb_Unfo_inCM_v31_officialMC_Reverse_WithResCorr_ak3PF_akPu3PF_noGplus_SecondHalfOnly_Converged_usedParameterizedUnfold0_jtpt35_bJets_clo0_chi100_v8_eta_-2.00To2.00_.root");
    TH1D *hfCheck3 = (TH1D*)ff2->Get("hGen_cent1");
    hfCheck2->Add(hfCheck3);
    //hfCheck2 = (TH1D*)hfCheck2->Rebin(16,hfCheck2->GetName(),xbinsXchk);
    
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
    h->SetYTitle("#frac{1}{T_{pA}} #frac{d#sigma}{dp_{T}} (1/GeV/c)");
    h->SetMaximum(1E-6);
    h->SetMinimum(1E-11);
    h->Draw();
    h_den->SetMarkerColor(4);
    h_den->Draw("same");
    hfCheck2->Scale(1./70.);
    hfCheck2->SetLineColor(kGreen+2);
    hf->Draw("same");
    he2->SetMarkerColor(kOrange+2);
    he2->Draw("same");
    //hfCheck2->Draw("same");
    //h_den = hfCheck2;
    //c1->SetLogy();

    for(int j=1; j<=h->GetNbinsX(); j++){
      cout << h->GetBinLowEdge(j) << " ";
    }
    cout << h->GetBinLowEdge(h->GetNbinsX())+h->GetBinWidth(h->GetNbinsX()) << endl;

    for(int j=1; j<=h_den->GetNbinsX(); j++){
      cout << h_den->GetBinLowEdge(j) << " ";
    }
    cout << h_den->GetBinLowEdge(h_den->GetNbinsX())+h_den->GetBinWidth(h_den->GetNbinsX()) << endl;

    TH1F *hpow, *hpow15;
    //if(drawPowheg){
      TFile *fpowheg50 = new TFile("Powheg/CMS-Dijet-Analysis-5tev-kt0-pt50Cut-mstwCL60PDF.root");
      hpow = (TH1F*)fpowheg50->Get("d05-x01-y01")->Clone("hpow"); //inclusive eta, ak3 jets
      TFile *fpowheg15 = new TFile("Powheg/CMS-Dijet-Analysis-5tev-kt0-pt20Cut-mstwCL60PDF.root");
      hpow15 = (TH1F*)fpowheg15->Get("d05-x01-y01")->Clone("hpow15");
      hpow15->Scale(25.4); //pthat20 relative (to pthat50) efficiency
      //hpow15->Scale(405.); //pthat 5 relative efficiency
      //hpow15->Scale(15.95);
      for(int i=1; i<=hpow->GetNbinsX(); i++){
	if(hpow->GetBinLowEdge(i)<100){
	  //cout << "edge: "<< hpow->GetBinLowEdge(i) << endl;
	  //cout << "hpow: "<< hpow->GetBinContent(i) << " hpow15: "<< hpow15->GetBinContent(i) << endl;
	  hpow->SetBinContent(i, hpow15->GetBinContent(i));
	  hpow->SetBinError(i, hpow15->GetBinError(i));
	}
      }
      hpow->SetMarkerStyle(20);
      hpow->SetMarkerColor(kOrange+2);
      hpow->SetLineColor(kOrange+2);
      hpow->Scale(2./10.); //fix eta selection
      hpow->Scale(1./70E-3); //correct for sigma_pp
      hpow->Scale(1.851e-3); //correct pthat50 efficiency (148086/80M)
      //hpow->Scale(4.7e-2);
      hpow->Scale(1e-9); //correct ncoll & translate to mb, not pb
      //hpow->Scale(1./50.);
      //hpow->Scale(0.85/(30.));
      if(drawPowheg) hpow->Draw("same");
      // }

    TLegend *t1 = new TLegend(0.5,0.697,0.879,0.906);
    t1->AddEntry(h,"pPb data b-jet spec.","p");
    t1->AddEntry(h2,"pPb MC b-jet spec.","p");
    t1->AddEntry(h_den,"pp MC b-jet spec.","p");
    if(drawPowheg) t1->AddEntry(hpow,"Powheg Spec.","p");
    t1->SetFillColor(0);
    t1->Draw("same");

    TCanvas *c2 = new TCanvas("c2","",600,600); //800,800
    //formatCanvas(c2);
    TH1D *hcln = (TH1D*)h->Clone("hcln");
    cout << "nbin: "<< hcln->GetNbinsX() << endl;
    //double xbinsrebin[17] = {0,5,10,15,20,25,30,35,40,55,70,90,110,140,170,220,400};
    //TH1D *hcln = new TH1D("hcln","",16,xbinsrebin);
    for(int ibin=1; ibin<=hcln->GetNbinsX(); ibin++){
      cout << "content bin " << ibin << " " << hcln->GetBinContent(ibin) << endl;
      cout << "powheg content bin " << ibin << " " << hpow->GetBinContent(ibin) << endl;
      //hcln->SetBinError(ibin,hcln_tmp->GetBinError(ibin));
      //hcln->SetBinContent(ibin,hcln_tmp->GetBinContent(ibin));
      }
    TH1D *hdencln = (TH1D*)h_den->Clone("hdencln");
    TH1D *hclnpow = (TH1D*)hcln->Clone("hclnpow");
    //cout << "WARNING! Scaling denominator for neutrinoless definition change!" << endl;
    //hdencln->Scale(1./1.243);
    if(!drawPowheg) hcln->Divide(hdencln);
    else{
      for(int ibin=1; ibin<=16; ibin++){
	int powbin=ibin;
	if(!useUnfolded) powbin+=8;
	if(!useUnfolded && ibin>8) continue;
	hcln->SetBinError(powbin,sqrt(pow(h_den->GetBinError(ibin)/h_den->GetBinContent(ibin),2)+pow(hpow->GetBinError(powbin)/hpow->GetBinContent(powbin),2)));
	hcln->SetBinContent(powbin,h_den->GetBinContent(ibin)/hpow->GetBinContent(powbin));
      }
    }
    //hcln->Divide(hcln);
    
    //Apply systematic errors to the inclusive eta plot
    TGraphErrors *systErr[7];
    TGraphErrors *unfoldErr[7];

    //set systematics
    const int nBins = 18;
    double xp[nBins], yp[nBins], xerr[nBins], yerr[nBins];
    //double yerrTot[nBins] = {0,0,0,0,0,0,0,0,0.13729, 0.1306, 0.1235, 0.146166, 0.16266, 0.1437,0.1816,0.23813}; //before recalibration
    double yerrTot[nBins] = {0,0,0,0,0,0,0,0,0,0.131694,0.0976, 0.0957, 0.1139, 0.1209, 0.1773, 0.2162,0,0}; //after recalibration
    double unfoldErrTot[nBins] = {0,0,0,0,0,0,0,0,0,0.05,0.05,0.05,0.05,0.05,0.05,0.06,0,0};
    double unfoldErrTotEta1[nBins] = {0,0,0,0,0,0,0,0,0,0.15,0.08,0.06,0.06,0.08,0.04,0.04,0,0};
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
    hcln->SetMarkerStyle(20);
    hcln->SetMarkerColor(1);
    hcln->SetMarkerSize(1.2);
    for(int i=1; i<hcln->GetNbinsX(); i++){
      if(hcln->GetBinLowEdge(i) < 55){
	hcln->SetBinContent(i,-1);
	hcln->SetBinError(i,0.0001);
      }
    }
    hcln->Draw();
    hclnpow->SetMarkerColor(kOrange+2);
    hclnpow->SetLineColor(kOrange+2);
    if(drawPowheg) hclnpow->Draw("same");
    TLine *at2 = new TLine(0,1,400,1);
    at2->SetLineStyle(7);
    at2->SetLineWidth(1);
    pythiaErr->SetFillColor(kRed-7);
    //pythiaErr->SetFillStyle(1); //3844
    if(!doInclusiveJet){
      for(int i=0; i<7; i++){
	unfoldErr[i]->SetFillColor(kCyan-7);
	systErr[i]->SetFillColor(kYellow-7);
	systErr[i]->SetMarkerStyle(20);
	systErr[i]->Draw("2,same");
	//unfoldErr[i]->Draw("2,same");
      }
    }
    pythiaErr->Draw("2,same");
    double lumix[1] = {10}; double lumiy[1] = {1}; double lumierrx[1] = {10}; double lumierry[1] = {0.035};
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

  TLegend *leg1 = new TLegend(0.41,0.729,0.88,0.886);
  leg1->SetFillColor(0);
  //leg1->AddEntry(inclJetSys,"Inclusive-Jet R_{pA}, |#eta_{CM}|<0.5","fp");
  //leg1->AddEntry(raaErr,"b-jet R_{AA}, (0-100%), |#eta|<2","fp");
  leg1->AddEntry(systErr[0],"pPb Data, -2 < #eta_{CM} < 2","fp");
  leg1->AddEntry(grepPbtheory,"Ref. [32]","lp");
  //leg1->AddEntry(pythiaErr,"Reference Unc.","f");
  //leg1->AddEntry(lumiErr,"Luminosity Unc.","f");
  //leg1->AddEntry(systErr[0],"Total Systematic Error","f");
  //leg1->AddEntry(unfoldErr[0],"Unfolding Uncertainty","f");
  leg1->Draw("same");

  TLegend *leg3 = new TLegend(0.2,0.207,0.567,0.328);
  leg3->SetFillColor(0);
  leg3->AddEntry(lumiErr,"pPb Luminosity Unc.","f");
  leg3->AddEntry(pythiaErr,"pPb Reference Unc.","f");
  leg3->Draw("same");

  grepPbtheory->Draw("C,same");

  c2->SetTopMargin(0.0734);
  c2->RedrawAxis();

  TLatex *cmsP = new TLatex(16.7,2.25,"CMS ");
  cmsP->SetTextFont(62);
  cmsP->SetTextSize(0.0558);
  cmsP->Draw("same");
  TLatex *l1 = new TLatex(-0.5,0.18,Form("0-100%%",etalo,etahi));
  //l1->Draw("same");
  TLatex *l4 = new TLatex(236.8,2.56,"35 nb^{-1} (5.02 TeV)");//; PbPb L = 150 #mub^{-1}");
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
