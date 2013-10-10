void evalTagEff(int RecOrGen=0, float discCut=0.6)
{
  
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  
  TFile *fB = new TFile("histos/ppMC_ppReco_ak3PF_BjetTrig_noIPupperCut_9-15.root");

  TNtuple *ntB = (TNtuple *) fB->Get("nt");

  TH1F *hPt1B=new TH1F("hPt1B","hPt1B",54,30.,300.);
  TH1F *hPt1BTag=new TH1F("hPt1BTag","hPt1BTag",54,30.,300.);
  TH1F *hPt2B=new TH1F("hPt2B","hPt2B",54,30.,300.);
  TH1F *hPt2BTag=new TH1F("hPt2BTag","hPt2BTag",54,30.,300.);
  TH1F *hPt3B=new TH1F("hPt3B","hPt3B",54,30.,300.);
  TH1F *hPt3BTag=new TH1F("hPt3BTag","hPt3BTag",54,30.,300.);

  TH1F *hPt2SubB=new TH1F("hPt2SubB","hPt2SubB",54,30.,300.);
  TH1F *hPt2SubBTag=new TH1F("hPt2SubBTag","hPt2SubBTag",54,30.,300.);

  TH1F *hPt1IncB=new TH1F("hPt1IncB","hPt1IncB",54,30.,300.);
  TH1F *hPt1IncBTag=new TH1F("hPt1IncBTag","hPt1IncBTag",54,30.,300.);
  TH1F *hPt2IncB=new TH1F("hPt2IncB","hPt2IncB",54,30.,300.);
  TH1F *hPt2IncBTag=new TH1F("hPt2IncBTag","hPt2IncBTag",54,30.,300.);
  TH1F *hPt3IncB=new TH1F("hPt3IncB","hPt3IncB",54,30.,300.);
  TH1F *hPt3IncBTag=new TH1F("hPt3IncBTag","hPt3IncBTag",54,30.,300.);

  hPt1B->Sumw2();
  hPt2B->Sumw2();
  hPt3B->Sumw2();
  hPt2SubB->Sumw2();
  hPt1IncB->Sumw2();
  hPt2IncB->Sumw2();
  hPt3IncB->Sumw2();

  hPt1BTag->Sumw2();
  hPt2BTag->Sumw2();
  hPt3BTag->Sumw2();
  hPt2SubBTag->Sumw2();
  hPt1IncBTag->Sumw2();
  hPt2IncBTag->Sumw2();
  hPt3IncBTag->Sumw2();
  
  if(RecOrGen){
    ntB->Draw("refpt1>>hPt1B","weight*(jtpt1>80&&abs(refparton_flavorForB1)==5)");
    ntB->Draw("refpt1>>hPt1BTag",Form("weight*(jtpt1>80&&abs(refparton_flavorForB1)==5&&discr_prob1>=%f)",discCut));
    ntB->Draw("refpt2>>hPt2B","weight*(jtpt1>80&&abs(refparton_flavorForB2)==5)");
    ntB->Draw("refpt2>>hPt2BTag",Form("weight*(jtpt1>80&&abs(refparton_flavorForB2)==5&&discr_prob2>=%f)",discCut));
    ntB->Draw("refpt3>>hPt3B","weight*(jtpt1>80&&abs(refparton_flavorForB3)==5)");
    ntB->Draw("refpt3>>hPt3BTag",Form("weight*(jtpt1>80&&abs(refparton_flavorForB3)==5&&discr_prob3>=%f)",discCut));
    ntB->Draw("refpt2>>hPt2SubB","weight*(abs(refparton_flavorForB2)==5&&jtpt1>80&&jtpt2>30&&acos(cos(jtphi1-jtphi2))>2./3.*acos(-1.))");
    ntB->Draw("refpt2>>hPt2SubBTag",Form("weight*(abs(refparton_flavorForB2)==5&&jtpt1>80&&jtpt2>30&&acos(cos(jtphi1-jtphi2))>2./3.*acos(-1.)&&discr_prob2>=%f)",discCut));

    ntB->Draw("refpt1>>hPt1IncB","weight*(abs(refparton_flavorForB1)==5)");
    ntB->Draw("refpt1>>hPt1IncBTag",Form("weight*(abs(refparton_flavorForB1)==5&&discr_prob1>=%f)",discCut));
    ntB->Draw("refpt2>>hPt2IncB","weight*(abs(refparton_flavorForB2)==5)");
    ntB->Draw("refpt2>>hPt2IncBTag",Form("weight*(abs(refparton_flavorForB2)==5&&discr_prob2>=%f)",discCut));
    ntB->Draw("refpt3>>hPt3IncB","weight*(abs(refparton_flavorForB3)==5)");
    ntB->Draw("refpt3>>hPt3IncBTag",Form("weight*(abs(refparton_flavorForB3)==5&&discr_prob3>=%f)",discCut));

  }
  else{
    ntB->Draw("jtpt1>>hPt1B","weight*(jtpt1>80&&abs(refparton_flavorForB1)==5)");
    ntB->Draw("jtpt1>>hPt1BTag",Form("weight*(jtpt1>80&&abs(refparton_flavorForB1)==5&&discr_prob1>=%f)",discCut));
    ntB->Draw("jtpt2>>hPt2B","weight*(jtpt1>80&&abs(refparton_flavorForB2)==5)");
    ntB->Draw("jtpt2>>hPt2BTag",Form("weight*(jtpt1>80&&abs(refparton_flavorForB2)==5&&discr_prob2>=%f)",discCut));
    ntB->Draw("jtpt3>>hPt3B","weight*(jtpt1>80&&abs(refparton_flavorForB3)==5)");
    ntB->Draw("jtpt3>>hPt3BTag",Form("weight*(jtpt1>80&&abs(refparton_flavorForB3)==5&&discr_prob3>=%f)",discCut));
    ntB->Draw("jtpt2>>hPt2SubB","weight*(jtpt1>80&&abs(refparton_flavorForB2)==5&&jtpt1>80&&jtpt2>30&&acos(cos(jtphi1-jtphi2))>2./3.*acos(-1.))");
    ntB->Draw("jtpt2>>hPt2SubBTag",Form("weight*(jtpt1>80&&abs(refparton_flavorForB2)==5&&jtpt1>80&&jtpt2>30&&acos(cos(jtphi1-jtphi2))>2./3.*acos(-1.)&&discr_prob2>=%f)",discCut));
 
    ntB->Draw("jtpt1>>hPt1IncB","weight*(abs(refparton_flavorForB1)==5)");
    ntB->Draw("jtpt1>>hPt1IncBTag",Form("weight*(abs(refparton_flavorForB1)==5&&discr_prob1>=%f)",discCut));
    ntB->Draw("jtpt2>>hPt2IncB","weight*(abs(refparton_flavorForB2)==5)");
    ntB->Draw("jtpt2>>hPt2IncBTag",Form("weight*(abs(refparton_flavorForB2)==5&&discr_prob2>=%f)",discCut));
    ntB->Draw("jtpt3>>hPt3IncB","weight*(abs(refparton_flavorForB3)==5)");
    ntB->Draw("jtpt3>>hPt3IncBTag",Form("weight*(abs(refparton_flavorForB3)==5&&discr_prob3>=%f)",discCut));

 }
  TH1F *hNum=(TH1F*)hPt1BTag->Clone("hNum");
  hNum->Add(hPt2BTag);
  hNum->Add(hPt3BTag);

  TH1F *hDen=(TH1F*)hPt1B->Clone("hDen");
  hDen->Add(hPt2B);
  hDen->Add(hPt3B);

 

  TCanvas *c1=new TCanvas("c1","c1",600,600);

  TH1F *hEff = (TH1F*)hNum->Clone("hEff");
  hEff->Reset();
  hEff->Divide(hNum,hDen,1,1,"B");
  hEff->SetMinimum(0.);
  hEff->SetMarkerStyle(8);
  if(RecOrGen)hEff->SetXTitle("Gen jet p_{T} (GeV/c)");
  else hEff->SetXTitle("Reco jet p_{T} (GeV/c)");
  hEff->SetYTitle("Efficiency");
  hEff->Draw();


 

  //TF1 *fitf = new TF1("fitf","[0]+[1]*log(x*[2])+[3]*log(x*[4])*log(x*[4])",30,300);
  //fitf->SetParameters(0.005,0.1,0.5,0.1,0.5);
  if(!RecOrGen){
    TF1 *fitf = new TF1("fitf","[0]+[1]*pow(log(x*[2]),2)",30,300);
    fitf->SetParameters(1.,0.1,0.008);
    fitf->SetLineColor(1);
    hEff->Fit(fitf,"","",30,300);
  }

 TH1F *hNumInc=(TH1F*)hPt1IncBTag->Clone("hNumInc");
  hNumInc->Add(hPt2IncBTag);
  hNumInc->Add(hPt3IncBTag);

  TH1F *hDenInc=(TH1F*)hPt1IncB->Clone("hDenInc");
  hDenInc->Add(hPt2IncB);
  hDenInc->Add(hPt3IncB);

  TH1F *hEffInc = (TH1F*)hNumInc->Clone("hEffInc");
  hEffInc->Reset();
  hEffInc->Divide(hNumInc,hDenInc,1,1,"B");
  hEffInc->SetMinimum(0.);
  hEffInc->SetMarkerStyle(4);
  hEffInc->SetMarkerColor(kGreen-2);
  hEffInc->SetLineColor(kGreen-2);
  hEffInc->Draw("same");

 if(!RecOrGen){
    TF1 *fitf2 = new TF1("fitf2","[0]+[1]*pow(log(x*[2]),2)",30,300);
    fitf2->SetParameters(1.,0.1,0.1);
    hEffInc->Fit(fitf2,"N","",30,300);
    fitf2->SetLineWidth(2);
    fitf2->SetLineColor(kGreen-2);
    fitf2->Draw("same");
  }


  TH1F *hNumSub=(TH1F*)hPt2SubBTag->Clone("hNumSub");
  TH1F *hDenSub=(TH1F*)hPt2SubB->Clone("hDenSub");

  TH1F *hEffSub = (TH1F*)hNumSub->Clone("hEffSub");
  hEffSub->Reset();
  hEffSub->Divide(hNumSub,hDenSub,1,1,"B");
  hEffSub->SetLineColor(2);
  hEffSub->SetMarkerColor(2);
  hEffSub->SetMarkerStyle(4);
  hEffSub->Draw("same");

  TLegend *leg=new TLegend(0.3,0.3,0.6,0.5);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(hEff,"Jet(40||60||80) && p_{T,1}>80 GeV/c","p");
  leg->AddEntry(hEffSub,"+ Sub-leading","p");
  leg->AddEntry(hEffInc,"Inclusive","p");
  leg->Draw();
}
