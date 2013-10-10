void formatHisto(TH1F *h, char *title){
  h->Sumw2();
  h->SetMarkerStyle(8);
  h->SetXTitle(title);
  h->SetYTitle("Entries");
}

void setFillColor(TH1F *h[3]){
  h[0]->SetFillColor(kRed);
  h[1]->SetFillColor(kGreen);
  h[2]->SetFillColor(kBlue);
}

void formatCanvas(TCanvas *c){
  c->Divide(1,2,0.01,0.01);
  c->cd(1);
  c->GetPad(1)->SetLogy();
  c->GetPad(1)->SetPad(0.,0.225,1.,1.);
  c->GetPad(2)->SetPad(0.,0.0,1.,0.3);
  c->GetPad(2)->SetBottomMargin(0.3);
  c->GetPad(2)->SetGridy(1);
}

void formatRatioHist(TH1F *h){
  h->GetYaxis()->SetRangeUser(0.4,1.6);
  h->GetYaxis()->SetNdivisions(505);
  h->SetYTitle("Data/MC");
  h->GetYaxis()->CenterTitle(1);
  h->GetYaxis()->SetLabelSize(14);
  h->GetYaxis()->SetLabelFont(43);
}

void stackHistos(TH1F *hArray[3]){
  hArray[0]->SetMarkerSize(0);
  hArray[1]->SetMarkerSize(0);
  hArray[2]->SetMarkerSize(0);
  hArray[1]->Add(hArray[0]);
  hArray[2]->Add(hArray[1]);
}

void scaleHistos(TH1F *hdata, TH1F *hArray[3]){
  double scale = hdata->Integral()/hArray[2]->Integral();
  for(int i=0;i<3;i++)hArray[i]->Scale(scale);
}


void compareDataMC(int ppPbPb=0, int plotSV=1, int plotTracks=1, int savePlots=1) {

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  
  gStyle->SetLabelFont(43,"xyz");
  gStyle->SetLabelSize(14,"xyz");
  gStyle->SetTitleFont(43,"xyz");
  gStyle->SetTitleSize(16,"xyz");
  
  gStyle->SetTitleOffset(3.5,"x");
  
  gROOT->ForceStyle(1);

  TFile *fMC[3], *fdata;
  if(ppPbPb){
    fMC[0] = new TFile("./histos/PbPbBMC_pt30by3_ipHICalibCentWeight.root");
    fMC[1] = new TFile("./histos/PbPbCMC_pt30by3_ipHICalibCentWeight.root");
    fMC[2] = new TFile("./histos/PbPbQCDMC_pt30by3_ipHICalibCentWeight.root");
    fdata = new TFile("./histos/PbPbdata_pt30by3_jpHICalibRepass_withDup.root");
  }
  else{
    fMC[0] = new TFile("./histos/ppMC_ppReco_ak3PF_BjetTrig_noIPupperCut_9-15.root");
    fMC[1] = new TFile("./histos/ppMC_ppReco_ak3PF_CjetTrig_noIPupperCut_9-15.root");
    fMC[2] = new TFile("./histos/ppMC_ppReco_ak3PF_QCDjetTrig_noIPupperCut_9-15.root");
    fdata = new TFile("./histos/ppdata_ppReco_ak3PF_jetTrig_noIPupperCut_9-15.root");
  }
  // declare histos
  TH1F *hjtpt, *hrawpt, *hjteta, *hjtphi, *hnsvtx, *hsvtxntrk, *hsvtxdl, *hsvtxdls, *hsvtxm, *hsvtxmSV3, *hsvtxpt, *hsvtxptSV3, *hnIPtrk, *hnselIPtrk, *hdiscr_csvSimple, *hdiscr_prob, *hdiscr_ssvHighEff, *hdiscr_ssvHighPur, *hmuptrel, *hmuptrelSV2, *hmuptrelSV3, *hipPt, *hipProb0, *hipProb1, *hip2d, *hip2dSig, *hip3d, *hip3dSig, *hip2d1, *hip2dSig, *hip3d1, *hip3dSig1, *hip2d2, *hip2dSig2, *hip3d2, *hip3dSig2, *hip2d3, *hip2dSig3, *hip3d3, *hip3dSig3, *hipDist2Jet, *hipDist2JetSig, *hipClosest2Jet;
  TH1F *hjtptMC[4], *hrawptMC[3], *hjtetaMC[3], *hjtphiMC[3], *hnsvtxMC[3], *hsvtxntrkMC[3], *hsvtxdlMC[3], *hsvtxdlsMC[3], *hsvtxmMC[3], *hsvtxmSV3MC[3], *hsvtxptMC[3], *hsvtxptSV3MC[3], *hnIPtrkMC[3], *hnselIPtrkMC[3], *hdiscr_csvSimpleMC[3], *hdiscr_probMC[3], *hdiscr_ssvHighEffMC[3], *hdiscr_ssvHighPurMC[3], *hmuptrelMC[3], *hmuptrelSV2MC[3], *hmuptrelSV3MC[3], *hipPtMC[3], *hipProb0MC[3], *hipProb1MC[3], *hip2dMC[3], *hip2dSigMC[3], *hip3dMC[3], *hip3dSigMC[3], *hip2d1MC[3], *hip2dSig1MC[3], *hip3d1MC[3], *hip3dSig1MC[3],*hip2d2MC[3], *hip2dSig2MC[3], *hip3d2MC[3], *hip3dSig2MC[3],*hip2d3MC[3], *hip2dSig3MC[3], *hip3d3MC[3], *hip3dSig3MC[3], *hipDist2JetMC[3], *hipDist2JetSigMC[3], *hipClosest2JetMC[3];

  // grab histos
  string suffix[4]={"B","C","L","U"};
  
  string name="hjtpt";
  hjtpt = (TH1F*) fdata->Get(name.c_str());
  for(int i=0;i<3;i++) hjtptMC[i] = (TH1F*) fMC[i]->Get((name+suffix[i]).c_str());

  string name="hrawpt";
  hrawpt = (TH1F*) fdata->Get(name.c_str());
  for(int i=0;i<3;i++) hrawptMC[i] = (TH1F*) fMC[i]->Get((name+suffix[i]).c_str());

  string name="hjteta";
  hjteta = (TH1F*) fdata->Get(name.c_str());
  for(int i=0;i<3;i++) hjtetaMC[i] = (TH1F*) fMC[i]->Get((name+suffix[i]).c_str());

  string name="hjtphi";
  hjtphi = (TH1F*) fdata->Get(name.c_str());
  for(int i=0;i<3;i++) hjtphiMC[i] = (TH1F*) fMC[i]->Get((name+suffix[i]).c_str());

  string name="hnsvtx";
  hnsvtx = (TH1F*) fdata->Get(name.c_str());
  for(int i=0;i<3;i++) hnsvtxMC[i] = (TH1F*) fMC[i]->Get((name+suffix[i]).c_str());

  name="hsvtxntrk";
  hsvtxntrk = (TH1F*) fdata->Get(name.c_str());
  for(int i=0;i<3;i++) hsvtxntrkMC[i] = (TH1F*) fMC[i]->Get((name+suffix[i]).c_str());

  name="hsvtxdl";
  hsvtxdl = (TH1F*) fdata->Get(name.c_str());
  for(int i=0;i<3;i++) hsvtxdlMC[i] = (TH1F*) fMC[i]->Get((name+suffix[i]).c_str());

  name="hsvtxdls";
  hsvtxdls = (TH1F*) fdata->Get(name.c_str());
  for(int i=0;i<3;i++) hsvtxdlsMC[i] = (TH1F*) fMC[i]->Get((name+suffix[i]).c_str());

  name="hsvtxm";
  hsvtxm = (TH1F*) fdata->Get(name.c_str());
  for(int i=0;i<3;i++) hsvtxmMC[i] = (TH1F*) fMC[i]->Get((name+suffix[i]).c_str());

  name="hsvtxmSV3";
  hsvtxmSV3 = (TH1F*) fdata->Get(name.c_str());
  for(int i=0;i<3;i++) hsvtxmSV3MC[i] = (TH1F*) fMC[i]->Get((name+suffix[i]).c_str());

  name="hsvtxpt";
  hsvtxpt = (TH1F*) fdata->Get(name.c_str());
  for(int i=0;i<3;i++) hsvtxptMC[i] = (TH1F*) fMC[i]->Get((name+suffix[i]).c_str());

  name="hsvtxptSV3";
  hsvtxptSV3 = (TH1F*) fdata->Get(name.c_str());
  for(int i=0;i<3;i++) hsvtxptSV3MC[i] = (TH1F*) fMC[i]->Get((name+suffix[i]).c_str());

  name="hnIPtrk";
  hnIPtrk = (TH1F*) fdata->Get(name.c_str());
  for(int i=0;i<3;i++) hnIPtrkMC[i] = (TH1F*) fMC[i]->Get((name+suffix[i]).c_str());

  name="hnselIPtrk";
  hnselIPtrk = (TH1F*) fdata->Get(name.c_str());
  for(int i=0;i<3;i++) hnselIPtrkMC[i] = (TH1F*) fMC[i]->Get((name+suffix[i]).c_str());

  name="hdiscr_csvSimple";
  hdiscr_csvSimple = (TH1F*) fdata->Get(name.c_str());
  for(int i=0;i<3;i++) hdiscr_csvSimpleMC[i] = (TH1F*) fMC[i]->Get((name+suffix[i]).c_str());

  name="hdiscr_prob";
  hdiscr_prob = (TH1F*) fdata->Get(name.c_str());
  for(int i=0;i<3;i++) hdiscr_probMC[i] = (TH1F*) fMC[i]->Get((name+suffix[i]).c_str());

  name="hdiscr_ssvHighEff";
  hdiscr_ssvHighEff = (TH1F*) fdata->Get(name.c_str());
  for(int i=0;i<3;i++) hdiscr_ssvHighEffMC[i] = (TH1F*) fMC[i]->Get((name+suffix[i]).c_str());

  name="hdiscr_ssvHighPur";
  hdiscr_ssvHighPur = (TH1F*) fdata->Get(name.c_str());
  for(int i=0;i<3;i++) hdiscr_ssvHighPurMC[i] = (TH1F*) fMC[i]->Get((name+suffix[i]).c_str());

  name="hmuptrel";
  hmuptrel = (TH1F*) fdata->Get(name.c_str());
  for(int i=0;i<3;i++) hmuptrelMC[i] = (TH1F*) fMC[i]->Get((name+suffix[i]).c_str());

  name="hmuptrelSV2";
  hmuptrelSV2 = (TH1F*) fdata->Get(name.c_str());
  for(int i=0;i<3;i++) hmuptrelSV2MC[i] = (TH1F*) fMC[i]->Get((name+suffix[i]).c_str());

  name="hmuptrelSV3";
  hmuptrelSV3 = (TH1F*) fdata->Get(name.c_str());
  for(int i=0;i<3;i++) hmuptrelSV3MC[i] = (TH1F*) fMC[i]->Get((name+suffix[i]).c_str());

  name="hipPt";
  hipPt = (TH1F*) fdata->Get(name.c_str());
  for(int i=0;i<3;i++) hipPtMC[i] = (TH1F*) fMC[i]->Get((name+suffix[i]).c_str());

  name="hipProb0";
  hipProb0 = (TH1F*) fdata->Get(name.c_str());
  for(int i=0;i<3;i++) hipProb0MC[i] = (TH1F*) fMC[i]->Get((name+suffix[i]).c_str());

  name="hipProb1";
  hipProb1 = (TH1F*) fdata->Get(name.c_str());
  for(int i=0;i<3;i++) hipProb1MC[i] = (TH1F*) fMC[i]->Get((name+suffix[i]).c_str());

  name="hip2d";
  hip2d = (TH1F*) fdata->Get(name.c_str());
  hip2d->Rebin(5);
  for(int i=0;i<3;i++){
    hip2dMC[i] = (TH1F*) fMC[i]->Get((name+suffix[i]).c_str());
    hip2dMC[i]->Rebin(5);
  }
  name="hip2dSig";
  hip2dSig = (TH1F*) fdata->Get(name.c_str());
  hip2dSig->Rebin(5);
  for(int i=0;i<3;i++){
    hip2dSigMC[i] = (TH1F*) fMC[i]->Get((name+suffix[i]).c_str());
    hip2dSigMC[i]->Rebin(5);
  }
  name="hip3d";
  hip3d = (TH1F*) fdata->Get(name.c_str());
  hip3d->Rebin(5);
  for(int i=0;i<3;i++) {
    hip3dMC[i] = (TH1F*) fMC[i]->Get((name+suffix[i]).c_str());
    hip3dMC[i]->Rebin(5);
  }
  name="hip3dSig";
  hip3dSig = (TH1F*) fdata->Get(name.c_str());
  hip3dSig->Rebin(5);
  for(int i=0;i<3;i++) {
    hip3dSigMC[i] = (TH1F*) fMC[i]->Get((name+suffix[i]).c_str());
    hip3dSigMC[i]->Rebin(5);
  }

  name="hip2d1";
  hip2d1 = (TH1F*) fdata->Get(name.c_str());
  hip2d1->Rebin(5);
  for(int i=0;i<3;i++){
    hip2d1MC[i] = (TH1F*) fMC[i]->Get((name+suffix[i]).c_str());
    hip2d1MC[i]->Rebin(5);
  }
  name="hip2dSig1";
  hip2dSig1 = (TH1F*) fdata->Get(name.c_str());
  hip2dSig1->Rebin(5);
  for(int i=0;i<3;i++){
    hip2dSig1MC[i] = (TH1F*) fMC[i]->Get((name+suffix[i]).c_str());
    hip2dSig1MC[i]->Rebin(5);
  }
  name="hip3d1";
  hip3d1 = (TH1F*) fdata->Get(name.c_str());
  hip3d1->Rebin(5);
  for(int i=0;i<3;i++) {
    hip3d1MC[i] = (TH1F*) fMC[i]->Get((name+suffix[i]).c_str());
    hip3d1MC[i]->Rebin(5);
  }
  name="hip3dSig1";
  hip3dSig1 = (TH1F*) fdata->Get(name.c_str());
  hip3dSig1->Rebin(5);
  for(int i=0;i<3;i++) {
    hip3dSig1MC[i] = (TH1F*) fMC[i]->Get((name+suffix[i]).c_str());
    hip3dSig1MC[i]->Rebin(5);
  }

  name="hip2d2";
  hip2d2 = (TH1F*) fdata->Get(name.c_str());
  hip2d2->Rebin(5);
  for(int i=0;i<3;i++){
    hip2d2MC[i] = (TH1F*) fMC[i]->Get((name+suffix[i]).c_str());
    hip2d2MC[i]->Rebin(5);
  }
  name="hip2dSig2";
  hip2dSig2 = (TH1F*) fdata->Get(name.c_str());
  hip2dSig2->Rebin(5);
  for(int i=0;i<3;i++){
    hip2dSig2MC[i] = (TH1F*) fMC[i]->Get((name+suffix[i]).c_str());
    hip2dSig2MC[i]->Rebin(5);
  }
  name="hip3d2";
  hip3d2 = (TH1F*) fdata->Get(name.c_str());
  hip3d2->Rebin(5);
  for(int i=0;i<3;i++) {
    hip3d2MC[i] = (TH1F*) fMC[i]->Get((name+suffix[i]).c_str());
    hip3d2MC[i]->Rebin(5);
  }
  name="hip3dSig2";
  hip3dSig2 = (TH1F*) fdata->Get(name.c_str());
  hip3dSig2->Rebin(5);
  for(int i=0;i<3;i++) {
    hip3dSig2MC[i] = (TH1F*) fMC[i]->Get((name+suffix[i]).c_str());
    hip3dSig2MC[i]->Rebin(5);
  }

  name="hip2d3";
  hip2d3 = (TH1F*) fdata->Get(name.c_str());
  hip2d3->Rebin(5);
  for(int i=0;i<3;i++){
    hip2d3MC[i] = (TH1F*) fMC[i]->Get((name+suffix[i]).c_str());
    hip2d3MC[i]->Rebin(5);
  }
  name="hip2dSig3";
  hip2dSig3 = (TH1F*) fdata->Get(name.c_str());
  hip2dSig3->Rebin(5);
  for(int i=0;i<3;i++){
    hip2dSig3MC[i] = (TH1F*) fMC[i]->Get((name+suffix[i]).c_str());
    hip2dSig3MC[i]->Rebin(5);
  }
  name="hip3d3";
  hip3d3 = (TH1F*) fdata->Get(name.c_str());
  hip3d3->Rebin(5);
  for(int i=0;i<3;i++) {
    hip3d3MC[i] = (TH1F*) fMC[i]->Get((name+suffix[i]).c_str());
    hip3d3MC[i]->Rebin(5);
  }
  name="hip3dSig3";
  hip3dSig3 = (TH1F*) fdata->Get(name.c_str());
  hip3dSig3->Rebin(5);
  for(int i=0;i<3;i++) {
    hip3dSig3MC[i] = (TH1F*) fMC[i]->Get((name+suffix[i]).c_str());
    hip3dSig3MC[i]->Rebin(5);
  }
  name="hipDist2Jet";
  hipDist2Jet = (TH1F*) fdata->Get(name.c_str());
  for(int i=0;i<3;i++) hipDist2JetMC[i] = (TH1F*) fMC[i]->Get((name+suffix[i]).c_str());

  name="hipDist2JetSig";
  hipDist2JetSig = (TH1F*) fdata->Get(name.c_str());
  for(int i=0;i<3;i++) hipDist2JetSigMC[i] = (TH1F*) fMC[i]->Get((name+suffix[i]).c_str());

  name="hipClosest2Jet";
  hipClosest2Jet = (TH1F*) fdata->Get(name.c_str());
  for(int i=0;i<3;i++) hipClosest2JetMC[i] = (TH1F*) fMC[i]->Get((name+suffix[i]).c_str());


  // turn off heavy flavor for testing
  /*
//  for(int i=0;i<2;i++) hjtptMC[i]->Scale(0.);
//  for(int i=0;i<2;i++) hnsvtxMC[i]->Scale(0.);
//  for(int i=0;i<2;i++) hsvtxntrkMC[i]->Scale(0.);
//  for(int i=0;i<2;i++) hsvtxdlMC[i]->Scale(0.);
//  for(int i=0;i<2;i++) hsvtxdlsMC[i]->Scale(0.);
//  for(int i=0;i<2;i++) hsvtxmMC[i]->Scale(0.);
//  for(int i=0;i<2;i++) hsvtxptMC[i]->Scale(0.);
//  for(int i=0;i<2;i++) hnIPtrkMC[i]->Scale(0.);
//  for(int i=0;i<2;i++) hnselIPtrkMC[i]->Scale(0.);
//  for(int i=0;i<2;i++) hdiscr_csvSimpleMC[i]->Scale(0.);
  */

  // stack em
  stackHistos(hjtptMC);
  stackHistos(hrawptMC);
  stackHistos(hjtetaMC);
  stackHistos(hjtphiMC);
  stackHistos(hnsvtxMC);
  stackHistos(hsvtxntrkMC);
  stackHistos(hsvtxdlMC);
  stackHistos(hsvtxdlsMC);
  stackHistos(hsvtxmMC);
  stackHistos(hsvtxmSV3MC);
  stackHistos(hsvtxptMC);
  stackHistos(hsvtxptSV3MC);
  stackHistos(hnIPtrkMC);
  stackHistos(hnselIPtrkMC);
  stackHistos(hdiscr_csvSimpleMC);
  stackHistos(hdiscr_probMC);
  stackHistos(hdiscr_ssvHighEffMC);
  stackHistos(hdiscr_ssvHighPurMC);
  stackHistos(hmuptrelMC);
  stackHistos(hmuptrelSV2MC);
  stackHistos(hmuptrelSV3MC);
  stackHistos(hipPtMC);
  stackHistos(hipProb0MC);
  stackHistos(hipProb1MC);

  stackHistos(hip2dMC);
  stackHistos(hip2dSigMC);
  stackHistos(hip3dMC);
  stackHistos(hip3dSigMC);

  stackHistos(hip2d1MC);
  stackHistos(hip2dSig1MC);
  stackHistos(hip3d1MC);
  stackHistos(hip3dSig1MC);

  stackHistos(hip2d2MC);
  stackHistos(hip2dSig2MC);
  stackHistos(hip3d2MC);
  stackHistos(hip3dSig2MC);

  stackHistos(hip2d3MC);
  stackHistos(hip2dSig3MC);
  stackHistos(hip3d3MC);
  stackHistos(hip3dSig3MC);

  stackHistos(hipDist2JetMC);
  stackHistos(hipDist2JetSigMC);
  stackHistos(hipClosest2JetMC);

  // scale em
  //double scale = hjtpt->Integral()/hjtptMC[2]->Integral();

  scaleHistos(hjtpt,hjtptMC);
  scaleHistos(hrawpt,hrawptMC);
  scaleHistos(hjteta,hjtetaMC);
  scaleHistos(hjtphi,hjtphiMC);
  scaleHistos(hnsvtx,hnsvtxMC);
  scaleHistos(hsvtxntrk,hsvtxntrkMC);
  scaleHistos(hsvtxdl,hsvtxdlMC);
  scaleHistos(hsvtxdls,hsvtxdlsMC);
  scaleHistos(hsvtxm,hsvtxmMC);
  scaleHistos(hsvtxmSV3,hsvtxmSV3MC);
  scaleHistos(hsvtxpt,hsvtxptMC);
  scaleHistos(hsvtxptSV3,hsvtxptSV3MC);
  scaleHistos(hnIPtrk,hnIPtrkMC);
  scaleHistos(hnselIPtrk,hnselIPtrkMC);
  scaleHistos(hdiscr_csvSimple,hdiscr_csvSimpleMC);
  scaleHistos(hdiscr_prob,hdiscr_probMC);
  scaleHistos(hdiscr_ssvHighEff,hdiscr_ssvHighEffMC);
  scaleHistos(hdiscr_ssvHighPur,hdiscr_ssvHighPurMC);
  scaleHistos(hmuptrel,hmuptrelMC);
  scaleHistos(hmuptrelSV2,hmuptrelSV2MC);
  scaleHistos(hmuptrelSV3,hmuptrelSV3MC);
  scaleHistos(hipPt,hipPtMC);
  scaleHistos(hipProb0,hipProb0MC);
  scaleHistos(hipProb1,hipProb1MC);

  scaleHistos(hip2d,hip2dMC);
  scaleHistos(hip2dSig,hip2dSigMC);
  scaleHistos(hip3d,hip3dMC);
  scaleHistos(hip3dSig,hip3dSigMC);

  scaleHistos(hip2d1,hip2d1MC);
  scaleHistos(hip2dSig1,hip2dSig1MC);
  scaleHistos(hip3d1,hip3d1MC);
  scaleHistos(hip3dSig1,hip3dSig1MC);

  scaleHistos(hip2d2,hip2d2MC);
  scaleHistos(hip2dSig2,hip2dSig2MC);
  scaleHistos(hip3d2,hip3d2MC);
  scaleHistos(hip3dSig2,hip3dSig2MC);

  scaleHistos(hip2d3,hip2d3MC);
  scaleHistos(hip2dSig3,hip2dSig3MC);
  scaleHistos(hip3d3,hip3d3MC);
  scaleHistos(hip3dSig3,hip3dSig3MC);
  scaleHistos(hipDist2Jet,hipDist2JetMC);
  scaleHistos(hipDist2JetSig,hipDist2JetSigMC);
  scaleHistos(hipClosest2Jet,hipClosest2JetMC);

  /*  
  for(int i=0;i<3;i++) hjtptMC[i]->Scale(scale);
  for(int i=0;i<3;i++) hnsvtxMC[i]->Scale(scale);
  for(int i=0;i<3;i++) hsvtxntrkMC[i]->Scale(scale);
  for(int i=0;i<3;i++) hsvtxdlMC[i]->Scale(scale);
  for(int i=0;i<3;i++) hsvtxdlsMC[i]->Scale(scale);
  for(int i=0;i<3;i++) hsvtxmMC[i]->Scale(scale);
  for(int i=0;i<3;i++) hsvtxmSV3MC[i]->Scale(scale);
  for(int i=0;i<3;i++) hsvtxptMC[i]->Scale(scale);
  for(int i=0;i<3;i++) hsvtxptSV3MC[i]->Scale(scale);
  for(int i=0;i<3;i++) hnIPtrkMC[i]->Scale(scale);
  for(int i=0;i<3;i++) hnselIPtrkMC[i]->Scale(scale);
  for(int i=0;i<3;i++) hdiscr_csvSimpleMC[i]->Scale(scale);
  for(int i=0;i<3;i++) hdiscr_probMC[i]->Scale(scale);
  for(int i=0;i<3;i++) hdiscr_ssvHighEffMC[i]->Scale(scale);
  for(int i=0;i<3;i++) hdiscr_ssvHighPurMC[i]->Scale(scale);
  for(int i=0;i<3;i++) hmuptrelMC[i]->Scale(scale);
  for(int i=0;i<3;i++) hmuptrelSV2MC[i]->Scale(scale);
  for(int i=0;i<3;i++) hmuptrelSV3MC[i]->Scale(scale);
  for(int i=0;i<3;i++) hipPtMC[i]->Scale(scale);
  for(int i=0;i<3;i++) hipProb0MC[i]->Scale(scale);
  for(int i=0;i<3;i++) hipProb1MC[i]->Scale(scale);
  for(int i=0;i<3;i++) hip2dMC[i]->Scale(scale);
  for(int i=0;i<3;i++) hip2dSigMC[i]->Scale(scale);
  for(int i=0;i<3;i++) hip3dMC[i]->Scale(scale);
  for(int i=0;i<3;i++) hip3dSigMC[i]->Scale(scale);
  for(int i=0;i<3;i++) hipDist2JetMC[i]->Scale(scale);
  for(int i=0;i<3;i++) hipDist2JetSigMC[i]->Scale(scale);
  for(int i=0;i<3;i++) hipClosest2JetMC[i]->Scale(scale);
  */
  // style em
  setFillColor(hjtptMC);
  setFillColor(hrawptMC);
  setFillColor(hjtetaMC);
  setFillColor(hjtphiMC);
  setFillColor(hnsvtxMC);
  setFillColor(hsvtxntrkMC);
  setFillColor(hsvtxdlMC);
  setFillColor(hsvtxdlsMC);
  setFillColor(hsvtxmMC);
  setFillColor(hsvtxmSV3MC);
  setFillColor(hsvtxptMC);
  setFillColor(hsvtxptSV3MC);
  setFillColor(hnIPtrkMC);
  setFillColor(hnselIPtrkMC);
  setFillColor(hdiscr_csvSimpleMC);
  setFillColor(hdiscr_probMC);
  setFillColor(hdiscr_ssvHighEffMC);
  setFillColor(hdiscr_ssvHighPurMC);
  setFillColor(hmuptrelMC);
  setFillColor(hmuptrelSV2MC);
  setFillColor(hmuptrelSV3MC);
  setFillColor(hipPtMC);
  setFillColor(hipProb0MC);
  setFillColor(hipProb1MC);

  setFillColor(hip2dMC);
  setFillColor(hip2dSigMC);
  setFillColor(hip3dMC);
  setFillColor(hip3dSigMC);

  setFillColor(hip2d1MC);
  setFillColor(hip2dSig1MC);
  setFillColor(hip3d1MC);
  setFillColor(hip3dSig1MC);

  setFillColor(hip2d2MC);
  setFillColor(hip2dSig2MC);
  setFillColor(hip3d2MC);
  setFillColor(hip3dSig2MC);

  setFillColor(hip2d3MC);
  setFillColor(hip2dSig3MC);
  setFillColor(hip3d3MC);
  setFillColor(hip3dSig3MC);

  setFillColor(hipDist2JetMC);
  setFillColor(hipDist2JetSigMC);
  setFillColor(hipClosest2JetMC);

  TLegend *leg = new TLegend(0.699,0.62,0.89,0.92);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  if(ppPbPb)leg->AddEntry(hjtpt,"PbPb data, 0-100%","p");
  else leg->AddEntry(hjtpt,"pp data","p");
  leg->AddEntry(hjtptMC[0],"b-jet","f");
  leg->AddEntry(hjtptMC[1],"c-jet","f");
  leg->AddEntry(hjtptMC[2],"udsg-jet","f");


  // draw em

  char *title = "Jet p_{T} (GeV/c)";
  TCanvas *c0=new TCanvas("c0",title,200,10,600,480);
  formatHisto(hjtpt,title);
  formatCanvas(c0);
  //hjtpt->Rebin(4);
  hjtpt->GetXaxis()->SetRangeUser(80,200);
  hjtpt->Draw();   
  for(int i=2;i>-1;i--){
    //hjtptMC[i]->Rebin(4);
    hjtptMC[i]->Draw("hist,e,same");
  }
  hjtpt->Draw("same");
  leg->Draw();
  c0->GetPad(1)->RedrawAxis();

  c0->cd(2);
  TH1F *hRjtpt = hjtpt->Clone("hRjtpt");
  formatRatioHist(hRjtpt);
  hRjtpt->Divide(hjtptMC[2]);
  hRjtpt->Draw();
  c0->SaveAs("plotsNewpp/jetPt.pdf");

  char *title = "Jet raw p_{T} (GeV/c)";
  TCanvas *c0bis=new TCanvas("c0bis",title,200,10,600,480);
  formatHisto(hrawpt,title);
  formatCanvas(c0bis);
  hrawpt->Draw();   
  for(int i=2;i>-1;i--) hrawptMC[i]->Draw("hist,e,same");
  hrawpt->Draw("same");
  leg->Draw();
  c0bis->GetPad(1)->RedrawAxis();

  c0bis->cd(2);
  TH1F *hRrawpt = hrawpt->Clone("hRrawpt");
  formatRatioHist(hRrawpt);
  hRrawpt->Divide(hrawptMC[2]);
  hRrawpt->Draw();
  c0bis->SaveAs("plotsNewpp/jetRawPt.pdf");

  char *title = "Jet #eta";
  TCanvas *c0a=new TCanvas("c0a",title,200,10,600,480);
  formatHisto(hjteta,title);
  formatCanvas(c0a);
  hjteta->SetMaximum(500.*hjteta->GetMaximum());
  hjteta->SetMinimum(10.);
  hjteta->Draw();   
  for(int i=2;i>-1;i--) hjtetaMC[i]->Draw("hist,e,same");
  hjteta->Draw("same");
  leg->Draw();
  c0a->GetPad(1)->RedrawAxis();

  c0a->cd(2);
  TH1F *hRjteta = hjteta->Clone("hRjteta");
  formatRatioHist(hRjteta);
  hRjteta->Divide(hjtetaMC[2]);
  hRjteta->Draw();
  c0a->SaveAs(Form("plotsNewpp/%s.pdf",title));

  char *title = "Jet #phi";
  TCanvas *c0b=new TCanvas("c0b",title,200,10,600,480);
  formatHisto(hjtphi,title);
  formatCanvas(c0b);
  hjtphi->SetMaximum(100.*hjteta->GetMaximum());
  hjtphi->SetMinimum(10.);
  hjtphi->Draw();   
  for(int i=2;i>-1;i--) hjtphiMC[i]->Draw("hist,e,same");
  hjtphi->Draw("same");
  leg->Draw();
  if(savePlots)c0b->GetPad(1)->RedrawAxis();

  c0b->cd(2);
  TH1F *hRjtphi = hjtphi->Clone("hRjtphi");
  formatRatioHist(hRjtphi);
  hRjtphi->Divide(hjtphiMC[2]);
  hRjtphi->Draw();
  if(savePlots)c0b->SaveAs(Form("plotsNewpp/%s.pdf",title));

  if(plotSV){
  title = "# of secondary vertices";
  TCanvas *c1=new TCanvas("c1",title,200,10,600,480);
  formatHisto(hnsvtx,title);
  formatCanvas(c1);
  hnsvtx->Draw();   
  for(int i=2;i>-1;i--) hnsvtxMC[i]->Draw("hist,e,same");
  hnsvtx->Draw("same");
  leg->Draw();
  if(savePlots)c1->GetPad(1)->RedrawAxis();

  c1->cd(2);
  TH1F *hRnsvtx = (TH1F*) hnsvtx->Clone("hRnsvtx");
  formatRatioHist(hRnsvtx);
  hRnsvtx->Divide(hnsvtxMC[2]);
  hRnsvtx->Draw();
  if(savePlots)c1->SaveAs(Form("plotsNewpp/%s.pdf",title));

  title = "tracks per SV";
  TCanvas *c2=new TCanvas("c2",title,200,10,600,480);
  formatHisto(hsvtxntrk,title);
  formatCanvas(c2);
  hsvtxntrk->Draw();
  for(int i=2;i>-1;i--) hsvtxntrkMC[i]->Draw("hist,e,same");
  hsvtxntrk->Draw("same");
  c2->GetPad(1)->RedrawAxis();
  leg->Draw();

  c2->cd(2);
  TH1F *hRsvtxntrk = hsvtxntrk->Clone("hRsvtxntrk");
  formatRatioHist(hRsvtxntrk);
  hRsvtxntrk->Divide(hsvtxntrkMC[2]);
  hRsvtxntrk->Draw();
  if(savePlots)c2->SaveAs(Form("plotsNewpp/%s.pdf",title));

  title = "flight distance";
  TCanvas *c3=new TCanvas("c3",title,200,10,600,480);
  formatHisto(hsvtxdl,title);
  formatCanvas(c3);
  hsvtxdl->Draw();
  for(int i=2;i>-1;i--) hsvtxdlMC[i]->Draw("hist,e,same");
  hsvtxdl->Draw("same");
  leg->Draw();
  c3->GetPad(1)->RedrawAxis();

  c3->cd(2);
  TH1F *hRsvtxdl = hsvtxdl->Clone("hRsvtxdl");
  formatRatioHist(hRsvtxdl);
  hRsvtxdl->Divide(hsvtxdlMC[2]);
  hRsvtxdl->Draw();
  if(savePlots)c3->SaveAs(Form("plotsNewpp/%s.pdf",title));

  title = "flight distance significance";
  TCanvas *c4=new TCanvas("c4",title,200,10,600,480);
  formatHisto(hsvtxdls,title);
  formatCanvas(c4);
  hsvtxdls->Draw();
  for(int i=2;i>-1;i--) hsvtxdlsMC[i]->Draw("hist,e,same");
  hsvtxdls->Draw("same");
  leg->Draw();
  c4->GetPad(1)->RedrawAxis();

  c4->cd(2);
  TH1F *hRsvtxdls = hsvtxdls->Clone("hRsvtxdls");
  formatRatioHist(hRsvtxdls);
  hRsvtxdls->Divide(hsvtxdlsMC[2]);
  hRsvtxdls->Draw();
  if(savePlots)c4->SaveAs(Form("plotsNewpp/%s.pdf",title));

  title = "SV mass";
  TCanvas *c5=new TCanvas("c5",title,200,10,600,480);
  formatHisto(hsvtxm,title);
  formatCanvas(c5);
  hsvtxm->Draw();
  for(int i=2;i>-1;i--) hsvtxmMC[i]->Draw("hist,e,same");
  hsvtxm->Draw("same");
  leg->Draw();
  c5->GetPad(1)->RedrawAxis();

  c5->cd(2);
  TH1F *hRsvtxm = hsvtxm->Clone("hRsvtxm");
  formatRatioHist(hRsvtxm);
  hRsvtxm->Divide(hsvtxmMC[2]);
  hRsvtxm->Draw();
  if(savePlots)c5->SaveAs(Form("plotsNewpp/%s.pdf",title));

  title = "3 track SV mass";
  TCanvas *c5b=new TCanvas("c5b",title,200,10,600,480);
  formatHisto(hsvtxmSV3,title);
  formatCanvas(c5b);
  hsvtxmSV3->Draw();
  for(int i=2;i>-1;i--) hsvtxmSV3MC[i]->Draw("hist,e,same");
  hsvtxmSV3->Draw("same");
  leg->Draw();
  c5b->GetPad(1)->RedrawAxis();

  c5b->cd(2);
  TH1F *hRsvtxmSV3 = hsvtxmSV3->Clone("hRsvtxmSV3");
  formatRatioHist(hRsvtxmSV3);
  hRsvtxmSV3->Divide(hsvtxmSV3MC[2]);
  hRsvtxmSV3->Draw();
  if(savePlots)c5b->SaveAs(Form("plotsNewpp/%s.pdf",title));

  title = "SV pT";
  TCanvas *c6=new TCanvas("c6",title,200,10,600,480);
  formatHisto(hsvtxpt,title);
  formatCanvas(c6);
  hsvtxpt->Draw();
  for(int i=2;i>-1;i--) hsvtxptMC[i]->Draw("hist,e,same");
  hsvtxpt->Draw("same");
  leg->Draw();
  c6->GetPad(1)->RedrawAxis();

  c6->cd(2);
  TH1F *hRsvtxpt = hsvtxpt->Clone("hRsvtxpt");
  formatRatioHist(hRsvtxpt);
  hRsvtxpt->Divide(hsvtxptMC[2]);
  hRsvtxpt->Draw();
  if(savePlots)c6->SaveAs(Form("plotsNewpp/%s.pdf",title));

  title = "3 track SV pT";
  TCanvas *c6b=new TCanvas("c6b",title,200,10,600,480);
  formatHisto(hsvtxptSV3,title);
  formatCanvas(c6b);
  hsvtxptSV3->Draw();
  for(int i=2;i>-1;i--) hsvtxptSV3MC[i]->Draw("hist,e,same");
  hsvtxptSV3->Draw("same");
  leg->Draw();
  c6b->GetPad(1)->RedrawAxis();

  c6b->cd(2);
  TH1F *hRsvtxptSV3 = hsvtxptSV3->Clone("hRsvtxptSV3");
  formatRatioHist(hRsvtxptSV3);
  hRsvtxptSV3->Divide(hsvtxptSV3MC[2]);
  hRsvtxptSV3->Draw();
  if(savePlots)c6b->SaveAs(Form("plotsNewpp/%s.pdf",title));
  }

  if(plotTracks){
  title = "# of IP tracks";
  TCanvas *c7=new TCanvas("c7",title,200,10,600,480);
  formatHisto(hnIPtrk,title);
  formatCanvas(c7);
  hnIPtrk->GetXaxis()->SetRangeUser(0,60);
  hnIPtrk->SetMaximum(10.*hnIPtrk->GetMaximum());
  hnIPtrk->Draw();
  for(int i=2;i>-1;i--) hnIPtrkMC[i]->Draw("hist,e,same");
  hnIPtrk->Draw("same");
  leg->Draw();
  c7->GetPad(1)->RedrawAxis();

  c7->cd(2);
  TH1F *hRnIPtrk = hnIPtrk->Clone("hRnIPtrk");
  formatRatioHist(hRnIPtrk);
  hRnIPtrk->Divide(hnIPtrkMC[2]);
  hRnIPtrk->GetXaxis()->SetRangeUser(0,60);
  hRnIPtrk->Draw();
  if(savePlots)c7->SaveAs(Form("plotsNewpp/%s.pdf",title));

  title = "# of selected IP tracks";
  TCanvas *c8=new TCanvas("c8",title,200,10,600,480);
  formatHisto(hnselIPtrk,title);
  formatCanvas(c8);
  hnselIPtrk->GetXaxis()->SetRangeUser(0,60);
  hnselIPtrk->SetMaximum(10.*hnselIPtrk->GetMaximum());
  hnselIPtrk->Draw();
  for(int i=2;i>-1;i--) hnselIPtrkMC[i]->Draw("hist,e,same");
  hnselIPtrk->Draw("same");
  leg->Draw();
  c8->GetPad(1)->RedrawAxis();

  c8->cd(2);
  TH1F *hRnselIPtrk = hnselIPtrk->Clone("hRnselIPtrk");
  formatRatioHist(hRnselIPtrk);
  hRnselIPtrk->Divide(hnselIPtrkMC[2]);
  hRnselIPtrk->GetXaxis()->SetRangeUser(0,60);
  hRnselIPtrk->Draw();
  if(savePlots)c8->SaveAs(Form("plotsNewpp/%s.pdf",title));
  }

  if(plotSV){
  title = "CSV Discriminator";
  TCanvas *c9=new TCanvas("c9",title,200,10,600,480);
  formatHisto(hdiscr_csvSimple,title);
  formatCanvas(c9);
  hdiscr_csvSimple->Draw();
  for(int i=2;i>-1;i--) hdiscr_csvSimpleMC[i]->Draw("hist,e,same");
  hdiscr_csvSimple->Draw("same");
  leg->Draw();
  c9->GetPad(1)->RedrawAxis();

  c9->cd(2);
  TH1F *hRdiscr_csvSimple = hdiscr_csvSimple->Clone("hRdiscr_csvSimple");
  formatRatioHist(hRdiscr_csvSimple);
  hRdiscr_csvSimple->Divide(hdiscr_csvSimpleMC[2]);
  hRdiscr_csvSimple->Draw();
  if(savePlots)c9->SaveAs(Form("plotsNewpp/%s.pdf",title));
  }

  if(plotTracks){
  title = "Jet Prob Discriminator";
  TCanvas *c10=new TCanvas("c10",title,200,10,600,480);
  formatHisto(hdiscr_prob,title);
  formatCanvas(c10);
  hdiscr_prob->Draw();
  for(int i=2;i>-1;i--) hdiscr_probMC[i]->Draw("hist,e,same");
  hdiscr_prob->Draw("same");
  leg->Draw();
  c10->GetPad(1)->RedrawAxis();

  c10->cd(2);
  TH1F *hRdiscr_prob = hdiscr_prob->Clone("hRdiscr_prob");
  formatRatioHist(hRdiscr_prob);
  hRdiscr_prob->Divide(hdiscr_probMC[2]);
  hRdiscr_prob->Draw();
  if(savePlots)c10->SaveAs(Form("plotsNewpp/%s.pdf",title));
  }

  if(plotSV){
  title = "SSVHE";
  TCanvas *c10a=new TCanvas("c10a",title,200,10,600,480);
  formatHisto(hdiscr_ssvHighEff,title);
  formatCanvas(c10a);
  hdiscr_ssvHighEff->Draw();
  for(int i=2;i>-1;i--) hdiscr_ssvHighEffMC[i]->Draw("hist,e,same");
  hdiscr_ssvHighEff->Draw("same");
  leg->Draw();
  c10a->GetPad(1)->RedrawAxis();

  c10a->cd(2);
  TH1F *hRdiscr_ssvHighEff = hdiscr_ssvHighEff->Clone("hRdiscr_ssvHighEff");
  formatRatioHist(hRdiscr_ssvHighEff);
  hRdiscr_ssvHighEff->Divide(hdiscr_ssvHighEffMC[2]);
  hRdiscr_ssvHighEff->Draw();
  if(savePlots)c10a->SaveAs(Form("plotsNewpp/%s.pdf",title));

  title = "SSVHP";
  TCanvas *c10b=new TCanvas("c10b",title,200,10,600,480);
  formatHisto(hdiscr_ssvHighPur,title);
  formatCanvas(c10b);
  hdiscr_ssvHighPur->Draw();
  for(int i=2;i>-1;i--) hdiscr_ssvHighPurMC[i]->Draw("hist,e,same");
  hdiscr_ssvHighPur->Draw("same");
  leg->Draw();
  c10b->GetPad(1)->RedrawAxis();

  c10b->cd(2);
  TH1F *hRdiscr_ssvHighPur = hdiscr_ssvHighPur->Clone("hRdiscr_ssvHighPur");
  formatRatioHist(hRdiscr_ssvHighPur);
  hRdiscr_ssvHighPur->Divide(hdiscr_ssvHighPurMC[2]);
  hRdiscr_ssvHighPur->Draw();
  if(savePlots)c10b->SaveAs(Form("plotsNewpp/%s.pdf",title));
  }

  if(plotTracks){
  title = "hipPt";
  TCanvas *c11=new TCanvas("c11",title,200,10,600,480);
  formatHisto(hipPt,title);
  formatCanvas(c11);
  hipPt->Draw();
  for(int i=2;i>-1;i--) hipPtMC[i]->Draw("hist,e,same");
  hipPt->Draw("same");
  leg->Draw();
  c11->GetPad(1)->RedrawAxis();

  c11->cd(2);
  TH1F *hRipPt = hipPt->Clone("hRipPt");
  formatRatioHist(hRipPt);
  hRipPt->Divide(hipPtMC[2]);
  hRipPt->Draw();
  if(savePlots)c11->SaveAs(Form("plotsNewpp/%s.pdf",title));

  title = "hipProb0";
  TCanvas *c12=new TCanvas("c12",title,200,10,600,480);
  formatHisto(hipProb0,title);
  formatCanvas(c12);
  hipProb0->Draw();
  for(int i=2;i>-1;i--) hipProb0MC[i]->Draw("hist,e,same");
  hipProb0->Draw("same");
  leg->Draw();
  c12->GetPad(1)->RedrawAxis();

  c12->cd(2);
  TH1F *hRipProb0 = hipProb0->Clone("hRipProb0");
  formatRatioHist(hRipProb0);
  hRipProb0->Divide(hipProb0MC[2]);
  hRipProb0->Draw();
  if(savePlots)c12->SaveAs(Form("plotsNewpp/%s.pdf",title));

  title = "hipProb1";
  TCanvas *c13=new TCanvas("c13",title,200,10,600,480);
  formatHisto(hipProb1,title);
  formatCanvas(c13);
  hipProb1->Draw();
  for(int i=2;i>-1;i--) hipProb1MC[i]->Draw("hist,e,same");
  hipProb1->Draw("same");
  leg->Draw();
  c13->GetPad(1)->RedrawAxis();

  c13->cd(2);
  TH1F *hRipProb1 = hipProb1->Clone("hRipProb1");
  formatRatioHist(hRipProb1);
  hRipProb1->Divide(hipProb1MC[2]);
  hRipProb1->Draw();
  if(savePlots)c13->SaveAs(Form("plotsNewpp/%s.pdf",title));
  /*  // doing it track by track instead
  title = "2D IP";
  TCanvas *c14=new TCanvas("c14",title,200,10,600,480);
  formatHisto(hip2d,title);
  formatCanvas(c14);
  hip2d->Draw();
  for(int i=2;i>-1;i--) hip2dMC[i]->Draw("hist,e,same");
  hip2d->Draw("same");
  leg->Draw();
  c14->GetPad(1)->RedrawAxis();

  c14->cd(2);
  TH1F *hRip2d = hip2d->Clone("hRip2d");
  formatRatioHist(hRip2d);
  hRip2d->Divide(hip2dMC[2]);
  hRip2d->Draw();
  if(savePlots)c14->SaveAs(Form("plotsNewpp/%s.pdf",title));

  title = "2D IP Significance";
  TCanvas *c15=new TCanvas("c15",title,200,10,600,480);
  formatHisto(hip2dSig,title);
  formatCanvas(c15);
  hip2dSig->Draw();
  for(int i=2;i>-1;i--) hip2dSigMC[i]->Draw("hist,e,same");
  hip2dSig->Draw("same");
  leg->Draw();
  c15->GetPad(1)->RedrawAxis();

  c15->cd(2);
  TH1F *hRip2dSig = hip2dSig->Clone("hRip2dSig");
  formatRatioHist(hRip2dSig);
  hRip2dSig->Divide(hip2dSigMC[2]);
  hRip2dSig->Draw();
  if(savePlots)c15->SaveAs(Form("plotsNewpp/%s.pdf",title));
  
  title = "3D IP";
  TCanvas *c16=new TCanvas("c16",title,200,10,600,480);
  formatHisto(hip3d,title);
  formatCanvas(c16);
  hip3d->Draw();
  for(int i=2;i>-1;i--) hip3dMC[i]->Draw("hist,e,same");
  hip3d->Draw("same");
  leg->Draw();
  c16->GetPad(1)->RedrawAxis();

  c16->cd(2);
  TH1F *hRip3d = hip3d->Clone("hRip3d");
  formatRatioHist(hRip3d);
  hRip3d->Divide(hip3dMC[2]);
  hRip3d->Draw();
  if(savePlots)c16->SaveAs(Form("plotsNewpp/%s.pdf",title));

  title = "3D IP Significance";
  TCanvas *c17=new TCanvas("c17",title,200,10,600,480);
  formatHisto(hip3dSig,title);
  formatCanvas(c17);
  hip3dSig->Draw();
  for(int i=2;i>-1;i--) hip3dSigMC[i]->Draw("hist,e,same");
  hip3dSig->Draw("same");
  leg->Draw();
  c17->GetPad(1)->RedrawAxis();

  c17->cd(2);
  TH1F *hRip3dSig = hip3dSig->Clone("hRip3dSig");
  formatRatioHist(hRip3dSig);
  hRip3dSig->Divide(hip3dSigMC[2]);
  hRip3dSig->Draw();
  if(savePlots)c17->SaveAs(Form("plotsNewpp/%s.pdf",title));
  */

  title = "2D IP 1st track";
  TCanvas *c14a=new TCanvas("c14a",title,200,10,600,480);
  formatHisto(hip2d1,title);
  formatCanvas(c14a);
  hip2d1->Draw();
  for(int i=2;i>-1;i--) hip2d1MC[i]->Draw("hist,e,same");
  hip2d1->Draw("same");
  leg->Draw();
  c14a->GetPad(1)->RedrawAxis();

  c14a->cd(2);
  TH1F *hRip2d1 = hip2d1->Clone("hRip2d1");
  formatRatioHist(hRip2d1);
  hRip2d1->Divide(hip2d1MC[2]);
  hRip2d1->Draw();
  if(savePlots)c14a->SaveAs(Form("plotsNewpp/%s.pdf",title));

  title = "2D IP Significance 1st track";
  TCanvas *c15a=new TCanvas("c15a",title,200,10,600,480);
  formatHisto(hip2dSig1,title);
  formatCanvas(c15a);
  hip2dSig1->Draw();
  for(int i=2;i>-1;i--) hip2dSig1MC[i]->Draw("hist,e,same");
  hip2dSig1->Draw("same");
  leg->Draw();
  c15a->GetPad(1)->RedrawAxis();

  c15a->cd(2);
  TH1F *hRip2dSig1 = hip2dSig1->Clone("hRip2dSig1");
  formatRatioHist(hRip2dSig1);
  hRip2dSig1->Divide(hip2dSig1MC[2]);
  hRip2dSig1->Draw();
  if(savePlots)c15a->SaveAs(Form("plotsNewpp/%s.pdf",title));

  title = "3D1 IP 1st track";
  TCanvas *c16a=new TCanvas("c16a",title,200,10,600,480);
  formatHisto(hip3d1,title);
  formatCanvas(c16a);
  hip3d1->Draw();
  for(int i=2;i>-1;i--) hip3d1MC[i]->Draw("hist,e,same");
  hip3d1->Draw("same");
  leg->Draw();
  c16a->GetPad(1)->RedrawAxis();

  c16a->cd(2);
  TH1F *hRip3d1 = hip3d1->Clone("hRip3d1");
  formatRatioHist(hRip3d1);
  hRip3d1->Divide(hip3d1MC[2]);
  hRip3d1->Draw();
  if(savePlots)c16a->SaveAs(Form("plotsNewpp/%s.pdf",title));

  title = "3D IP Significance 1st track";
  TCanvas *c17a=new TCanvas("c17a",title,200,10,600,480);
  formatHisto(hip3dSig1,title);
  formatCanvas(c17a);
  hip3dSig1->Draw();
  for(int i=2;i>-1;i--) hip3dSig1MC[i]->Draw("hist,e,same");
  hip3dSig1->Draw("same");
  leg->Draw();
  c17a->GetPad(1)->RedrawAxis();

  c17a->cd(2);
  TH1F *hRip3dSig1 = hip3dSig1->Clone("hRip3dSig1");
  formatRatioHist(hRip3dSig1);
  hRip3dSig1->Divide(hip3dSig1MC[2]);
  hRip3dSig1->Draw();
  if(savePlots)c17a->SaveAs(Form("plotsNewpp/%s.pdf",title));



  title = "2D IP 2nd track";
  TCanvas *c14b=new TCanvas("c14b",title,200,10,600,480);
  formatHisto(hip2d2,title);
  formatCanvas(c14b);
  hip2d2->Draw();
  for(int i=2;i>-1;i--) hip2d2MC[i]->Draw("hist,e,same");
  hip2d2->Draw("same");
  leg->Draw();
  c14b->GetPad(1)->RedrawAxis();

  c14b->cd(2);
  TH1F *hRip2d2 = hip2d2->Clone("hRip2d2");
  formatRatioHist(hRip2d2);
  hRip2d2->Divide(hip2d2MC[2]);
  hRip2d2->Draw();
  if(savePlots)c14b->SaveAs(Form("plotsNewpp/%s.pdf",title));

  title = "2D IP Significance 2nd track";
  TCanvas *c15b=new TCanvas("c15b",title,200,10,600,480);
  formatHisto(hip2dSig2,title);
  formatCanvas(c15b);
  hip2dSig2->Draw();
  for(int i=2;i>-1;i--) hip2dSig2MC[i]->Draw("hist,e,same");
  hip2dSig2->Draw("same");
  leg->Draw();
  c15b->GetPad(1)->RedrawAxis();

  c15b->cd(2);
  TH1F *hRip2dSig2 = hip2dSig2->Clone("hRip2dSig2");
  formatRatioHist(hRip2dSig2);
  hRip2dSig2->Divide(hip2dSig2MC[2]);
  hRip2dSig2->Draw();
  if(savePlots)c15b->SaveAs(Form("plotsNewpp/%s.pdf",title));

  title = "3D IP 2nd track";
  TCanvas *c16b=new TCanvas("c16b",title,200,10,600,480);
  formatHisto(hip3d2,title);
  formatCanvas(c16b);
  hip3d2->Draw();
  for(int i=2;i>-1;i--) hip3d2MC[i]->Draw("hist,e,same");
  hip3d2->Draw("same");
  leg->Draw();
  c16b->GetPad(1)->RedrawAxis();

  c16b->cd(2);
  TH1F *hRip3d2 = hip3d2->Clone("hRip3d2");
  formatRatioHist(hRip3d2);
  hRip3d2->Divide(hip3d2MC[2]);
  hRip3d2->Draw();
  if(savePlots)c16b->SaveAs(Form("plotsNewpp/%s.pdf",title));

  title = "3D IP Significance 2nd track";
  TCanvas *c17b=new TCanvas("c17b",title,200,10,600,480);
  formatHisto(hip3dSig2,title);
  formatCanvas(c17b);
  hip3dSig2->Draw();
  for(int i=2;i>-1;i--) hip3dSig2MC[i]->Draw("hist,e,same");
  hip3dSig2->Draw("same");
  leg->Draw();
  c17b->GetPad(1)->RedrawAxis();

  c17b->cd(2);
  TH1F *hRip3dSig2 = hip3dSig2->Clone("hRip3dSig2");
  formatRatioHist(hRip3dSig2);
  hRip3dSig2->Divide(hip3dSig2MC[2]);
  hRip3dSig2->Draw();
  if(savePlots)c17b->SaveAs(Form("plotsNewpp/%s.pdf",title));



  title = "2D IP 3rd track";
  TCanvas *c14c=new TCanvas("c14c",title,200,10,600,480);
  formatHisto(hip2d3,title);
  formatCanvas(c14c);
  hip2d3->Draw();
  for(int i=2;i>-1;i--) hip2d3MC[i]->Draw("hist,e,same");
  hip2d3->Draw("same");
  leg->Draw();
  c14c->GetPad(1)->RedrawAxis();

  c14c->cd(2);
  TH1F *hRip2d3 = hip2d3->Clone("hRip2d3");
  formatRatioHist(hRip2d3);
  hRip2d3->Divide(hip2d3MC[2]);
  hRip2d3->Draw();
  if(savePlots)c14c->SaveAs(Form("plotsNewpp/%s.pdf",title));

  title = "2D IP Significance 3rd track";
  TCanvas *c15c=new TCanvas("c15c",title,200,10,600,480);
  formatHisto(hip2dSig3,title);
  formatCanvas(c15c);
  hip2dSig3->Draw();
  for(int i=2;i>-1;i--) hip2dSig3MC[i]->Draw("hist,e,same");
  hip2dSig3->Draw("same");
  leg->Draw();
  c15c->GetPad(1)->RedrawAxis();

  c15c->cd(2);
  TH1F *hRip2dSig3 = hip2dSig3->Clone("hRip2dSig3");
  formatRatioHist(hRip2dSig3);
  hRip2dSig3->Divide(hip2dSig3MC[2]);
  hRip2dSig3->Draw();
  if(savePlots)c15c->SaveAs(Form("plotsNewpp/%s.pdf",title));

  title = "3D IP 3rd track";
  TCanvas *c16c=new TCanvas("c16c",title,200,10,600,480);
  formatHisto(hip3d3,title);
  formatCanvas(c16c);
  hip3d3->Draw();
  for(int i=2;i>-1;i--) hip3d3MC[i]->Draw("hist,e,same");
  hip3d3->Draw("same");
  leg->Draw();
  c16c->GetPad(1)->RedrawAxis();

  c16c->cd(2);
  TH1F *hRip3d3 = hip3d3->Clone("hRip3d3");
  formatRatioHist(hRip3d3);
  hRip3d3->Divide(hip3d3MC[2]);
  hRip3d3->Draw();
  if(savePlots)c16c->SaveAs(Form("plotsNewpp/%s.pdf",title));

  title = "3D IP Significance 3rd track";
  TCanvas *c17c=new TCanvas("c17c",title,200,10,600,480);
  formatHisto(hip3dSig3,title);
  formatCanvas(c17c);
  hip3dSig3->Draw();
  for(int i=2;i>-1;i--) hip3dSig3MC[i]->Draw("hist,e,same");
  hip3dSig3->Draw("same");
  leg->Draw();
  c17c->GetPad(1)->RedrawAxis();

  c17c->cd(2);
  TH1F *hRip3dSig3 = hip3dSig3->Clone("hRip3dSig3");
  formatRatioHist(hRip3dSig3);
  hRip3dSig3->Divide(hip3dSig3MC[2]);
  hRip3dSig3->Draw();
  if(savePlots)c17c->SaveAs(Form("plotsNewpp/%s.pdf",title));



  title = "Distance to jet axis";
  TCanvas *c18=new TCanvas("c18",title,200,10,600,480);
  formatHisto(hipDist2Jet,title);
  formatCanvas(c18);
  hipDist2Jet->Draw();
  for(int i=2;i>-1;i--) hipDist2JetMC[i]->Draw("hist,e,same");
  hipDist2Jet->Draw("same");
  leg->Draw();
  c18->GetPad(1)->RedrawAxis();

  c18->cd(2);
  TH1F *hRipDist2Jet = hipDist2Jet->Clone("hRipDist2Jet");
  formatRatioHist(hRipDist2Jet);
  hRipDist2Jet->Divide(hipDist2JetMC[2]);
  hRipDist2Jet->Draw();
  if(savePlots)c18->SaveAs(Form("plotsNewpp/%s.pdf",title));

  title = "hipDist2JetSig";
  TCanvas *c19=new TCanvas("c19",title,200,10,600,480);
  formatHisto(hipDist2JetSig,title);
  formatCanvas(c19);
  hipDist2JetSig->Draw();
  for(int i=2;i>-1;i--) hipDist2JetSigMC[i]->Draw("hist,e,same");
  hipDist2JetSig->Draw("same");
  leg->Draw();
  c19->GetPad(1)->RedrawAxis();

  c19->cd(2);
  TH1F *hRipDist2JetSig = hipDist2JetSig->Clone("hRipDist2JetSig");
  formatRatioHist(hRipDist2JetSig);
  hRipDist2JetSig->Divide(hipDist2JetSigMC[2]);
  hRipDist2JetSig->Draw();
  if(savePlots)c19->SaveAs(Form("plotsNewpp/%s.pdf",title));

  title = "Decay length";
  TCanvas *c20=new TCanvas("c20",title,200,10,600,480);
  formatHisto(hipClosest2Jet,title);
  formatCanvas(c20);
  hipClosest2Jet->Draw();
  for(int i=2;i>-1;i--) hipClosest2JetMC[i]->Draw("hist,e,same");
  hipClosest2Jet->Draw("same");
  leg->Draw();
  c20->GetPad(1)->RedrawAxis();

  c20->cd(2);
  TH1F *hRipClosest2Jet = hipClosest2Jet->Clone("hRipClosest2Jet");
  formatRatioHist(hRipClosest2Jet);
  hRipClosest2Jet->Divide(hipClosest2JetMC[2]);
  hRipClosest2Jet->Draw();
  if(savePlots)c20->SaveAs(Form("plotsNewpp/%s.pdf",title));
  }
  /*
  title = "muon pTrel";
  TCanvas *c21=new TCanvas("c21",title,200,10,600,480);
  formatHisto(hmuptrel,title);
  formatCanvas(c21);
  hmuptrel->Draw();
  for(int i=2;i>-1;i--) hmuptrelMC[i]->Draw("hist,e,same");
  hmuptrel->Draw("same");
  leg->Draw();
  c21->GetPad(1)->RedrawAxis();

  c21->cd(2);
  TH1F *hRmuptrel = hmuptrel->Clone("hmuptrel");
  formatRatioHist(hRmuptrel);
  hRmuptrel->Divide(hmuptrelMC[2]);
  hRmuptrel->Draw();

  title = "muon pTrel (requiring 2-track SV)";
  TCanvas *c21a=new TCanvas("c21a",title,200,10,600,480);
  formatHisto(hmuptrelSV2,title);
  formatCanvas(c21a);
  hmuptrelSV2->Draw();
  for(int i=2;i>-1;i--) hmuptrelSV2MC[i]->Draw("hist,e,same");
  hmuptrelSV2->Draw("same");
  leg->Draw();
  c21a->GetPad(1)->RedrawAxis();

  c21a->cd(2);
  TH1F *hRmuptrelSV2 = hmuptrelSV2->Clone("hmuptrelSV2");
  formatRatioHist(hRmuptrelSV2);
  hRmuptrelSV2->Divide(hmuptrelSV2MC[2]);
  hRmuptrelSV2->Draw();

  title = "muon pTrel (requiring 3-track SV)";
  TCanvas *c21b=new TCanvas("c21b",title,200,10,600,480);
  formatHisto(hmuptrelSV3,title);
  formatCanvas(c21b);
  hmuptrelSV3->Draw();
  for(int i=2;i>-1;i--) hmuptrelSV3MC[i]->Draw("hist,e,same");
  hmuptrelSV3->Draw("same");
  leg->Draw();
  c21b->GetPad(1)->RedrawAxis();

  c21b->cd(2);
  TH1F *hRmuptrelSV3 = hmuptrelSV3->Clone("hmuptrelSV3");
  formatRatioHist(hRmuptrelSV3);
  hRmuptrelSV3->Divide(hmuptrelSV3MC[2]);
  hRmuptrelSV3->Draw();
  */


}
