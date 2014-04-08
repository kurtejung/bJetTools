#include <iostream>
#include <stdio.h>

#include <TRandom2.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>

#include "utilities.h"
#include "bayesianUnfold.h"
#include "prior.h"
#include "generateSmearingMatrix.h"

using namespace std;


//==============================================================================
// Unfolding Ying Lu 08 07 11
// Update Yen-Jie Lee 06.22.12
//==============================================================================

void Unfold2(int algo= 3,bool useSpectraFromFile=0, bool useMatrixFromFile=0, int doToy = 0, int isMC = 0,char *spectraFileName = (char*)"pbpb_spectra_ak3PF.root",double recoJetPtCut = 30., int nBayesianIter = 6, int doBjets=1, int doTrigEffCorr=0, int priorGenOrBin=0, int doSkew=0, int doJESshift=0, bool doResSmear=false, int centBin=0) 
{

  bool doParameterizedMatrix = 1;
  double genJetPtCut=20;
  
  gStyle->SetErrorX(0.);
  gStyle->SetPaintTextFormat("3.2f");
  gStyle->SetOptLogz(1);
  gStyle->SetPadRightMargin(0.13);	
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetNdivisions(505,"x");

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  
  // input files
  char *fileNamePP_mc = NULL;
  char *fileNamePbPb_mc = NULL;
  char *fileNamePP_data = NULL;
  char *fileNamePbPb_data = NULL;
  
  int cBinLo=0;
  int cBinHi=100;

  float etalo=-2.;
  float etahi=2.;

  float ppetalo = etalo+0.465;
  float ppetahi = etahi+0.465;

  if(centBin==1) cBinHi=4;
  if(centBin==2) {
    cBinLo=4;
    cBinHi=12;
  }
  if(centBin==3) cBinLo=12;
  if(centBin==4){
    cBinLo=12;
    cBinHi=20;
  }
  if(centBin==5) cBinLo=20;

  // pp file needs replacing

  if(doBjets)fileNamePP_data = (char*)Form("output/ppMC_akPu3PF_bJet_TEST_matchBin0_gsp0_discrSSVHEGT2_eta%.1fTo%.1f.root",ppetalo,ppetahi);
  else fileNamePP_data = (char*)"output/ppMC_ak3PF_matchBinToYaxian0_gsp0_discrGT2_eta-0.5To1.5.root";

  if(doBjets) fileNamePbPb_data = (char*)Form("output/NewFormatV16_FixedTrgMergingFinal_akPu3PF_consistentEta_FirstDataHalf_fixBin_bFractionMCTemplate_pPbpp1_gsp0_jetptcut30_SSVHEat2.0FixCL0_bin_0_100_eta_%.1f_%.1f.root",etalo,etahi);
  else fileNamePbPb_data = (char*)"histos/PUForest_akpu3pf_etaLT1_FirstHalf_noGplus_forUnfolding.root";
 
  if(doBjets) fileNamePP_mc = (char*)"histos/ppMC_ppReco_akPu3PF_QCDjetTrig_etashift_Fix2Sample_MCWeightFinalWithVz_noTrgSelection_Full.root";
  else fileNamePP_mc = (char*)"histos/ppMC_ppReco_akPu3PF_QCDjetTrig_etashift_Fix2Sample_MCWeightFinalWithVz_noTrgSelection_Full.root";
  
  if(doBjets)fileNamePbPb_mc = (char*) "histos/pPbMC_ppReco_akPu3PF_QCDjetTrig_etashift_noTrgSelection_MCWeightFinalWithSmoothedVzCentWeight_vZcut_noHBHE_8.root";
  else fileNamePbPb_mc = (char*) "histos/pPbMC_ppReco_akPu3PF_QCDjetTrig_etashift_noTrgSelection_MCWeightFinalWithSmoothedVzCentWeight_vZcut_noHBHE_8.root";


  // grab ntuples
  TFile *infPbPb_mc = new TFile(fileNamePbPb_mc);
  TFile *infPP_mc = new TFile(fileNamePP_mc);
  

  string bJetString = "Inc";
  if(doBjets) bJetString = "bJets";
  
  string jetType = "akPu3PF";

  string skewString = "";
  if(doSkew==1)skewString = "_skewPos";
  else if(doSkew==2)skewString = "_skewNeg";

  string resSmearString = "";
  if(doResSmear) resSmearString = "_resSmear";

  string shiftString = "";
  if(doJESshift==1) shiftString = "_shiftJESdown";
  if(doJESshift==2) shiftString = "_shiftJESup";

  // Output file
  TFile *pbpb_Unfo;
  if (isMC) pbpb_Unfo = new TFile(Form("pPb_Unfo_%s_%s_eta_%.1fTo%.1f_MC_v2.root",jetType.c_str(),algoName[algo],etalo,etahi),"RECREATE");
  else pbpb_Unfo  = new TFile(Form("pPb_Unfo_%s_%s_noGplus_FirstHalf_Converged_jtpt%.0f_%s_clo%d_chi%d_v8%s%s%s_eta_%.1fTo%.1f_.root",jetType.c_str(),algoName[algo],recoJetPtCut,bJetString.c_str(),cBinLo,cBinHi,skewString.c_str(),resSmearString.c_str(),shiftString.c_str(),etalo,etahi),"RECREATE");

  // Histograms used by RooUnfold
  UnfoldingHistos *uhist[nbins_cent+1];
		
  // Initialize Histograms   
	
  for (int i=0;i<=nbins_cent;i++) uhist[i] = new UnfoldingHistos(i);
	

  if (isMC) cout<<"This is a MC closure test"<<endl;
  else cout<< "This is a data analysis"<<endl;    		     
	     	
  // Setup jet data branches, basically the jet tree branches are assigned to this object when we loop over the events
	
  JetDatapPb *dataPbPb   = new JetDatapPb(fileNamePbPb_mc,(char*)"nt"); // pPb data	
  JetDataPP *dataPP = new JetDataPP(fileNamePP_mc,(char*)"nt");	// pp data
	
  TFile *fSpectra(0);		
  if (useSpectraFromFile||useMatrixFromFile){
    fSpectra = new TFile(spectraFileName,"read");
  }
  
  // Come back to the output file dir
  pbpb_Unfo->cd();

  cout <<"MC..."<<endl;	
  
  //Now with added trigger efficiency corrections for both pPb and pp!

  TRandom *rand =  new TRandom();

  TF1 *fResPbPb =  new TF1("fResPbPb","[0]+[1]/x+[2]/x/x",0,1000);
  if(centBin==0)fResPbPb->SetParameters( 3.09252e-02,1.40151e+01,-3.34065e+02);
  else if(centBin==1)fResPbPb->SetParameters(  3.57669e-02, 1.38110e+01, -3.43024e+02);
  else if(centBin==2)fResPbPb->SetParameters(   4.74219e-02, 9.44555e+00, -1.72944e+02);
  else if(centBin==3)fResPbPb->SetParameters(  5.00618e-02, 6.67202e+00, -5.57270e+01);
  else if(centBin==4)fResPbPb->SetParameters(  5.00618e-02, 6.67202e+00, -5.57270e+01);
  else if(centBin==5)fResPbPb->SetParameters(  5.00618e-02, 6.67202e+00, -5.57270e+01);
  if(centBin>3 && doResSmear) cout<<" DONT FORGET TO FIX RESOLUTION "<<endl; /// not actually used anymore

		
  // Fill PbPb MC   
  if (!useMatrixFromFile) {
    for (Long64_t jentry2=0; jentry2<dataPbPb->tJet->GetEntries();jentry2++) {
      dataPbPb->tJet->GetEntry(jentry2);
      
      // change when we switch to centrality binning
      int cBin = 0;
      
      double jtpt = dataPbPb->jtpt;      
      if(doBjets) jtpt = dataPbPb->jtpt;
      
	
      if ( dataPbPb->refpt  < 0 ) continue;
      if ( dataPbPb->jteta  > etahi || dataPbPb->jteta < etalo ) continue;
      if ( dataPbPb->refpt < 20. ) continue;
      if (dataPbPb->jtpt < recoJetPtCut) continue;
      //if ( dataPbPb->refpt < genJetPtCut) continue;
      //if ( !dataPbPb->gPlus ) continue;
      if (doBjets && fabs(dataPbPb->refparton_flavorForB)!=5) continue;
      if(dataPbPb->refpt/dataPbPb->pthat > 1.0) continue;
      //if (doBjets&& dataPbPb->discr_ssvHighEff<2) continue;
      //if (!doTrigEffCorr && dataPbPb->isTrig <1) continue;
      if(dataPbPb->subid != 0) continue;
      //if ( dataPbPb->isTrig <1) continue;
      //if( dataPbPb->bin < cBinLo || dataPbPb->bin >= cBinHi) continue;
      //if( dataPbPb->jtpt/dataPbPb->rawpt>1.5) continue;

      if(doResSmear){
	double sigSmear = fResPbPb->Eval(dataPbPb->refpt) * sqrt(0.1*0.1 + 2*0.1);
	jtpt *= (1.+rand->Gaus(0,sigSmear));
	if(jtpt< 0) cout<<" LESS THAN ZERO "<<endl;
      }

      if(doJESshift==1) jtpt*=(1.-sqrt(8.)*0.01);
      else if(doJESshift==2) jtpt*=(1.+sqrt(8.)*0.01);


      //if(jtpt < recoJetPtCut) continue;

      //if(!doBjets)if(dataPbPb->refpt < 50 && dataPbPb->jtptA>120) continue;
      //if(doBjets)if(dataPbPb->refpt < 50 && dataPbPb->jtptB>120) continue;

      float weight = dataPbPb->weight;

      if(doTrigEffCorr){
	for(int iBin = 0; iBin<nbins_rec; iBin++){
	  if(jtpt > boundaries_rec[iBin] && jtpt < boundaries_rec[iBin+1]){
	    if(doBjets) weight/= trigEffbJet[iBin];
	    else weight/= trigEffInc[iBin];
	  }							  
	}
      }
      if (!isMC||jentry2 % 2 == 1) {
	uhist[cBin]-> hMatrix->Fill(dataPbPb->refpt,jtpt,dataPbPb->weight);
      }	  
      // if (jentry2 % 1 == 0) {
	uhist[cBin]-> hGen->Fill(dataPbPb->refpt,dataPbPb->weight);   
	//uhist[cBin]-> hMeas->Fill(jtpt,dataPbPb->weight);
	// } 
    }
    for(int jj=1; jj<uhist[0]->hGen->GetNbinsX(); jj++){
      // cout << "bin " << jj << " content: "<< uhist[0]->hGen->GetBinContent(jj) << endl;
    }

    //pp will just fill the last index of the centrality array
    TF1 *fResPP =  new TF1("fResPP","[0]+[1]/x+[2]/x/x",0,1000);
    fResPP->SetParameters(4.62986e-02,6.78514e+00,-1.44082e+02);


    // fill pp MC
    for (Long64_t jentry2=0; jentry2<dataPP->tJet->GetEntries();jentry2++) {
      dataPP->tJet->GetEntry(jentry2);

      //uhist[nbins_cent]-> hMatrix->Fill(dataPP->refpt,dataPP->jtpt,dataPP->weight);

      if(dataPP->refpt>0 && dataPP->jteta<etahi+0.465 && dataPP->jteta>etalo+0.465 && ((doBjets && fabs(dataPP->refparton_flavorForB)!=5) || !doBjets)){
	
      }
	 if ( dataPP->refpt<genJetPtCut) continue;
	 //if ( dataPP->rawpt<20) continue;
      if ( dataPP->jteta  > etahi+0.465 || dataPP->jteta < etalo+0.465 ) continue;
      //if ( !dataPP->gPlus ) continue;
      ///if ( dataPP->refpt<0) dataPP->refpt=0;
      if ( doBjets && fabs(dataPP->refparton_flavorForB)!=5) continue;
      //if ( doBjets && dataPP->discr_ssvHighEff<2) continue;
      //if( dataPP->jtpt/dataPP->rawpt>1.5) continue;  // protect against bad JEC fit

      if(doResSmear){
	double sigSmear = fResPP->Eval(dataPP->refpt) * sqrt(0.1*0.1 + 2*0.1);
	dataPP->jtpt *= (1.+rand->Gaus(0,sigSmear));
	if(dataPP->jtpt< 0) cout<<" LESS THAN ZERO "<<endl;
      }

      if(doJESshift==1) dataPP->jtpt*=(1.-0.02);
      else if(doJESshift==2) dataPP->jtpt*=(1.+0.02);

      if(dataPP->jtpt < recoJetPtCut) continue;

      float jtpt = dataPP->jtpt;
      float weight = dataPP->weight;
      if(doTrigEffCorr){
	for(int iBin = 0; iBin<nbins_rec; iBin++){
	  if(jtpt > boundaries_rec[iBin] && jtpt < boundaries_rec[iBin+1]){
	    if(doBjets) weight/= trigEffppbJet[iBin];
	    else weight/= trigEffppInc[iBin];
	  }							  
	}
      }
      
      if (!isMC||jentry2 % 2 == 1) {
       	uhist[nbins_cent]-> hMatrix->Fill(dataPP->refpt,dataPP->jtpt,dataPP->weight);
      }	  
      //if (jentry2 % 2 == 0) {
	uhist[nbins_cent]-> hGen->Fill(dataPP->refpt,dataPP->weight);   
	uhist[nbins_cent]-> hMeas->Fill(dataPP->jtpt,dataPP->weight); 
	// }        
    }
  }


  if (isMC==0) {
    // Use measured histogram from Matt & Kurt's file
	   
    // PbPb file:

    TFile *infMatt = new TFile(fileNamePbPb_data);
    TH1F *hMattPbPb = NULL;
    TH1F *hTagEffPbPb = NULL;

    /*TFile *backup = new TFile("histos/pPbdata_ppreco_akPu3PF_jetTrig_etashift_FinalTC_WxieForest4_AddSvtxInfo_NoEvtSelection_NoResCorr_FullMerge.root");
    TTree *nt = (TTree*)backup->Get("nt");
    hMattPbPb = new TH1F("hIncJetsData","",nbins_rec,boundaries_rec);
    hMattPbPb->Sumw2();
    nt->Draw("jtpt>>hIncJetsData","weight*(abs(jteta)<1 && evtSelection && abs(vz)<15 && run<211300 && rawpt>=22)");
    // for(int i=1; i<=hMattPbPb->GetNbinsX(); i++){
    //  hMattPbPb->SetBinError(i,hMattPbPb->GetBinError(i)/hMattPbPb->GetBinWidth(i));
    //  hMattPbPb->SetBinContent(i,hMattPbPb->GetBinContent(i)/hMattPbPb->GetBinWidth(i));
   //}
   hMattPbPb->Scale(1./(0.85*20.7*1.e6));*/
    
    if(doBjets){
      hMattPbPb = (TH1F*) infMatt->Get("hRawBData");
      hTagEffPbPb = (TH1F*) infMatt->Get("hBEfficiencyMC");
    }
    

    else hMattPbPb = (TH1F*) infMatt->Get("hIncJetsData");
    //divideBinWidth(hMattPbPb);
    

           
    // Need to match the binning carefully, please double check whenever you change the binning
    for (int i=1;i<=hMattPbPb->GetNbinsX();i++)
      {
     	float binContent = hMattPbPb->GetBinContent(i);  
	float binError = hMattPbPb->GetBinError(i); 

	if(doBjets){
	  float tagEff =hTagEffPbPb->GetBinContent(i);
	  float tagEffErr =     hTagEffPbPb->GetBinError(i);   

	  cout << "tageff, bin " << i << " : "<< tagEff << endl;
	  
	  if(tagEff>0){
	    // careful of the order here!
	    binError=binContent/tagEff *sqrt(tagEffErr*tagEffErr/tagEff/tagEff + binError*binError/binContent/binContent);
	    binContent /= tagEff;
	  }
	  else cout<<"pPb TAGEFF = 0"<<endl;	  
	}

	cout << "pPb bin content " << i << " "<< binContent/(20.7E9*2110E-3*6.9) << endl;

	float binCenter = hMattPbPb->GetBinCenter(i);  
	//if(binCenter - hMattPbPb->GetBinWidth(i)/2.  < 22.) continue;
	
	int ibin=uhist[0]->hMeas->FindBin(binCenter);
	
	/*for(ibin=1;ibin<=uhist[0]->hMeas->GetNbinsX();ibin++){
	 
	  float testLowEdge = uhist[0]->hMeas->GetBinLowEdge(ibin);
	  float testBinWidth = uhist[0]->hMeas->GetBinWidth(ibin);
	  if(binCenter>testLowEdge && binCenter < testLowEdge+testBinWidth) break;
	  	}*/
	uhist[0]->hMeas->SetBinContent(ibin,0); //reset to make sure nothing from MC is inside
	uhist[0]->hMeas->SetBinError(ibin,0);
	
	if(doTrigEffCorr){
	  float trigEff = 0;
	  if(doBjets) trigEff = trigEffbJet[ibin-1];
	  else  trigEff = trigEffInc[ibin-1];

	  cout<<" binCenter = "<<binCenter<<" trigEff "<<trigEff<<endl;

	  if(trigEff>0){
	    // careful of the order here!
	    binContent /= trigEff;
	    binError /= trigEff;
	  }
	  else cout<<"pPb TRIGEFF = 0"<<endl;	  
	}   
	
	//  cout <<"get bins " << testLowEdge<<" content =" <<binContent<<endl ;
        uhist[0]->hMeas->SetBinContent(ibin,binContent);  
        uhist[0]->hMeas->SetBinError(ibin,binError);  
      }

    // pp file:
    TFile *infMattPP = new TFile(fileNamePP_data);
    TFile *MCtruth = new TFile(fileNamePP_data);
    TH1F *hMattPP = NULL;
    TH1F *hTagEffPP = NULL;
    if(doBjets){
      hMattPP = (TH1F*) MCtruth->Get("hRawBMC"); //MODIFIED BECAUSE pA pp reference is MC!
      hTagEffPP = (TH1F*) infMattPP->Get("hBEfficiencyMC");
    }
    //else hMattPP = (TH1F*) infMattPP->Get("hpp");
    else hMattPP = (TH1F*) infMattPP->Get("hIncJetsMC");
    //divideBinWidth(hMattPP);	   
     for (int i=1;i<=hMattPP->GetNbinsX();i++)
      {
      	float binContent = hMattPP->GetBinContent(i);  
	float binError = hMattPP->GetBinError(i);  
	
	/*if(doBjets){
	  float tagEff =hTagEffPP->GetBinContent(i);
	  float tagEffErr =     hTagEffPP->GetBinError(i);   
	  if(tagEff>0){
	    // careful of the order here!
	    binError=binContent/tagEff *sqrt(tagEffErr*tagEffErr/tagEff/tagEff + binError*binError/binContent/binContent);
	    binContent /= tagEff;
	  }
	  else cout<<"pp TAGEFF = 0"<<endl;	  
	  }*/
	
     	float binCenter = hMattPP->GetBinCenter(i);  
	if(binCenter - hMattPP->GetBinWidth(i)/2.  < recoJetPtCut) continue;

	int ibin=0;

	for(ibin=1;ibin<=uhist[nbins_cent]->hMeas->GetNbinsX();ibin++){
	  float testLowEdge = uhist[nbins_cent]->hMeas->GetBinLowEdge(ibin);
	  float testBinWidth = uhist[nbins_cent]->hMeas->GetBinWidth(ibin);
	  if(binCenter>testLowEdge && binCenter < testLowEdge+testBinWidth) break;
	}
	//	uhist[nbins_cent]->hMeas->SetBinContent(ibin,binContent);  
	//	uhist[nbins_cent]->hMeas->SetBinError(ibin,binError);  

	cout << "pp bin content " << i << " "<< uhist[nbins_cent]->hGen->GetBinContent(ibin)/70. << endl;
	cout << "pp bin content from file " << i << " "<< binContent/70. << endl;
	/*
	cout<<" i "<<i<<endl;
	int newBin = i+uhist[nbins_cent]->hMeas->FindBin(61)-1;

	cout<<" newBin "<<newBin<<endl;
	cout<<"hMattPP->GetBinCenter(i) "<<hMattPP->GetBinCenter(i)<<endl;
	cout<<"uhist[nbins_cent]->hMeas->GetBinCenter(newBin) "<<uhist[nbins_cent]->hMeas->GetBinCenter(newBin)<<endl;
	*/
      }
  }



    TLatex *tpp =  new TLatex(0.5,0.9,"pp");
    tpp->SetTextFont(42);
    tpp->SetTextSize(0.042);
    tpp->SetNDC();

    TLatex *tpbpb =  NULL;
    if(centBin==0)tpbpb = new TLatex(0.5,0.9,"pPb, 0-100%");
    if(centBin==1)tpbpb = new TLatex(0.5,0.9,"pPb, 0-10%");
    if(centBin==2)tpbpb = new TLatex(0.5,0.9,"pPb, 10-30%");
    if(centBin==3)tpbpb = new TLatex(0.5,0.9,"pPb, 30-100%");
    if(centBin==4)tpbpb = new TLatex(0.5,0.9,"pPb, 30-50%");
    if(centBin==5)tpbpb = new TLatex(0.5,0.9,"pPb, 50-100%");
    tpbpb->SetTextFont(42);
    tpbpb->SetTextSize(0.042);
    tpbpb->SetNDC();

   //=========================================Response Matrix========================================================= 

  cout <<"Response Matrix..."<<endl;
	
  TCanvas * cMatrix = new TCanvas("cMatrix","Matrix",1200,600);
  cMatrix->Divide(2,1);
  TH2F * hResponse[nbins_cent+1];
  for (int i=0;i<=nbins_cent;i++){
    cMatrix->cd(i+1);
    if(doParameterizedMatrix){
      uhist[i]->hMatrix = (TH2F*)generateSmearingMatrix(i,uhist[i]->hGen,uhist[i]->hMeas,genJetPtCut,recoJetPtCut,etalo,1)->Clone(Form("hMatrix_%d",i));
      //cMatrix->cd(i); 
      uhist[i]->hMatrix->SetMaximum(1);
      uhist[i]->hMatrix->SetMinimum(1E-7);
      uhist[i]->hMatrix->Draw("colz");
    }
    else if (!useMatrixFromFile) {
     
      TF1 *f = new TF1("f","[0]*pow(x+[2],[1])");
      f->SetParameters(1e10,-8.8,40);
      for (int y=1;y<=uhist[i]->hMatrix->GetNbinsY();y++) {
	double sum=0;
	for (int x=1;x<=uhist[i]->hMatrix->GetNbinsX();x++) {
	  if (uhist[i]->hMatrix->GetBinContent(x,y)<=1*uhist[i]->hMatrix->GetBinError(x,y)) {
	    uhist[i]->hMatrix->SetBinContent(x,y,0);
	    uhist[i]->hMatrix->SetBinError(x,y,0);
	  }
	  sum+=uhist[i]->hMatrix->GetBinContent(x,y);
	}

	for (int x=1;x<=uhist[i]->hMatrix->GetNbinsX();x++) {	   
	  double ratio = 1;
	  uhist[i]->hMatrix->SetBinContent(x,y,uhist[i]->hMatrix->GetBinContent(x,y)*ratio);
	  uhist[i]->hMatrix->SetBinError(x,y,uhist[i]->hMatrix->GetBinError(x,y)*ratio);
	}
      }
    }
    uhist[i]->hResponse = (TH2F*)uhist[i]->hMatrix->Clone(Form("hResponse_cent%d",i));
    TH1F *hProj = (TH1F*)uhist[i]->hResponse->ProjectionY(Form("hProj_cent%d",i));


    for (int y=1;y<=uhist[i]->hResponse->GetNbinsY();y++) {
      double sum=0;
      for (int x=1;x<=uhist[i]->hResponse->GetNbinsX();x++) {
	if (uhist[i]->hResponse->GetBinContent(x,y)<=1*uhist[i]->hResponse->GetBinError(x,y)) {
	  uhist[i]->hResponse->SetBinContent(x,y,0);
	  uhist[i]->hResponse->SetBinError(x,y,0);
	}
	sum+=uhist[i]->hResponse->GetBinContent(x,y);
      }

      for (int x=1;x<=uhist[i]->hResponse->GetNbinsX();x++) {  	
	if (sum==0) continue;
	double ratio =1.;
	//	if (uhist[i]->hMeas->GetBinContent(y)==0) ratio = 1e-100/sum;
	//         else ratio = uhist[i]->hMeas->GetBinContent(y)/sum;
	ratio = 1./sum ;
	if (hProj->GetBinContent(y)==0) ratio = 1e-100/sum;
	else ratio = hProj->GetBinContent(y)/sum;
	uhist[i]->hResponse->SetBinContent(x,y,uhist[i]->hResponse->GetBinContent(x,y)*ratio);
        uhist[i]->hResponse->SetBinError(x,y,uhist[i]->hResponse->GetBinError(x,y)*ratio);
      }
    }

    uhist[i]->hResponseNorm = (TH2F*)uhist[i]->hMatrix->Clone(Form("hResponseNorm_cent%d",i));
    for (int x=1;x<=uhist[i]->hResponseNorm->GetNbinsX();x++) {
      double sum=0;
      for (int y=1;y<=uhist[i]->hResponseNorm->GetNbinsY();y++) {
	if (uhist[i]->hResponseNorm->GetBinContent(x,y)<=0*uhist[i]->hResponseNorm->GetBinError(x,y)) {
	  uhist[i]->hResponseNorm->SetBinContent(x,y,0);
	  uhist[i]->hResponseNorm->SetBinError(x,y,0);
	}
	sum+=uhist[i]->hResponseNorm->GetBinContent(x,y);
      }

      for (int y=1;y<=uhist[i]->hResponseNorm->GetNbinsY();y++) {  	
	if (sum==0) continue;
	double ratio = 1./sum;
	uhist[i]->hResponseNorm->SetBinContent(x,y,uhist[i]->hResponseNorm->GetBinContent(x,y)*ratio);
	uhist[i]->hResponseNorm->SetBinError(x,y,uhist[i]->hResponseNorm->GetBinError(x,y)*ratio);
      }
    }

    if (!useMatrixFromFile) uhist[i]->hMatrixFit = uhist[i]->hMatrix;
    uhist[i]->hMatrixFit->SetName(Form("hMatrixFit_cent%d",i));
  } 
  /*
    TCanvas * cMatrix = new TCanvas("cMatrix","Matrix",800,400);
    cMatrix->Divide(2,1);
    //  cMatrix->cd(1);
    */    for (int i=0;i<=nbins_cent;i++){
    cMatrix->cd(i+1); 
    //     cMatrix->cd(i+1);		
    uhist[i]->hResponse->SetTitleOffset(1.4, "Y");
    uhist[i]->hResponse->SetTitleOffset(1.2, "X");
    //  if(isMC)uhist[i]->hResponse->SetMinimum(1.e-8);
    //  else
    uhist[i]->hResponse->SetMinimum(1.e-10);
    //uhist[i]->hResponse->DrawCopy("colz");
    /*    uhist[i]->hMatrix->SetTitleOffset(1.4, "Y");
	  uhist[i]->hMatrix->SetTitleOffset(1.2, "X");
	  uhist[i]->hMatrix->Draw("colz");
		
	  cMatrix->cd(i+1)->Modified();
	  cMatrix->cd(i+1);
    */  }

  // cMatrix->SetSelected(cMatrix);
  cMatrix->Update();

  pbpb_Unfo->cd();
	
  cout << "==================================== UNFOLD ===================================" << endl;
	
  //char chmet[100]; 
	
  // ======================= Reconstructed pp and PbPb spectra =========================================================
  TCanvas * cPbPb = new TCanvas("cPbPb","Comparison",1200,600);
  cPbPb->Divide(2,1); 
  cPbPb->cd(1);
	
  TH1F *hRecoBW[nbins_cent+1], *hRecoBinByBinBW[nbins_cent+1], *hMeasBW[nbins_cent+1], *hGenBW[nbins_cent+1], *hReproducedBW[nbins_cent+1];


  TCanvas * cPbPbMeas = new TCanvas("cPbPbMeas","Measurement",1200,600);
  cPbPbMeas->Divide(2,1); 
  cPbPbMeas->cd(1);
	

  for (int i=0;i<=nbins_cent;i++) {
    cPbPb->cd(i+1)->SetLogy();   
    // Do Bin-by-bin
    TH1F *hBinByBinCorRaw = (TH1F*)uhist[i]->hResponse->ProjectionY(); hBinByBinCorRaw->Sumw2();
    TH1F *hMCGen           = (TH1F*)uhist[i]->hResponse->ProjectionX(); hMCGen->Sumw2(); // gen 

    for(int jj=1; jj<=hBinByBinCorRaw->GetNbinsX(); jj++){
      //cout << "Y, bin " << jj << " " << hBinByBinCorRaw->GetBinContent(jj) << endl;
    }
     for(int jj=1; jj<=hMCGen->GetNbinsX(); jj++){
       //cout << "X, bin " << jj << " " << hMCGen->GetBinContent(jj) << endl;
     }
    hBinByBinCorRaw->Divide(hMCGen);
    TF1 *f = new TF1("f","[0]+[1]*x");

    hBinByBinCorRaw->Fit("f","LL ","",recoJetPtCut,850);
    TH1F* hBinByBinCor = (TH1F*)hBinByBinCorRaw->Clone();//functionHist(f,hBinByBinCorRaw,Form("hBinByBinCor_cent%d",i));
    uhist[i]->hRecoBinByBin = (TH1F*) uhist[i]->hMeas->Clone(Form("hRecoBinByBin_cent%d",i));
    uhist[i]->hRecoBinByBin->Divide(hBinByBinCor);
		
    // Do unfolding
    //if (isMC) uhist[i]->hMeas = (TH1F*)uhist[i]->hMatrix->ProjectionY()->Clone(Form("hMeas_cent%d",i));
    TH1F *hPrior;//=(TH1F*) functionHist(fPow,uhist[i]->hMeas,Form("hPrior_cent%d",i));
    hPrior = (TH1F*)uhist[i]->hGen->Clone(Form("hPrior_cent%d",i));
    //if(priorGenOrBin) hPrior=(TH1F*)uhist[i]->hRecoBinByBin->Clone(Form("hPrior_cent%d",i));
    //else hPrior=(TH1F*)uhist[i]->hGen->Clone(Form("hPrior_cent%d",i));
    //hPrior = (TH1F*)uhist[i]->hMeas->Clone(Form("hPrior_cent%d",i));

    //TF1 *fSkew = new TF1("fShift","0.2/100.*x+0.7",50,600);
    if(doSkew){
      TF1 *fSkew = NULL;
      if(doSkew==1)fSkew = new TF1("fShift","0.5/100.*x+0.25",40,400);
      else fSkew = new TF1("fShift","-0.5/100.*x+1.75",40,400);
      for(int ib=0;ib<hPrior->GetNbinsX();ib++){
	hPrior->SetBinContent(ib+1,hPrior->GetBinContent(ib+1)*fSkew->Eval(hPrior->GetBinCenter(ib+1)));
      }
    }


    //hPrior=(TH1F*)uhist[i]->hRecoBinByBin->Clone("hPrior");
    removeZero(hPrior);
    for(int jj=1; jj<=hPrior->GetNbinsX(); jj++){
      //cout << "hprior bin : " << jj << " bin center: " << hPrior->GetBinCenter(jj) << " content: " << hPrior->GetBinContent(jj) << endl;
    }

    // what's in yen-jie code
    //prior myPrior(uhist[i]->hMatrixFit,uhist[i]->hMeas,0);
    prior myPrior(uhist[i]->hMatrixFit,hPrior,0);
    myPrior.unfold(uhist[i]->hMeas,nBayesianIter);


    // TH1F *hReweighted = (TH1F*)(TH1F*)uhist[i]->hResponse->ProjectionY(Form("hReweighted_cent%d",i));
		
    bayesianUnfold myUnfoldingJECSys(uhist[i]->hMatrixFit,hPrior,0);
    myUnfoldingJECSys.unfold(uhist[i]->hMeasJECSys,nBayesianIter);
    bayesianUnfold myUnfoldingSmearSys(uhist[i]->hMatrixFit,hPrior,0);
    myUnfoldingSmearSys.unfold(uhist[i]->hMeasSmearSys,nBayesianIter);
    bayesianUnfold myUnfolding(uhist[i]->hMatrixFit,hPrior,0); //changed to myPrior.hPrior from hPrior
    myUnfolding.unfold(uhist[i]->hMeas,nBayesianIter);
    cout <<"Unfolding bin "<<i<<endl;

    delete hBinByBinCorRaw;
    delete hMCGen;

    // Iteration Systematics
    for (int j=2;j<=nBayesianIter;j++)
      {

	bayesianUnfold myUnfoldingSys(uhist[i]->hMatrixFit,hPrior,0);
	myUnfoldingSys.unfold(uhist[i]->hMeas,j);
	uhist[i]->hRecoIterSys[j]  = (TH1F*) myUnfoldingSys.hPrior->Clone(Form("hRecoRAA_IterSys%d_cent%d",j,i));
      }
    
    
    uhist[i]->hReco         = (TH1F*) uhist[i]->hRecoIterSys[nBayesianIter]->Clone(Form("Unfolded_cent%i",i));
    uhist[i]->hRecoJECSys   = (TH1F*) myUnfoldingJECSys.hPrior->Clone(Form("UnfoldedJeCSys_cent%i",i));
    uhist[i]->hRecoSmearSys   = (TH1F*) myUnfoldingSmearSys.hPrior->Clone(Form("UnfoldedSmearSys_cent%i",i));
    uhist[i]->hRecoBinByBin->SetName(Form("UnfoldedBinByBin_cent%i",i));
    
    if (doToy) {
      TCanvas *cToy = new TCanvas("cToy","toy",600,600);
      cToy->cd();
      int nExp=1000;
      TH1F *hTmp[nbins_truth+1];
      TH1F *hTmp2[nbins_truth+1];
      for (int j=1;j<=nbins_truth;j++) {
	hTmp[j] = new TH1F(Form("hTmp%d",j),"",200,0,10.+uhist[i]->hReco->GetBinContent(j)*2);
	hTmp2[j] = new TH1F(Form("hTmp2%d",j),"",200,0,10.+uhist[i]->hRecoBinByBin->GetBinContent(j)*2);
      }
      for (int exp =0; exp<nExp; exp++) {
	TH1F *hToy = (TH1F*)uhist[i]->hMeas->Clone();   
	TH2F *hMatrixToy = (TH2F*)uhist[i]->hMatrixFit->Clone();
	hToy->SetName("hToy");
	if (exp%100==0) cout <<"Pseudo-experiment "<<exp<<endl;
	for (int j=1;j<=hToy->GetNbinsX();j++) {
	  double value = gRandom->Poisson(uhist[i]->hMeas->GetBinContent(j));
	  hToy->SetBinContent(j,value);
	}
				
	for (int j=1;j<=hMatrixToy->GetNbinsX();j++) {
	  for (int k=1;k<=hMatrixToy->GetNbinsY();k++) {
	    double value = gRandom->Gaus(uhist[i]->hMatrixFit->GetBinContent(j,k),uhist[i]->hMatrixFit->GetBinError(j,k));
	    hMatrixToy->SetBinContent(j,k,value);
	  }
	}

	prior myPriorToy(hMatrixToy,hToy,0.0);
	myPriorToy.unfold(hToy,1);
	bayesianUnfold myUnfoldingToy(hMatrixToy,myPriorToy.hPrior,0.0);
	myUnfoldingToy.unfold(hToy,nBayesianIter);
	TH1F *hRecoTmp = (TH1F*) myUnfoldingToy.hPrior->Clone();
				
	for (int j=1;j<=hRecoTmp->GetNbinsX();j++) {
	  hTmp[j]->Fill(hRecoTmp->GetBinContent(j));
	}
	delete hToy;
	delete hRecoTmp;
	delete hMatrixToy;
      }
      TF1 *fGaus = new TF1("fGaus","[0]*TMath::Gaus(x,[1],[2])");
      for (int j=1;j<=nbins_truth;j++)
	{

	  f->SetParameters(hTmp[j]->GetMaximum(),hTmp[j]->GetMean(),hTmp[j]->GetRMS());
				
	  if (hTmp[j]->GetMean()>0) {
	    hTmp[j]->Fit(fGaus,"LL Q ");
	    hTmp[j]->Fit(fGaus,"LL Q ");
	    uhist[i]->hReco->SetBinError(j,f->GetParameter(2));
	  }	       
	  f->SetParameters(hTmp2[j]->GetMaximum(),hTmp2[j]->GetMean(),hTmp2[j]->GetRMS());
	  if (hTmp2[j]->GetMean()>0) {
	    hTmp2[j]->Fit(fGaus,"LL Q ");
	    hTmp2[j]->Fit(fGaus,"LL Q ");
	    uhist[i]->hRecoBinByBin->SetBinError(j,f->GetParameter(2));
	  }	       
	  delete hTmp[j];
	  delete hTmp2[j];
	}
      cPbPb->cd(i+1);
    }

    uhist[i]->hMeas->SetMarkerStyle(20);
    uhist[i]->hMeas->SetMarkerColor(2);
    uhist[i]->hReco->SetMarkerStyle(25);
    uhist[i]->hReco->SetName(Form("hReco_cent%d",i));
    
    uhist[i]->hReco->SetXTitle("p_{T} (GeV/c)");    
    uhist[i]->hReco->SetYTitle("dN/dp_{T} (GeV/c)^{-1}");    
    uhist[i]->hReco->GetXaxis()->SetNdivisions(505);
    //uhist[i]->hReco->Draw("");    
    uhist[i]->hReco->SetAxisRange(40.,400.,"X");
    uhist[i]->hReco->Draw("");   
     
    uhist[i]->hGen->SetLineWidth(2);
    uhist[i]->hGen->SetLineColor(2);
    //if(isMC)uhist[i]->hGen->Draw("hist same");
    //uhist[i]->hReco->Draw("same");    
    uhist[i]->hRecoBinByBin->SetMarkerStyle(28);
    uhist[i]->hRecoBinByBin->Draw("same");    

    uhist[i]->hReco->SetAxisRange(40.,400.);
    TH1F *hReproduced = (TH1F*)myUnfolding.hReproduced->Clone(Form("hReproduced_cent%d",i));
    hReproduced->SetMarkerColor(4);
    hReproduced->SetMarkerStyle(24);
    //uhist[i]->hMeas->Draw("same");    

    hRecoBW[i] = (TH1F*)uhist[i]->hReco->Clone(Form("hReco%d",i));
    hRecoBinByBinBW[i] = (TH1F*)uhist[i]->hRecoBinByBin->Clone(Form("hRecoBinByBin%d",i));
    hMeasBW[i] = (TH1F*)uhist[i]->hMeas->Clone(Form("hMeas%d",i));
    if(isMC)hGenBW[i] = (TH1F*)uhist[i]->hGen->Clone(Form("hGen%d",i));
    hReproducedBW[i] = (TH1F*) hReproduced->Clone(Form("hReproduced%d",i));
    
    divideBinWidth(hRecoBW[i]);    
    if(isMC)divideBinWidth(hGenBW[i]);    
    divideBinWidth(hRecoBinByBinBW[i]);    
    divideBinWidth(hMeasBW[i]);    
    divideBinWidth(hReproducedBW[i]);

    hRecoBW[i]->Draw();
    if(isMC)hGenBW[i]->Draw("hist,same");
    hRecoBinByBinBW[i]->Draw("same");
    hMeasBW[i]->Draw("same");

    if(i==0) tpbpb->Draw();
    else tpp->Draw();


    uhist[i]->hReco->SetTitle("Baysian Unfolded");
    uhist[i]->hRecoBinByBin->SetTitle("Bin-by-bin Unfolded");

    TLegend *leg = new TLegend(0.45,0.65,0.85,0.90);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->AddEntry(uhist[i]->hMeas,"Measured","pl");
    leg->AddEntry(uhist[i]->hReco,"Bayesian unfolded","pl");
    leg->AddEntry(uhist[i]->hRecoBinByBin,"Bin-by-bin unfolded","pl");
    if(isMC)leg->AddEntry(uhist[i]->hGen,"Generator level truth","l");
    leg->Draw();

    cPbPbMeas->cd(i+1)->SetLogy();   
    //uhist[i]->hMeas->SetAxisRange(0,600,"X");
    //uhist[i]->hMeas->Draw();
    //hReproduced->Draw("same");
    hMeasBW[i]->SetAxisRange(40,600,"x");
    hMeasBW[i]->SetXTitle("p_{T} (GeV/c)");
    hMeasBW[i]->SetYTitle("dN/p_{T} (GeV/c)^{-1}");
    hMeasBW[i]->Draw();
    hReproducedBW[i]->Draw("same");
    if(i==0) tpbpb->Draw();
    else tpp->Draw();
    
    TLegend *leg2 = new TLegend(0.5,0.6,0.85,0.9);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    leg2->SetTextFont(42);
    leg2->AddEntry(uhist[i]->hMeas,"Measured","pl");
    leg2->AddEntry(hReproduced,"Reproduced","pl");

    leg2->Draw();


  }	     
  
 
  pbpb_Unfo->Write();

  SysData systematics;
  TLine *line = new TLine(40.,1,400.,1);

  // Iteration systematics
  TCanvas *cIterSys = new TCanvas("cIterSys","cIterSys",1200,600);
  cIterSys->Divide(2,1);
  cIterSys->cd(2);
  cIterSys->GetPad(2)->SetGridy();
  TH1F *hRecoIterSysPP[100];
  TH1F *hRebinPP_tmp         = rebin(uhist[nbins_cent]->hReco, (char*)"hRebinPP_tmp");
  TLegend *legBayesianIterPP = myLegend(0.4,0.7,0.9,0.9);
  legBayesianIterPP->SetTextFont(42);
  //legBayesianIterPP->AddEntry("","pp","");
         
  for (int j=2;j<nBayesianIter;j++) {
    hRecoIterSysPP[j] = rebin(uhist[nbins_cent]->hRecoIterSys[j],Form("hRecoIterSysPP_IterSys%d",j));
    hRecoIterSysPP[j]->SetLineColor(colorCode[j-2]);
    hRecoIterSysPP[j]->SetMarkerColor(colorCode[j-2]);
    hRecoIterSysPP[j]->Divide(hRebinPP_tmp);
    if (j==2){
      makeHistTitle(hRecoIterSysPP[j],(char*)"",(char*)"Jet p_{T} (GeV/c)",(char*)"Ratio (Nth Iteration / Nominal)");
      hRecoIterSysPP[j]->SetTitleOffset(1.4,"Y");
      hRecoIterSysPP[j]->SetTitleOffset(1.2,"X");
      hRecoIterSysPP[j]->SetAxisRange(40.,400.,"X");
      hRecoIterSysPP[j]->SetAxisRange(0.8,1.2,"Y");
      hRecoIterSysPP[j]->Draw("hist"); 
    } else {
      hRecoIterSysPP[j]->Draw("hist,same");
    }
         
    checkMaximumSys(systematics.hSysIter[nbins_cent],hRecoIterSysPP[j],0,1.1);
    if(j!=4)legBayesianIterPP->AddEntry(hRecoIterSysPP[j],Form("Iteration %d",j),"l");     
  }

  legBayesianIterPP->Draw();
  line->Draw();
  //drawEnvelope(systematics.hSysIter[nbins_cent],(char*)"hist same");

  tpp->Draw();


  cIterSys->cd(1);
  cIterSys->GetPad(1)->SetGridy();
  TH1F *hRecoIterSysPbPb[100];
  TH1F *hRebinPbPb_tmp         = rebin(uhist[0]->hReco, (char*)"hRebinPbPb_tmp");
  TLegend *legBayesianIterPbPb = myLegend(0.4,0.7,0.9,0.9);
  legBayesianIterPbPb->SetTextFont(42);
  //legBayesianIterPbPb->AddEntry("","PbPb","");
  for (int j=2;j<nBayesianIter;j++) {
    hRecoIterSysPbPb[j] = rebin(uhist[0]->hRecoIterSys[j],Form("hRecoIterSysPbPb_IterSys%d",j));
    hRecoIterSysPbPb[j]->SetLineColor(colorCode[j-2]);
    hRecoIterSysPbPb[j]->SetMarkerColor(colorCode[j-2]);
    hRecoIterSysPbPb[j]->Divide(hRebinPbPb_tmp);
    if (j==2){
      makeHistTitle(hRecoIterSysPbPb[j],(char*)"",(char*)"Jet p_{T} (GeV/c)",(char*)"Ratio (Nth Iteration / Nominal)");
      hRecoIterSysPbPb[j]->SetTitleOffset(1.4,"Y");
      hRecoIterSysPbPb[j]->SetTitleOffset(1.2,"X");
      hRecoIterSysPbPb[j]->SetAxisRange(40.,400.,"X");
      hRecoIterSysPbPb[j]->SetAxisRange(0.8,1.2,"Y");
      hRecoIterSysPbPb[j]->Draw("hist"); 
    } else {
      hRecoIterSysPbPb[j]->Draw("hist,same");
    }
         
    checkMaximumSys(systematics.hSysIter[0],hRecoIterSysPbPb[j],0,1.1);
    if(j!=4)legBayesianIterPbPb->AddEntry(hRecoIterSysPbPb[j],Form("Iteration %d",j),"l");     
  }
  legBayesianIterPbPb->Draw();
  line->Draw();
  tpbpb->Draw();
  //drawEnvelope(systematics.hSysIter[0],(char*)"hist same");

  if(isMC){
    TCanvas *cRatios = new TCanvas("cRatios","",1200,600);
    cRatios->Divide(2,1);
    
    TH1D *hRatioGen[2], *hRatioBbyB[2], *hRatioReco[2], *hRatioMeas[2];
    
    for(int i=0; i<2; i++){
      cRatios->cd(i+1);
      
      hRatioGen[i] = (TH1D*)uhist[i]->hGen->Clone(Form("hRatioGen_%d",i));
      hRatioBbyB[i] = (TH1D*)uhist[i]->hRecoBinByBin->Clone(Form("hRatioBbyB_%d",i));
      hRatioReco[i] = (TH1D*)uhist[i]->hReco->Clone(Form("hRatioReco_%d",i));
      hRatioMeas[i] = (TH1D*)uhist[i]->hMeas->Clone(Form("hRatioMeas_%d",i));
      hRatioMeas[i]->Divide(uhist[i]->hMeas,uhist[i]->hGen,1,1);
      hRatioBbyB[i]->Divide(uhist[i]->hRecoBinByBin,uhist[i]->hGen,1,1);
      hRatioReco[i]->Divide(uhist[i]->hReco,uhist[i]->hGen,1,1);
      
      hRatioBbyB[i]->SetMaximum(1.3);
      hRatioBbyB[i]->SetMinimum(0.6);
      hRatioBbyB[i]->Draw(); 
      hRatioMeas[i]->Draw("same");
      hRatioReco[i]->Draw("same");
      
      TLegend *legr = new TLegend(0.45,0.65,0.85,0.90);
      legr->SetBorderSize(0);
      legr->SetFillStyle(0);
      legr->SetTextFont(42);
      legr->AddEntry(hRatioReco[i],"Bayesian unfolded/Gen Truth","pl");
      legr->AddEntry(hRatioBbyB[i],"Bin-by-bin unfolded/Gen Truth","pl");
      legr->AddEntry(hRatioMeas[i],"Measured/Gen Truth","pl");
      legr->Draw();
    }
  }
}





