#ifndef __generate_smear_matrix__H
#define __generate_smear_matrix__H

#include "TH1F.h"
#include "TH2F.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TCanvas.h"
#include "TLine.h"

double _epsilon = 1e-6; // speed up calculations with acceptable loss of precision
const int nk = 4; // number of kernel parameters (excluding pt, eta)
//TF1 *parf = new TF1("parf","[0]*TMath::Sqrt(pow([1],2)/pow(x,2) + pow([2],2)/x + pow([3],2))",35,500); //need to make these global to save time
TF1 *parf = new TF1("parf","[0]/pow(x,[1])",50,300);
TF1 *_kernel = 0;

double ptresolution(double pt, double eta){

  //method taken from SMP-13-002
  // if(abs(eta)>1) parf->SetParameters(1.6205,74.597,2.170,0.3979); //calculated from pPb MC for 1<abs(eta CM)<2
  //else if(abs(eta)>0.5) parf->SetParameters(1.6163,22.67,3.33,0.395); //calculated from pPb MC for 0.5<abs(eta CM)<1
  //else parf->SetParameters(1.6163,13.6607,2.844,0.3953); //calculated from pPb MC for abs(eta CM) < 0.5

  parf->SetParameters(1.052,0.5261);

  if(pt<30) return 1.;
  else return parf->Eval(pt);

}

// Ansatz Kernel

Double_t smearedAnsatzKernel(Double_t *x, Double_t *p) {
  int cnt_a = 0;
  if (++cnt_a % 1000000==0) {
    cout << "." << flush;
  }

  const double pt = x[0]; // true pT
  const double ptmeas = p[0]; // measured pT
  const double eta = p[1]; // rapidity

  double res = 0.;
  res = ptresolution(pt, eta+1e-3) * pt;
  const double s = TMath::Gaus(ptmeas, pt, res, kTRUE);
  const double f = p[2] * exp(p[3]/pt) * pow(pt, p[4])
    * pow(1 - pt*cosh(eta) / 3500., p[5]);

  return (f * s);
}

Double_t smearedAnsatz(Double_t *x, Double_t *p) {

  if(!_kernel) _kernel = new TF1("_kernel",smearedAnsatzKernel,1.,1000.,nk+2);
  
  const double pt = x[0];
  const double eta = p[0];
  double res = 0.;
  res = ptresolution(pt, eta+1e-3) * pt;

  const double sigma = TMath::Min(res, 0.30);
  double ptmin = pt / (1. + 4.*sigma); // xmin*(1+4*sigma)=x
  ptmin = TMath::Max(1.,ptmin); // safety check
  double ptmax = pt / (1. - 2.*sigma); // xmax*(1-2*sigma)=x
  ptmax = TMath::Min(3500./cosh(eta),ptmax); // safety check
  const double par[nk+2] = {pt, eta, p[1], p[2], p[3], p[4]};
  _kernel->SetParameters(&par[0]);

  // Set pT bin limits needed in smearing matrix generation
  if (p[5]>0 && p[5]<3500./cosh(eta)) ptmin = p[5];
  if (p[6]>0 && p[6]<3500./cosh(eta)) ptmax = p[6];

  return ( _kernel->Integral(ptmin, ptmax, &par[0], _epsilon) );
}


TH2F *generateSmearingMatrix(int iter, TH1F *hgen, TH1F *hpt, double genPtmin=20., double recoPtmin=30., double etalo=-2., bool binBbin_or_bayes=1, int _debug=1){

  // initial fit of the NLO curve to a histogram
  TF1 *fnlo = new TF1(Form("fus%d",iter),
                      "[0]*exp([1]/x)*pow(x,[2])"
                      "*pow(1-x*cosh([4])/3500.,[3])", 10., 1000.);
  fnlo->SetParameters(2e6,-35.,-5.2,2.,etalo);
  fnlo->FixParameter(4,etalo);

  hgen->Fit(fnlo,"QRN");
  //  cout << "fnlo integral: "<< fnlo->Integral(20,500) << endl;
  //cout << "fnlo eval 50 & 100: "<< fnlo->Eval(50) << " " << fnlo->Eval(100) << endl;

  // Create smeared theory curve
  double maxpt = 3450./cosh(etalo);
  TF1 *fnlos = new TF1(Form("fs%d",iter),smearedAnsatz,genPtmin,maxpt,nk+3);
  
  fnlos->SetParameters(etalo, fnlo->GetParameter(0), fnlo->GetParameter(1),
                       fnlo->GetParameter(2), fnlo->GetParameter(3), 0, 0);

  //cout << "fnlos integral: "<< fnlos->Integral(20,500) << endl;
  //cout << "fnlos eval 20 50 and 100: "<< fnlos->Eval(20) << " " << fnlos->Eval(50) << " " << fnlos->Eval(100) << endl;

 if (_debug)
    cout << "Calculate forward smearing and unfold hpt" << endl << flush;

    // Calculate smearing matrix
  if (_debug) 
    cout << "Generating smearing matrix T..." << flush;
  double tmp_eps = _epsilon;

  // NB: GetArray only works if custom x binning
  //outdir->cd();

  // Deduce range and binning for true and measured spectra
  vector<double> vx; // true
  vector<double> vy; // measured
  for (int i = 1; i != hgen->GetNbinsX()+1; ++i) {

    double x = hgen->GetBinCenter(i);
    double x1 = hgen->GetBinLowEdge(i);
    double x2 = hgen->GetBinLowEdge(i+1);
    double y = hgen->GetBinContent(i);

    // if (x1>=recoPtmin && y>0) {
      if (vx.size()==0){ vx.push_back(x1); }
       vx.push_back(x2);
      //}
  }
  for (int i = 1; i != hpt->GetNbinsX()+1; ++i) {
    double x = hpt->GetBinCenter(i);
    double x1 = hpt->GetBinLowEdge(i);
    double x2 = hpt->GetBinLowEdge(i+1);
    double y = hpt->GetBinContent(i);

    //if (x1>=recoPtmin && y>0) {
      if (vy.size()==0){ vy.push_back(x1); }
       vy.push_back(x2);
      // }
  } // for i
  
  // copy over relevant part of hpt
  TH1F *hreco = new TH1F(Form("hreco_jet_%d",iter),";p_{T,reco} (GeV)",
			 vy.size()-1,&vy[0]);
  for (int i = 1; i != hreco->GetNbinsX()+1; ++i) {
    int j = hpt->FindBin(hreco->GetBinCenter(i));
    double dpt = hpt->GetBinWidth(j);
    hreco->SetBinContent(i, hpt->GetBinContent(j)/dpt);
    hreco->SetBinError(i, hpt->GetBinError(j)/dpt);
  }
  for (int i = 1; i != hgen->GetNbinsX()+1; ++i) {
    int j = hpt->FindBin(hgen->GetBinCenter(i));
    double dpt = hpt->GetBinWidth(j);
    hgen->SetBinContent(i, hgen->GetBinContent(j)/dpt);
    hgen->SetBinError(i, hgen->GetBinError(j)/dpt);
  }

  TH2F *mt = new TH2F(Form("mt_jet_%d",iter),"mt;p_{T,reco};p_{T,gen}",
		      vy.size()-1, &vy[0], vx.size()-1, &vx[0]);
  TH1F *mx = new TH1F(Form("mx_jet_%d",iter),"mx;p_{T,gen};#sigma/dp_{T}",
		      vx.size()-1, &vx[0]);
  TH1F *my = new TH1F(Form("my_jet_%d",iter),"my;p_{T,reco};#sigma/dp_{T}",
		      vy.size()-1, &vy[0]);

  // From http://hepunx.rl.ac.uk/~adye/software/unfold/RooUnfold.html
  // For 1-dimensional true and measured distribution bins Tj and Mi,
  // the response matrix element Rij gives the fraction of events
  // from bin Tj that end up measured in bin Mi.

  for (int i = 1; i != mt->GetNbinsX()+1; ++i) {

    double ptreco1 = mt->GetXaxis()->GetBinLowEdge(i);
    double ptreco2 = mt->GetXaxis()->GetBinLowEdge(i+1);
    double yreco = fnlo->Integral(ptreco1, ptreco2) / (ptreco2 - ptreco1);
    double ptreco = fnlo->GetX(yreco, ptreco1, ptreco2);
    cout << "ptreco1: "<< ptreco1 << " ptreco2: "<< ptreco2 << " yreco: "<< yreco << " ptreco: "<< ptreco << endl;

    for (int j = 1; j != mt->GetNbinsY()+1; ++j) {

      double ptgen1 = min(3500./cosh(etalo), mt->GetYaxis()->GetBinLowEdge(j));
      double ptgen2 = min(3500./cosh(etalo), mt->GetYaxis()->GetBinLowEdge(j+1));

      if (ptgen1>genPtmin && ptreco>recoPtmin && ptgen1*cosh(etalo) < 3500.) {
        fnlos->SetParameter(5, ptgen1);
        fnlos->SetParameter(6, ptgen2);
        // 2D integration over pTreco, pTgen simplified to 1D over pTgen
        mt->SetBinContent(i, j, fnlos->Eval(ptreco) / (ptreco2 - ptreco1));
        fnlos->SetParameter(5, 0);
        fnlos->SetParameter(6, 0);
      }
    } // for j
  } // for i

  for (int j = 1; j != mt->GetNbinsY()+1; ++j) {

    double ptgen1 = min(3500./cosh(etalo), mt->GetYaxis()->GetBinLowEdge(j));
    double ptgen2 = min(3500./cosh(etalo), mt->GetYaxis()->GetBinLowEdge(j+1));
    double ygen = fnlo->Integral(ptgen1, ptgen2);
    mx->SetBinContent(j, ygen);
  }
  
  for (int i = 1; i != mt->GetNbinsX()+1; ++i) {

    double yreco(0);
    for (int j = 1; j != mt->GetNbinsY()+1; ++j) {
      yreco += mt->GetBinContent(i, j);
    }
    my->SetBinContent(i, yreco);
  } // for i
  
  TH2F *mtu = (TH2F*)mt->Clone(Form("mtu_jet_%d",iter));
  for (int i = 1; i != mt->GetNbinsX()+1; ++i) {
    for (int j = 1; j != mt->GetNbinsY()+1; ++j) {
      if (mx->GetBinContent(i)!=0) {
	mtu->SetBinContent(i, j, mt->GetBinContent(i,j) / mx->GetBinContent(j));
      }
    } // for j
  } // for i


  // For BinByBin and SVD, need square matrix
  TH2F *mts(0);
  TH1F *mxs(0);


  mts = new TH2F(Form("mts_jet_%d",iter),"mts;p_{T,reco};p_{T,gen}",
		 vy.size()-1, &vy[0], vy.size()-1, &vy[0]);
  mxs = new TH1F(Form("mxs_jet_%d",iter),"mxs;p_{T,gen};#sigma/dp_{T}",
		 vy.size()-1, &vy[0]);
  for (int i = 1; i != mts->GetNbinsX()+1; ++i) {
    for (int j = 1; j != mts->GetNbinsY()+1; ++j) {
      double x = mts->GetBinCenter(i);
	double y = mts->GetBinCenter(j);
	int i2 = mt->GetXaxis()->FindBin(x);
	int j2 = mt->GetYaxis()->FindBin(y);
	mts->SetBinContent(i, j, mt->GetBinContent(i2, j2));
	mts->SetBinError(i, j, mt->GetBinError(i2, j2));
      }
    }
    for (int i = 1; i != mxs->GetNbinsX()+1; ++i) {
      double x = mxs->GetBinCenter(i);
      int i2 = mx->FindBin(x);
      mxs->SetBinContent(i, mx->GetBinContent(i2));
      mxs->SetBinError(i, mx->GetBinError(i2));
    }


//RooUnfoldResponse *uResp;
//if(binBbin_or_bayes) uResp = new RooUnfoldResponse(my, mx, mt);
// else uResp = uResp = new RooUnfoldResponse(my, mxs, mts);

_epsilon = tmp_eps;
if (_debug) 
  cout << "done." << endl << flush;

if(binBbin_or_bayes) return mtu;
 else return mts;
   
}

#endif
