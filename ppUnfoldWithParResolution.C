// d'Agostini ("Bayesian" or Richardson-Lucy) unfolding and
// response matrix generation based on input NLO spectrum + parameterized JER
// call from mk_dagostini.C
#include "TFile.h"
#include "TDirectory.h"
#include "TList.h"
#include "TObject.h"
#include "TKey.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TCAnvas.h"
#include "TLine.h"
#include "TStyle.h"
#include "TLegend.h"

#include "RooUnfold-1.1.1/src/RooUnfold.h"
#include "RooUnfold-1.1.1/src/RooUnfoldBayes.h"
#include "RooUnfold-1.1.1/src/RooUnfoldBinByBin.h"
#include "RooUnfold-1.1.1/src/RooUnfoldSvd.h"
#include "RooUnfold-1.1.1/src/RooUnfoldResponse.h"
//#include "RooUnfold.h"
//#include "RooUnfoldBayes.h"
//#include "RooUnfoldResponse.h"

#include "tdrstyle_mod.C"
#include "ptresolution.h"
#include "settings.h"
#include "tools.h"

#include <iostream>

using namespace std;

// Resolution function
bool _ispf = false; // global variable
bool _iscalo = false; // global variable
int _jk = 0; // global variable
bool _jet = false; // global variable
                                  
Double_t fPtRes(Double_t *x, Double_t *p) {

  if (_ispf) return ptresolution(x[0], p[0]);
  if (_iscalo) return ptresolution_jpt(x[0], p[0]);
  return 0.1;
}

// Ansatz Kernel
int cnt_a = 0;
const int nk = 4; // number of kernel parameters (excluding pt, eta)
Double_t smearedAnsatzKernel(Double_t *x, Double_t *p) {

  if (++cnt_a%1000000==0) {
    cout << "." << flush;
  }

  const double pt = x[0]; // true pT
  const double ptmeas = p[0]; // measured pT
  const double eta = p[1]; // rapidity

  double res = 0.;
  if (_iscalo) res = ptresolution_calo(pt, eta+1e-3) * pt;
  if (_ispf) res = ptresolution(pt, eta+1e-3) * pt;
  const double s = TMath::Gaus(ptmeas, pt, res, kTRUE);
  const double f = p[2] * exp(p[3]/pt) * pow(pt, p[4])
    * pow(1 - pt*cosh(eta) / 3500., p[5]);

  return (f * s);
}

// Smeared Ansatz
double _epsilon = 1e-12;
TF1 *_kernel = 0; // global variable, not pretty but works
Double_t smearedAnsatz(Double_t *x, Double_t *p) {

  if (!_kernel) _kernel = new TF1("_kernel",smearedAnsatzKernel,1.,1000.,nk+2);
  
  const double pt = x[0];
  const double eta = p[0];
  double res = 0.;
  if (_iscalo) res = ptresolution_calo(pt, eta+1e-3) * pt;
  if (_ispf) res = ptresolution(pt, eta+1e-3) * pt;
  const double sigma = min(res, 0.30);
  double ptmin = pt / (1. + 4.*sigma); // xmin*(1+4*sigma)=x
  ptmin = max(1.,ptmin); // safety check
  double ptmax = pt / (1. - 2.*sigma); // xmax*(1-2*sigma)=x
  ptmax = min(3500./cosh(eta),ptmax); // safety check
  const double par[nk+2] = {pt, eta, p[1], p[2], p[3], p[4]};
  _kernel->SetParameters(&par[0]);

  // Set pT bin limits needed in smearing matrix generation
  if (p[5]>0 && p[5]<3500./cosh(eta)) ptmin = p[5];
  if (p[6]>0 && p[6]<3500./cosh(eta)) ptmax = p[6];

  return ( _kernel->Integral(ptmin, ptmax, &par[0], _epsilon) );
}

void recurseFile(TDirectory *indir, TDirectory *indir2, TDirectory *outdir,
                 bool ismc);
void dagostiniUnfold_histo(TH1D *hpt, TH1D *hpt2, TDirectory *outdir,
			   bool ismc, bool kscale = false, string id = "");


void dagostiniUnfold(string type) {

  TFile *fin = new TFile(Form("output-%s-2b.root",type.c_str()),"READ");
  assert(fin && !fin->IsZombie());

  TFile *fin2 = new TFile(Form("output-%s-5.root",type.c_str()),"READ");
  assert(fin2 && !fin2->IsZombie());

  TFile *fout = new TFile(Form("output-%s-3c.root",type.c_str()),"RECREATE");
  assert(fout && !fout->IsZombie());

  _ak7 = (_algo=="AK7");
  if (_ak7) cout << "Using AK7 JER" << endl << flush;
  bool ismc = (type=="MC"||type=="HW");

  recurseFile(fin, fin2, fout, ismc);

  cout << "Output stored in " << fout->GetName() << endl;
  fout->Close();
  fout->Delete();

  fin->Close();
  fin->Delete();
}


void recurseFile(TDirectory *indir, TDirectory *indir2, TDirectory *outdir,
                 bool ismc) {

  TDirectory *curdir = gDirectory;

  // Automatically go through the list of keys (directories)
  TList *keys = indir->GetListOfKeys();
  TListIter itkey(keys);
  TObject *key, *obj;

  while ( (key = itkey.Next()) ) {

    obj = ((TKey*)key)->ReadObj(); assert(obj);

    // Found a subdirectory: copy it to output and go deeper
    if (obj->InheritsFrom("TDirectory")) {

      if (_debug) cout << key->GetName() << endl;

      assert(outdir->mkdir(obj->GetName()));
      assert(outdir->cd(obj->GetName()));
      TDirectory *outdir2 = gDirectory;

      assert(indir->cd(obj->GetName()));
      TDirectory *indir2a = gDirectory;

      if (indir2->cd(obj->GetName())) {
        TDirectory *indir2b = gDirectory;

        recurseFile(indir2a, indir2b, outdir2, ismc);
      }
    } // inherits from TDirectory

    // Found hpt plot: call unfolding routine
    if (obj->InheritsFrom("TH1") &&
        (string(obj->GetName())=="hpt" ||
	 string(obj->GetName())=="hpt_jet" ||
	 string(obj->GetName())=="hpt_jk1" ||
	 string(obj->GetName())=="hpt_jk2" ||
	 string(obj->GetName())=="hpt_jk3" ||
	 string(obj->GetName())=="hpt_jk4" ||
	 string(obj->GetName())=="hpt_jk5" ||
	 string(obj->GetName())=="hpt_jk6" ||
	 string(obj->GetName())=="hpt_jk7" ||
	 string(obj->GetName())=="hpt_jk8" ||
	 string(obj->GetName())=="hpt_jk9" ||
	 string(obj->GetName())=="hpt_jk10"
	 )) {
      
      cout << "+" << flush;

      _jk = 0;
      if (TString(obj->GetName()).Contains("hpt_jk")) {
	assert( sscanf(obj->GetName(), "hpt_jk%d", &_jk) == 1);
      }
      _jet = TString(obj->GetName()).Contains("hpt_jet");
      
      TH1D *hpt = (TH1D*)obj;
      TH1D *hpt2 = (TH1D*)indir2->Get("hnlo"); assert(hpt2);
      if (hpt2)
        dagostiniUnfold_histo(hpt, hpt2, outdir, ismc);
    } // hpt                  

    // Try to process friends similarly
    if (obj->InheritsFrom("TH1") &&
        (string(obj->GetName())=="hpt_ak5calo")) {

      cout << "-" << flush;

      _jk = 0; _jet = false;
      TH1D *hpt = (TH1D*)obj;
      TH1D *hpt2 = (TH1D*)indir2->Get("hnlo"); assert(hpt2);
      if (hpt2)
        dagostiniUnfold_histo(hpt, hpt2, outdir, ismc, false, "_ak5calo");
    } // hpt

  } // while key
  
  curdir->cd();
} // recurseFile    


void dagostiniUnfold_histo(TH1D *hpt, TH1D *hnlo, TDirectory *outdir,
			   bool ismc, bool kscale, string id) {

  float y1, y2;
  assert(sscanf(outdir->GetName(),"Eta_%f-%f",&y1,&y2)==2);
  const char *c = id.c_str();
  if (_jk) c = Form("_jk%d",_jk);
  if (_jet) c = "_jet";
  
  _ismcjer = ismc;
  if (id=="") { _ispf=true; _iscalo=false; }
  else        { _ispf=false; _iscalo=true; }

  // initial fit of the NLO curve to a histogram
  TF1 *fnlo = new TF1(Form("fus%s",c),
                      "[0]*exp([1]/x)*pow(x,[2])"
                      "*pow(1-x*cosh([4])/3500.,[3])", 10., 1000.);
  fnlo->SetParameters(2e14,-18,-5.2,8.9,y1);
  fnlo->FixParameter(4,y1);

  hnlo->Fit(fnlo,"QRN");

  // Graph of theory points with centered bins
  const double minerr = 0.02;
  TGraphErrors *gnlo = new TGraphErrors(0);
  TGraphErrors *gnlo2 = new TGraphErrors(0); // above + minerr
  gnlo->SetName("gnlo");
  gnlo2->SetName("gnlo2");
  for (int i = 1; i != hnlo->GetNbinsX()+1; ++i) {

    double y = hnlo->GetBinContent(i);
    double dy = hnlo->GetBinError(i);
    double ptmin = hnlo->GetBinLowEdge(i);
    double ptmax = hnlo->GetBinLowEdge(i+1);

    double y0 = fnlo->Integral(ptmin, ptmax) / (ptmax - ptmin);
    double x = fnlo->GetX(y0, ptmin, ptmax);

    int n = gnlo->GetN();
    tools::SetPoint(gnlo, n, x, y, 0, dy);
    tools::SetPoint(gnlo2, n, x, y, 0, tools::oplus(dy, minerr*y));
  }

  // Second fit to properly centered graph
  gnlo2->Fit(fnlo,"QRN");
  
  // Bin-centered data points
  TGraphErrors *gpt = new TGraphErrors(0);
  gpt->SetName(Form("gpt%s",c));
  for (int i = 1; i != hpt->GetNbinsX()+1; ++i) {

    double ptmin = hpt->GetBinLowEdge(i);
    double ptmax = hpt->GetBinLowEdge(i+1);
    double y = fnlo->Integral(ptmin, ptmax) / (ptmax - ptmin);
    double x = fnlo->GetX(y, ptmin, ptmax);
    double ym = hpt->GetBinContent(i);
    double ym_err = hpt->GetBinError(i);
    if (ym>0) {
      tools::SetPoint(gpt, gpt->GetN(), x, ym, 0., ym_err);
    }
  } // for i                                                 

  // Create smeared theory curve
  double maxpt = 3450./cosh(y1);
  TF1 *fnlos = new TF1(Form("fs%s",c),smearedAnsatz,xmin,maxpt,nk+3);
  fnlos->SetParameters(y1, fnlo->GetParameter(0), fnlo->GetParameter(1),
                       fnlo->GetParameter(2), fnlo->GetParameter(3), 0, 0);

 if (_debug)
    cout << "Calculate forward smearing and unfold hpt" << endl << flush;

  TGraphErrors *gfold_fwd = new TGraphErrors(0);
  gfold_fwd->SetName(Form("gfold_fwd%s",c));
  TGraphErrors *gcorrpt_fwd = new TGraphErrors(0);
  gcorrpt_fwd->SetName(Form("gcorrpt_fwd%s",c));
  TH1D *hcorrpt_fwd = (TH1D*)hpt->Clone(Form("hcorrpt_fwd%s",c));

  for (int i = 0; i != gpt->GetN(); ++i) {
    double x, y, ex, ey;
    tools::GetPoint(gpt, i, x, y, ex, ey);
    double k = fnlo->Eval(x) / fnlos->Eval(x);
    if (!TMath::IsNaN(k)) {

      tools::SetPoint(gfold_fwd, gfold_fwd->GetN(), x, k, ex, 0.);
      tools::SetPoint(gcorrpt_fwd, gcorrpt_fwd->GetN(), x, k*y, ex, k*ey);

      int j = hpt->FindBin(x);
      hcorrpt_fwd->SetBinContent(j, k*hpt->GetBinContent(j));
      hcorrpt_fwd->SetBinError(j, k*hpt->GetBinError(j));
    }
  }


  // Calculate smearing matrix
  if (_debug) 
    cout << "Generating smearing matrix T..." << flush;
  double tmp_eps = _epsilon;
  _epsilon = 1e-6; // speed up calculations with acceptable loss of precision

  // NB: GetArray only works if custom x binning
  outdir->cd();

  // Deduce range and binning for true and measured spectra
  vector<double> vx; // true
  vector<double> vy; // measured
  for (int i = 1; i != hpt->GetNbinsX()+1; ++i) {

    double x = hpt->GetBinCenter(i);
    double x1 = hpt->GetBinLowEdge(i);
    double x2 = hpt->GetBinLowEdge(i+1);
    double y = hpt->GetBinContent(i);

    if (x>=_recopt && y>0) {
      if (vx.size()==0) vx.push_back(x1);
      vx.push_back(x2);
    }

    if (x>=fitptmin && y>0) {
      if (vy.size()==0) vy.push_back(x1);
      vy.push_back(x2);
    }
  } // for i
  
  // copy over relevant part of hpt
  TH1D *hreco = new TH1D(Form("hreco%s",c),";p_{T,reco} (GeV)",
			 vy.size()-1,&vy[0]);
  for (int i = 1; i != hreco->GetNbinsX()+1; ++i) {
    int j = hpt->FindBin(hreco->GetBinCenter(i));
    double dpt = hpt->GetBinWidth(j);
    hreco->SetBinContent(i, hpt->GetBinContent(j)*dpt);
    hreco->SetBinError(i, hpt->GetBinError(j)*dpt);
  }

  TH2D *mt = new TH2D(Form("mt%s",c),"mt;p_{T,reco};p_{T,gen}",
		      vy.size()-1, &vy[0], vx.size()-1, &vx[0]);
  TH1D *mx = new TH1D(Form("mx%s",c),"mx;p_{T,gen};#sigma/dp_{T}",
		      vx.size()-1, &vx[0]);
  TH1D *my = new TH1D(Form("my%s",c),"my;p_{T,reco};#sigma/dp_{T}",
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
    
    for (int j = 1; j != mt->GetNbinsY()+1; ++j) {

      double ptgen1 = min(3500./cosh(y1), mt->GetYaxis()->GetBinLowEdge(j));
      double ptgen2 = min(3500./cosh(y1), mt->GetYaxis()->GetBinLowEdge(j+1));

      if (ptgen1>21. && ptreco>21. && ptgen1*cosh(y1) < 3500.) {
        fnlos->SetParameter(5, ptgen1);
        fnlos->SetParameter(6, ptgen2);
        // 2D integration over pTreco, pTgen simplified to 1D over pTgen
        mt->SetBinContent(i, j, fnlos->Eval(ptreco) * (ptreco2 - ptreco1));
        fnlos->SetParameter(5, 0);
        fnlos->SetParameter(6, 0);
      }
    } // for j
  } // for i

  for (int j = 1; j != mt->GetNbinsY()+1; ++j) {

    double ptgen1 = min(3500./cosh(y1), mt->GetYaxis()->GetBinLowEdge(j));
    double ptgen2 = min(3500./cosh(y1), mt->GetYaxis()->GetBinLowEdge(j+1));
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
  
  TH2D *mtu = (TH2D*)mt->Clone(Form("mtu%s",c));
  for (int i = 1; i != mt->GetNbinsX()+1; ++i) {
    for (int j = 1; j != mt->GetNbinsY()+1; ++j) {
      if (mx->GetBinContent(i)!=0) {
	mtu->SetBinContent(i, j, mt->GetBinContent(i,j) / mx->GetBinContent(j));
      }
    } // for j
  } // for i


  // For BinByBin and SVD, need square matrix
  TH2D *mts(0);
  TH1D *mxs(0);
  if (!_jk && !_jet) {

    mts = new TH2D(Form("mts%s",c),"mts;p_{T,reco};p_{T,gen}",
		   vy.size()-1, &vy[0], vy.size()-1, &vy[0]);
    mxs = new TH1D(Form("mxs%s",c),"mxs;p_{T,gen};#sigma/dp_{T}",
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
  } // !_jk


  _epsilon = tmp_eps;
  if (_debug) 
    cout << "done." << endl << flush;

  outdir->cd();
  if (!_jk) {
    hreco->Write();
    mx->Write();
    my->Write();
    mt->Write();
    mtu->Write();
  }

  // Now to actual unfolding business with the d'Agostini method
  if (_debug)
    cout << "Unfolding..." << flush;

  // RooUnfoldResponse constructor - create from already-filled histograms
  // "response" (mt) gives the response matrix, measured X truth.
  // "measured" (my) and "truth" (mx) give the projections of "response"
  // onto the X-axis and Y-axis respectively,
  // but with additional entries in "measured" (my) for measurements with
  // no corresponding truth (fakes/background) [not implemented] and
  // in "truth" (mx) for unmeasured events (inefficiency) [is implemented].
  // "measured" and/or "truth" can be specified as 0 (1D case only)
  // or an empty histograms (no entries) as a shortcut
  // to indicate, respectively, no fakes and/or no inefficiency.
  // RooUnfoldResponse(const TH1* measured,
  //                   const TH1* truth, const TH2* response,
  //                   const char* name, const char* title)
  RooUnfoldResponse *uResp = new RooUnfoldResponse(my, mx, mt);

  // RooUnfoldBayes (const RooUnfoldResponse* res, const TH1* meas,
  //                 Int_t niter= 4, Bool_t smoothit= false,
  //                 const char* name= 0, const char* title= 0);
  RooUnfoldBayes *uBayes = new RooUnfoldBayes(uResp, hreco, 4);

  if (_debug)
    uBayes->Print();

  //hreco = my;
  TH1D *hTrueBayes = (TH1D*)uBayes->Hreco(RooUnfold::kCovariance);
  assert(hTrueBayes);
  TMatrixD *mCov = new TMatrixD(uBayes->Ereco());
  assert(mCov);
  //TH1D *hcorrpt = (TH1D*)hTrueBayes->Clone(Form("hcorrpt%s",c));
  //delete uBayes; // ensure static members get destroyed before next instance

  TH2D *hCov = new TH2D(Form("hCov%s",c), Form("hCov%s;p_{T};p_{T}",c),
			vx.size()-1, &vx[0], vx.size()-1, &vx[0]);
  assert(hCov->GetNbinsX()==mCov->GetNrows());
  assert(hCov->GetNbinsY()==mCov->GetNcols());
  for (int i = 1; i != hCov->GetNbinsX()+1; ++i) {
    for (int j = 1; j != hCov->GetNbinsY()+1; ++j) {
      hCov->SetBinContent(i, j, (*mCov)[i-1][j-1]);
    } // for j
  } // for i

  // copy into original binning for plotting macros to work (e.g. drawClosure)
  TH1D *hcorrpt = (TH1D*)hpt->Clone(Form("hcorrpt%s",c));
  hcorrpt->Reset();
  for (int i = 1; i != hTrueBayes->GetNbinsX()+1; ++i) {

    int j = hpt->FindBin(hTrueBayes->GetBinCenter(i));
    hcorrpt->SetBinContent(j, hTrueBayes->GetBinContent(i));
    hcorrpt->SetBinError(j, hTrueBayes->GetBinError(i));
  }

  // RooUnfoldBinByBin (const RooUnfoldResponse* res, const TH1* meas,
  //                    const char* name=0, const char* title=0);
  //hreco = my;
  // BinByBin and SVD can only handle square matrix
  TH1D *hcorrpt_bin(0), *hcorrpt_svd(0);
  if (!_jk && !_jet) {

    RooUnfoldResponse *uResps = new RooUnfoldResponse(my, mxs, mts);
    RooUnfoldBinByBin *uBin = new RooUnfoldBinByBin(uResps, hreco);
    TH1D *hTrueBin = (TH1D*)uBin->Hreco(RooUnfold::kCovariance);
    assert(hTrueBin);

    hcorrpt_bin = (TH1D*)hpt->Clone(Form("hcorrpt_bin%s",c));
    hcorrpt_bin->Reset();
    for (int i = 1; i != hTrueBin->GetNbinsX()+1; ++i) {
      
      int j = hpt->FindBin(hTrueBin->GetBinCenter(i));
      hcorrpt_bin->SetBinContent(j, hTrueBin->GetBinContent(i));
      hcorrpt_bin->SetBinError(j, hTrueBin->GetBinError(i));
    }
    
    bool _svd = true;// && id=="";// && y1<1.5;
    //TH1D *hcorrpt_svd(0);
    if (_svd) {
      int kreg = int(vy.size()/2); //4
      int ntoys = 0;//10; // default 1000
      RooUnfoldSvd *uSVD = new RooUnfoldSvd(uResps, hreco, kreg, ntoys);
      TH1D *hTrueSVD = (TH1D*)uSVD->Hreco();//RooUnfold::kCovariance);
      assert(hTrueSVD);
      hcorrpt_svd = (TH1D*)hTrueSVD->Clone(Form("hcorrpt_svd%s",c));
      delete uSVD; // ensure static members get destroyed before next instance
      
      hcorrpt_svd = (TH1D*)hpt->Clone(Form("hcorrpt_svd%s",c));
      hcorrpt_svd->Reset();
      for (int i = 1; i != hTrueSVD->GetNbinsX()+1; ++i) {
	
	int j = hpt->FindBin(hTrueSVD->GetBinCenter(i));
	hcorrpt_svd->SetBinContent(j, hTrueSVD->GetBinContent(i));
	hcorrpt_svd->SetBinError(j, hTrueSVD->GetBinError(i));
      }
    }
    else  {
      hcorrpt_svd = (TH1D*)hcorrpt_bin->Clone(Form("hcorrpt_svd%s",c));
      hcorrpt_svd->Reset();
    }
  } // !_jk
  
  if (_debug)
    cout << "done." << endl << flush;


  // Store "unfolding correction" (unfolded / original) and corrected graph
  //cout << "Store graphs..." << endl << flush;
  outdir->cd();
  TGraphErrors *gfold = new TGraphErrors(0);
  gfold->SetName(Form("gfold%s",c));
  TGraphErrors *gcorrpt = new TGraphErrors(0);
  gcorrpt->SetName(Form("gcorrpt%s",c));
  //
  TGraphErrors *gfold_bin(0), *gcorrpt_bin(0);
  TGraphErrors *gfold_svd(0), *gcorrpt_svd(0);
  if (!_jk && !_jet) {

    gfold_bin = new TGraphErrors(0);
    gfold_bin->SetName(Form("gfold_bin%s",c));
    gcorrpt_bin = new TGraphErrors(0);
    gcorrpt_bin->SetName(Form("gcorrpt_bin%s",c));
    gfold_svd = new TGraphErrors(0);
    gfold_svd->SetName(Form("gfold_svd%s",c));
    gcorrpt_svd = new TGraphErrors(0);
    gcorrpt_svd->SetName(Form("gcorrpt_svd%s",c));
  } // !_jk

  // Normalize hcorrpt
  for (int i = 1; i != hcorrpt->GetNbinsX()+1; ++i) {
    double dpt = hcorrpt->GetBinWidth(i);
    hcorrpt->SetBinContent(i, hcorrpt->GetBinContent(i) / dpt);
    hcorrpt->SetBinError(i, hcorrpt->GetBinError(i) / dpt);
  }
  if (!_jk && !_jet) {

    for (int i = 1; i != hcorrpt_bin->GetNbinsX()+1; ++i) {
      double dpt = hcorrpt_bin->GetBinWidth(i);
      hcorrpt_bin->SetBinContent(i, hcorrpt_bin->GetBinContent(i) / dpt);
      hcorrpt_bin->SetBinError(i, hcorrpt_bin->GetBinError(i) / dpt);
    }
    for (int i = 1; i != hcorrpt_bin->GetNbinsX()+1; ++i) {
      double dpt = hcorrpt_svd->GetBinWidth(i);
      hcorrpt_svd->SetBinContent(i, hcorrpt_svd->GetBinContent(i) / dpt);
      hcorrpt_svd->SetBinError(i, hcorrpt_svd->GetBinError(i) / dpt);
    }
  } // !_jk

  for (int i = 0; i != gpt->GetN(); ++i) {

    double x, y, ex, ey;
    tools::GetPoint(gpt, i, x, y, ex, ey);

    int j = hcorrpt->FindBin(x);
    double ycorr = hcorrpt->GetBinContent(j);
    double eycorr = hcorrpt->GetBinError(j);

    double k = (ycorr && y ? ycorr / y : 1.);

    if (!TMath::IsNaN(k)) {

      double dk = 0;
      //if (eycorr/ycorr > ey/y)
      dk = eycorr/ycorr*k;
      tools::SetPoint(gfold, gfold->GetN(), x, k, ex, dk);
      tools::SetPoint(gcorrpt, gcorrpt->GetN(), x, ycorr, ex, eycorr);
    }
  } // for i

  if (!_jk && !_jet) {

    for (int i = 0; i != gpt->GetN(); ++i) {
      
      double x, y, ex, ey;
      tools::GetPoint(gpt, i, x, y, ex, ey);
      
      int j = hcorrpt_bin->FindBin(x);
      double ycorr_bin = hcorrpt_bin->GetBinContent(j);
      double eycorr_bin = hcorrpt_bin->GetBinError(j);
      double ycorr_svd = hcorrpt_svd->GetBinContent(j);
      double eycorr_svd = hcorrpt_svd->GetBinError(j);
      
      double k_bin = (ycorr_bin && y ? ycorr_bin / y : 1.);
      double k_svd = (ycorr_svd && y ? ycorr_svd / y : 1.);
      
      if (!TMath::IsNaN(k_bin)) {
	
	double dk = 0;
	//if (eycorr_bin/ycorr_bin > ey/y)
	dk = eycorr_bin/ycorr_bin*k_bin;
	tools::SetPoint(gfold_bin, gfold_bin->GetN(), x, k_bin, ex, dk);
	tools::SetPoint(gcorrpt_bin, gcorrpt_bin->GetN(),
			x, ycorr_bin, ex, eycorr_bin);
      }

      if (!TMath::IsNaN(k_svd)) {
	
	double dk = 0;
	//if (eycorr_svd/ycorr_svd > ey/y)
	dk = eycorr_svd/ycorr_svd*k_svd;
	tools::SetPoint(gfold_svd, gfold_svd->GetN(), x, k_svd, ex, dk);
	tools::SetPoint(gcorrpt_svd, gcorrpt_svd->GetN(),
			x, ycorr_svd, ex, eycorr_svd);
      }
      
    } // for i
  } // !_jk

  if (!_jk && !_jet) {

    //uResp->Write();
    fnlo->Write();
    // Calculating points for fs is taking significant time,
    // and even got stuck at some point when too far out of the range
    fnlos->SetRange(xmin, min(xmax, 3000./cosh(y1)));
    fnlos->SetNpx(1000); // otherwise ugly on log x-axis after write
    fnlos->Write();
    gpt->Write();
    gfold->Write();
    gfold_fwd->Write();
    gfold_bin->Write();
    gfold_svd->Write();
    gcorrpt->Write();
    gcorrpt_fwd->Write();
    gcorrpt_bin->Write();
    gcorrpt_svd->Write();
    hpt->Write();
    hcorrpt->Write();
    hcorrpt_fwd->Write();
    hcorrpt_bin->Write();
    hcorrpt_svd->Write();
    
    hCov->Write();
  }
  else if (!_jk) {
    hpt->Write();
    hcorrpt->Write();
    hCov->Write();
  }
  else {
    //gcorrpt->Write();
    hcorrpt->Write();
  }
  
} // dagostiniUnfold_histo


void drawDagostini(string type) {

  TDirectory *curdir = gDirectory;
  setTDRStyle();

  TFile *f = new TFile(Form("output-%s-3c.root",type.c_str()),"READ");
  //TFile *f = new TFile(Form("GRV23_AK7_42x_v2/output-%s-3c.root",type.c_str()),"READ");
  //TFile *f = new TFile(Form("GRV23_AK5_42x_v2/output-%s-3c.root",type.c_str()),"READ");
  assert(f && !f->IsZombie());

  assert(f->cd("Standard"));
  TDirectory *din = gDirectory;
  curdir->cd();
    
  TCanvas *c1 = new TCanvas("c1","c1",1200,800);
  //c1->SetTopMargin(0.20);
  c1->Divide(3,2);//,0.1,0.00);

  TCanvas *c1b = new TCanvas("c1b","c1b",600,600);

  TCanvas *c2 = new TCanvas("c2","c2",1200,800);
  c2->SetTopMargin(0.10);
  c2->Divide(3,2,-1,-1);

  TCanvas *c3 = new TCanvas("c3","c3",1200,800);
  c3->SetTopMargin(0.10);
  c3->Divide(3,2,-1,-1);

  TH1D *h = new TH1D("h",";p_{T} (GeV);Unfolding correction",
		     int(xmax-xmin),xmin,xmax);
  h->SetMinimum(0.45);
  h->SetMaximum(1.15);
  h->GetXaxis()->SetMoreLogLabels();
  h->GetXaxis()->SetNoExponent();

  TLegend *leg = new TLegend(0.25,0.74,0.45,0.97,"","brNDC");
  leg->SetFillStyle(kNone);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.045);

  TLatex *tex = new TLatex();
  tex->SetTextSize(0.045);

  const int ny = 6;
  for (int iy = 0; iy != ny; ++iy) {

    double y1 = 0.5*iy; double y2 = 0.5*(iy+1);
    assert(din->cd(Form("Eta_%1.1f-%1.1f",y1,y2)));
    TDirectory *d = gDirectory;

    TH1D *hreco = (TH1D*)d->Get("hreco"); assert(hreco);
    TH2D *h2resp = (TH2D*)d->Get("mtu"); assert(h2resp);
    TGraphErrors *gfwd = (TGraphErrors*)d->Get("gfold_fwd"); assert(gfwd);
    TGraphErrors *gbin = (TGraphErrors*)d->Get("gfold_bin"); assert(gbin);
    TGraphErrors *gsvd = (TGraphErrors*)d->Get("gfold_svd"); assert(gsvd);
    TGraphErrors *gbayes = (TGraphErrors*)d->Get("gfold"); assert(gbayes);
    curdir->cd();

    // Find the low end of reconstructed data used in unfolding
    double ptmin(0);
    for (int i = 1; i != hreco->GetNbinsX()+1 && !ptmin; ++i) {
      if (hreco->GetBinContent(i)>0) ptmin = hreco->GetBinLowEdge(i);
    }

    c1->cd(iy+1);
    gPad->SetLogx();
    gPad->SetLogy();
    gStyle->SetPalette(1);
    gPad->SetRightMargin(0.10);
    h2resp->GetXaxis()->SetMoreLogLabels();
    h2resp->GetXaxis()->SetNoExponent();
    h2resp->GetXaxis()->SetTitle("Measured p_{T,reco} (GeV)");
    h2resp->GetXaxis()->SetTitleOffset(1.5);
    h2resp->GetYaxis()->SetMoreLogLabels();
    h2resp->GetYaxis()->SetNoExponent();
    h2resp->GetYaxis()->SetTitle("True p_{T,gen} (GeV)");
    h2resp->GetYaxis()->SetTitleOffset(2.0);
    h2resp->GetZaxis()->SetRangeUser(1e-3,0.9999);//1.05);
    h2resp->DrawClone("COLZ");
    //h2resp->Draw("SAMEBOX");

    tex->SetTextSize(0.045);
    tex->SetTextAlign(31); // align right
    tex->SetNDC(kTRUE);
    tex->DrawLatex(0.8, 0.2, y1==0 ? "|y|<0.5" : Form("%1.1f<|y|<%1.1f",y1,y2));

    if (y1==0) {
      c1b->cd();
      gPad->SetLogx();
      gPad->SetLogy();
      gPad->SetRightMargin(0.13);//0.11);
      gPad->SetBottomMargin(0.14);//0.10);
      gPad->SetLeftMargin(0.18);//0.10);
      //h2resp->GetXaxis()->SetTitleOffset(1.4);
      //TH2D *h2 = (TH2D*)h2resp->DrawClone("COLZ");

      // copy TH2D over to a fresh one to get default graphics style
      TH2D *h2 = new TH2D("h2",";Measured Jet p_{T} (GeV);"
			  "True Jet p_{T} (GeV)",
			  h2resp->GetNbinsX(),
			  h2resp->GetXaxis()->GetXbins()->GetArray(),
			  h2resp->GetNbinsY(),
			  h2resp->GetYaxis()->GetXbins()->GetArray());
      for (int i = 1; i != h2->GetNbinsX()+1; ++i) {
	for (int j = 1; j != h2->GetNbinsY()+1; ++j) {
	  h2->SetBinContent(i, j, h2resp->GetBinContent(i, j));
	  //h2->SetBinError(i, j, h2resp->GetBinError(i, j));
	} // for j
      } // for i
      h2->Draw("COLZ");
      h2->GetXaxis()->SetMoreLogLabels();
      h2->GetXaxis()->SetNoExponent();
      //h2->GetYaxis()->SetMoreLogLabels();
      h2->GetYaxis()->SetNoExponent();
      h2->GetYaxis()->SetTitleOffset(1.5);
      h2->GetZaxis()->SetRangeUser(1e-3,0.9999);

      h2->GetXaxis()->SetRangeUser(fitptmin,1890.);//1327.);
      h2->GetYaxis()->SetRangeUser(_recopt,1890.);
      //h2resp->Draw("SAME COLZ");

      //tex->DrawLatex(0.8, 0.15, _algo=="AK7" ?
      tex->DrawLatex(0.8, 0.20, _algo=="AK7" ?
		     "Anti-k_{T} R=0.7, |y|<0.5" :
		     "Anti-k_{T} R=0.5, |y|<0.5");
      //cmsPrel(_lumi);
      cmsPrel(0); // simulation
    }

    c2->cd(iy+1);
    gPad->SetLogx();

    h->SetMinimum(0.45);
    h->SetMaximum(1.15);
    h->SetXTitle(iy==ny-1 ? "p_{T} (GeV)" : "");
    h->SetYTitle(iy==0 ? "Unfolding correction" : "");
    h->DrawClone("AXIS");
    
    TLine *l = new TLine();
    l->SetLineStyle(kDotted);
    l->DrawLine(ptmin,0.45,ptmin,1.15);
    l->SetLineStyle(kDashed);
    //l->DrawLine(56,0.45,56,1.15); // v2
    l->DrawLine(xminpas,0.45,xminpas,1.15);

    gfwd->SetName("gfwd");
    gfwd->SetLineWidth(2);
    gfwd->SetLineColor(kRed);
    gfwd->Draw("SAMEL");

    gbin->SetName("gbin");
    gbin->SetMarkerStyle(kOpenSquare);
    gbin->SetMarkerColor(kGreen+1);
    gbin->SetLineColor(kGreen+1);
    gbin->Draw("SAMEPz");

    gbayes->SetName("gbayes");
    gbayes->SetMarkerStyle(kFullCircle);
    gbayes->Draw("SAMEP");

    gsvd->SetName("gsvd");
    gsvd->SetMarkerStyle(kOpenDiamond);
    gsvd->SetMarkerColor(kCyan+1);
    gsvd->SetLineColor(kCyan+1);
    gsvd->Draw("SAMEP");

    tex->SetTextSize(iy<3 ? 0.053 : 0.045);
    tex->SetTextAlign(21); // align middle
    tex->SetNDC(kFALSE);
    tex->DrawLatex(150, 0.5, y1==0 ? "|y|<0.5" : Form("%1.1f<|y|<%1.1f",y1,y2));

    if (iy==2) {
      leg->AddEntry(gbayes,"RooUnfoldBayes","P");
      leg->AddEntry(gbin,"RooUnfoldBinByBin","P");
      leg->AddEntry(gsvd,"RooUnfoldSvd","P");
      leg->AddEntry(gfwd,"Forward smearing","L");
      leg->Draw();
    }

    c3->cd(iy+1);
    gPad->SetLogx();

    h->SetMinimum(0.90+0.0001);
    h->SetMaximum(1.10-0.0001);
    h->SetXTitle(iy==ny-1 ? "p_{T} (GeV)" : "");
    h->SetYTitle(iy==0 ? "Correction / Forward smearing" : "");
    h->DrawClone("AXIS");

    l->SetLineStyle(kDotted);
    l->DrawLine(ptmin,0.9,ptmin,1.1);
    l->SetLineStyle(kDashed);
    //l->DrawLine(56,0.9,56,1.1); // v2
    l->DrawLine(xminpas,0.9,xminpas,1.1);

    TGraphErrors *grfwd = tools::ratioGraphs(gfwd, gfwd);
    grfwd->Draw("SAMEL");
    TGraphErrors *grbin = tools::ratioGraphs(gbin, gfwd);
    grbin->Draw("SAMEPz");
    TGraphErrors *grbayes = tools::ratioGraphs(gbayes, gfwd);
    grbayes->Draw("SAMEPz");
    TGraphErrors *grsvd = tools::ratioGraphs(gsvd, gfwd);
    grsvd->Draw("SAMEPz");

    tex->DrawLatex(150,0.91, y1==0 ? "|y|<0.5" : Form("%1.1f<|y|<%1.1f",y1,y2));

    if (iy==2) leg->Draw();

  }

  const char *a = _algo.c_str();
  const char *t = type.c_str();

  c1->cd(3);
  tex->SetTextSize(0.045);
  tex->SetNDC(kTRUE);
  tex->DrawLatex(0.35, 0.85, Form("%s %s",t,_algo.c_str()));
  c1->cd(0);
  //cmsPrel(type=="DATA" ? _lumi : 0, true);
  c1->SaveAs(Form("pdf/roounfold_matrix_%s_%s.pdf",a,t));
  c1b->SaveAs(Form("pdf/roounfold_matrix0_%s_%s.pdf",a,t));

  c2->cd(2);
  tex->SetTextSize(0.053);
  tex->DrawLatex(0.50, 0.85, Form("%s %s",t,_algo.c_str()));
  c2->cd(0);
  cmsPrel(type=="DATA" ? _lumi : 0, true);
  c2->SaveAs(Form("pdf/roounfold_comparison_%s_%s.pdf",a,t));

  c3->cd(2);
  tex->DrawLatex(0.50, 0.85, Form("%s %s",t,_algo.c_str()));
  c3->cd(0);
  cmsPrel(type=="DATA" ? _lumi : 0, true);
  c3->SaveAs(Form("pdf/roounfold_ratiotofwd_%s_%s.pdf",a,t));
  

} // drawDagostini
