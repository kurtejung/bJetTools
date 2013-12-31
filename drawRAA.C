{
  
  TH1D *h = new TH1D();
  TH1D *h_den = new TH1D();
  h->Sumw2();
  h_den->Sumw2();
  
  TFile *f1 = new TFile("raa_inclJet_pPb_numberator.root","OLD");
  f1->cd();
  //h = (TH1D*)(f1->Get("hRawBData"));
  h = (TH1D*)(f1->Get("hIncJetsData"));
  h->Scale(1E-4./0.096); //TpA from http://arxiv.org/pdf/nucl-ex/0302016v3.pdf
  
  TFile *f2 = new TFile("raa_inclJet_pp_denomForpA.root","OLD");
  f2->cd();
  //h_den = (TH1D*)(f2->Get("hRawBData"));
  h_den = (TH1D*)(f2->Get("hIncJetsData"));
  
  h->Draw();
  h->SetTitle("");
  //h_den->SetMarkerColor(2);
  //h_den->Draw("same");
  h->Divide(h_den);
  h->SetMaximum(3);
  h->SetMinimum(0);
  h->SetStats(0);
  h->SetYTitle("#frac{1}{<N_{coll}>} #frac{d#sigma_{pA}/dp_{T}}{d#sigma_{pp}/dp_{T}}");
  h->SetTitle("Inclusive Jet R_{pA}");
  h->Draw("p");
  
}
