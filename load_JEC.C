{
  //  gROOT->ProcessLine(".exception"); //  let exceptions pass CINT, easier to read
  
  // Need to add this - matt
  gROOT->ProcessLine(".L FWCore/Utilities/src/typelookup.cc+");
  // For JEC residual (and pile-up)
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/Utilities.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/JetCorrectorParameters.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/SimpleJetCorrector.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/FactorizedJetCorrector.cc+");
  // For JEC uncertainty
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/SimpleJetCorrectionUncertainty.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/JetCorrectionUncertainty.cc+");

  // looks like Mikko's private stuff
  //  gROOT->ProcessLine(".L fillHistosInclB.C++g"); // g to compile in debug symbols 


}
