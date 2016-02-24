#ifdef __CINT__

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ nestedtypedefs;

#pragma link C++ class TKDInterpolatorBase+;
#pragma link C++ class TKDNodeInfo+;
#pragma link C++ class TKDNodeInfo::TKDNodeDraw+;
#pragma link C++ class TKDPDF+;
#pragma link C++ class TKDInterpolator+;
#pragma link C++ class TKDSpline+;

#pragma link C++ class AliTMinuitToolkit+;
#pragma link C++ class AliNDLocalRegression+;

#pragma link C++ class AliExternalInfo+;

#pragma link C++ namespace TStatToolkit;
#pragma link C++ function  TStatToolkit::LTMHisto(TH1 * , TVectorF &, Float_t);
#pragma link C++ function  TStatToolkit::LTMHisto(TH1 * , TVectorD &, Float_t);

/*
// RS At the moment is not recognized by the CINT...
#pragma link C++ function TStatToolkit::GetMinMax(const Double_t*, Long64_t, double &, double &);
#pragma link C++ function TStatToolkit::GetMinMax(const Float_t*, Long64_t, double &, double &);
#pragma link C++ function TStatToolkit::GetMinMax(const Int_t*, Long64_t, double &, double &);

#pragma link C++ function TStatToolkit::GetMinMaxMean(const Double_t*, Long64_t, double &, double &, double &);
#pragma link C++ function TStatToolkit::GetMinMaxMean(const Float_t*, Long64_t, double &, double &, double &);
#pragma link C++ function TStatToolkit::GetMinMaxMean(const Int_t*, Long64_t, double &, double &, double&);

#pragma link C++ function TStatToolkit::TruncatedMean(const TH1*, TVectorF *,Float_t,Float_t,Bool_t);
#pragma link C++ function TStatToolkit::TruncatedMean(const TH1*, TVectorD *,Float_t,Float_t,Bool_t);

#pragma link C++ function TStatToolkit::LTMHisto(TH1*,TVectorF &,Float_t);
#pragma link C++ function TStatToolkit::LTMHisto(TH1*,TVectorD &,Float_t);

#pragma link C++ function TStatToolkit::LTM(TH1*,TVectorF*,Float_t,Bool_t);
#pragma link C++ function TStatToolkit::LTM(TH1*,TVectorD*,Float_t,Bool_t);

#pragma link C++ function TStatToolkit::FitGaus(TH1*,TVectorF*,TMatrixF*,Float_t,Float_t,Bool_t);
#pragma link C++ function TStatToolkit::FitGaus(TH1*,TVectorD*,TMatrixD*,Float_t,Float_t,Bool_t);

#pragma link C++ function TStatToolkit::FitGaus(Float_t*,Int_t,Float_t,Float_t,TVectorF*,TMatrixF*,Bool_t);
#pragma link C++ function TStatToolkit::FitGaus(Float_t*,Int_t,Float_t,Float_t,TVectorD*,TMatrixD*,Bool_t);

*/

#endif
