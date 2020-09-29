#ifndef ALIHFDHADRONCORRSYSTUNC_H
#define ALIHFDHADRONCORRSYSTUNC_H
/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: $ */

/////////////////////////////////////////////////////////////
// class for systematic uncertainties on D meson -hadron correlation distribution
//
// Author: A. Rossi, andrea.rossi@cern.ch
//
// Responsible of the values set for the different mesons/datasets
//      D0:  in pp (2010 min bias) Fabio Colamaria, fabio.colamaria@ba.infn.it ;  p-Pb (2013 min bias): Fabio Colamaria, fabio.colamaria@ba.infn.it,  Somnath Kar, somnath.kar@cern.ch
//      D*+: in pp 2010 min. bias and p-Pb 2013 min. bias  Sandro Bjelogrlic, sandro.bjelogrlic@cern.ch
//      D+:  in pp 2010 min. bias and p-Pb 2013 min. bias  Jitendra Kumar, jitendra.kumar@cern.ch
//
/////////////////////////////////////////////////////////////
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <TString.h>
class AliHFDhadronCorrSystUnc : public TNamed{
  
 public:
  AliHFDhadronCorrSystUnc();
  AliHFDhadronCorrSystUnc(const char* name);
  ~AliHFDhadronCorrSystUnc();
  void InitEmptyHistosFromTemplate();
  void InitStandardUncertaintiesPP2010(Int_t meson,Double_t ptD,Double_t minptAss,Double_t maxptAss);
  void InitStandardUncertaintiesPPb2013(Int_t meson,Double_t ptD,Double_t minptAss,Double_t maxptAss);
  void InitStandardUncertaintiesPPb2016(Int_t meson,Double_t ptD,Double_t minptAss,Double_t maxptAss);
  void InitStandardUncertaintiesPPb2016in020(Int_t meson,Double_t ptD,Double_t minptAss,Double_t maxptAss);
  void InitStandardUncertaintiesPPb2016in2060(Int_t meson,Double_t ptD,Double_t minptAss,Double_t maxptAss);
  void InitStandardUncertaintiesPPb2016in60100(Int_t meson,Double_t ptD,Double_t minptAss,Double_t maxptAss);
  void InitStandardUncertaintiesPP2017(Int_t meson,Double_t ptD,Double_t minptAss,Double_t maxptAss);
  void InitStandardUncertaintiesPP13TeV(Int_t meson,Double_t ptD,Double_t minptAss,Double_t maxptAss);
  void InitStandardUncertaintiesPP13TeVin001(Int_t meson,Double_t ptD,Double_t minptAss,Double_t maxptAss);
  void InitStandardUncertaintiesPP13TeVin0110(Int_t meson,Double_t ptD,Double_t minptAss,Double_t maxptAss);
  void InitStandardUncertaintiesPP13TeVin1030(Int_t meson,Double_t ptD,Double_t minptAss,Double_t maxptAss);
  void InitStandardUncertaintiesPP13TeVin30100(Int_t meson,Double_t ptD,Double_t minptAss,Double_t maxptAss);  
  
  // Method with uncertainties for pp 2010, Dzero and D*+ and pt assoc> 0.3 GeV/c, with values used for HP2013
  void InitStandardUncertaintiesPP2010DzeroLowPtAss03();
  void InitStandardUncertaintiesPP2010DzeroMidPtAss03();
  void InitStandardUncertaintiesPP2010DzeroHighPtAss03();
  
  void InitStandardUncertaintiesPP2010DstarLowPtAss03();
  void InitStandardUncertaintiesPP2010DstarMidPtAss03();
  void InitStandardUncertaintiesPP2010DstarHighPtAss03();


  // Method with uncertainties for pp 2010, all kinematic cases but those approved for HP2013
  void InitStandardUncertaintiesPP2010DplusLowPtAss03();
  void InitStandardUncertaintiesPP2010DplusMidPtAss03();
  void InitStandardUncertaintiesPP2010DplusHighPtAss03();

  void InitStandardUncertaintiesPP2010DzeroLowPtAss03to1();
  void InitStandardUncertaintiesPP2010DzeroMidPtAss03to1();
  void InitStandardUncertaintiesPP2010DzeroHighPtAss03to1();

  void InitStandardUncertaintiesPP2010DstarLowPtAss03to1();
  void InitStandardUncertaintiesPP2010DstarMidPtAss03to1();
  void InitStandardUncertaintiesPP2010DstarHighPtAss03to1();

  void InitStandardUncertaintiesPP2010DplusLowPtAss03to1();
  void InitStandardUncertaintiesPP2010DplusMidPtAss03to1();
  void InitStandardUncertaintiesPP2010DplusHighPtAss03to1();


  void InitStandardUncertaintiesPP2010DzeroLowPtAss1();
  void InitStandardUncertaintiesPP2010DzeroMidPtAss1();
  void InitStandardUncertaintiesPP2010DzeroHighPtAss1();

  void InitStandardUncertaintiesPP2010DstarLowPtAss1();
  void InitStandardUncertaintiesPP2010DstarMidPtAss1();
  void InitStandardUncertaintiesPP2010DstarHighPtAss1();

  void InitStandardUncertaintiesPP2010DplusLowPtAss1();
  void InitStandardUncertaintiesPP2010DplusMidPtAss1();
  void InitStandardUncertaintiesPP2010DplusHighPtAss1();


  // Method with uncertainties for pPb 2013
  void InitStandardUncertaintiesPPb2013DzeroLowPtAss03();
  void InitStandardUncertaintiesPPb2013DzeroMidPtAss03();
  void InitStandardUncertaintiesPPb2013DzeroHighPtAss03();

  void InitStandardUncertaintiesPPb2013DstarLowPtAss03();
  void InitStandardUncertaintiesPPb2013DstarMidPtAss03();
  void InitStandardUncertaintiesPPb2013DstarHighPtAss03();

  void InitStandardUncertaintiesPPb2013DplusLowPtAss03();
  void InitStandardUncertaintiesPPb2013DplusMidPtAss03();
  void InitStandardUncertaintiesPPb2013DplusHighPtAss03();

  void InitStandardUncertaintiesPPb2013DzeroLowPtAss03to1();
  void InitStandardUncertaintiesPPb2013DzeroMidPtAss03to1();
  void InitStandardUncertaintiesPPb2013DzeroHighPtAss03to1();

  void InitStandardUncertaintiesPPb2013DstarLowPtAss03to1();
  void InitStandardUncertaintiesPPb2013DstarMidPtAss03to1();
  void InitStandardUncertaintiesPPb2013DstarHighPtAss03to1();

  void InitStandardUncertaintiesPPb2013DplusLowPtAss03to1();
  void InitStandardUncertaintiesPPb2013DplusMidPtAss03to1();
  void InitStandardUncertaintiesPPb2013DplusHighPtAss03to1();


  void InitStandardUncertaintiesPPb2013DzeroLowPtAss1();
  void InitStandardUncertaintiesPPb2013DzeroMidPtAss1();
  void InitStandardUncertaintiesPPb2013DzeroHighPtAss1();

  void InitStandardUncertaintiesPPb2013DstarLowPtAss1();
  void InitStandardUncertaintiesPPb2013DstarMidPtAss1();
  void InitStandardUncertaintiesPPb2013DstarHighPtAss1();

  void InitStandardUncertaintiesPPb2013DplusLowPtAss1();
  void InitStandardUncertaintiesPPb2013DplusMidPtAss1();
  void InitStandardUncertaintiesPPb2013DplusHighPtAss1();


  // Method with uncertainties for pPb 2016
  void InitStandardUncertaintiesPPb2016DzeroLowPtAss03to99();
  void InitStandardUncertaintiesPPb2016DzeroMidPtAss03to99();
  void InitStandardUncertaintiesPPb2016DzeroHighPtAss03to99();
  void InitStandardUncertaintiesPPb2016DzeroVeryHighPtAss03to99();

  void InitStandardUncertaintiesPPb2016DstarLowPtAss03to99();
  void InitStandardUncertaintiesPPb2016DstarMidPtAss03to99();
  void InitStandardUncertaintiesPPb2016DstarHighPtAss03to99();
  void InitStandardUncertaintiesPPb2016DstarVeryHighPtAss03to99();

  void InitStandardUncertaintiesPPb2016DplusLowPtAss03to99();
  void InitStandardUncertaintiesPPb2016DplusMidPtAss03to99();
  void InitStandardUncertaintiesPPb2016DplusHighPtAss03to99();
  void InitStandardUncertaintiesPPb2016DplusVeryHighPtAss03to99();

  void InitStandardUncertaintiesPPb2016DzeroLowPtAss03to1();
  void InitStandardUncertaintiesPPb2016DzeroMidPtAss03to1();
  void InitStandardUncertaintiesPPb2016DzeroHighPtAss03to1();
  void InitStandardUncertaintiesPPb2016DzeroVeryHighPtAss03to1();

  void InitStandardUncertaintiesPPb2016DstarLowPtAss03to1();
  void InitStandardUncertaintiesPPb2016DstarMidPtAss03to1();
  void InitStandardUncertaintiesPPb2016DstarHighPtAss03to1();
  void InitStandardUncertaintiesPPb2016DstarVeryHighPtAss03to1();

  void InitStandardUncertaintiesPPb2016DplusLowPtAss03to1();
  void InitStandardUncertaintiesPPb2016DplusMidPtAss03to1();
  void InitStandardUncertaintiesPPb2016DplusHighPtAss03to1();
  void InitStandardUncertaintiesPPb2016DplusVeryHighPtAss03to1();

  void InitStandardUncertaintiesPPb2016DzeroLowPtAss1to99();
  void InitStandardUncertaintiesPPb2016DzeroMidPtAss1to99();
  void InitStandardUncertaintiesPPb2016DzeroHighPtAss1to99();
  void InitStandardUncertaintiesPPb2016DzeroVeryHighPtAss1to99();

  void InitStandardUncertaintiesPPb2016DstarLowPtAss1to99();
  void InitStandardUncertaintiesPPb2016DstarMidPtAss1to99();
  void InitStandardUncertaintiesPPb2016DstarHighPtAss1to99();
  void InitStandardUncertaintiesPPb2016DstarVeryHighPtAss1to99();

  void InitStandardUncertaintiesPPb2016DplusLowPtAss1to99();
  void InitStandardUncertaintiesPPb2016DplusMidPtAss1to99();
  void InitStandardUncertaintiesPPb2016DplusHighPtAss1to99();
  void InitStandardUncertaintiesPPb2016DplusVeryHighPtAss1to99();

  void InitStandardUncertaintiesPPb2016DzeroLowPtAss2to99();
  void InitStandardUncertaintiesPPb2016DzeroMidPtAss2to99();
  void InitStandardUncertaintiesPPb2016DzeroHighPtAss2to99();
  void InitStandardUncertaintiesPPb2016DzeroVeryHighPtAss2to99();

  void InitStandardUncertaintiesPPb2016DstarLowPtAss2to99();
  void InitStandardUncertaintiesPPb2016DstarMidPtAss2to99();
  void InitStandardUncertaintiesPPb2016DstarHighPtAss2to99();
  void InitStandardUncertaintiesPPb2016DstarVeryHighPtAss2to99();

  void InitStandardUncertaintiesPPb2016DplusLowPtAss2to99();
  void InitStandardUncertaintiesPPb2016DplusMidPtAss2to99();
  void InitStandardUncertaintiesPPb2016DplusHighPtAss2to99();
  void InitStandardUncertaintiesPPb2016DplusVeryHighPtAss2to99();

  void InitStandardUncertaintiesPPb2016DzeroLowPtAss3to99();
  void InitStandardUncertaintiesPPb2016DzeroMidPtAss3to99();
  void InitStandardUncertaintiesPPb2016DzeroHighPtAss3to99();
  void InitStandardUncertaintiesPPb2016DzeroVeryHighPtAss3to99();

  void InitStandardUncertaintiesPPb2016DstarLowPtAss3to99();
  void InitStandardUncertaintiesPPb2016DstarMidPtAss3to99();
  void InitStandardUncertaintiesPPb2016DstarHighPtAss3to99();
  void InitStandardUncertaintiesPPb2016DstarVeryHighPtAss3to99();

  void InitStandardUncertaintiesPPb2016DplusLowPtAss3to99();
  void InitStandardUncertaintiesPPb2016DplusMidPtAss3to99();
  void InitStandardUncertaintiesPPb2016DplusHighPtAss3to99();
  void InitStandardUncertaintiesPPb2016DplusVeryHighPtAss3to99();

  void InitStandardUncertaintiesPPb2016DzeroLowPtAss1to2();
  void InitStandardUncertaintiesPPb2016DzeroMidPtAss1to2();
  void InitStandardUncertaintiesPPb2016DzeroHighPtAss1to2();
  void InitStandardUncertaintiesPPb2016DzeroVeryHighPtAss1to2();

  void InitStandardUncertaintiesPPb2016DstarLowPtAss1to2();
  void InitStandardUncertaintiesPPb2016DstarMidPtAss1to2();
  void InitStandardUncertaintiesPPb2016DstarHighPtAss1to2();
  void InitStandardUncertaintiesPPb2016DstarVeryHighPtAss1to2();

  void InitStandardUncertaintiesPPb2016DplusLowPtAss1to2();
  void InitStandardUncertaintiesPPb2016DplusMidPtAss1to2();
  void InitStandardUncertaintiesPPb2016DplusHighPtAss1to2();
  void InitStandardUncertaintiesPPb2016DplusVeryHighPtAss1to2();

  void InitStandardUncertaintiesPPb2016DzeroLowPtAss2to3();
  void InitStandardUncertaintiesPPb2016DzeroMidPtAss2to3();
  void InitStandardUncertaintiesPPb2016DzeroHighPtAss2to3();
  void InitStandardUncertaintiesPPb2016DzeroVeryHighPtAss2to3();

  void InitStandardUncertaintiesPPb2016DstarLowPtAss2to3();
  void InitStandardUncertaintiesPPb2016DstarMidPtAss2to3();
  void InitStandardUncertaintiesPPb2016DstarHighPtAss2to3();
  void InitStandardUncertaintiesPPb2016DstarVeryHighPtAss2to3();

  void InitStandardUncertaintiesPPb2016DplusLowPtAss2to3();
  void InitStandardUncertaintiesPPb2016DplusMidPtAss2to3();
  void InitStandardUncertaintiesPPb2016DplusHighPtAss2to3();
  void InitStandardUncertaintiesPPb2016DplusVeryHighPtAss2to3();  

  // Method with uncertainties for pPb 2016 0-20%
  void InitStandardUncertaintiesPPb2016DzeroLowPtAss03to99in020();
  void InitStandardUncertaintiesPPb2016DzeroMidPtAss03to99in020();
  void InitStandardUncertaintiesPPb2016DzeroHighPtAss03to99in020();
  void InitStandardUncertaintiesPPb2016DzeroVeryHighPtAss03to99in020();

  void InitStandardUncertaintiesPPb2016DstarLowPtAss03to99in020();
  void InitStandardUncertaintiesPPb2016DstarMidPtAss03to99in020();
  void InitStandardUncertaintiesPPb2016DstarHighPtAss03to99in020();
  void InitStandardUncertaintiesPPb2016DstarVeryHighPtAss03to99in020();

  void InitStandardUncertaintiesPPb2016DplusLowPtAss03to99in020();
  void InitStandardUncertaintiesPPb2016DplusMidPtAss03to99in020();
  void InitStandardUncertaintiesPPb2016DplusHighPtAss03to99in020();
  void InitStandardUncertaintiesPPb2016DplusVeryHighPtAss03to99in020();

  void InitStandardUncertaintiesPPb2016DzeroLowPtAss03to1in020();
  void InitStandardUncertaintiesPPb2016DzeroMidPtAss03to1in020();
  void InitStandardUncertaintiesPPb2016DzeroHighPtAss03to1in020();
  void InitStandardUncertaintiesPPb2016DzeroVeryHighPtAss03to1in020();

  void InitStandardUncertaintiesPPb2016DstarLowPtAss03to1in020();
  void InitStandardUncertaintiesPPb2016DstarMidPtAss03to1in020();
  void InitStandardUncertaintiesPPb2016DstarHighPtAss03to1in020();
  void InitStandardUncertaintiesPPb2016DstarVeryHighPtAss03to1in020();

  void InitStandardUncertaintiesPPb2016DplusLowPtAss03to1in020();
  void InitStandardUncertaintiesPPb2016DplusMidPtAss03to1in020();
  void InitStandardUncertaintiesPPb2016DplusHighPtAss03to1in020();
  void InitStandardUncertaintiesPPb2016DplusVeryHighPtAss03to1in020();

  void InitStandardUncertaintiesPPb2016DzeroLowPtAss1to99in020();
  void InitStandardUncertaintiesPPb2016DzeroMidPtAss1to99in020();
  void InitStandardUncertaintiesPPb2016DzeroHighPtAss1to99in020();
  void InitStandardUncertaintiesPPb2016DzeroVeryHighPtAss1to99in020();

  void InitStandardUncertaintiesPPb2016DstarLowPtAss1to99in020();
  void InitStandardUncertaintiesPPb2016DstarMidPtAss1to99in020();
  void InitStandardUncertaintiesPPb2016DstarHighPtAss1to99in020();
  void InitStandardUncertaintiesPPb2016DstarVeryHighPtAss1to99in020();

  void InitStandardUncertaintiesPPb2016DplusLowPtAss1to99in020();
  void InitStandardUncertaintiesPPb2016DplusMidPtAss1to99in020();
  void InitStandardUncertaintiesPPb2016DplusHighPtAss1to99in020();
  void InitStandardUncertaintiesPPb2016DplusVeryHighPtAss1to99in020();

    // Method with uncertainties for pPb 2016 20-60%
  void InitStandardUncertaintiesPPb2016DzeroLowPtAss03to99in2060();
  void InitStandardUncertaintiesPPb2016DzeroMidPtAss03to99in2060();
  void InitStandardUncertaintiesPPb2016DzeroHighPtAss03to99in2060();
  void InitStandardUncertaintiesPPb2016DzeroVeryHighPtAss03to99in2060();

  void InitStandardUncertaintiesPPb2016DstarLowPtAss03to99in2060();
  void InitStandardUncertaintiesPPb2016DstarMidPtAss03to99in2060();
  void InitStandardUncertaintiesPPb2016DstarHighPtAss03to99in2060();
  void InitStandardUncertaintiesPPb2016DstarVeryHighPtAss03to99in2060();

  void InitStandardUncertaintiesPPb2016DplusLowPtAss03to99in2060();
  void InitStandardUncertaintiesPPb2016DplusMidPtAss03to99in2060();
  void InitStandardUncertaintiesPPb2016DplusHighPtAss03to99in2060();
  void InitStandardUncertaintiesPPb2016DplusVeryHighPtAss03to99in2060();

  void InitStandardUncertaintiesPPb2016DzeroLowPtAss03to1in2060();
  void InitStandardUncertaintiesPPb2016DzeroMidPtAss03to1in2060();
  void InitStandardUncertaintiesPPb2016DzeroHighPtAss03to1in2060();
  void InitStandardUncertaintiesPPb2016DzeroVeryHighPtAss03to1in2060();

  void InitStandardUncertaintiesPPb2016DstarLowPtAss03to1in2060();
  void InitStandardUncertaintiesPPb2016DstarMidPtAss03to1in2060();
  void InitStandardUncertaintiesPPb2016DstarHighPtAss03to1in2060();
  void InitStandardUncertaintiesPPb2016DstarVeryHighPtAss03to1in2060();

  void InitStandardUncertaintiesPPb2016DplusLowPtAss03to1in2060();
  void InitStandardUncertaintiesPPb2016DplusMidPtAss03to1in2060();
  void InitStandardUncertaintiesPPb2016DplusHighPtAss03to1in2060();
  void InitStandardUncertaintiesPPb2016DplusVeryHighPtAss03to1in2060();

  void InitStandardUncertaintiesPPb2016DzeroLowPtAss1to99in2060();
  void InitStandardUncertaintiesPPb2016DzeroMidPtAss1to99in2060();
  void InitStandardUncertaintiesPPb2016DzeroHighPtAss1to99in2060();
  void InitStandardUncertaintiesPPb2016DzeroVeryHighPtAss1to99in2060();

  void InitStandardUncertaintiesPPb2016DstarLowPtAss1to99in2060();
  void InitStandardUncertaintiesPPb2016DstarMidPtAss1to99in2060();
  void InitStandardUncertaintiesPPb2016DstarHighPtAss1to99in2060();
  void InitStandardUncertaintiesPPb2016DstarVeryHighPtAss1to99in2060();

  void InitStandardUncertaintiesPPb2016DplusLowPtAss1to99in2060();
  void InitStandardUncertaintiesPPb2016DplusMidPtAss1to99in2060();
  void InitStandardUncertaintiesPPb2016DplusHighPtAss1to99in2060();
  void InitStandardUncertaintiesPPb2016DplusVeryHighPtAss1to99in2060();

    // Method with uncertainties for pPb 2016 60-100%
  void InitStandardUncertaintiesPPb2016DzeroLowPtAss03to99in60100();
  void InitStandardUncertaintiesPPb2016DzeroMidPtAss03to99in60100();
  void InitStandardUncertaintiesPPb2016DzeroHighPtAss03to99in60100();
  void InitStandardUncertaintiesPPb2016DzeroVeryHighPtAss03to99in60100();

  void InitStandardUncertaintiesPPb2016DstarLowPtAss03to99in60100();
  void InitStandardUncertaintiesPPb2016DstarMidPtAss03to99in60100();
  void InitStandardUncertaintiesPPb2016DstarHighPtAss03to99in60100();
  void InitStandardUncertaintiesPPb2016DstarVeryHighPtAss03to99in60100();

  void InitStandardUncertaintiesPPb2016DplusLowPtAss03to99in60100();
  void InitStandardUncertaintiesPPb2016DplusMidPtAss03to99in60100();
  void InitStandardUncertaintiesPPb2016DplusHighPtAss03to99in60100();
  void InitStandardUncertaintiesPPb2016DplusVeryHighPtAss03to99in60100();

  void InitStandardUncertaintiesPPb2016DzeroLowPtAss03to1in60100();
  void InitStandardUncertaintiesPPb2016DzeroMidPtAss03to1in60100();
  void InitStandardUncertaintiesPPb2016DzeroHighPtAss03to1in60100();
  void InitStandardUncertaintiesPPb2016DzeroVeryHighPtAss03to1in60100();

  void InitStandardUncertaintiesPPb2016DstarLowPtAss03to1in60100();
  void InitStandardUncertaintiesPPb2016DstarMidPtAss03to1in60100();
  void InitStandardUncertaintiesPPb2016DstarHighPtAss03to1in60100();
  void InitStandardUncertaintiesPPb2016DstarVeryHighPtAss03to1in60100();

  void InitStandardUncertaintiesPPb2016DplusLowPtAss03to1in60100();
  void InitStandardUncertaintiesPPb2016DplusMidPtAss03to1in60100();
  void InitStandardUncertaintiesPPb2016DplusHighPtAss03to1in60100();
  void InitStandardUncertaintiesPPb2016DplusVeryHighPtAss03to1in60100();

  void InitStandardUncertaintiesPPb2016DzeroLowPtAss1to99in60100();
  void InitStandardUncertaintiesPPb2016DzeroMidPtAss1to99in60100();
  void InitStandardUncertaintiesPPb2016DzeroHighPtAss1to99in60100();
  void InitStandardUncertaintiesPPb2016DzeroVeryHighPtAss1to99in60100();

  void InitStandardUncertaintiesPPb2016DstarLowPtAss1to99in60100();
  void InitStandardUncertaintiesPPb2016DstarMidPtAss1to99in60100();
  void InitStandardUncertaintiesPPb2016DstarHighPtAss1to99in60100();
  void InitStandardUncertaintiesPPb2016DstarVeryHighPtAss1to99in60100();

  void InitStandardUncertaintiesPPb2016DplusLowPtAss1to99in60100();
  void InitStandardUncertaintiesPPb2016DplusMidPtAss1to99in60100();
  void InitStandardUncertaintiesPPb2016DplusHighPtAss1to99in60100();
  void InitStandardUncertaintiesPPb2016DplusVeryHighPtAss1to99in60100();

  void InitStandardUncertaintiesPPb2016DummyValues();

  // Method with uncertainties for pp 2017
  void InitStandardUncertaintiesPP2017DzeroVeryLowPtAss03to99();
  void InitStandardUncertaintiesPP2017DzeroLowPtAss03to99();
  void InitStandardUncertaintiesPP2017DzeroMidPtAss03to99();
  void InitStandardUncertaintiesPP2017DzeroHighPtAss03to99();
  void InitStandardUncertaintiesPP2017DzeroVeryHighPtAss03to99();

  void InitStandardUncertaintiesPP2017DstarVeryLowPtAss03to99();
  void InitStandardUncertaintiesPP2017DstarLowPtAss03to99();
  void InitStandardUncertaintiesPP2017DstarMidPtAss03to99();
  void InitStandardUncertaintiesPP2017DstarHighPtAss03to99();
  void InitStandardUncertaintiesPP2017DstarVeryHighPtAss03to99();

  void InitStandardUncertaintiesPP2017DplusVeryLowPtAss03to99();
  void InitStandardUncertaintiesPP2017DplusLowPtAss03to99();
  void InitStandardUncertaintiesPP2017DplusMidPtAss03to99();
  void InitStandardUncertaintiesPP2017DplusHighPtAss03to99();
  void InitStandardUncertaintiesPP2017DplusVeryHighPtAss03to99();

  void InitStandardUncertaintiesPP2017DzeroVeryLowPtAss03to1();
  void InitStandardUncertaintiesPP2017DzeroLowPtAss03to1();
  void InitStandardUncertaintiesPP2017DzeroMidPtAss03to1();
  void InitStandardUncertaintiesPP2017DzeroHighPtAss03to1();
  void InitStandardUncertaintiesPP2017DzeroVeryHighPtAss03to1();

  void InitStandardUncertaintiesPP2017DstarVeryLowPtAss03to1();
  void InitStandardUncertaintiesPP2017DstarLowPtAss03to1();
  void InitStandardUncertaintiesPP2017DstarMidPtAss03to1();
  void InitStandardUncertaintiesPP2017DstarHighPtAss03to1();
  void InitStandardUncertaintiesPP2017DstarVeryHighPtAss03to1();

  void InitStandardUncertaintiesPP2017DplusVeryLowPtAss03to1();
  void InitStandardUncertaintiesPP2017DplusLowPtAss03to1();
  void InitStandardUncertaintiesPP2017DplusMidPtAss03to1();
  void InitStandardUncertaintiesPP2017DplusHighPtAss03to1();
  void InitStandardUncertaintiesPP2017DplusVeryHighPtAss03to1();

  void InitStandardUncertaintiesPP2017DzeroVeryLowPtAss1to99();
  void InitStandardUncertaintiesPP2017DzeroLowPtAss1to99();
  void InitStandardUncertaintiesPP2017DzeroMidPtAss1to99();
  void InitStandardUncertaintiesPP2017DzeroHighPtAss1to99();
  void InitStandardUncertaintiesPP2017DzeroVeryHighPtAss1to99();

  void InitStandardUncertaintiesPP2017DstarVeryLowPtAss1to99();
  void InitStandardUncertaintiesPP2017DstarLowPtAss1to99();
  void InitStandardUncertaintiesPP2017DstarMidPtAss1to99();
  void InitStandardUncertaintiesPP2017DstarHighPtAss1to99();
  void InitStandardUncertaintiesPP2017DstarVeryHighPtAss1to99();

  void InitStandardUncertaintiesPP2017DplusVeryLowPtAss1to99();
  void InitStandardUncertaintiesPP2017DplusLowPtAss1to99();
  void InitStandardUncertaintiesPP2017DplusMidPtAss1to99();
  void InitStandardUncertaintiesPP2017DplusHighPtAss1to99();
  void InitStandardUncertaintiesPP2017DplusVeryHighPtAss1to99();

  void InitStandardUncertaintiesPP2017DzeroVeryLowPtAss2to99();
  void InitStandardUncertaintiesPP2017DzeroLowPtAss2to99();
  void InitStandardUncertaintiesPP2017DzeroMidPtAss2to99();
  void InitStandardUncertaintiesPP2017DzeroHighPtAss2to99();
  void InitStandardUncertaintiesPP2017DzeroVeryHighPtAss2to99();

  void InitStandardUncertaintiesPP2017DstarVeryLowPtAss2to99();
  void InitStandardUncertaintiesPP2017DstarLowPtAss2to99();
  void InitStandardUncertaintiesPP2017DstarMidPtAss2to99();
  void InitStandardUncertaintiesPP2017DstarHighPtAss2to99();
  void InitStandardUncertaintiesPP2017DstarVeryHighPtAss2to99();

  void InitStandardUncertaintiesPP2017DplusVeryLowPtAss2to99();
  void InitStandardUncertaintiesPP2017DplusLowPtAss2to99();
  void InitStandardUncertaintiesPP2017DplusMidPtAss2to99();
  void InitStandardUncertaintiesPP2017DplusHighPtAss2to99();
  void InitStandardUncertaintiesPP2017DplusVeryHighPtAss2to99();

  void InitStandardUncertaintiesPP2017DzeroVeryLowPtAss3to99();
  void InitStandardUncertaintiesPP2017DzeroLowPtAss3to99();
  void InitStandardUncertaintiesPP2017DzeroMidPtAss3to99();
  void InitStandardUncertaintiesPP2017DzeroHighPtAss3to99();
  void InitStandardUncertaintiesPP2017DzeroVeryHighPtAss3to99();

  void InitStandardUncertaintiesPP2017DstarVeryLowPtAss3to99();
  void InitStandardUncertaintiesPP2017DstarLowPtAss3to99();
  void InitStandardUncertaintiesPP2017DstarMidPtAss3to99();
  void InitStandardUncertaintiesPP2017DstarHighPtAss3to99();
  void InitStandardUncertaintiesPP2017DstarVeryHighPtAss3to99();

  void InitStandardUncertaintiesPP2017DplusVeryLowPtAss3to99();
  void InitStandardUncertaintiesPP2017DplusLowPtAss3to99();
  void InitStandardUncertaintiesPP2017DplusMidPtAss3to99();
  void InitStandardUncertaintiesPP2017DplusHighPtAss3to99();
  void InitStandardUncertaintiesPP2017DplusVeryHighPtAss3to99();

  void InitStandardUncertaintiesPP2017DzeroVeryLowPtAss1to2();
  void InitStandardUncertaintiesPP2017DzeroLowPtAss1to2();
  void InitStandardUncertaintiesPP2017DzeroMidPtAss1to2();
  void InitStandardUncertaintiesPP2017DzeroHighPtAss1to2();
  void InitStandardUncertaintiesPP2017DzeroVeryHighPtAss1to2();

  void InitStandardUncertaintiesPP2017DstarVeryLowPtAss1to2();
  void InitStandardUncertaintiesPP2017DstarLowPtAss1to2();
  void InitStandardUncertaintiesPP2017DstarMidPtAss1to2();
  void InitStandardUncertaintiesPP2017DstarHighPtAss1to2();
  void InitStandardUncertaintiesPP2017DstarVeryHighPtAss1to2();

  void InitStandardUncertaintiesPP2017DplusVeryLowPtAss1to2();
  void InitStandardUncertaintiesPP2017DplusLowPtAss1to2();
  void InitStandardUncertaintiesPP2017DplusMidPtAss1to2();
  void InitStandardUncertaintiesPP2017DplusHighPtAss1to2();
  void InitStandardUncertaintiesPP2017DplusVeryHighPtAss1to2();

  void InitStandardUncertaintiesPP2017DzeroVeryLowPtAss2to3();
  void InitStandardUncertaintiesPP2017DzeroLowPtAss2to3();
  void InitStandardUncertaintiesPP2017DzeroMidPtAss2to3();
  void InitStandardUncertaintiesPP2017DzeroHighPtAss2to3();
  void InitStandardUncertaintiesPP2017DzeroVeryHighPtAss2to3();

  void InitStandardUncertaintiesPP2017DstarVeryLowPtAss2to3();
  void InitStandardUncertaintiesPP2017DstarLowPtAss2to3();
  void InitStandardUncertaintiesPP2017DstarMidPtAss2to3();
  void InitStandardUncertaintiesPP2017DstarHighPtAss2to3();
  void InitStandardUncertaintiesPP2017DstarVeryHighPtAss2to3();

  void InitStandardUncertaintiesPP2017DplusVeryLowPtAss2to3();
  void InitStandardUncertaintiesPP2017DplusLowPtAss2to3();
  void InitStandardUncertaintiesPP2017DplusMidPtAss2to3();
  void InitStandardUncertaintiesPP2017DplusHighPtAss2to3();
  void InitStandardUncertaintiesPP2017DplusVeryHighPtAss2to3();  

  // Method with uncertainties for pp 13 TeV
  void InitStandardUncertaintiesPP13TeVDzeroLowPtAss03to99();
  void InitStandardUncertaintiesPP13TeVDzeroMidPtAss03to99();
  void InitStandardUncertaintiesPP13TeVDzeroHighPtAss03to99();
  void InitStandardUncertaintiesPP13TeVDzeroVeryHighPtAss03to99();
  void InitStandardUncertaintiesPP13TeVDzeroExtremelyHighPtAss03to99();

  void InitStandardUncertaintiesPP13TeVDzeroLowPtAss03to1();
  void InitStandardUncertaintiesPP13TeVDzeroMidPtAss03to1();
  void InitStandardUncertaintiesPP13TeVDzeroHighPtAss03to1();
  void InitStandardUncertaintiesPP13TeVDzeroVeryHighPtAss03to1();
  void InitStandardUncertaintiesPP13TeVDzeroExtremelyHighPtAss03to1();

  void InitStandardUncertaintiesPP13TeVDzeroLowPtAss1to2();
  void InitStandardUncertaintiesPP13TeVDzeroMidPtAss1to2();
  void InitStandardUncertaintiesPP13TeVDzeroHighPtAss1to2();
  void InitStandardUncertaintiesPP13TeVDzeroVeryHighPtAss1to2();
  void InitStandardUncertaintiesPP13TeVDzeroExtremelyHighPtAss1to2();

  void InitStandardUncertaintiesPP13TeVDzeroLowPtAss2to3();
  void InitStandardUncertaintiesPP13TeVDzeroMidPtAss2to3();
  void InitStandardUncertaintiesPP13TeVDzeroHighPtAss2to3();
  void InitStandardUncertaintiesPP13TeVDzeroVeryHighPtAss2to3();
  void InitStandardUncertaintiesPP13TeVDzeroExtremelyHighPtAss2to3();

  void InitStandardUncertaintiesPP13TeVDzeroLowPtAss3to99();
  void InitStandardUncertaintiesPP13TeVDzeroMidPtAss3to99();
  void InitStandardUncertaintiesPP13TeVDzeroHighPtAss3to99();
  void InitStandardUncertaintiesPP13TeVDzeroVeryHighPtAss3to99();
  void InitStandardUncertaintiesPP13TeVDzeroExtremelyHighPtAss3to99();        

  // Method with uncertainties for pp 13 TeV - 0-0.1%
  void InitStandardUncertaintiesPP13TeVin001DzeroLowPtAss03to99();
  void InitStandardUncertaintiesPP13TeVin001DzeroMidPtAss03to99();
  void InitStandardUncertaintiesPP13TeVin001DzeroHighPtAss03to99();
  void InitStandardUncertaintiesPP13TeVin001DzeroVeryHighPtAss03to99();

  void InitStandardUncertaintiesPP13TeVin001DzeroLowPtAss03to1();
  void InitStandardUncertaintiesPP13TeVin001DzeroMidPtAss03to1();
  void InitStandardUncertaintiesPP13TeVin001DzeroHighPtAss03to1();
  void InitStandardUncertaintiesPP13TeVin001DzeroVeryHighPtAss03to1();

  void InitStandardUncertaintiesPP13TeVin001DzeroLowPtAss1to2();
  void InitStandardUncertaintiesPP13TeVin001DzeroMidPtAss1to2();
  void InitStandardUncertaintiesPP13TeVin001DzeroHighPtAss1to2();
  void InitStandardUncertaintiesPP13TeVin001DzeroVeryHighPtAss1to2();

  void InitStandardUncertaintiesPP13TeVin001DzeroLowPtAss2to3();
  void InitStandardUncertaintiesPP13TeVin001DzeroMidPtAss2to3();
  void InitStandardUncertaintiesPP13TeVin001DzeroHighPtAss2to3();
  void InitStandardUncertaintiesPP13TeVin001DzeroVeryHighPtAss2to3();

  void InitStandardUncertaintiesPP13TeVin001DzeroLowPtAss3to99();
  void InitStandardUncertaintiesPP13TeVin001DzeroMidPtAss3to99();
  void InitStandardUncertaintiesPP13TeVin001DzeroHighPtAss3to99();
  void InitStandardUncertaintiesPP13TeVin001DzeroVeryHighPtAss3to99();     

  // Method with uncertainties for pp 13 TeV - 0.1-10%
  void InitStandardUncertaintiesPP13TeVin0110DzeroLowPtAss03to99();
  void InitStandardUncertaintiesPP13TeVin0110DzeroMidPtAss03to99();
  void InitStandardUncertaintiesPP13TeVin0110DzeroHighPtAss03to99();
  void InitStandardUncertaintiesPP13TeVin0110DzeroVeryHighPtAss03to99();

  void InitStandardUncertaintiesPP13TeVin0110DzeroLowPtAss03to1();
  void InitStandardUncertaintiesPP13TeVin0110DzeroMidPtAss03to1();
  void InitStandardUncertaintiesPP13TeVin0110DzeroHighPtAss03to1();
  void InitStandardUncertaintiesPP13TeVin0110DzeroVeryHighPtAss03to1();

  void InitStandardUncertaintiesPP13TeVin0110DzeroLowPtAss1to2();
  void InitStandardUncertaintiesPP13TeVin0110DzeroMidPtAss1to2();
  void InitStandardUncertaintiesPP13TeVin0110DzeroHighPtAss1to2();
  void InitStandardUncertaintiesPP13TeVin0110DzeroVeryHighPtAss1to2();

  void InitStandardUncertaintiesPP13TeVin0110DzeroLowPtAss2to3();
  void InitStandardUncertaintiesPP13TeVin0110DzeroMidPtAss2to3();
  void InitStandardUncertaintiesPP13TeVin0110DzeroHighPtAss2to3();
  void InitStandardUncertaintiesPP13TeVin0110DzeroVeryHighPtAss2to3();

  void InitStandardUncertaintiesPP13TeVin0110DzeroLowPtAss3to99();
  void InitStandardUncertaintiesPP13TeVin0110DzeroMidPtAss3to99();
  void InitStandardUncertaintiesPP13TeVin0110DzeroHighPtAss3to99();
  void InitStandardUncertaintiesPP13TeVin0110DzeroVeryHighPtAss3to99();  

  // Method with uncertainties for pp 13 TeV - 10-30%
  void InitStandardUncertaintiesPP13TeVin1030DzeroLowPtAss03to99();
  void InitStandardUncertaintiesPP13TeVin1030DzeroMidPtAss03to99();
  void InitStandardUncertaintiesPP13TeVin1030DzeroHighPtAss03to99();
  void InitStandardUncertaintiesPP13TeVin1030DzeroVeryHighPtAss03to99();

  void InitStandardUncertaintiesPP13TeVin1030DzeroLowPtAss03to1();
  void InitStandardUncertaintiesPP13TeVin1030DzeroMidPtAss03to1();
  void InitStandardUncertaintiesPP13TeVin1030DzeroHighPtAss03to1();
  void InitStandardUncertaintiesPP13TeVin1030DzeroVeryHighPtAss03to1();

  void InitStandardUncertaintiesPP13TeVin1030DzeroLowPtAss1to2();
  void InitStandardUncertaintiesPP13TeVin1030DzeroMidPtAss1to2();
  void InitStandardUncertaintiesPP13TeVin1030DzeroHighPtAss1to2();
  void InitStandardUncertaintiesPP13TeVin1030DzeroVeryHighPtAss1to2();

  void InitStandardUncertaintiesPP13TeVin1030DzeroLowPtAss2to3();
  void InitStandardUncertaintiesPP13TeVin1030DzeroMidPtAss2to3();
  void InitStandardUncertaintiesPP13TeVin1030DzeroHighPtAss2to3();
  void InitStandardUncertaintiesPP13TeVin1030DzeroVeryHighPtAss2to3();

  void InitStandardUncertaintiesPP13TeVin1030DzeroLowPtAss3to99();
  void InitStandardUncertaintiesPP13TeVin1030DzeroMidPtAss3to99();
  void InitStandardUncertaintiesPP13TeVin1030DzeroHighPtAss3to99();
  void InitStandardUncertaintiesPP13TeVin1030DzeroVeryHighPtAss3to99();  

  // Method with uncertainties for pp 13 TeV - 30-100%
  void InitStandardUncertaintiesPP13TeVin30100DzeroLowPtAss03to99();
  void InitStandardUncertaintiesPP13TeVin30100DzeroMidPtAss03to99();
  void InitStandardUncertaintiesPP13TeVin30100DzeroHighPtAss03to99();
  void InitStandardUncertaintiesPP13TeVin30100DzeroVeryHighPtAss03to99();

  void InitStandardUncertaintiesPP13TeVin30100DzeroLowPtAss03to1();
  void InitStandardUncertaintiesPP13TeVin30100DzeroMidPtAss03to1();
  void InitStandardUncertaintiesPP13TeVin30100DzeroHighPtAss03to1();
  void InitStandardUncertaintiesPP13TeVin30100DzeroVeryHighPtAss03to1();

  void InitStandardUncertaintiesPP13TeVin30100DzeroLowPtAss1to2();
  void InitStandardUncertaintiesPP13TeVin30100DzeroMidPtAss1to2();
  void InitStandardUncertaintiesPP13TeVin30100DzeroHighPtAss1to2();
  void InitStandardUncertaintiesPP13TeVin30100DzeroVeryHighPtAss1to2();

  void InitStandardUncertaintiesPP13TeVin30100DzeroLowPtAss2to3();
  void InitStandardUncertaintiesPP13TeVin30100DzeroMidPtAss2to3();
  void InitStandardUncertaintiesPP13TeVin30100DzeroHighPtAss2to3();
  void InitStandardUncertaintiesPP13TeVin30100DzeroVeryHighPtAss2to3();

  void InitStandardUncertaintiesPP13TeVin30100DzeroLowPtAss3to99();
  void InitStandardUncertaintiesPP13TeVin30100DzeroMidPtAss3to99();
  void InitStandardUncertaintiesPP13TeVin30100DzeroHighPtAss3to99();
  void InitStandardUncertaintiesPP13TeVin30100DzeroVeryHighPtAss3to99();  

  /////////////

  TGraphAsymmErrors* GetUncGraphFromHistos(TH1D *hRef,TH1D *hMin,TH1D *hMax);
  void BuildGraphsRelUnc();
  void BuildGraphsUnc(TH1D *hRef);
  TCanvas* BuildSystUncertaintyPlotVsDeltaPhi(TH1D *hCorrPlot,Int_t doInit);
  void BuildTotalNonFDUncHisto();
  void BuildTotalUncHisto();
  void BuildTotalNonFlatUncHisto();  
  TH1D *GetVariedHisto(const TH1D *hIn,const TGraphAsymmErrors *gr,Int_t minmax);
  TH1D *GetHistoTotFlatMin(){return fhtotFlatMin;}
  TH1D *GetHistoTotFlatMax(){return fhtotFlatMax;}

  TH1D *GetHistoYieldUnc(){
    return fhYieldExtraction;
  }

  TH1D *GetHistoBackSubUncMin(){
    return fhBackSubtractionMin;
  }

  TH1D *GetHistoBackSubUncMax(){
    return fhBackSubtractionMax;
  }
  
  TH1D *GetHistoTemplate(){
    return fhDeltaPhiTemplate;
  }
  
  TH1D *GetHistoMCclosureTestMin(){
    return fhMCclosureTestMin;
  }
  TH1D *GetHistoMCclosureTestMax(){
    return fhMCclosureTestMax;
  }

  TH1D *GetHistoMCcorrectionsMin(){
    return fhMCcorrectionsMin;
  }
  TH1D *GetHistoMCcorrectionsMax(){
    return fhMCcorrectionsMax;
  }

  TH1D *GetHistoMCDefficiencyMin(){
    return fhMCDefficiencyMin;
  }
  TH1D *GetHistoMCDefficiencyMax(){
    return fhMCDefficiencyMax;
  }

  TH1D *GetHistoSecContaminationMin(){
    return fhSecContaminationMin;
  }
  TH1D *GetHistoSecContaminationMax(){
    return fhSecContaminationMax;
  }

  TH1D *GetHistoFDmin(){
    return fhBeautyFDmin;
  }

  TH1D *GetHistoFDmax(){
    return fhBeautyFDmax;
  }
  
  void SetHistoTemplate(TH1D *h,TString strname="",Bool_t clone=kTRUE);
  void SetHistoMCclosureTestMin(TH1D *h,TString strname="",Bool_t clone=kTRUE);
  void SetHistoMCclosureTestMax(TH1D *h,TString strname="",Bool_t clone=kTRUE);
  void SetHistoMCcorrectionsMin(TH1D *h,TString strname="",Bool_t clone=kTRUE);
  void SetHistoMCcorrectionsMax(TH1D *h,TString strname="",Bool_t clone=kTRUE);
  void SetHistoMCDefficiencyMin(TH1D *h,TString strname="",Bool_t clone=kTRUE);
  void SetHistoMCDefficiencyMax(TH1D *h,TString strname="",Bool_t clone=kTRUE);
  void SetHistoSecContaminationMin(TH1D *h,TString strname="",Bool_t clone=kTRUE);
  void SetHistoSecContaminationMax(TH1D *h,TString strname="",Bool_t clone=kTRUE);
  void SetHistoBeautyFDmin(TH1D *h,TString strname="",Bool_t clone=kTRUE);
  void SetHistoBeautyFDmax(TH1D *h,TString strname="",Bool_t clone=kTRUE);
  void SetHistoYieldExtraction(TH1D *h,TString strname="",Bool_t clone=kTRUE);
  void SetHistoBackSubtraction(TH1D *hMax,TString strname="",Bool_t clone=kTRUE,TH1D *hMin=0x0);

  
  TGraphAsymmErrors* GetTotUncGraph(){return fgrTotal;}
  TGraphAsymmErrors* GetTotNonFlatUncGraph(){return fgrTotalNonFlatDPhi;}
  TGraphAsymmErrors* GetTotFlatUncGraph(){return fgrTotalFlatDPhi;}
  
 private:
  Int_t fmeson;                       // 0=D0, 1=D*, 2=D+
  TString fstrmeson;                  // meson name
  TString fstrptAss;                  // string with pt range associated tracks
  TString fstrptD;                  // string with pt range D meson
  TH1D *fhDeltaPhiTemplate;            // histo do define the binning in DeltaPhi
  TH1D *fhYieldExtraction;            // yield extr unc
  TH1D *fhBackSubtractionMin;            // uncertainty from variation of SB range, etc.
  TH1D *fhBackSubtractionMax;            // uncertainty from variation of SB range, etc.
  TH1D *fhBeautyFDmin;                   // feed down uncertainty
  TH1D *fhBeautyFDmax;                   // feed down uncertainty
  TH1D *fhMCclosureTestMin;              // mc closure
  TH1D *fhMCclosureTestMax;              // mc closure
  TH1D *fhMCcorrectionsMin;              // mc corrections ( associated track selection variation)  
  TH1D *fhMCcorrectionsMax;              // mc corrections ( associated track selection variation)  
  TH1D *fhMCDefficiencyMin;              // mc corrections (D cut variation )
  TH1D *fhMCDefficiencyMax;              // mc corrections (D cut variation ) 
  TH1D *fhSecContaminationMin;           // contamination from secondaries
  TH1D *fhSecContaminationMax;           // contamination from secondaries
  TH1D *fhTotalMin;                      //
  TH1D *fhTotalMax;                      //
  TH1D *fhTotalNonFDMin;                //
  TH1D *fhTotalNonFDMax;                 //
  TH1D *fhTotalNonFlatDPhiMin;           //
  TH1D *fhTotalNonFlatDPhiMax;           //
  TH1D *fhtotFlatMin;                     //
  TH1D *fhtotFlatMax;                     //
  TGraphAsymmErrors *fgrYieldUnc;        //  
  TGraphAsymmErrors *fgrBackSubUnc;        //  
  TGraphAsymmErrors *fgrMCcorrections;   //
  TGraphAsymmErrors *fgrMCDefficiency;   //
  TGraphAsymmErrors *fgrSecContamination;   //
  TGraphAsymmErrors *fgrMCclosureTest;   //
  TGraphAsymmErrors *fgrBeautyFD;        //
  TGraphAsymmErrors *fgrYieldUncRel;        // 
  TGraphAsymmErrors *fgrBackSubUncRel;        //   
  TGraphAsymmErrors *fgrMCcorrectionsRel;   //
  TGraphAsymmErrors *fgrMCDefficiencyRel;   //
  TGraphAsymmErrors *fgrSecContaminationRel;   //
  TGraphAsymmErrors *fgrMCclosureTestRel;   //
  TGraphAsymmErrors *fgrBeautyFDRel;        //
  TGraphAsymmErrors  *fgrTotal;         //
  TGraphAsymmErrors  *fgrTotalRel;         //
  TGraphAsymmErrors  *fgrTotalNonFD;         //
  TGraphAsymmErrors  *fgrTotalNonFlatDPhi;         //
  TGraphAsymmErrors  *fgrTotalNonFlatDPhiRel;         //
  TGraphAsymmErrors  *fgrTotalFlatDPhi;         //
  TGraphAsymmErrors  *fgrTotalFlatDPhiRel;         //


  ClassDef(AliHFDhadronCorrSystUnc,1);
};



#endif
