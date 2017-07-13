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
  void InitStandardUncertaintiesPP2010(Int_t meson,Double_t ptD,Double_t minptAss);
  void InitStandardUncertaintiesPPb2013(Int_t meson,Double_t ptD,Double_t minptAss);
  void InitStandardUncertaintiesPPb2016(Int_t meson,Double_t ptD,Double_t minptAss,Double_t maxptAss);
  // Method with uncertainties for pp 2010, Dzero and D*+ and pt assoc> 0.3 GeV/c, with values used for HP2013
  void InitStandardUncertaintiesPP2010DzeroLowPtAss03HP();
  void InitStandardUncertaintiesPP2010DzeroMidPtAss03HP();
  void InitStandardUncertaintiesPP2010DzeroHighPtAss03HP();
  
  void InitStandardUncertaintiesPP2010DstarLowPtAss03HP();
  void InitStandardUncertaintiesPP2010DstarMidPtAss03HP();
  void InitStandardUncertaintiesPP2010DstarHighPtAss03HP();


  // Method with uncertainties for pp 2010, all kinematic cases but those approved for HP2013
  void InitStandardUncertaintiesPP2010DplusLowPtAss03();
  void InitStandardUncertaintiesPP2010DplusMidPtAss03();
  void InitStandardUncertaintiesPP2010DplusHighPtAss03();

  void InitStandardUncertaintiesPP2010DzeroLowPtAss05();
  void InitStandardUncertaintiesPP2010DzeroMidPtAss05();
  void InitStandardUncertaintiesPP2010DzeroHighPtAss05();

  void InitStandardUncertaintiesPP2010DstarLowPtAss05();
  void InitStandardUncertaintiesPP2010DstarMidPtAss05();
  void InitStandardUncertaintiesPP2010DstarHighPtAss05();

  void InitStandardUncertaintiesPP2010DplusLowPtAss05();
  void InitStandardUncertaintiesPP2010DplusMidPtAss05();
  void InitStandardUncertaintiesPP2010DplusHighPtAss05();


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

  void InitStandardUncertaintiesPPb2013DzeroLowPtAss05();
  void InitStandardUncertaintiesPPb2013DzeroMidPtAss05();
  void InitStandardUncertaintiesPPb2013DzeroHighPtAss05();

  void InitStandardUncertaintiesPPb2013DstarLowPtAss05();
  void InitStandardUncertaintiesPPb2013DstarMidPtAss05();
  void InitStandardUncertaintiesPPb2013DstarHighPtAss05();

  void InitStandardUncertaintiesPPb2013DplusLowPtAss05();
  void InitStandardUncertaintiesPPb2013DplusMidPtAss05();
  void InitStandardUncertaintiesPPb2013DplusHighPtAss05();


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
