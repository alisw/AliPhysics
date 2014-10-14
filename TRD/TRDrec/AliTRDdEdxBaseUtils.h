/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
//
//
//
//  Xianguo Lu 
//  lu@physi.uni-heidelberg.de
//  Xianguo.Lu@cern.ch
//
/*
grep " AliTRDdEdxBaseUtils::" AliTRDdEdxBaseUtils.cxx | grep "=" -v  | grep -v "[6]" | grep -v printf  |wc
grep "(" AliTRDdEdxBaseUtils.h | grep ";" | grep -v grep | grep -v ClassDef | grep -v "{" | grep -v typedef | wc
*/


#ifndef ALITRDDEDXBASEUTILS_H
#define ALITRDDEDXBASEUTILS_H

#ifndef TVECTORD_H
#include "TVectorD.h"
#endif

#ifndef THNSPARSE_H
#include "THnBase.h"
#endif

#ifndef TTREESTREAM_H
#include "TTreeStream.h"
#endif 

class TH1D;
class TH2D;
class TObjArray;

class AliESDEvent;
class AliESDtrack;
class AliTRDcluster;
class AliTRDtrackV1;
class AliTRDseedV1;

class AliTRDdEdxBaseUtils
{
 public:
  //===================================================================================
  //                                   Math and Histogram
  //===================================================================================
  static void BinLogX(TAxis *axis);
  static void GetCDFCuts(const TH1D *hh, Int_t ncut, Double_t cuts[], const Double_t cdfs[], Double_t thres);
  static Double_t GetMeanRMS(Double_t nn, Double_t sum, Double_t w2s, Double_t * grms=0x0, Double_t * gerr=0x0);
  static Double_t TruncatedMean(Int_t nx, const Double_t xdata[], Double_t lowfrac, Double_t highfrac, Double_t * grms=0x0, Double_t * gerr=0x0, Double_t *wws=0x0);
  static Double_t TruncatedMean(const TH1 *hh, Double_t lowfrac, Double_t highfrac, Double_t * grms=0x0, Double_t * gerr=0x0);
  static void FitSlicesY(const TH2D *hh, TH1D *&hnor, TH1D *&hmpv, TH1D *&hwid, TH1D *&hres, Double_t thres, Double_t lowfrac, Double_t highfrac);

  //===================================================================================
  //                                TRD Analysis Fast Tool
  //===================================================================================
  static Int_t GetNtracklet(const AliESDEvent *esd);
  static AliTRDtrackV1 * GetTRDtrackV1(const AliESDtrack * esdtrack);
  static Bool_t IsInSameStack(const AliTRDtrackV1 *trdtrack);
  static AliTRDseedV1 * GetFirstTracklet(const AliTRDtrackV1 *trdtrack);
  static AliTRDseedV1 * GetLastTracklet(const AliTRDtrackV1 *trdtrack);
  static void GetFirstSectorStackMomentum(const AliTRDtrackV1 *trdtrack, Int_t & isec, Int_t & istk, Double_t & mom);
  static Double_t GetRedefinedPhi(Double_t phi);
  static Double_t Getdydx(const AliTRDseedV1 *tracklet);
  static Double_t Getdzdx(const AliTRDseedV1 *tracklet);
  static Double_t Getdldx(const AliTRDseedV1 *tracklet);

  //===================================================================================
  //                                 Detector, Data and Control Constant
  //===================================================================================
  
  static Int_t NTRDchamber(){return 18*5*6;} //540
  static Int_t NTRDtimebin(){return NTRDchamber()*31;} //16740
  static Int_t ToDetector(Int_t gtb);
  static Int_t ToTimeBin(Int_t gtb);
  static Int_t ToSector(Int_t gtb);
  static Int_t ToStack(Int_t gtb);
  static Int_t ToLayer(Int_t gtb);

  static void CheckRunB(TString listrun1kg, Int_t run, TString & type);
  static TString GetRunType(Int_t run);

  static void SetQ0Frac(Double_t q0){ fgQ0Frac = q0; }
  static void SetQ1Frac(Double_t q1){ fgQ1Frac = q1; }
  static void SetTimeBinCountCut(Double_t tbc){ fgTimeBinCountCut = tbc; }
  static void SetCalibTPCnclsCut(Int_t tpc){ fgCalibTPCnclsCut = tpc; }
  static void SetExBOn(Bool_t kon){ fgExBOn = kon; }
  static void SetPadGainOn(Bool_t kon){ fgPadGainOn = kon; }
  static void SetQScale(Double_t scale){ fgQScale = scale; }
 
  static Double_t Q0Frac(){return fgQ0Frac;}
  static Double_t Q1Frac(){return fgQ1Frac;}
  static Double_t TimeBinCountCut(){return fgTimeBinCountCut;}
  static Int_t CalibTPCnclsCut(){return fgCalibTPCnclsCut;}
  static Bool_t IsExBOn(){return fgExBOn;}
  static Bool_t IsPadGainOn(){return fgPadGainOn;}
  static Double_t QScale(){return fgQScale;}
  
  static void PrintControl();
  
  //===================================================================================
  //                                 dEdx Parameterization
  //===================================================================================
  static void FastFitdEdxTR(TH1 * hh);

  static Double_t ALEPH(const Double_t * xx,  const Double_t * par);
  static Double_t MeandEdx(const Double_t * xx, const Double_t * par);
  static Double_t MeanTR(const Double_t * xx, const Double_t * par);
  static Double_t MeandEdxTR(const Double_t * xx, const Double_t * par);

  static Double_t QMeanTPC(Double_t bg);
  static Double_t Q0MeanTRDpp(Double_t bg);
  static Double_t Q1MeanTRDpp(Double_t bg);
  static Double_t Q0MeanTRDPbPb(Double_t bg);
  static Double_t Q1MeanTRDPbPb(Double_t bg);

  typedef Double_t (*FFunc)(const Double_t *xx, const Double_t *par);
  
  static Double_t MeandEdxLogx(const Double_t * xx, const Double_t * par){return ToLogx(MeandEdx, xx, par);}
  static Double_t MeanTRLogx(const Double_t * xx, const Double_t * par){return ToLogx(MeanTR, xx, par);}
  static Double_t MeandEdxTRLogx(const Double_t * xx, const Double_t * par){return ToLogx(MeandEdxTR, xx, par);}

 private:
  //dEdx Parameterization
  static Double_t ToLogx(FFunc func, const Double_t * xx, const Double_t * par);

  //Control Constant
  static Double_t fgQ0Frac; //q0frac
  static Double_t fgQ1Frac; //q1frac
  static Double_t fgTimeBinCountCut; //tbcut
  static Int_t    fgCalibTPCnclsCut; //tpccut
  static Bool_t   fgExBOn;    //exbon
  static Bool_t fgPadGainOn; //pad gain
  static Double_t fgQScale; //Qscale

};

#endif  
