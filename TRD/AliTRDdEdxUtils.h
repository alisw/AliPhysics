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
grep " AliTRDdEdxUtils::" AliTRDdEdxUtils.cxx | grep "=" -v  | grep -v "[6]" | grep -v printf  |wc
grep "(" AliTRDdEdxUtils.h | grep ";" | grep -v grep | grep -v ClassDef | grep -v "{" | grep -v typedef | wc
*/


#ifndef ALITRDDEDXUTILS_H
#define ALITRDDEDXUTILS_H

#ifndef TVECTORD_H
#include "TVectorD.h"
#endif

#ifndef THNSPARSE_H
#include "THn.h"
#include "THnSparse.h"
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

class AliTRDdEdxUtils
{
 public:

  //===================================================================================
  //                                   Math and Histogram
  //===================================================================================
  static void GetCDFCuts(const TH1D *hh, const Int_t ncut, Double_t cuts[], const Double_t cdfs[], const Double_t thres);
  static Double_t GetMeanRMS(const Double_t nn, const Double_t sum, const Double_t w2s, Double_t * grms=0x0, Double_t * gerr=0x0);
  static Double_t TruncatedMean(const Int_t nx, const Double_t xdata[], const Double_t lowfrac, const Double_t highfrac, Double_t * grms=0x0, Double_t * gerr=0x0, Double_t *wws=0x0);
  static Double_t TruncatedMean(const TH1 *hh, const Double_t lowfrac, const Double_t highfrac, Double_t * grms=0x0, Double_t * gerr=0x0);
  static void FitSlicesY(const TH2D *hh, TH1D *&hnor, TH1D *&hmpv, TH1D *&hwid, TH1D *&hres, const Double_t thres, const Double_t lowfrac, const Double_t highfrac);

  //===================================================================================
  //                                TRD Analysis Fast Tool
  //===================================================================================
  static Int_t GetNtracklet(const AliESDEvent *esd);
  static AliTRDtrackV1 * GetTRDtrackV1(const AliESDtrack * esdtrack);
  static Bool_t IsInSameStack(const AliTRDtrackV1 *trdtrack);
  static Bool_t GetFirstSectorStackMomentum(const AliTRDtrackV1 *trdtrack, Int_t & isec, Int_t & istk, Double_t & mom);

  //===================================================================================
  //                                Calibration
  //===================================================================================
  static void SetCalibFile(const TString file){fgCalibFile = file;}
  static void DeleteCalibObj();
  static void IniCalibObj();

  static void SetObjPHQ(TObjArray * obj){fgObjPHQ = obj;}
  static Bool_t GenerateDefaultPHQOCDB(const TString path="local://./");

  static Double_t GetCalibTPCscale(const Int_t tpcncls, const Double_t tpcsig);

  static void DeleteCalibHist();
  static void IniCalibHist(TList *list, const Bool_t kPHQonly=kFALSE);
  static Bool_t ReadCalibHist(const TString filename, const TString listname);

  static TObjArray * GetObjPHQ();
  static TObjArray * GetHistPHQ(){return fgHistPHQ;}
  static TObjArray * GetObjPHQ(const Bool_t kinvq, const Double_t mag, const Int_t charge);
  static THnF * GetHistPHQ(const Bool_t kinvq, const Double_t mag, const Int_t charge);

  static void FillCalibHist(const Int_t ncls, const TVectorD *arrayQ, const TVectorD *arrayX, THnF * hcalib, const Double_t scale);
  static void FillCalibHist(const AliTRDtrackV1 *trdv1, const Bool_t kinvq, const Double_t mag, const Int_t charge, const Double_t scale);

  static void CalibOutput(const TList *l0, Int_t run);
  static TObjArray* GetCalibObj(const THnF *hh, Int_t run=-999, TList *lout=0x0, TTreeSRedirector *calibStream=0x0);

  static THnF * GetHistGain(){return fgHistGain;}
  static THnF * GetHistT0(){return fgHistT0;}
  static THnF * GetHistVd(){return fgHistVd;}

  //===================================================================================
  //                                   dEdx calculation
  //===================================================================================
  static Double_t ToyCook(const Bool_t kinvq, Int_t &ncluster, TVectorD *arrayQ, TVectorD *arrayX, const TObjArray *cobj=0x0);
  static Double_t CombineddEdx(const Bool_t kinvq, Int_t &concls, TVectorD *coarrayQ, TVectorD *coarrayX, const Int_t tpcncls, const TVectorD *tpcarrayQ, const TVectorD *tpcarrayX, const Int_t trdncls, const TVectorD *trdarrayQ, const TVectorD *trdarrayX);

  //===================================================================================
  //                                   dEdx Getter and Setter
  //===================================================================================
  static Int_t GetArrayClusterQ(const Bool_t kinvq, TVectorD *arrayQ, TVectorD *arrayX, const AliTRDtrackV1 *trdtrack, Int_t timeBin0=-1, Int_t timeBin1=1000, Int_t tstep=1);
  static Int_t UpdateArrayX(const Int_t ncls, TVectorD* arrayX);
  static void SetChamberQT(const AliTRDtrackV1 *trdtrack, const Int_t kcalib, THnF * hgain=0x0, THnF * ht0=0x0, THnF * hvd=0x0);

  static Int_t GetNchamber() {return fgNchamber;}
  static Double_t GetChamberQ(const Int_t ich)  {return fgChamberQ[ich];}
  static Double_t GetChamberTmean(const Int_t ich) {return fgChamberTmean[ich];}
  static Double_t GetTrackTmean() {return fgTrackTmean;}
  
  //===================================================================================
  //                                 dEdx Parameterization
  //===================================================================================
  
  static Double_t MeandEdx(const Double_t * xx, const Double_t * par);
  static Double_t MeanTR(const Double_t * xx, const Double_t * par);
  static Double_t MeandEdxTR(const Double_t * xx, const Double_t * par);

  static Double_t QMeanTPC(const Double_t bg);
  static Double_t Q0MeanTRDpp(const Double_t bg);
  static Double_t Q1MeanTRDpp(const Double_t bg);
  static Double_t Q0MeanTRDPbPb(const Double_t bg);
  static Double_t Q1MeanTRDPbPb(const Double_t bg);

  typedef Double_t (*FFunc)(const Double_t *xx, const Double_t *par);
  
  static Double_t MeandEdxLogx(const Double_t * xx, const Double_t * par){return ToLogx(MeandEdx, xx, par);}
  static Double_t MeanTRLogx(const Double_t * xx, const Double_t * par){return ToLogx(MeanTR, xx, par);}
  static Double_t MeandEdxTRLogx(const Double_t * xx, const Double_t * par){return ToLogx(MeandEdxTR, xx, par);}

  //===================================================================================
  //                                 Detector, Data and Control Constant
  //===================================================================================
  static Int_t NTRDchamber(){return 18*5*6;} //540
  static Int_t NTRDtimebin(){return NTRDchamber()*31;} //16740
  static Int_t ToDetector(const Int_t gtb);
  static Int_t ToTimeBin(const Int_t gtb);
  static Int_t ToSector(const Int_t gtb);
  static Int_t ToStack(const Int_t gtb);
  static Int_t ToLayer(const Int_t gtb);

  static TString GetRunType(const Int_t run);

  static void SetPadGainOn(const Bool_t kon){ fgPadGainOn = kon; }
  static void SetExBOn(const Bool_t kon){ fgExBOn = kon; }
  static void SetQScale(const Double_t scale){ fgQScale = scale; }
  static void SetQ0Frac(const Double_t q0){ fgQ0Frac = q0; }
  static void SetQ1Frac(const Double_t q1){ fgQ1Frac = q1; }
  static void SetTimeBinCountCut(const Double_t tbc){ fgTimeBinCountCut = tbc; }
  static void SetCalibTPCnclsCut(const Int_t tpc){ fgCalibTPCnclsCut = tpc; }

  static Bool_t IsPadGainOn(){return fgPadGainOn;}
  static Bool_t IsExBOn(){return fgExBOn;}
  static Double_t QScale(){return fgQScale;}
  static Double_t Q0Frac(){return fgQ0Frac;}
  static Double_t Q1Frac(){return fgQ1Frac;}
  static Double_t TimeBinCountCut(){return fgTimeBinCountCut;}
  static Int_t CalibTPCnclsCut(){return fgCalibTPCnclsCut;}

  static void PrintControl();
  //===================================================================================
  //===================================================================================

  private:

  //dEdx Getter and Setter
  static Double_t GetAngularCorrection(const AliTRDseedV1 *seed);
  static Double_t GetPadGain(const Int_t det, const Int_t icol, const Int_t irow);
  static Double_t GetRNDClusterQ(AliTRDcluster *cl);
  static Double_t GetClusterQ(const Bool_t kinvq, const AliTRDseedV1 * seed, const Int_t itb);
  
  //dEdx Parameterization
  static Double_t ToLogx(FFunc func, const Double_t * xx, const Double_t * par);

  //Calibration
  static Int_t ApplyCalib(const Int_t nc0, TVectorD *arrayQ, TVectorD *arrayX, const TObjArray *cobj);
  static void GetPHCountMeanRMS(const TH1D *hnor, TH1D *&hmean);
  static Int_t GetPHQIterator(const Bool_t kinvq, const Double_t mag, const Int_t charge);
  static TString GetPHQName(const Bool_t kobj, const Int_t iter);

  //Detector, Data and Control Constant
  static THnF *fgHistGain;//PH hist
  static THnF *fgHistT0;//PH hist
  static THnF *fgHistVd;//PH hist
  static TObjArray * fgHistPHQ;//array containing 8 THnF!

  static TString fgCalibFile; //private path for calibration object

  static TObjArray * fgObjGain;//chamber gain obj
  static TObjArray * fgObjT0;//t0 obj
  static TObjArray * fgObjVd;//Vd obj
  static TObjArray * fgObjPHQ;//array containing 8 TObjArray!

  static Int_t fgNchamber; //number of chamber used in cookdedx
  static Double_t fgChamberQ[6];  //dqdl in chamber [i]
  static Double_t fgChamberTmean[6]; //Q-weighted timebin \sum Q*T / \sum Q

  static Double_t fgTrackTmean; //mean timebin over track

  static Bool_t fgPadGainOn; //pad gain
  static Bool_t   fgExBOn;    //exbon
  static Double_t fgQScale; //Qscale
  static Double_t fgQ0Frac; //q0frac
  static Double_t fgQ1Frac; //q1frac
  static Double_t fgTimeBinCountCut; //tbcut
  static Int_t    fgCalibTPCnclsCut; //tpccut
};

#endif
