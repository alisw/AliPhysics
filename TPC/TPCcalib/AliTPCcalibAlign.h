#ifndef ALITPCCALIBALIGN_H
#define ALITPCCALIBALIGN_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////
////
////

class TFile;
class TGraphErrors;
class TH1;
class THnSparse;
class THnBase;
class THn;
#include <TLinearFitter.h>
#include <TMatrixDfwd.h>
class TObjArray;
class TTree;

#include "AliTPCcalibBase.h"
class AliExternalTrackParam;
class AliTPCPointCorrection;
class AliTPCseed;
class AliVEvent;

class AliTPCcalibAlign:public AliTPCcalibBase {
public:
  enum HistoType {kY=0, kZ =1, kPhi=2, kTheta=3, 
		  kYPhi=4, kZTheta=5, 
		  kYz=6,kZz=7,kPhiZ=8,kThetaZ=9};
  enum FitType{ k6=0, k9=1, k12=2};
  AliTPCcalibAlign();
  AliTPCcalibAlign(const Text_t *name, const Text_t *title);
  AliTPCcalibAlign(const AliTPCcalibAlign &align);
  //
  virtual ~AliTPCcalibAlign();
  void     Process(AliVEvent *event);
  virtual void ProcessSeed(AliTPCseed *track);
  virtual void Process(AliTPCseed */*track*/){ return ;}
  virtual void Analyze();
  virtual void Terminate();  
  virtual Long64_t Merge(TCollection* const list);
  void ExportTrackPoints(AliVEvent *event);
  //
  //
  void MakeReportDy(TFile *output); 
  void MakeReportDyPhi(TFile *output);
  //
  void UpdatePointCorrection(AliTPCPointCorrection * correction);
  //
  virtual void EvalFitters(Int_t minPoints=20);
  TH1 * GetHisto(HistoType type, Int_t s1, Int_t s2, Bool_t force=kFALSE);
  void  MakeTree(const char *fname="alignTree.root", Int_t minPoints=20);
  TGraphErrors * MakeGraph(Int_t sec0, Int_t sec1, Int_t dsec, 
			   Int_t i0, Int_t i1, FitType type); 
  Int_t  RefitLinear(const AliTPCseed * seed, Int_t isec, Double_t *fitParam, Int_t refSector, TMatrixD &param, TMatrixD&covar, Double_t xRef, Bool_t both=kFALSE);
  
  void ProcessTracklets(const AliExternalTrackParam &t1,
			const AliExternalTrackParam &t2,
			const AliTPCseed * seed,
			Int_t s1,Int_t s2);
  
  void UpdateClusterDeltaField(const AliTPCseed * seed);
  void UpdateAlignSector(const AliTPCseed * seed,Int_t isec); 
  Int_t GetIndex(Int_t s1,Int_t s2) const {return 72*s1+s2;}
  //
  inline const TMatrixD     * GetTransformation(Int_t s1,Int_t s2, Int_t fitType);
  //
  inline TLinearFitter* GetFitter12(Int_t s1,Int_t s2);
  inline TLinearFitter* GetFitter9(Int_t s1,Int_t s2);
  inline TLinearFitter* GetFitter6(Int_t s1,Int_t s2);
  //
  Bool_t GetTransformation12(Int_t s1,Int_t s2,TMatrixD &a);
  Bool_t GetTransformation9(Int_t s1,Int_t s2,TMatrixD &a);
  Bool_t GetTransformation6(Int_t s1,Int_t s2,TMatrixD &a);
  Int_t  AcceptTracklet(const AliExternalTrackParam &tp1,
			const AliExternalTrackParam &tp2) const;
  Int_t  AcceptTracklet(const Double_t *t1,
			const Double_t *t2) const;

  void ProcessDiff(const AliExternalTrackParam &t1,
		   const AliExternalTrackParam &t2,
		   const AliTPCseed *seed,
		   Int_t s1,Int_t s2);
  void ProcessAlign(Double_t * t1, Double_t * t2, Int_t s1,Int_t s2);

//   Bool_t GetTransformationCovar12(Int_t s1,Int_t s2,TMatrixD &a, Bool_t norm=kFALSE);
//   Bool_t GetTransformationCovar9(Int_t s1,Int_t s2,TMatrixD &a, Bool_t norm=kFALSE);
//   Bool_t GetTransformationCovar6(Int_t s1,Int_t s2,TMatrixD &a, Bool_t norm=kFALSE);
  void Add(AliTPCcalibAlign * align);
  const Int_t *GetPoints() const {return fPoints;}
  //void     Process(AliESDtrack *const track, Int_t runNo=-1){AliTPCcalibBase::Process(track,runNo);}
  void     Process(AliVTrack *const track, Int_t runNo=-1){AliTPCcalibBase::Process(track,runNo);}
  TLinearFitter* GetOrMakeFitter12(Int_t s1,Int_t s2);
  TLinearFitter* GetOrMakeFitter9(Int_t s1,Int_t s2);
  TLinearFitter* GetOrMakeFitter6(Int_t s1,Int_t s2);
  void Process12(const Double_t *t1, const Double_t *t2,
		 TLinearFitter *fitter) const;
  void Process9(const Double_t *const t1, const Double_t *const t2, TLinearFitter *fitter) const;
  void Process6(const Double_t *const t1, const Double_t *const t2, TLinearFitter *fitter) const;
  void GlobalAlign6(Int_t minPoints, Float_t sysError, Int_t niter);
  //
  // Cluster comparison Part
  //
  //
  // For visualization and test purposes
  //
  Double_t Correct(Int_t type, Int_t value, Int_t s1, Int_t s2, Double_t x, Double_t y, Double_t z, Double_t phi,Double_t theta); 
  static Double_t SCorrect(Int_t type, Int_t value, Int_t s1, Int_t s2, Double_t x, Double_t y, Double_t z, Double_t phi,Double_t theta){return Instance()->Correct(type,value,s1,s2,x,y,z,phi,theta);}
  static AliTPCcalibAlign* Instance();
  void SetInstance(AliTPCcalibAlign* const param){fgInstance = param;}
  static void Constrain1Pt(AliExternalTrackParam &t1, const AliExternalTrackParam &t2, Bool_t noField);
  void SetNoField(Bool_t noField){ fNoField=noField;}

  //
  // Kalman fileter for sectors
  //
  void MakeSectorKalman();
  void UpdateSectorKalman(Int_t sector, Int_t quadrant0, Int_t quadrant1,  TMatrixD *const p0, TMatrixD *const c0, TMatrixD *const p1, TMatrixD *const c1);
  void UpdateSectorKalman(TMatrixD &par0, TMatrixD &cov0, TMatrixD &para1, TMatrixD &cov1);
  Double_t GetCorrectionSector(Int_t coord, Int_t sector, Double_t lx, Double_t ly, Double_t lz); 
  static Double_t SGetCorrectionSector(Int_t coord, Int_t sector, Double_t lx, Double_t ly, Double_t lz); 

  //
  // Kalman filter for full TPC
  //
  void MakeKalman();
  void UpdateKalman(Int_t sector0, Int_t sector1,  TMatrixD &p0, TMatrixD &c0, TMatrixD &p1, TMatrixD &c1);
  void UpdateKalman(TMatrixD &par0, TMatrixD &cov0, TMatrixD &para1, TMatrixD &cov1);
  //
  //private:
  static Int_t CheckCovariance(TMatrixD &covar);
  //
  //
  void MakeResidualHistos();
  void MakeResidualHistosTracklet();
  THn * GetClusterDelta(Int_t index) const  { return fClusterDelta[index];}
  THnSparse * GetTrackletDelta(Int_t index) const  { return fTrackletDelta[index];}
public:
  
  void FillHisto(const Double_t *t1,
		 const Double_t *t2,
		 Int_t s1,Int_t s2);
  void FillHisto(AliExternalTrackParam *tp1,
		 AliExternalTrackParam *tp2,
		 Int_t s1,Int_t s2);

  static void SetMergeEntriesCut(Double_t entriesCut){fgkMergeEntriesCut = entriesCut;}
protected:
  THn     *fClusterDelta[2];  //clusters residuals
  THnSparse     *fTrackletDelta[4]; //track residuals

  TObjArray fDphiHistArray;    // array of residual histograms  phi      -kPhi
  TObjArray fDthetaHistArray;  // array of residual histograms  theta    -kTheta
  TObjArray fDyHistArray;      // array of residual histograms  y        -kY
  TObjArray fDzHistArray;      // array of residual histograms  z        -kZ
  //
  TObjArray fDyPhiHistArray;      // array of residual histograms  y     -kYPhi
  TObjArray fDzThetaHistArray;    // array of residual histograms  z-z   -kZTheta
  //
  TObjArray fDphiZHistArray;      // array of residual histograms  phi   -kPhiz
  TObjArray fDthetaZHistArray;    // array of residual histograms  theta -kThetaz
  TObjArray fDyZHistArray;        // array of residual histograms  y     -kYz
  TObjArray fDzZHistArray;        // array of residual histograms  z     -kZz
  //
  //
  TObjArray fFitterArray12;    // array of fitters
  TObjArray fFitterArray9;     // array of fitters
  TObjArray fFitterArray6;     // array of fitters
  //
  TObjArray fMatrixArray12;    // array of transnformtation matrix
  TObjArray fMatrixArray9;     // array of transnformtation matrix
  TObjArray fMatrixArray6;     // array of transnformtation matrix 
  //
  //
  //
  //
  TObjArray fCombinedMatrixArray6;      // array  combeined transformation matrix
  //
  //
  Int_t fPoints[72*72];        // number of points in the fitter 
  Bool_t fNoField;             // flag - no field data
  // refernce x
  Double_t fXIO;               // OROC-IROC boundary
  Double_t fXmiddle;           // center of the TPC sector local X
  Double_t fXquadrant;         // x quadrant
  //
  // Kalman filter for sector internal  alignemnt
  //
  TObjArray fArraySectorIntParam; // array of sector alignment parameters
  TObjArray fArraySectorIntCovar; // array of sector alignment covariances 
  //
  // Kalman filter for global alignment
  //
  TMatrixD  *fSectorParamA;     // Kalman parameter   for A side
  TMatrixD  *fSectorCovarA;     // Kalman covariance  for A side 
  TMatrixD  *fSectorParamC;     // Kalman parameter   for A side
  TMatrixD  *fSectorCovarC;     // Kalman covariance  for A side 
  //
  //
  //
  Bool_t    fUseInnerOuter;         // flag- use Inner Outer sector for left righ alignment
  
  static AliTPCcalibAlign*   fgInstance; //! Instance of this class (singleton implementation)
  static Double_t            fgkMergeEntriesCut;  //maximal number of entries for merging  -can be modified via setter
private:
  AliTPCcalibAlign&  operator=(const AliTPCcalibAlign&);// not implemented

  // IMPORTANT: If you change the data members, 
  // please do not forget to increment the ClassDef and to update the Streamer in AliTPCcalibAlign.cxx
  ClassDef(AliTPCcalibAlign,7)
};


TLinearFitter* AliTPCcalibAlign::GetFitter12(Int_t s1,Int_t s2) {
  return static_cast<TLinearFitter*>(fFitterArray12[GetIndex(s1,s2)]);
}
TLinearFitter* AliTPCcalibAlign::GetFitter9(Int_t s1,Int_t s2) {
  return static_cast<TLinearFitter*>(fFitterArray9[GetIndex(s1,s2)]);
}
TLinearFitter* AliTPCcalibAlign::GetFitter6(Int_t s1,Int_t s2) {
  return static_cast<TLinearFitter*>(fFitterArray6[GetIndex(s1,s2)]);
}

const TMatrixD * AliTPCcalibAlign::GetTransformation(Int_t s1,Int_t s2, Int_t fitType){
  if (fitType==0) return static_cast<TMatrixD*>(fMatrixArray6[GetIndex(s1,s2)]);
  if (fitType==1) return static_cast<TMatrixD*>(fMatrixArray9[GetIndex(s1,s2)]);
  if (fitType==2) return static_cast<TMatrixD*>(fMatrixArray12[GetIndex(s1,s2)]);
  return 0;
}



#endif
