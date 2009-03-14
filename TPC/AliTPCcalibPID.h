#ifndef ALITPCCALIBCPID_H
#define ALITPCCALIBCPID_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliTPCcalibBase.h"
#include "AliTPCCalPad.h"
#include "TH2F.h"
#include "TF1.h"
#include "THnSparse.h"
class TH1F;
class TList;
class AliESDEvent;
class AliESDtrack;
class AliTPCseed;

#include "TTreeStream.h"


class AliTPCcalibPID:public AliTPCcalibBase {
public:
  AliTPCcalibPID(); 
  AliTPCcalibPID(const Text_t *name, const Text_t *title);
  virtual ~AliTPCcalibPID();
  
  virtual void           Process(AliESDEvent *event);
  virtual Long64_t       Merge(TCollection *li);
  virtual void           Analyze();
  void                   MakeReport();
  //
  //
  TH1F   *          GetHistNTracks(){return fHistNTracks;};
  TH1F   *          GetHistClusters(){return fClusters;};
  TH2F   *          GetHistPileUp(){return fPileUp;};
  TH2F   *          GetHistLandau(){return fLandau;};
  //
  THnSparseS *      GetHistQmax(){return fDeDxQmax;};
  THnSparseS *      GetHistQtot(){return fDeDxQtot;};
  THnSparseS *      GetHistRatio(){return fDeDxRatio;};
  THnSparseS *      GetHistShortMediumRatio(){return fDeDxShortMediumRatio;};
  THnSparseS *      GetHistLongMediumRatio(){return fDeDxLongMediumRatio;};
  //
  void SetMIPvalue(Float_t mip){fMIP = mip;};
  void SetLowerTrunc(Float_t lowerTrunc){fLowerTrunc = lowerTrunc;};
  void SetUpperTrunc(Float_t upperTrunc){fUpperTrunc = upperTrunc;};
  void SetUseShapeNorm(Bool_t useShapeNorm){fUseShapeNorm = useShapeNorm;};
  void SetUsePosNorm(Bool_t usePosNorm){fUsePosNorm = usePosNorm;};
  void SetPadNorm(Int_t padNorm){fUsePadNorm = padNorm;};
  void SetIsCosmic(Bool_t isCosmic){fIsCosmic = isCosmic;};
  //
  //
  static void       BinLogX(THnSparse * h, Int_t axisDim);   // method for correct histogram binning


private:
  //
  // parameter specifications
  //
  Float_t fMIP;
  Float_t fLowerTrunc;
  Float_t fUpperTrunc;
  Bool_t  fUseShapeNorm;
  Bool_t  fUsePosNorm;
  Int_t   fUsePadNorm;
  //
  Bool_t  fIsCosmic;
  //
  // histograms
  //
  TH1F  *fHistNTracks;            //  histogram showing number of ESD tracks per event
  TH1F  *fClusters;               //  histogram showing the number of clusters per track
  TH2F  *fPileUp;                 //  histogram which shows correlation between time mismatch and dEdx signal
  TH2F  *fLandau;                 //  histogran which shows Landau distribution for the three pad geometries
  //
  THnSparseS * fDeDxQmax;               //  histogram which shows dEdx (Qmax) as a function of z,sin(phi),tan(theta),p,betaGamma
  THnSparseS * fDeDxQtot;               //  histogram which shows dEdx (Qtot) as a function of z,sin(phi),tan(theta),p,betaGamma
  THnSparseS * fDeDxRatio;              //  histogram which shows dEdx ratio (Qmax/Qtot) as a function of z,sin(phi),tan(theta),p,betaGamma
  THnSparseS * fDeDxShortMediumRatio;   //  histogram which shows dEdx ratio (QmaxShort/QmaxMedium) as a function of z,sin(phi),tan(theta),p,betaGamma
  THnSparseS * fDeDxLongMediumRatio;    //  histogram which shows dEdx ratio (QmaxLong/QmaxMedium) as a function of z,sin(phi),tan(theta),p,betaGamma
  //
  AliTPCcalibPID(const AliTPCcalibPID&); 
  AliTPCcalibPID& operator=(const AliTPCcalibPID&); 

  ClassDef(AliTPCcalibPID, 1); 
};

#endif


