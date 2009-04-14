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
  void                   MakeReport(const char * outputpath);
  void                   DrawRatioTot(Int_t ipad, const char* outputpath);
  void                   DrawRatioMax(Int_t ipad, const char* outputpath);
  void                   DrawRatiodEdx(Float_t demin, Float_t demax, const char* outputpath);
  void                   DrawResolBGQtot(Int_t minClusters, Int_t maxClusters, Float_t minp, Float_t maxp,  const char *outputpath, Bool_t resol=kTRUE); //
  void                   DrawResolBGQmax(Int_t minClusters, Int_t maxClusters, Float_t minp, Float_t maxp, const char *outputpath, Bool_t resol=kTRUE);
  //
  TH1F   *          GetHistNTracks(){return fHistNTracks;};
  TH1F   *          GetHistClusters(){return fClusters;};
  TH2F   *          GetHistPileUp(){return fPileUp;};
  TH2F   *          GetHistLandau(){return fLandau;};
  //
  THnSparseS *      GetHistQmax(){return fDeDxQmax;};
  THnSparseS *      GetHistQtot(){return fDeDxQtot;};
  THnSparseS *      GetHistRatioMaxTot(){return fDeDxRatioMaxTot;};
  THnSparseS *      GetHistRatioQmax(){return fDeDxRatioQmax;};
  THnSparseS *      GetHistRatioQtot(){return fDeDxRatioQtot;};
  THnSparseS *      GetHistRatioTruncQmax(){return fDeDxRatioTruncQmax;};
  THnSparseS *      GetHistRatioTruncQtot(){return fDeDxRatioTruncQtot;};
  //
  void SetMIPvalue(Float_t mip){fMIP = mip;};
  void SetLowerTrunc(Float_t lowerTrunc){fLowerTrunc = lowerTrunc;};
  void SetUpperTrunc(Float_t upperTrunc){fUpperTrunc = upperTrunc;};
  void SetUseShapeNorm(Bool_t useShapeNorm){fUseShapeNorm = useShapeNorm;};
  void SetUsePosNorm(Int_t usePosNorm){fUsePosNorm = usePosNorm;};
  void SetPadNorm(Int_t padNorm){fUsePadNorm = padNorm;};
  void SetIsCosmic(Bool_t isCosmic){fIsCosmic = isCosmic;};
  //
  //
  static void       BinLogX(THnSparse * h, Int_t axisDim);   // method for correct histogram binning
  void DumpTree(THnSparse * hndim, const char * outname);
  void DumpTrees();
  void     Process(AliESDtrack *track, Int_t runNo=-1){AliTPCcalibBase::Process(track,runNo);};
  void     Process(AliTPCseed *track){return AliTPCcalibBase::Process(track);}
private:
  //
  // parameter specifications
  //
  Float_t fMIP;
  Float_t fLowerTrunc;
  Float_t fUpperTrunc;
  Bool_t  fUseShapeNorm;
  Int_t  fUsePosNorm;
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
  //
  // ratio histograms
  //
  THnSparseS * fDeDxRatioMaxTot;              //  histogram which shows dEdx ratio (Qmax/Qtot) as a function of z,sin(phi),tan(theta),dEdx,dEdx*dl
  THnSparseS * fDeDxRatioQmax;   // dEdx ratio (tracklet/track) as a function of z,sin(phi),tan(theta),dEdx,dEdx*dl
  THnSparseS * fDeDxRatioQtot;   // dEdx ratio (tracklet/track) as a function of z,sin(phi),tan(theta),dEdx,dEdx*dl
  THnSparseS * fDeDxRatioTruncQtot;   // dEdx ratio (tracklet/track) as a function of z,sin(phi),tan(theta),dEdx,dEdx*dl
  THnSparseS * fDeDxRatioTruncQmax;   // dEdx ratio (tracklet/track) as a function of z,sin(phi),tan(theta),dEdx,dEdx*dl

  //
  AliTPCcalibPID(const AliTPCcalibPID&); 
  AliTPCcalibPID& operator=(const AliTPCcalibPID&); 

  ClassDef(AliTPCcalibPID, 1); 
};

#endif


