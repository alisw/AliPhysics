#ifndef ALITPCCALIBGAINMULT_H
#define ALITPCCALIBGAINMULT_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliTPCcalibBase.h"
#include "AliTPCCalPad.h"
#include "TH3F.h"
#include "TF1.h"
#include "THnSparse.h"
class TH1F;
class TList;
class AliESDEvent;
class AliESDtrack;
class AliTPCseed;

#include "TTreeStream.h"


class AliTPCcalibGainMult:public AliTPCcalibBase {
public:
  AliTPCcalibGainMult(); 
  AliTPCcalibGainMult(const Text_t *name, const Text_t *title);
  virtual ~AliTPCcalibGainMult();
  //
  virtual void           Process(AliESDEvent *event);
  virtual Long64_t       Merge(TCollection *li);
  virtual void           Analyze();
  //
  TH1F   *          GetHistNTracks() const {return fHistNTracks;};
  TH1F   *          GetHistClusterShape() const {return fHistClusterShape;};
  TH3F   *          GetHistQA() const {return fHistQA;};
  //
  THnSparseF *      GetHistGainSector() const {return fHistGainSector;};
  THnSparseF *      GetHistPadEqual() const {return fHistPadEqual;};
  THnSparseF *      GetHistGainMult() const {return fHistGainMult;};
  //
  void SetMIPvalue(Float_t mip){fMIP = mip;};
  void SetLowerTrunc(Float_t lowerTrunc){fLowerTrunc = lowerTrunc;};
  void SetUpperTrunc(Float_t upperTrunc){fUpperTrunc = upperTrunc;};
  void SetUseMax(Bool_t useMax){fUseMax = useMax;};
  //
  //
  static void       BinLogX(const TH1 *h);   // method for correct histogram binning
  void     Process(AliESDtrack *track, Int_t runNo=-1){AliTPCcalibBase::Process(track,runNo);};
  void     Process(AliTPCseed *track){return AliTPCcalibBase::Process(track);}
  //
  void     MakeLookup(THnSparse * hist, Char_t * outputFile);
  //
  void     UpdateGainMap();
  void     UpdateClusterParam();


private:
  //
  // parameter specifications
  //
  Float_t fMIP;                  // MIP position to be in fMIP
  Float_t fLowerTrunc;           // lower truncation for dEdx
  Float_t fUpperTrunc;           // upper truncation for dEdx
  //
  Bool_t fUseMax;                 // flag if Qmax or Qtot should be used
  //
  // histograms
  //
  TH1F  *fHistNTracks;            //  histogram showing number of ESD tracks per event
  TH1F  *fHistClusterShape;       //  histogram to check the cluster shape
  TH3F  *fHistQA;                 //  dE/dx histogram showing the final spectrum
  //
  THnSparseF * fHistGainSector;   //  histogram which shows MIP peak for each of the 3x36 sectors (pad region)
  THnSparseF * fHistPadEqual;     //  histogram for the equalization of the gain in the different pad regions -> pass0
  THnSparseF * fHistGainMult;     //  histogram which shows decrease of MIP signal as a function
  //
  AliTPCcalibGainMult(const AliTPCcalibGainMult&); 
  AliTPCcalibGainMult& operator=(const AliTPCcalibGainMult&); 

  ClassDef(AliTPCcalibGainMult, 1); 
};

#endif


