#ifndef ALIEXTERNALCOMPARISON_H
#define ALIEXTERNALCOMPARISON_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


#include "TNamed.h"
#include "TMatrixD.h"
class THnSparse;
class TCollection;
class AliExternalTrackParam;
class TParticle;
class AliTrackReference;
class TObjArray;
 
class AliExternalComparison:public TNamed {
public:
  AliExternalComparison(); 
  AliExternalComparison(const Text_t *name, const Text_t *title);
  virtual ~AliExternalComparison();
  virtual Long64_t       Merge(TCollection *li);
  virtual void           Add(AliExternalComparison*comp);
  virtual void           Analyze();
  virtual void           Process(const AliExternalTrackParam *param0, const AliExternalTrackParam *param1);
  virtual void           Process(const AliExternalTrackParam *param0, TParticle *part);
  virtual void           Process(const AliExternalTrackParam *param0, TParticle *part, const AliTrackReference *ref);
  //
  virtual THnSparse      * GetHisto(Int_t ivar, Int_t type);
  void MakeComparisonTree(const char * outname);
  static AliExternalTrackParam *MakeExternalParam(TParticle *part);
  static AliExternalTrackParam *MakeExternalParam(TParticle *part, const AliTrackReference *ref);
  Bool_t                 AcceptPair(const AliExternalTrackParam *param0, const AliExternalTrackParam *param1);

  virtual void           SetDefaultRange(Float_t scale=0.3, Float_t arm=160, Int_t nbins=200);  
  virtual void           SetDefaultCuts();  

  //
  void SetDistCut(Float_t dP0, Float_t dP1,Float_t dP2,Float_t dP3, Float_t dP4);
  void SetPullDistCut(Float_t dnP0, Float_t dnP1,Float_t dnP2,Float_t dnP3, Float_t dnP4);
  void SetResolRange(Int_t param, Float_t min, Float_t max, Int_t nbins);

protected:
  void    MakeHistos();
public:
  TObjArray * fResolHistos;             // resolution histogram
  TObjArray * fPullHistos;              // pull       histogram
  TMatrixD  * fRangeMatrix;             // range matrix
  TMatrixD  * fCutMatrix;               // cut matrix
  //
  AliExternalComparison(const AliExternalComparison&); 
  AliExternalComparison& operator=(const AliExternalComparison&); 

  ClassDef(AliExternalComparison, 1); 
};

#endif


