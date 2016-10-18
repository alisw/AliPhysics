#ifndef ALIGENPARAM_H
#define ALIGENPARAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//
// Class to generate particles from using parametrized pT and y distributions.
// Distributions are obtained from pointer to object of type AliGenLib.
// (For example AliGenMUONlib)
//
// andreas.morsch@cern.ch
//

#include "AliGenMC.h"

class AliPythia;
class TParticle;
class AliGenLib;
class TF1;

typedef enum { kAnalog, kNonAnalog} Weighting_t;
//-------------------------------------------------------------
class AliGenParam : public AliGenMC
{
public:
  AliGenParam();
  AliGenParam(Int_t npart, const AliGenLib * Library, Int_t param,   const char*  tname = 0);
  AliGenParam(Int_t npart, Int_t param, const char* tname = 0, const char*  name  = 0);
  AliGenParam(Int_t npart, Int_t param,
              Double_t (*PtPara)(const Double_t*, const Double_t*),
              Double_t (*YPara )(const Double_t*, const Double_t*),
              Double_t (*V2Para)(const Double_t*, const Double_t*),
              Int_t    (*IpPara)(TRandom*)           );

  virtual ~AliGenParam();
  virtual void GenerateN(Int_t ntimes);
  virtual void Generate();
  virtual void Init();
  // select particle type
  virtual void SetParam(Int_t param) {fParam = param;}
  //Setting the flag for Background transportation while using SetForceDecay()
  void SetSelectAll(Bool_t selectall) {fSelectAll = selectall;}
  // force decay type
  virtual void SetWeighting(Weighting_t flag = kAnalog) {fAnalog = flag;}
  virtual void SetDeltaPt(Float_t delta=0.01) {fDeltaPt = delta;}
  virtual void SetDecayer(AliDecayer* decayer) {fDecayer = decayer;}
  virtual void SetForceGammaConversion(Bool_t force=kTRUE) {fForceConv = force;}
  virtual void SetKeepParent(Bool_t keep=kTRUE){fKeepParent= keep;} //Store parent even if it does not have childs within cuts
  virtual void SetKeepIfOneChildSelected(Bool_t keep=kTRUE){fKeepIfOneChildSelected = keep;} //Accept parent and child even if other children are not within cut.
  virtual void SetPreserveFullDecayChain(Int_t preserve = kFALSE) {fPreserveFullDecayChain = preserve;} //Prevent flagging(/skipping) of decay daughter particles; preserves complete forced decay chain

  virtual void Draw(const char * opt);
  TF1 *  GetPt() { return fPtPara;}
  TF1 *  GetY() {return fYPara;}
  Float_t GetRelativeArea(Float_t ptMin, Float_t ptMax, Float_t yMin, Float_t yMax, Float_t phiMin, Float_t phiMax);

  static TVector3 OrthogonalVector(TVector3 &inVec);
  static void RotateVector( Double_t *pin, Double_t *pout, Double_t costheta, Double_t sintheta,
                            Double_t cosphi, Double_t sinphi);
  static double ScreenFunction1(double d);
  static double ScreenFunction2(double d);
  double RandomEnergyFraction(double Z, double E);
  double RandomPolarAngle();
  double RandomMass(Double_t mh);
  Int_t VirtualGammaPairProduction(TClonesArray *particles, Int_t nPart);
  Int_t ForceGammaConversion(TClonesArray *particles, Int_t nPart);
  virtual void SetSeed(UInt_t /*seed*/) {;}

  // allow explicit setting of functions in case of streaming
  void SetParamsExplicitly ( const AliGenLib * Library, Int_t param, const char* tname) {
    fPtParaFunc = Library->GetPt(param, tname);
    fYParaFunc  = Library->GetY (param, tname);
    fIpParaFunc = Library->GetIp(param, tname);
    fV2ParaFunc = Library->GetV2(param, tname);
  }

  // retrive particle type
  Int_t GetParam() {return fParam; }

protected:
  Double_t (*fPtParaFunc)(const Double_t*, const Double_t*); //! Pointer to Pt parametrisation function
  Double_t (*fYParaFunc )(const Double_t*, const Double_t*); //! Pointer to Y parametrisation function
  Int_t    (*fIpParaFunc )(TRandom*);    //! Pointer to particle type parametrisation function
  Double_t (*fV2ParaFunc )(const Double_t*, const Double_t*);//! Pointer to V2 parametrisation function
  TF1* fPtPara;              // Transverse momentum parameterisation
  TF1* fYPara;               // Rapidity parameterisation
  TF1*        fV2Para;       // v2 parametrization
  TF1*        fdNdPhi;       // Phi distribution depending on v2
  Int_t       fParam;        // Parameterisation type
  Float_t     fdNdy0;        // central multiplicity per event
  Float_t     fYWgt;         // Y-weight
  Float_t     fPtWgt;        // Pt-weight
  Float_t     fBias;         // Biasing factor
  Int_t       fTrials;       // Number of trials
  Float_t     fDeltaPt;      // pT sampling in steps of fDeltaPt
  Bool_t      fSelectAll;    // Flag for transportation of Background while using SetForceDecay()
  AliDecayer  *fDecayer;     // ! Pointer to pythia object for decays
  Bool_t      fForceConv;    //
  Bool_t      fKeepParent;   //  Store parent even if it does not have childs within cuts
  Bool_t      fKeepIfOneChildSelected; //Accept parent and child even if other children are not within cut.
  Bool_t      fPreserveFullDecayChain; //Prevent flagging(/skipping) of decay daughter particles; preserves complete forced decay chain

private:
  AliGenParam(const AliGenParam &Param);
  AliGenParam & operator=(const AliGenParam & rhs);

  ClassDef(AliGenParam, 4) // Generator using parameterised pt- and y-distribution
};
#endif
