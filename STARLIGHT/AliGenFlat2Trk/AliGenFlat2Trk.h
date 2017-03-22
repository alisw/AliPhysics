// -*- C++ -*-
// $Id$

#ifndef ALIGENFLAT2TRK_H_cm021218
#define ALIGENFLAT2TRK_H_cm021218

#include <TLorentzVector.h>
#include <TString.h>

#include <AliGenMC.h>

class TF1;
class AliGenEventHeader;

//-------------------------------------------------------------
class AliGenFlat2Trk : public AliGenMC
{
public:
  AliGenFlat2Trk();
  AliGenFlat2Trk(Int_t    pid,       // =13,
		 Double_t mInvMax,   // =1.5,
		 Double_t ptPairMin, // =0.90,
		 Double_t ptPairMax, // =0.150,
		 Double_t yPairMin,  // =-1
		 Double_t yPairMax,  // =+1
		 Int_t    tries,     // =1 no weigths, =10*1000 weights
		 Bool_t   runAsMCGenerator = kTRUE,
		 TString  fTheta = ""); // default: spin1->pi+pi-

  virtual ~AliGenFlat2Trk();

  void Generate();
  void Init();

  // set
  void SetDebug(Int_t debug) { fDebug = debug; }

  void SetPid(Int_t pid);

  Double_t SetMinvMin(Double_t mMin); // returns max(mMin, threshold)
  void     SetMinvMax(Double_t mMax) { fMinvMax = mMax; }

  void SetPtPairMin(Double_t ptPairMin) { fPtPairMin = ptPairMin; }
  void SetPtPairMax(Double_t ptPairMax) { fPtPairMax = ptPairMax; }

  void SetYPairMin(Double_t yMin) { fYPairMin = yMin; }
  void SetYPairMax(Double_t yMax) { fYPairMax = yMax; }

  // get
  Double_t GetMinvMin() const { return fMinvMin; }
  Double_t GetMinvMax() const { return fMinvMax; }

  Double_t GetPtPairMin() const { return fPtPairMin; }
  Double_t GetPtPairMax() const { return fPtPairMax; }

  Double_t GetYPairMin() const { return fYPairMin; }
  Double_t GetYPairMax() const { return fYPairMax; }

private:
  AliGenFlat2Trk(const AliGenFlat2Trk& gen);
  AliGenFlat2Trk& operator=(const AliGenFlat2Trk& gen);

public:
  struct PPair {
    PPair()
      : fV0()
      , fV1()
      , fWeight(1.0) {}

    PPair(const TLorentzVector &v0,
	  const TLorentzVector &v1,
	  Double_t weight)
      : fV0(v0)
      , fV1(v1)
      , fWeight(weight) {}

    Double_t M()    const { return (fV0+fV1).M(); }
    Double_t Perp() const { return (fV0+fV1).Perp(); }
    Double_t Y()    const { return (fV0+fV1).Y(); }

    TLorentzVector fV0;
    TLorentzVector fV1;
    Double_t       fWeight;
  } ;

  PPair GeneratePair() const;

private:
  PPair GeneratePair(Double_t rapidity,
 		     Double_t pt,
 		     Double_t mInv) const;

  Int_t    fPid;       // particle ID
  Double_t fMass;      // particle mass

  Double_t fMinvMin;   // invariant-mass lower bound
  Double_t fMinvMax;   // invariant-mass upper bound

  Double_t fPtPairMin; // pair-pT upper bound
  Double_t fPtPairMax; // pair-pT upper bound

  Double_t fYPairMin; // pair-rapidity upper bound
  Double_t fYPairMax; // pair-rapidity upper bound

  Int_t    fTries;
  Bool_t   fRunAsMCGenerator;

  Int_t    fDebug;      //
  Int_t    fEvent;      // internal event number
  Int_t    fNpartProd;  // total number of particles produced
  TF1*     fTheta;      // decay angle distribution

  AliGenEventHeader *fHeader; //!

  ClassDef(AliGenFlat2Trk,2);
} ;

#endif // ALIGENFLAT2TRK_H_cm021218
