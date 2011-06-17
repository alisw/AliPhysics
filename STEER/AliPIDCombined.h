#ifndef ALIPIDCOMBINED_H
#define ALIPIDCOMBINED_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//---------------------------------------------------------------//
//        Base class for combining PID response of               //
//        of different detectors                                 //
//        and compute Bayesian probabilities                     //
//                                                               //
//   Origin: Pietro Antonioli, INFN-BO, pietro.antonioli@cern.ch //
//                                                               //
//---------------------------------------------------------------//



#include <TNamed.h>
#include <AliPID.h>
#include <AliPIDResponse.h>
#include <TH1F.h>

//class TH1;
class AliPIDResponse;

class AliPIDCombined : public TNamed {
public:
  AliPIDCombined();
  AliPIDCombined(const TString& name, const TString& title);
  virtual ~AliPIDCombined();

  void SetDetectorMask(Int_t mask) {fDetectorMask=mask;}
  Int_t GetDetectorMask() const {return fDetectorMask;}
  void SetRejectMismatchMask(Int_t mask) {fRejectMismatchMask=mask;}
  Int_t GetRejectMismatchMask() const {return fRejectMismatchMask;}
  void SetEnablePriors(Bool_t flag) {fEnablePriors=flag;}
  Bool_t GetEnablePriors() const {return fEnablePriors;}
  void SetPriorDistribution(AliPID::EParticleType type,TH1F *prior);
  //  const TH1* GetPriorDistribution(AliPID::EParticleType type) const {return (TH1*)fPriorsDistributions[type];}
  TH1* GetPriorDistribution(AliPID::EParticleType type)  const {return (TH1*)fPriorsDistributions[type];}
	
  UInt_t ComputeProbabilities(const AliVTrack *track, const AliPIDResponse *response, Double_t* bayesProbabilities) const;
  void SetSelectedSpecies(Int_t selectedSpecies) {fSelectedSpecies = selectedSpecies;}
  Int_t GetSelectedSpecies() const {return fSelectedSpecies;}

protected:
  void GetPriors(const AliVTrack *track,Double_t* priors) const;
  void ComputeBayesProbabilities(Double_t* bayesProbabilities,const Double_t* probDensity, const Double_t* priors) const;
  void SetCombinedStatus(const AliPIDResponse::EDetPidStatus status,UInt_t *mask, const AliPIDResponse::EDetCode bit, Double_t* p) const;

private:
  AliPIDCombined(const AliPIDCombined&);
  AliPIDCombined &operator=(const AliPIDCombined&);

  Int_t fDetectorMask;       // Detectors included in combined pid
  Int_t fRejectMismatchMask; // Detectors set return flat prob. if mismatch detected 
  Bool_t fEnablePriors;      // Enable bayesian PID (if kFALSE priors set flat)
  Int_t fSelectedSpecies;    // Number of selected species to study
  TH1F *fPriorsDistributions[AliPID::kSPECIES+AliPID::kSPECIESLN]; // priors

  ClassDef(AliPIDCombined,1);
};

#endif
