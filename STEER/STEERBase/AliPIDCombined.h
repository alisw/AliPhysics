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
#include "AliPID.h"
#include "AliPIDResponse.h"
#include <TH1F.h>
#include <TH2F.h>

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
  
  void GetPriors(const AliVTrack *track,Double_t* p,const AliPIDResponse *response,UInt_t detUsed) const;
  
  void SetDefaultTPCPriors();
	
  UInt_t ComputeProbabilities(const AliVTrack *track, const AliPIDResponse *response, Double_t* bayesProbabilities,Double_t* priorsOwn=NULL) const;
  void SetSelectedSpecies(Int_t selectedSpecies) {fSelectedSpecies = selectedSpecies;}
  Int_t GetSelectedSpecies() const {return fSelectedSpecies;}

  static Float_t GetTOFmismatchProb() {return fTOFmismProb;}

protected:
  void GetPriors(const AliVTrack *track,Double_t* priors,Float_t centrality=-1,Bool_t isPPB=kFALSE) const;
  void ComputeBayesProbabilities(Double_t* bayesProbabilities,const Double_t* probDensity, const Double_t* priors, Double_t* probDensityMism=NULL) const;
  void SetCombinedStatus(const AliPIDResponse::EDetPidStatus status,UInt_t *mask, const AliPIDResponse::EDetCode bit) const;
  void SetCombinedStatus(const AliPIDResponse::EDetPidStatus status,UInt_t *mask, const AliPIDResponse::EDetCode bit, Double_t* p,const Float_t probMis) const;

private:
  AliPIDCombined(const AliPIDCombined&);
  AliPIDCombined &operator=(const AliPIDCombined&);

  Int_t fDetectorMask;       // Detectors included in combined pid
  Int_t fRejectMismatchMask; // Detectors set return flat prob. if mismatch detected 
  Bool_t fEnablePriors;      // Enable bayesian PID (if kFALSE priors set flat)
  Int_t fSelectedSpecies;    // Number of selected species to study
  TH1F *fPriorsDistributions[AliPID::kSPECIESC]; // priors
  Bool_t fUseDefaultTPCPriors; // switch to use Defaul TPC Priors
  static TH2F *fDefaultPriorsTPC[AliPID::kSPECIESC]; // Default priors for TPC tracks
  static TH2F *fDefaultPriorsTPCpPb[AliPID::kSPECIESC]; // Default priors for TPC tracks in pPb
  static Float_t fTOFmismProb; //TOF mismatch probability

  ClassDef(AliPIDCombined, 4);   // Combined PID using priors
};

#endif
