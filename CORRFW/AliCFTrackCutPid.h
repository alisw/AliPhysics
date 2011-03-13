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
// author : A. Mastroserio
//

#ifndef ALICFTRACKCUTPID_H
#define ALICFTRACKCUTPID_H

#include "AliCFCutBase.h"
#include "AliPID.h"
#include "AliESDtrack.h"
#include "AliAODTrack.h"
#include <TString.h>
#include <TObject.h>
#include <TH1F.h>
#include <TF1.h>
//__________________________________________________________________________________
// CUT ON TRACK PID
//__________________________________________________________________________________
class AliESDtrack;

class AliCFTrackCutPid : public AliCFCutBase
{
  public :
    AliCFTrackCutPid() ;
  AliCFTrackCutPid(const Char_t* name, const Char_t* title) ;
  AliCFTrackCutPid(const AliCFTrackCutPid& c) ;
  AliCFTrackCutPid& operator=(const AliCFTrackCutPid& c) ;
  
  virtual ~AliCFTrackCutPid();
  
  enum EDetType {kITS = 0, kTPC, kTRD,  kTOF, kHMPID=4, kNoDet=-11};
  enum EDetNum {kNdets=5};
  enum ENoId {kCheckProb = -10, kCheckResp = -11, kDetRestr = -12};
  
  // Setters 
  
  void SetDetectors(TString dets);                            
  void SetPriors(Double_t r[AliPID::kSPECIES]);                    
  void SetProbabilityCut(Double_t cut) {fCut=cut;}                  
  void SetParticleType(Int_t iType, Bool_t tocombine) {fgParticleType=iType; fgIsComb=tocombine;} 
  void SetMinDiffResp(Bool_t check, Double_t mindiff) {fCheckResponse=check; fMinDiffResponse=mindiff;}  
  void SetMinDiffProb(Bool_t check, Double_t mindiff) {fCheckSelection=check; fMinDiffProbability=mindiff;} 
  void SetPriorFunctions(TF1 *func[AliPID::kSPECIES]);
  void SetANDstatus(TString dets);
  void SetDetectorProbabilityRestriction(TString det, Int_t iPart, Double_t upperprob); 
  void SetHistogramAxis(Int_t nbins, Double_t xmin, Double_t xmax) {fNbins=nbins; fXmin = xmin; fXmax = xmax;}
  void SetAODmode(Bool_t isaod = kFALSE) {fgIsAOD=isaod;}  
  void SetProbThreshold(Double_t value) {fProbThreshold=value;}
 
  //returns the track identification number  
  Int_t GetID(ULong_t status[kNdets+1], Double_t pid[kNdets+1][AliPID::kSPECIES]) const;  
  //returns the track identification number in caso of an AliAODTrack
  Int_t GetAODID(AliAODTrack *aodtrack) const;
  
 
  //main 
  virtual Bool_t IsSelected(TObject *track); 
  Bool_t IsSelected(TList* /*list*/) {return kTRUE;}  
  //histos are added to a list
  void AddQAHistograms(TList *qalist);
  
  
 private:

  //loads the track detector responses and the track status
  void TrackInfo(const AliESDtrack *pTrk,ULong_t status[kNdets+1], Double_t pid[kNdets+1][AliPID::kSPECIES]) const;

  //identifies the track
  Int_t Identify(Double_t pid[AliPID::kSPECIES]) const;

  //identifies the track filling the QA histograms
  Int_t IdentifyQA(const Double_t pid[AliPID::kSPECIES],Int_t idets) const;

  void SetPPriors(AliESDtrack *pTrk);                          
  ULong_t StatusForAND(ULong_t status[kNdets+1]) const; 
  void InitialiseHisto();
  void DefineHistograms();                                 // histo booking  
  Bool_t Check(const Double_t *p, Int_t iPsel, Double_t minDiff) const;
  void CombPID(ULong_t status[kNdets+1],Double_t pid[kNdets+1][AliPID::kSPECIES],Double_t *combpid) const;
  
  Double_t fCut;                                            // probability cut
  Double_t fMinDiffResponse;                                // minimum difference between detector resposes
  Double_t fMinDiffProbability;                             // minimum difference between probability values
  Int_t fgParticleType;                                     // requested particle type
  Bool_t fgIsComb;                                          // flag for the combined pid
  Bool_t fgIsAOD;                                           // flag for AOD QA histograms
  Bool_t fCheckResponse;                                    // flag to check the minimum difference of det responsess
  Bool_t fCheckSelection;                                   // flag to check the minimum difference of probabilities
  Bool_t fIsPpriors;                                        // flag for momentum dependent priors
  Bool_t fIsDetAND;                                         // flag for AND with multiple detectors
  Double_t fXmin;                                           // x min QA histo
  Double_t fXmax;                                           // x max QA histo
  Int_t fNbins;                                             // n bins QA histo 
  Int_t fDetRestr;                                          // id of the detector for the restriction
  Int_t fiPartRestr;                                        // id of the particle for the restriction
  Double_t fDetProbRestr;                                   // probability restriction value
  Double_t fProbThreshold;                                  // if different from 0, the assigned PID will be set to 
                                                            // fgParticleType if the probability is larger than this threshold,
                                                            // regardless it is the highest or not (!)

  Double_t fPriors[AliPID::kSPECIESN];                       // a priori concentrations
  TF1 *fPriorsFunc[AliPID::kSPECIES];                       // momentum dependent priors
  Bool_t fDets[kNdets];                                     // boolean(s) corresponding to the chosen detector(s) 
  Bool_t fDetsInAnd[kNdets];                                // detector to be in AND for the combined PID
  TH1F *fhResp[kNdets][AliPID::kSPECIES];                   // QA histo
  TH1F *fhProb[kNdets][AliPID::kSPECIES];                   // QA histo
  TH1F *fhCombResp[AliPID::kSPECIES];                       // QA histo
  TH1F *fhCombProb[AliPID::kSPECIES];                       // QA histo
  
  ClassDef(AliCFTrackCutPid,1);
};
#endif
    
