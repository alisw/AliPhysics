#ifndef ALIAODPIDHF_H
#define ALIAODPIDHF_H

/* Copyright(c) 1998-2006, ALICE Experiment at CERN, All rights reserved. *
 *  * See cxx source for full Copyright notice                               */

//***********************************************************
//// Class AliAODPidHF
//// class for PID with AliAODRecoDecayHF 
//// Authors: D. Caffarri caffarri@bo.infn.it, A.Dainese andrea.dainese@pd.infn.it, S. Dash dash@to.infn.it, F. Prino prino@to.infn.it, R. Romita r.romita@gsi.de, Y. Wang yifei@pi0.physi.uni-heidelberg.de
////***********************************************************

#include "AliAODPid.h"
#include "AliAODTrack.h"

class AliAODPidHF : public AliAODPid{

 public:

 AliAODPidHF();
 AliAODPidHF(const AliAODPidHF& pid);
 AliAODPidHF& operator=(const AliAODPidHF& pid);
 virtual ~AliAODPidHF();

 //Setters
 void SetSigma(Double_t sigma){fSigma=sigma;return;}
 void SetPriors(Double_t *priors){fPriors=priors;return;}
 

 Int_t RawSignalPID (AliAODTrack *track, TString detector);
 Bool_t IsKaonRaw (AliAODTrack *track, TString detector);
 Bool_t IsPionRaw (AliAODTrack *track, TString detector);
 Bool_t IsProtonRaw (AliAODTrack *track, TString detector);
 Bool_t IsElectronRaw (AliAODTrack *track, TString detector);
 void BayesianProbability(AliAODTrack *track,TString detectors,Double_t *pid);
 void CombinedProbability(AliAODTrack *track,Bool_t *type); //0 = pion, 1 = kaon, 2 = proton
 Bool_t CheckStatus(AliAODTrack *track,TString detectors);

 protected:

 Int_t ApplyPidTPCRaw(AliAODTrack *track,Int_t specie);
 Int_t ApplyPidTOFRaw(AliAODTrack *track,Int_t specie);
 Int_t ApplyPidITSRaw(AliAODTrack *track,Int_t specie);
 void BayesianProbabilityITS(AliAODTrack *track,Double_t *prob);
 void BayesianProbabilityTPC(AliAODTrack *track,Double_t *prob);
 void BayesianProbabilityTOF(AliAODTrack *track,Double_t *prob);
 void BayesianProbabilityTRD(AliAODTrack *track,Double_t *prob);

 private:
 Double_t fSigma; // sigma for the raw signal PID 
 Double_t *fPriors; // set of priors


 ClassDef(AliAODPidHF,3)

};

#endif
