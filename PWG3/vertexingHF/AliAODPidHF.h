#ifndef ALIAODPIDHF_H
#define ALIAODPIDHF_H

/* Copyright(c) 1998-2006, ALICE Experiment at CERN, All rights reserved. *
 *  * See cxx source for full Copyright notice                               */

//***********************************************************
//// Class AliAODPidHF
//// class for PID with AliAODRecoDecayHF 
//// Authors: D. Caffarri caffarri@pd.infn.it, A.Dainese andrea.dainese@pd.infn.it, S. Dash dash@to.infn.it, F. Prino prino@to.infn.it, R. Romita r.romita@gsi.de, Y. Wang yifei@pi0.physi.uni-heidelberg.de
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
 void SetSigma(Double_t *sigma){fnSigma=sigma;return;}
 void SetSigma(Int_t idet,Double_t sigma){fnSigma[idet]=sigma;return;}
 void SetTofSigma(Double_t sigma){fTOFSigma=sigma;return;}
 void SetPriors(Double_t *priors){fPriors=priors;return;}
// void SetPLimit(Double_t *plim){fPLimit=plim;return;}
 void SetPLimit(Double_t *plim){for(Int_t i=0;i<2;i++){fPLimit[i]=plim[i];}return;}
 void SetAsym(Bool_t asym){fAsym=asym;return;}
 void SetTPC(Bool_t tpc){fTPC=tpc;return;}
 void SetTOF(Bool_t tof){fTOF=tof;return;}
 void SetITS(Bool_t its){fITS=its;return;}
 void SetTRD(Bool_t trd){fTRD=trd;return;}
 void SetMatch(Bool_t match){fMatch=match;return;}
 void SetCompat(Bool_t comp){fCompat=comp;return;}
 

 Int_t RawSignalPID (AliAODTrack *track, TString detector) const;
 Bool_t IsKaonRaw (AliAODTrack *track, TString detector) const;
 Bool_t IsPionRaw (AliAODTrack *track, TString detector) const;
 Bool_t IsProtonRaw (AliAODTrack *track, TString detector) const;
 Bool_t IsElectronRaw (AliAODTrack *track, TString detector) const;
 void BayesianProbability(AliAODTrack *track,Double_t *pid) const;
 void CombinedProbability(AliAODTrack *track,Bool_t *type) const; //0 = pion, 1 = kaon, 2 = proton
 Bool_t CheckStatus(AliAODTrack *track,TString detectors) const;

 Bool_t TPCRawAsym(AliAODTrack* track,Int_t specie) const;
 Int_t MatchTPCTOF(AliAODTrack *track,Int_t mode,Int_t specie,Bool_t compat);

 Int_t MakeRawPid(AliAODTrack *track,Int_t specie); //general method to perform PID using raw signals



 protected:

 Int_t ApplyPidTPCRaw(AliAODTrack *track,Int_t specie) const;
 Int_t ApplyPidTOFRaw(AliAODTrack *track,Int_t specie) const;
 Int_t ApplyPidITSRaw(AliAODTrack *track,Int_t specie) const;
 void BayesianProbabilityITS(AliAODTrack *track,Double_t *prob) const;
 void BayesianProbabilityTPC(AliAODTrack *track,Double_t *prob) const;
 void BayesianProbabilityTOF(AliAODTrack *track,Double_t *prob) const;
 void BayesianProbabilityTRD(AliAODTrack *track,Double_t *prob) const;

 private:
 Int_t fnNSigma; // number of sigmas
 Double_t *fnSigma; // [fnNSigma], sigma for the raw signal PID: 0-2 for TPC, 3 for TOF, 4 for ITS 
 Double_t fTOFSigma; // TOF precision 
 Int_t fnPriors; //number of priors
 Double_t *fPriors; // [fnPriors], set of priors
 Int_t fnPLimit; //number of Plimit
 Double_t *fPLimit; // [fnPLimit], limit of p intervals for asimmetric PID: fPLimit<p[0], fPLimit[0]<p<fPLimit[1], p>fPLimit[1]
 Bool_t fAsym; // asimmetric PID required
 Bool_t fTPC; //switch to include or exclude TPC 
 Bool_t fTOF; // switch to include or exclude TOF
 Bool_t fITS; //switch to include or exclude ITS
 Bool_t fTRD; // switch to include or exclude TRD
 Bool_t fMatch; //switch to combine the info from more detectors: 1 = || , 2 = &, 3 = p region
 Bool_t fCompat; // compatibility region : useful only if fMatch=1
 


 ClassDef(AliAODPidHF,3) // AliAODPid for heavy flavor PID

};

#endif
