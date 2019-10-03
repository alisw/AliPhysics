//
// This cut implements all the checks done to accept a track as a proton or pion
// for the pPb analysis using 2013 runs. 
// It is based on standard cuts on track quality and nsigma cuts
// with respect to the TPC and TOF signals for the PID.
//

#include <Riostream.h>

#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliRsnCutDelta.h"

ClassImp(AliRsnCutDelta)

//__________________________________________________________________________________________________
AliRsnCutDelta::AliRsnCutDelta(const char *name,AliPID::EParticleType pid, Bool_t TPCMethod, Bool_t usePPb) :
AliRsnCut(name, AliRsnTarget::kDaughter),
  fPID(pid),
  fNSigmaTPC(3.0),
  fLimitTPC(0.350),
  fNSigmaTOF(3.0),
  fTPCMethod(TPCMethod),
  fNSigmaTPCProton(3.0),
  fNSigmaTPCPion(3.0),
  fNSigmaTOFProton(3.0),
  fNSigmaTOFPion(3.0),
  fTPCMomProton(1.1),
  fTOFMomProton(2.5),
  fEtaRange(0.8),
  fTPCNCluster(70),
  fPtDepDCASigma(1.0),
  fCutQuality(Form("%sQuality", name))
{
  //
  // Constructor
  // Initialize the contained cuts and sets defaults
  //    
  if(usePPb) {
    Double_t a=0; a=fPtDepDCASigma*0.0105; 
    Double_t b=0; b=fPtDepDCASigma*0.0350;   
    // track quality 
    fCutQuality.SetPtRange(0.15, 1E+20);
    fCutQuality.SetEtaRange(-fEtaRange, fEtaRange);
    fCutQuality.SetDCARPtFormula(Form("%f+%f/pt^1.1",a,b));    
    fCutQuality.SetDCAZmax(2.0);
    fCutQuality.SetSPDminNClusters(1);
    fCutQuality.SetMinNCrossedRowsTPC(70, kTRUE);
    fCutQuality.SetMinNCrossedRowsOverFindableClsTPC(0.8, kTRUE);
    fCutQuality.SetMaxChi2TPCConstrainedGlobal(36);
    fCutQuality.SetTPCmaxChi2(4.0);
    fCutQuality.SetRejectKinkDaughters();
    fCutQuality.SetAODTestFilterBit(5);
  }  else {
    Double_t a=0; a=fPtDepDCASigma*0.0026*7; 
    Double_t b=0; b=fPtDepDCASigma*0.05*7;   
    // track quality 
    fCutQuality.SetPtRange(0.15, 1E+20);
    fCutQuality.SetEtaRange(-fEtaRange, fEtaRange);
    fCutQuality.SetDCARPtFormula(Form("%f+%f/pt^1.01",a,b));    //("0.0182+0.0350/pt^1.01");
    fCutQuality.SetDCAZmax(2.0);
    fCutQuality.SetSPDminNClusters(1);
    fCutQuality.SetITSminNClusters(0);
    fCutQuality.SetITSmaxChi2(1E+20);
    fCutQuality.SetTPCminNClusters(fTPCNCluster);
    fCutQuality.SetTPCmaxChi2(4.0);
    fCutQuality.SetRejectKinkDaughters();
    fCutQuality.SetAODTestFilterBit(5);
  }
   
   
}
//
//__________________________________________________________________________________________________
Bool_t AliRsnCutDelta::IsSelected(TObject *obj)
{
  //
  // Global check
  //

  // coherence check
  if (!TargetOK(obj)) return kFALSE;
   
  // check track
  AliVTrack *track = fDaughter->Ref2Vtrack();
  if (!track) {
    if (!fDaughter->GetRef()) AliWarning("NULL ref");
    return kFALSE;
  }
   
  // check flags
  if ((track->GetStatus() & AliESDtrack::kTPCin   ) == 0) return kFALSE;
  if ((track->GetStatus() & AliESDtrack::kTPCrefit) == 0) return kFALSE;
  if ((track->GetStatus() & AliESDtrack::kITSrefit) == 0) return kFALSE;
   
  // quality
  if (!fCutQuality.IsSelected(obj)) return kFALSE;


   
  // check initialization of PID object
  AliPIDResponse *pid = fEvent->GetPIDResponse();
  if (!pid) {
    AliFatal("NULL PID response");
    return kFALSE;
  }
   
  // check if TOF is matched
  // and computes all values used in the PID cut
  Bool_t   isTOF        	= MatchTOF(track);
  Double_t pTPC         	= track->GetTPCmomentum();
  Double_t pTOF            	= track->P();
  Double_t nsTPC        	= TMath::Abs(pid->NumberOfSigmasTPC(track, fPID));
  Double_t nsTPCProton  	= TMath::Abs(pid->NumberOfSigmasTPC(track, AliPID::kProton));
  Double_t nsTPCPion	  	= TMath::Abs(pid->NumberOfSigmasTPC(track, AliPID::kPion));
  Double_t nsTOFProton  	= TMath::Abs(pid->NumberOfSigmasTOF(track, AliPID::kProton));
  Double_t nsTOF        	= TMath::Abs(pid->NumberOfSigmasTOF(track, fPID));
  Double_t maxTPC       	= 1E20;
  Double_t maxTOF       	= 1E20;

  if (fTPCMethod) {
    if(fPID == AliPID::kPion) { 
 
      if(nsTPCPion>fNSigmaTPCPion) return kFALSE;
      maxTPC = fNSigmaTPCPion; 
      if(nsTPC > maxTPC) return kFALSE;
      return kTRUE;

    }
    else if (fPID == AliPID::kProton) { 

      if(nsTPCPion<=fNSigmaTPCPion) return kFALSE;
      if(pTPC>=fTPCMomProton) return kFALSE; 
      if(nsTPCProton>fNSigmaTPCProton) return kFALSE;
      maxTPC = fNSigmaTPCProton; 
      if(nsTPC > maxTPC) return kFALSE;
      return kTRUE;

    }

  } else {

    if(fPID == AliPID::kProton) {  

      if(!isTOF) return kFALSE;
      if(nsTPCProton>fNSigmaTPCProton) return kFALSE;
      if(pTOF >= fTOFMomProton) return kFALSE;
      if(nsTOFProton>fNSigmaTOFProton) return kFALSE;
      maxTOF = fNSigmaTOFProton;
      if (nsTOF > maxTOF) return kFALSE; 
      return kTRUE;

    }

    else if (fPID == AliPID::kPion) { 

      if(isTOF && nsTPCProton<=fNSigmaTPCProton && pTOF < fTOFMomProton && nsTOFProton<=fNSigmaTOFProton) return kFALSE; 

      if(nsTPCPion>fNSigmaTPCPion) return kFALSE;
      maxTPC = fNSigmaTPCPion;
      if (nsTPC > maxTPC) return kFALSE; 
      return kTRUE;

    }//pion pid

  } //else tpcmethod
  return kTRUE;
}

//__________________________________________________________________________________________________
Bool_t AliRsnCutDelta::MatchTOF(const AliVTrack *vtrack)
{
//
// Checks if the track has matched the TOF detector
//

   if (!vtrack) {
      AliWarning("NULL argument: impossible to check status");
      return kFALSE;
   }

   Bool_t isTOFout = ((vtrack->GetStatus() & AliESDtrack::kTOFout) != 0);
   Bool_t isTIME   = ((vtrack->GetStatus() & AliESDtrack::kTIME) != 0);

   return (isTOFout && isTIME);
}
