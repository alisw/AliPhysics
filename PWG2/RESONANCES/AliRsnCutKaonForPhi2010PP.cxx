//
// This cut implements all the checks done to accept a track as a Kaon
// for the pp analysis using 2010 runs. 
// It is based on standard cuts on track quality and nsigma cuts
// with respect to the TPC and TOF signals for the PID.
//

#include <Riostream.h>

#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliESDpid.h"
#include "AliAODpidUtil.h"

#include "AliRsnCutKaonForPhi2010PP.h"

ClassImp(AliRsnCutKaonForPhi2010PP)

//__________________________________________________________________________________________________
AliRsnCutKaonForPhi2010PP::AliRsnCutKaonForPhi2010PP(const char *name) :
   AliRsnCut(name, AliRsnTarget::kDaughter, -3.0, 3.0),
   fNSigmaTPCLow(5.0),
   fNSigmaTPCHigh(3.0),
   fLimitTPC(0.350),
   fNSigmaTOF(3.0),
   fMyPID(0x0),
   fCutQuality(Form("%sQuality", name))
{
//
// Constructor
// Initialize the contained cuts and sets defaults
//

   // track quality
   //fCutQuality.AddStatusFlag(AliESDtrack::kTPCin   , kTRUE);
   //fCutQuality.AddStatusFlag(AliESDtrack::kTPCrefit, kTRUE);
   //fCutQuality.AddStatusFlag(AliESDtrack::kITSrefit, kTRUE);
   fCutQuality.SetPtRange(0.15, 1E+20);
   fCutQuality.SetEtaRange(-0.8, 0.8);
   fCutQuality.SetDCARPtFormula("0.0182+0.0350/pt^1.01");
   fCutQuality.SetDCAZmax(2.0);
   fCutQuality.SetSPDminNClusters(1);
   fCutQuality.SetITSminNClusters(0);
   fCutQuality.SetITSmaxChi2(1E+20);
   fCutQuality.SetTPCminNClusters(70);
   fCutQuality.SetTPCmaxChi2(4.0);
   fCutQuality.SetRejectKinkDaughters();
   fCutQuality.SetAODTestFilterBit(5);
}

//__________________________________________________________________________________________________
AliRsnCutKaonForPhi2010PP::AliRsnCutKaonForPhi2010PP(const AliRsnCutKaonForPhi2010PP &copy) :
   AliRsnCut(copy),
   fNSigmaTPCLow(copy.fNSigmaTPCLow),
   fNSigmaTPCHigh(copy.fNSigmaTPCHigh),
   fLimitTPC(copy.fLimitTPC),
   fNSigmaTOF(copy.fNSigmaTOF),
   fMyPID(copy.fMyPID),
   fCutQuality(copy.fCutQuality)
{
//
// Copy constructor
//
}

//__________________________________________________________________________________________________
AliRsnCutKaonForPhi2010PP& AliRsnCutKaonForPhi2010PP::operator=(const AliRsnCutKaonForPhi2010PP &copy)
{
//
// Assignment operator
//

   AliRsnCut::operator=(copy);
   
   fNSigmaTPCLow = copy.fNSigmaTPCLow;
   fNSigmaTPCHigh = copy.fNSigmaTPCHigh;
   fLimitTPC = copy.fLimitTPC;
   fNSigmaTOF = copy.fNSigmaTOF;
   fMyPID = copy.fMyPID;
   fCutQuality = copy.fCutQuality;
   
   return *this;
}

//__________________________________________________________________________________________________
void AliRsnCutKaonForPhi2010PP::InitMyPID(Bool_t isMC, Bool_t isESD)
{
//
// Initialize manual PID object
//

   if (isESD) 
      fMyPID = new AliESDpid(isMC);
   else
      fMyPID = new AliAODpidUtil(isMC);
}

//__________________________________________________________________________________________________
Bool_t AliRsnCutKaonForPhi2010PP::IsSelected(TObject *obj)
{
//
// Global check
//

   // coherence check
   if (!TargetOK(obj)) return kFALSE;
   
   // check track
   AliVTrack *track = dynamic_cast<AliVTrack*>(fDaughter->GetRef());
   if (!track) return kFALSE;
   
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
   
   // PID TPC :
   // depends on momentum
   // and if local PID object is initialized, it is used instead of that got from manager
   Double_t mom = (Double_t)TMath::Abs(track->GetTPCmomentum());
   if (mom < fLimitTPC) 
      SetRangeD(0.0, fNSigmaTPCLow);
   else
      SetRangeD(0.0, fNSigmaTPCHigh);
   if (fMyPID) 
      fCutValueD = TMath::Abs(fMyPID->NumberOfSigmasTPC(track, AliPID::kKaon));
   else
      fCutValueD = TMath::Abs(pid->NumberOfSigmasTPC(track, AliPID::kKaon));
   if (!OkRangeD()) return kFALSE;
   
   // if TOF is not matched, end here
   // otherwise check TOF
   if (!MatchTOF(track)) 
      return kTRUE;
   else {
      SetRangeD(0.0, fNSigmaTOF);
      fCutValueD = TMath::Abs(pid->NumberOfSigmasTOF(track, AliPID::kKaon));
      return OkRangeD();
   }
}
