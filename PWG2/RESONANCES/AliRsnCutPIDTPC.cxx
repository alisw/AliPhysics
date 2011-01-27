//
// Class AliRsnCutPIDTPC
//
// General implementation of a single cut strategy, which can be:
// - a value contained in a given interval  [--> IsBetween()   ]
// - a value equal to a given reference     [--> MatchesValue()]
//
// In all cases, the reference value(s) is (are) given as data members
// and each kind of cut requires a given value type (Int, UInt, Double),
// but the cut check procedure is then automatized and chosen thanks to
// an enumeration of the implemented cut types.
// At the end, the user (or any other point which uses this object) has
// to use the method IsSelected() to check if this cut has been passed.
//
// authors: Martin Vala (martin.vala@cern.ch)
//          Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//

#include <Riostream.h>

#include "AliESDpid.h"
#include "AliTOFT0maker.h"
#include "AliTOFcalib.h"
#include "AliCDBManager.h"
#include "AliITSPIDResponse.h"

#include "AliRsnEvent.h"
#include "AliRsnDaughter.h"
#include "AliRsnCutPIDTPC.h"

ClassImp(AliRsnCutPIDTPC)

//_________________________________________________________________________________________________
AliRsnCutPIDTPC::AliRsnCutPIDTPC
(const char *name, AliPID::EParticleType type, Double_t momLimit, Double_t cut1, Double_t cut2) :
  AliRsnCut(name, AliRsnCut::kDaughter, 0.0, 0.0),
  fPIDtype(type),
  fMomentumLimit(momLimit),
  fLargeCut(TMath::Max(cut1, cut2)),
  fSmallCut(TMath::Min(cut1, cut2)),
  fESDpid(),
  fAODpid()
{
//
// Main constructor.
//
}

//_________________________________________________________________________________________________
AliRsnCutPIDTPC::AliRsnCutPIDTPC
(const AliRsnCutPIDTPC& copy) :
  AliRsnCut(copy),
  fPIDtype(copy.fPIDtype),
  fMomentumLimit(copy.fMomentumLimit),
  fLargeCut(copy.fLargeCut),
  fSmallCut(copy.fSmallCut),
  fESDpid(copy.fESDpid),
  fAODpid(copy.fAODpid)
{
//
// Copy constructor.
//
}

//_________________________________________________________________________________________________
AliRsnCutPIDTPC& AliRsnCutPIDTPC::operator=(const AliRsnCutPIDTPC& copy)
{
//
// Assignment operator
//

  AliRsnCut::operator=(copy);

  fPIDtype = copy.fPIDtype;
  fMomentumLimit = copy.fMomentumLimit;
  fLargeCut = copy.fLargeCut;
  fSmallCut = copy.fSmallCut;
  fESDpid = copy.fESDpid;
  fAODpid = copy.fAODpid;
  
  return (*this);
}

//_________________________________________________________________________________________________
void AliRsnCutPIDTPC::SetBBParam(Double_t p0, Double_t p1, Double_t p2, Double_t p3, Double_t p4)
{
//
// Properly set the Bethe-Bloch parameters in all places where it is needed.
//

  fESDpid.GetTPCResponse().SetBetheBlochParameters(p0, p1, p2, p3, p4);
  fAODpid.GetTPCResponse().SetBetheBlochParameters(p0, p1, p2, p3, p4);
}

//_________________________________________________________________________________________________
Bool_t AliRsnCutPIDTPC::IsSelected(TObject *object)
{
//
// Cut checker.
//

  // coherence check
  if (!TargetOK(object)) return kFALSE;
  
  // reject not TPC tracks
  // status is checked in the same way for all tracks
  AliVTrack *vtrack = dynamic_cast<AliVTrack*>(fDaughter->GetRef());
  if (!vtrack)
  {
    AliDebug(AliLog::kDebug + 2, Form("This object is not either an ESD nor AOD track, it is an %s", fDaughter->GetRef()->ClassName()));
    return kFALSE;
  }
  ULong_t status = (ULong_t)vtrack->GetStatus();
  if ((status & AliESDtrack::kTPCin) == 0)
  {
    AliDebug(AliLog::kDebug + 2, "Track is not found in TPC");
    return kFALSE;
  }
  
  // common evaluation variables
  Double_t mom;

  // retrieve real object type
  AliESDtrack *esdTrack = fDaughter->GetRefESDtrack();
  AliAODTrack *aodTrack = fDaughter->GetRefAODtrack();
  if (esdTrack) 
  {
    AliDebug(AliLog::kDebug + 2, "Checking an ESD track");
    
    mom = esdTrack->GetInnerParam()->P();
    
    AliTPCPIDResponse &tpcrsp = fESDpid.GetTPCResponse();
    fCutValueD = tpcrsp.GetNumberOfSigmas(mom, esdTrack->GetTPCsignal(), esdTrack->GetTPCsignalN(), fPIDtype);
  }
  else if (aodTrack)
  {
    AliDebug(AliLog::kDebug + 2, "Checking an AOD track");
    
    AliAODPid *pidObj = aodTrack->GetDetPid();
    mom = pidObj->GetTPCmomentum();
  
    fCutValueD = fAODpid.NumberOfSigmasTPC(aodTrack, fPIDtype);
  }
  else
  {
    AliDebug(AliLog::kDebug + 2, Form("This object is not either an ESD nor AOD track, it is an %s", fDaughter->GetRef()->ClassName()));
    return kFALSE;
  }
  
  // determine cut range from the momentum
  if (mom < fMomentumLimit)
  {
    fMinD = -fLargeCut;
    fMaxD =  fLargeCut;
  }
  else
  {
    fMinD = -fSmallCut;
    fMaxD =  fSmallCut;
  }
  
  // check the cut using the standard AliRsnCut facilities
  return OkRangeD();
}

//_________________________________________________________________________________________________
void AliRsnCutPIDTPC::Print(const Option_t *) const
{
//
// Print information on this cut
//

  AliInfo(Form("Cut name, type                  : %s %s", GetName(), ClassName()));
  AliInfo(Form("TPC PID cut: limit, large, small: %.3f %.3f %.3f", fMomentumLimit, fLargeCut, fSmallCut));
}
