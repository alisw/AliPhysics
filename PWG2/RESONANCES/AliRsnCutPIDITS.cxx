//
// Class AliRsnCutPIDITS
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
#include "AliRsnCutPIDITS.h"

ClassImp(AliRsnCutPIDITS)

//_________________________________________________________________________________________________
AliRsnCutPIDITS::AliRsnCutPIDITS
(const char *name, AliPID::EParticleType type, Bool_t isMC, Double_t momLimit, Double_t cut1, Double_t cut2) :
  AliRsnCut(name, AliRsnCut::kDaughter, 0.0, 0.0),
  fPIDtype(type),
  fIsMC(isMC),
  fMomentumLimit(momLimit),
  fLargeCut(TMath::Max(cut1, cut2)),
  fSmallCut(TMath::Min(cut1, cut2)),
  fESDpid(),
  fAODpid()
{
//
// Main constructor.
//

  SetMC(isMC);
}

//_________________________________________________________________________________________________
AliRsnCutPIDITS::AliRsnCutPIDITS
(const AliRsnCutPIDITS& copy) :
  AliRsnCut(copy),
  fPIDtype(copy.fPIDtype),
  fIsMC(copy.fIsMC),
  fMomentumLimit(copy.fMomentumLimit),
  fLargeCut(copy.fLargeCut),
  fSmallCut(copy.fSmallCut),
  fESDpid(copy.fESDpid),
  fAODpid(copy.fAODpid)
{
//
// Copy constructor.
//

  SetMC(copy.fIsMC);
}

//_________________________________________________________________________________________________
AliRsnCutPIDITS& AliRsnCutPIDITS::operator=(const AliRsnCutPIDITS& copy)
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
  
  SetMC(copy.fIsMC);
  
  return (*this);
}

//_________________________________________________________________________________________________
void AliRsnCutPIDITS::SetMC(Bool_t yn)
{
//
// Properly set the PID response
//

  fIsMC = yn;
  AliITSPIDResponse itsresponse(fIsMC);
  fESDpid.GetITSResponse() = itsresponse;
  fAODpid.GetITSResponse() = itsresponse;
}

//_________________________________________________________________________________________________
Bool_t AliRsnCutPIDITS::IsSelected(TObject *object)
{
//
// Cut checker.
//

  // coherence check
  if (!TargetOK(object)) return kFALSE;
  
  // reject not ITS tracks
  // status is checked in the same way for all tracks
  AliVTrack *vtrack = dynamic_cast<AliVTrack*>(fDaughter->GetRef());
  if (!vtrack)
  {
    AliDebug(AliLog::kDebug + 2, Form("This object is not either an ESD nor AOD track, it is an %s", fDaughter->GetRef()->ClassName()));
    return kFALSE;
  }
  
  // check status, to know it track is an ITS+TPC or ITS standalone
  // and reject it if it is of none of the allowed types
  Bool_t isSA = kFALSE;
  if (IsITSTPC(vtrack)) isSA = kFALSE;
  else if (IsITSSA(vtrack)) isSA = kTRUE;
  else
  {
    AliWarning("Track is neither ITS+TPC nor ITS standalone");
    return kFALSE;
  }
  
  // common evaluation variables
  Double_t mom = vtrack->P();
  Int_t    k, nITSpidLayers;

  // retrieve real object type
  AliESDtrack *esdTrack = fDaughter->GetRefESDtrack();
  AliAODTrack *aodTrack = fDaughter->GetRefAODtrack();
  if (esdTrack) 
  {
    AliDebug(AliLog::kDebug + 2, "Checking an ESD track");
    
    // count PID layers and reject if they are too few
    nITSpidLayers = 0;
    UChar_t itsCluMap = esdTrack->GetITSClusterMap();
    for(k = 2; k < 6; k++) if(itsCluMap & (1 << k)) ++nITSpidLayers;
    if (nITSpidLayers < 3)
    {
      AliDebug(AliLog::kDebug+2, "Rejecting track with too few ITS pid layers");
      return kFALSE;
    }
  
    // create the PID response object and compute nsigma
    AliITSPIDResponse &itsrsp = fESDpid.GetITSResponse();
    fCutValueD = itsrsp.GetNumberOfSigmas(mom, esdTrack->GetITSsignal(), fPIDtype, nITSpidLayers, isSA);
  }
  else if (aodTrack)
  {
    AliDebug(AliLog::kDebug + 2, "Checking an AOD track");
    
    // count PID layers and reject if they are too few
    nITSpidLayers = 0;
    for(k = 2; k < 6; k++) if (TESTBIT(aodTrack->GetITSClusterMap(), k)) ++nITSpidLayers;
    if (nITSpidLayers < 3)
    {
      AliDebug(AliLog::kDebug+2, "Rejecting track with too few ITS pid layers");
      return kFALSE;
    }
    
    // compute nsigma
    fCutValueD = fAODpid.NumberOfSigmasITS(aodTrack, fPIDtype);
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
void AliRsnCutPIDITS::Print(const Option_t *) const
{
//
// Print information on this cut
//

  AliInfo(Form("Cut name, type                  : %s %s", GetName(), ClassName()));
  AliInfo(Form("ITS PID cut: limit, large, small: %.3ff %.3ff %.3ff", fMomentumLimit, fLargeCut, fSmallCut));
}
