//
// Class AliRsnCutAOD2010
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

#include <TBits.h>

#include "AliESDtrack.h"
#include "AliAODTrack.h"

#include "AliRsnEvent.h"
#include "AliRsnDaughter.h"
#include "AliRsnCutAOD2010.h"

ClassImp(AliRsnCutAOD2010)

//_________________________________________________________________________________________________
AliRsnCutAOD2010::AliRsnCutAOD2010() :
  AliRsnCut(AliRsnCut::kDaughter),
  
  fIsMC(kFALSE),
  fCheckITS(kTRUE),
  fCheckTPC(kTRUE),
  fCheckTOF(kTRUE),
  fUseGlobal(kTRUE),
  fUseITSSA(kTRUE),
    
  fPIDtype(AliPID::kKaon),
  
  fTPCminNclusters(70),
  fTPCmaxChi2(4.0),
  fTPCmaxNSigmaDCA(7.0),
  fTPClowBand(5.0),
  fTPChighBand(3.0),
  fTPClowLimit(0.35),
  
  fITSminNclusters(4),
  fITSmaxChi2(2.5),
  fITSmaxNSigmaDCA(7.0),
  fITSband(3.0),
  
  fTOFlowLimit(-2.5),
  fTOFhighLimit(3.5),
  fPID()
{
//
// Default constructor.
//
  
  Int_t i = 0;
  for (i = 0; i < 3; i++) fTPCparamDCA[i] = fITSparamDCA[i] = 0.0;
  for (i = 0; i < 5; i++) fTPCparamBB[i] = 0.0;
}

//_________________________________________________________________________________________________
AliRsnCutAOD2010::AliRsnCutAOD2010
(const char *name) :
  AliRsnCut(name, AliRsnCut::kDaughter, 0.0, 0.0),

  fIsMC(kFALSE),
  fCheckITS(kTRUE),
  fCheckTPC(kTRUE),
  fCheckTOF(kTRUE),
  fUseGlobal(kTRUE),
  fUseITSSA(kTRUE),
    
  fPIDtype(AliPID::kKaon),
  
  fTPCminNclusters(70),
  fTPCmaxChi2(4.0),
  fTPCmaxNSigmaDCA(7.0),
  fTPClowBand(5.0),
  fTPChighBand(3.0),
  fTPClowLimit(0.35),
  
  fITSminNclusters(4),
  fITSmaxChi2(2.5),
  fITSmaxNSigmaDCA(7.0),
  fITSband(3.0),
  
  fTOFlowLimit(-2.5),
  fTOFhighLimit(3.5),
  fPID()
{
//
// Main constructor.
//

  Int_t i = 0;
  for (i = 0; i < 3; i++) fTPCparamDCA[i] = fITSparamDCA[i] = 0.0;
  for (i = 0; i < 5; i++) fTPCparamBB[i] = 0.0;
}

//_________________________________________________________________________________________________
AliRsnCutAOD2010::AliRsnCutAOD2010
(const AliRsnCutAOD2010& copy) :
  AliRsnCut(copy),
  
  fIsMC(copy.fIsMC),
  fCheckITS(copy.fCheckITS),
  fCheckTPC(copy.fCheckTPC),
  fCheckTOF(copy.fCheckTOF),
  fUseGlobal(copy.fUseGlobal),
  fUseITSSA(copy.fUseITSSA),
  
  fPIDtype(copy.fPIDtype),
  
  fTPCminNclusters(copy.fTPCminNclusters),
  fTPCmaxChi2(copy.fTPCmaxChi2),
  fTPCmaxNSigmaDCA(copy.fTPCmaxNSigmaDCA),
  fTPClowBand(copy.fTPClowBand),
  fTPChighBand(copy.fTPChighBand),
  fTPClowLimit(copy.fTPClowLimit),
  
  fITSminNclusters(copy.fITSminNclusters),
  fITSmaxChi2(copy.fITSmaxChi2),
  fITSmaxNSigmaDCA(copy.fITSmaxNSigmaDCA),
  fITSband(copy.fITSband),
  
  fTOFlowLimit(copy.fTOFlowLimit),
  fTOFhighLimit(copy.fTOFhighLimit),
  fPID(copy.fPID)
{
//
// Copy constructor
//

  Int_t i = 0;
  for (i = 0; i < 3; i++) 
  {
    fTPCparamDCA[i] = copy.fTPCparamDCA[i];
    fITSparamDCA[i] = copy.fITSparamDCA[i];
  }
  for (i = 0; i < 5; i++) fTPCparamBB[i] = copy.fTPCparamBB[i];
}

//_________________________________________________________________________________________________
Bool_t AliRsnCutAOD2010::IsSelected(TObject *obj1, TObject* /*obj2*/)
{
//
// Cut checker.
//

  // coherence check: require an AOD track
  AliRsnDaughter *daughter = dynamic_cast<AliRsnDaughter*>(obj1);
  if (!daughter) return kFALSE;
  AliAODTrack *track = daughter->GetRefAODtrack();
  if (!track) return kFALSE;
  
  // step #0: check presence of an SPD cluster
  Int_t nITS = 0, nSPD = 0;
  nSPD  = TESTBIT(track->GetITSClusterMap(), 0);
  nSPD += TESTBIT(track->GetITSClusterMap(), 1);
  nITS  = track->GetITSNcls() - nSPD;
  if (nSPD <= 0)
  {
    AliDebug(AliLog::kDebug + 2, "No SPD clusters in this track. Rejected");
    return kFALSE;
  }

  // step #1: check status flags and reject track if it does not match any possibility
  Bool_t  isTPC   = track->IsOn(AliESDtrack::kTPCin) && track->IsOn(AliESDtrack::kTPCrefit) && track->IsOn(AliESDtrack::kITSrefit);
  Bool_t  isITSSA = !track->IsOn(AliESDtrack::kTPCin) && track->IsOn(AliESDtrack::kITSrefit) && (!track->IsOn(AliESDtrack::kITSpureSA)) && track->IsOn(AliESDtrack::kITSpid);
  Bool_t  isTOF   = track->IsOn(AliESDtrack::kTOFout) && track->IsOn(AliESDtrack::kTIME);
  if (!isTPC && !isITSSA) 
  {
    AliDebug(AliLog::kDebug + 2, "Track is not either a TPC track or a ITS standalone. Rejected");
    return kFALSE;
  }
  else if (isTPC && !fUseGlobal)
  {
    AliDebug(AliLog::kDebug + 2, "Global tracks not used. Rejected");
    return kFALSE;
  }
  else if (isITSSA && !fUseITSSA)
  {
    AliDebug(AliLog::kDebug + 2, "ITS standalone not used. Rejected");
    return kFALSE;
  }
  
  // step #2: check number of clusters
  Int_t count = 0;
  if (isTPC)
  {
    count = track->GetTPCNcls();
    if (count < fTPCminNclusters)
    {
      AliDebug(AliLog::kDebug + 2, "Too few clusters. Rejected");
      return kFALSE;
    }
  }
  else // then is ITS-SA
  {
    count = track->GetITSNcls();
    if (count < fITSminNclusters)
    {
      AliDebug(AliLog::kDebug + 2, "Too few clusters. Rejected");
      return kFALSE;
    }
  }
  
  // step #3: check chi square
  if (isTPC)
  {
    if (track->Chi2perNDF() > fTPCmaxChi2)
    {
      AliDebug(AliLog::kDebug + 2, "Bad chi2. Rejected");
      return kFALSE;
    }
  }
  else
  {
    if (track->Chi2perNDF() > fITSmaxChi2)
    {
      AliDebug(AliLog::kDebug + 2, "Bad chi2. Rejected");
      return kFALSE;
    }
  }
  
  // step #4: reject kink daughters
  AliAODVertex *vertex = track->GetProdVertex();
  if (isTPC && vertex != 0x0)
  {
    if (vertex->GetType() == AliAODVertex::kKink && vertex->HasDaughter(track))
    {
      AliDebug(AliLog::kDebug + 2, "Kink daughter. Rejected");
      return kFALSE;
    }
  }
  
  // step #5: DCA cut (transverse)
  Double_t sigmaDCA = 0.0, nsigma = 0.0;
  if (isTPC)
  {
    sigmaDCA = fTPCparamDCA[0] + fTPCparamDCA[1] / TMath::Power(track->Pt(), fTPCparamDCA[2]);
    nsigma = fTPCmaxNSigmaDCA;
  }
  else 
  {
    sigmaDCA = fITSparamDCA[0] + fITSparamDCA[1] / TMath::Power(track->Pt(), fITSparamDCA[2]);
    nsigma = fITSmaxNSigmaDCA;
  }
  if (track->DCA() > nsigma * sigmaDCA)
  {
    AliDebug(AliLog::kDebug + 2, "Excceeded cut in DCA. Rejected");
    return kFALSE;
  }
  
  // step #6 PID cuts
  Double_t bandTPC   = 0.0;
  Double_t nsigmaTPC = 0.0;
  Double_t nsigmaITS = 0.0;
  Double_t nsigmaTOF = 0.0;
  if (isTPC)   nsigmaTPC = fPID.NumberOfSigmasTPC(track, fPIDtype);
  if (isITSSA) nsigmaITS = fPID.NumberOfSigmasITS(track, fPIDtype);
  if (isTOF)   nsigmaTOF = fPID.NumberOfSigmasTOF(track, fPIDtype);
  if (isITSSA && fCheckITS)
  {
    if (nITS < 3) return kFALSE;
    if (TMath::Abs(nsigmaITS) > fITSband)
    {
      AliDebug(AliLog::kDebug + 2, "Bad ITS PID. Rejected");
      return kFALSE;
    }
    else
    {
      AliDebug(AliLog::kDebug + 2, "Good ITS PID. Accepted");
      return kFALSE;
    }
  }
  else
  {
    if (fCheckTPC)
    {
      if (track->P() > fTPClowLimit) bandTPC = fTPChighBand; else bandTPC = fTPClowBand;
      if (TMath::Abs(nsigmaTPC) > bandTPC) 
      {
        AliDebug(AliLog::kDebug + 2, "Bad TPC PID. Rejected");
        return kFALSE;
      }
      else
      {
        AliDebug(AliLog::kDebug + 2, "Good TPC PID");
        if (fCheckTOF && isTOF)
        {
          if (nsigmaTOF < fTOFlowLimit || nsigmaTOF > fTOFhighLimit)
          {
            AliDebug(AliLog::kDebug + 2, "Bad TOF PID. Rejected");
            return kFALSE;
          }
          else
          {
            AliDebug(AliLog::kDebug + 2, "Good TOF PID. Accepted");
            return kFALSE;
          }
        }
        else
        {
          AliDebug(AliLog::kDebug + 2, "TOF not checked. Accepted");
          return kTRUE;
        }
      }
    }
    else
    {
      if (fCheckTOF && isTOF)
      {
        if (nsigmaTOF < fTOFlowLimit || nsigmaTOF > fTOFhighLimit)
        {
          AliDebug(AliLog::kDebug + 2, "Bad TOF PID. Rejected");
          return kFALSE;
        }
        else
        {
          AliDebug(AliLog::kDebug + 2, "Good TOF PID. Accepted");
          return kFALSE;
        }
      }
      else
      {
        AliDebug(AliLog::kDebug + 2, "No PID checked. Accepted");
        return kTRUE;
      }
    }
  }
  
  return kTRUE;
}
