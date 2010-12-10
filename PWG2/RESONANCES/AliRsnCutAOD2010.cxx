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
AliRsnCutAOD2010::AliRsnCutAOD2010(const char *name, Bool_t isMC) :
  AliRsnCut(name, AliRsnCut::kDaughter, 0.0, 0.0),

  fIsMC(kFALSE),
  fCheckITS(kTRUE),
  fCheckTPC(kTRUE),
  fCheckTOF(kTRUE),
  fUseGlobal(kTRUE),
  fUseITSSA(kTRUE),
  
  fMaxEta(1E6),
    
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
// Sets all parameters to the currently used values, and requires
// to know if we are running on data or MonteCarlo to set some others.
//

  fTPCparamDCA[0] = 0.0050;
  fTPCparamDCA[1] = 0.0070;
  fTPCparamDCA[2] = 1.0000;
  
  fITSparamDCA[0] = 0.0085;
  fITSparamDCA[1] = 0.0026;
  fITSparamDCA[2] = 1.5500;
  
  SetMC(isMC);
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
  
  fMaxEta(copy.fMaxEta),
  
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
}

//_________________________________________________________________________________________________
void AliRsnCutAOD2010::SetMC(Bool_t yn)
{
//
// Sets some parameters depending on MC or dataanalysis
//
  
  fIsMC = yn;
  
  if (fIsMC)
  {
    AliDebug(AliLog::kDebug + 2, "Setting for MC");
    fPID.GetTPCResponse().SetBetheBlochParameters(2.15898 / 50.0, 1.75295E1, 3.40030E-9, 1.96178, 3.91720);
  }
  else
  {
    AliDebug(AliLog::kDebug + 2, "Setting for DATA");
    fPID.GetTPCResponse().SetBetheBlochParameters(1.41543 / 50.0, 2.63394E1, 5.0411E-11, 2.12543, 4.88663);
  }
  
  AliITSPIDResponse itsrsp(fIsMC);
  fPID.GetITSResponse() = itsrsp;
}

//_________________________________________________________________________________________________
Bool_t AliRsnCutAOD2010::IsSelected(TObject *object)
{
//
// Cut checker.
//

  // coherence check: require an AOD track
  AliRsnDaughter *daughter = dynamic_cast<AliRsnDaughter*>(object);
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
  Bool_t isTPC, isTOF;
  if (!track->IsOn(AliESDtrack::kTPCin) && track->IsOn(AliESDtrack::kITSrefit) && (!track->IsOn(AliESDtrack::kITSpureSA)))
  {
    isTPC = kFALSE;
    isTOF = kFALSE;
    if (!fUseITSSA)
    {
      AliDebug(AliLog::kDebug + 2, "ITS standalone not used. Rejected");
      return kFALSE;
    }
  }
  else if (track->IsOn(AliESDtrack::kTPCin) && track->IsOn(AliESDtrack::kTPCrefit) && track->IsOn(AliESDtrack::kITSrefit))
  {
    isTPC = kTRUE;
    if (!fUseGlobal)
    {
      AliDebug(AliLog::kDebug + 2, "ITS standalone not used. Rejected");
      return kFALSE;
    }
    if (track->IsOn(AliESDtrack::kTOFout) && track->IsOn(AliESDtrack::kTIME))
      isTOF = kTRUE;
    else
      isTOF = kFALSE;
  }
  else
  {
    AliDebug(AliLog::kDebug + 2, "Track is not either a TPC track or a ITS standalone. Rejected");
    return kFALSE;
  }
  
  // step #2: check number of clusters
  if (isTPC)
  {
    if (track->GetTPCNcls() < fTPCminNclusters)
    {
      AliDebug(AliLog::kDebug + 2, "Too few clusters. Rejected");
      return kFALSE;
    }
  }
  else
  {
    if (track->GetITSNcls() < fITSminNclusters)
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
    if (vertex->GetType() == AliAODVertex::kKink)
    {
      AliDebug(AliLog::kDebug + 2, "Kink daughter. Rejected");
      return kFALSE;
    }
  }
  
  // step #5: DCA cut (transverse)
  Double_t dz[2], cov[3], sigmaDCA = 0.0, nsigma = 0.0;
  vertex = AliRsnTarget::GetCurrentEvent()->GetRefAOD()->GetPrimaryVertex();
  if (!vertex)
  {
    AliDebug(AliLog::kDebug + 2, "NULL vertex");
    return kFALSE;
  }
  if (!track->PropagateToDCA(vertex, AliRsnTarget::GetCurrentEvent()->GetRefAOD()->GetMagneticField(), kVeryBig, dz, cov))
  {
    AliDebug(AliLog::kDebug + 2, "Failed propagation to vertex");
    return kFALSE;
  }
  // compute the pt-dependent sigma
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
  // check the DCA
  if (dz[0] > nsigma * sigmaDCA)
  {
    AliDebug(AliLog::kDebug + 2, "Excceeded cut in DCA. Rejected");
    return kFALSE;
  }
  
  // step #6: check eta range
  if (TMath::Abs(track->Eta()) >= fMaxEta)
  {
    AliDebug(AliLog::kDebug + 2, "Outside ETA acceptance");
    return kFALSE;
  }
  
  // step #7: PID cuts
  if (isTPC)
  {
    if (fCheckTPC)
    {
      AliAODPid *pidObj    = track->GetDetPid();
      Double_t   mom       = pidObj->GetTPCmomentum();
      Double_t   nsigmaTPC = fPID.NumberOfSigmasTPC(track, fPIDtype);
      Double_t   bandTPC   = fTPChighBand;
      if (mom <= fTPClowLimit) bandTPC = fTPClowBand;
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
          Double_t nsigmaTOF = (Double_t)fPID.NumberOfSigmasTOF(track, fPIDtype);
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
        Double_t nsigmaTOF = (Double_t)fPID.NumberOfSigmasTOF(track, fPIDtype);
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
  else
  {
    if (fCheckITS)
    {
      if (nITS < 3 || !track->IsOn(AliESDtrack::kITSpid)) return kFALSE;
      Double_t nsigmaITS = (Double_t)fPID.NumberOfSigmasITS(track, fPIDtype);
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
      AliDebug(AliLog::kDebug + 2, "No PID checked. Accepted");
      return kTRUE;
    }
  }
}
