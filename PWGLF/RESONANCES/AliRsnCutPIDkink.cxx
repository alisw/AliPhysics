//
// Class AliRsnCutPIDkink
//
// This class is based on the class AliRsnCutPIDNSigma with the aim
// of selecting kaons via their weak decays to pions and muons with
// the neutral daughter escaping detection (kink topology).

// Modifications from the original code by
// Martin Vala (martin.vala@cern.ch) and Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
// for the Kaon kink PID are made by Paraskevi Ganoti (Paraskevi.Ganoti@cern.ch)
//
#include "TF1.h"
#include "TMath.h"
#include "AliPIDResponse.h"
#include "AliESDpid.h"
#include "AliAODpidUtil.h"

#include "AliRsnCutPIDkink.h"
#include "AliESDkink.h"

ClassImp(AliRsnCutPIDkink)

//_________________________________________________________________________________________________
AliRsnCutPIDkink::AliRsnCutPIDkink() :
   AliRsnCut("cut", AliRsnTarget::kDaughter),
   fSpecies(AliPID::kUnknown),
   fDetector(kDetectors),
   fTrackNSigma(0.0),
   fTrackMom(0.0),
   fMyPID(0x0),
   fRanges("AliRsnPIDRange", 0)

{
//
  f1=new TF1("f1","((atan([0]*[1]*(1.0/(sqrt((x^2)*(1.0-([1]^2))-([0]^2)*([1]^2))))))*180.)/[2]",1.1,10.0);
  f1->SetParameter(0,0.493677);
  f1->SetParameter(1,0.9127037);
  f1->SetParameter(2,TMath::Pi());
  
  
  f2=new TF1("f2","((atan([0]*[1]*(1.0/(sqrt((x^2)*(1.0-([1]^2))-([0]^2)*([1]^2))))))*180.)/[2]",0.1,10.0);
  f2->SetParameter(0,0.13957018);
  f2->SetParameter(1,0.2731374);
  f2->SetParameter(2,TMath::Pi());
   // Main constructor.
//
}

//_________________________________________________________________________________________________
AliRsnCutPIDkink::AliRsnCutPIDkink
(const char *name, AliPID::EParticleType species, EDetector det) :
   AliRsnCut(name, AliRsnTarget::kDaughter),
   fSpecies(species),
   fDetector(det),
   fTrackNSigma(0.0),
   fTrackMom(0.0),
   fMyPID(0x0),
   fRanges("AliRsnPIDRange", 0)

{
//
  f1=new TF1("f1","((atan([0]*[1]*(1.0/(sqrt((x^2)*(1.0-([1]^2))-([0]^2)*([1]^2))))))*180.)/[2]",1.1,10.0);
  f1->SetParameter(0,0.493677);
  f1->SetParameter(1,0.9127037);
  f1->SetParameter(2,TMath::Pi());
  
  
  f2=new TF1("f2","((atan([0]*[1]*(1.0/(sqrt((x^2)*(1.0-([1]^2))-([0]^2)*([1]^2))))))*180.)/[2]",0.1,10.0);
  f2->SetParameter(0,0.13957018);
  f2->SetParameter(1,0.2731374);
  f2->SetParameter(2,TMath::Pi());
// Main constructor.
//
}

//_________________________________________________________________________________________________
AliRsnCutPIDkink::AliRsnCutPIDkink
(const AliRsnCutPIDkink &copy) :
   AliRsnCut(copy),
   fSpecies(copy.fSpecies),
   fDetector(copy.fDetector),
   fTrackNSigma(0.0),
   fTrackMom(0.0),
   fMyPID(copy.fMyPID),
   fRanges(copy.fRanges),
   f1(copy.f1),
   f2(copy.f2)

{
//
// Copy constructor.
//
}

//_________________________________________________________________________________________________
AliRsnCutPIDkink &AliRsnCutPIDkink::operator=(const AliRsnCutPIDkink &copy)
{
//
// Assignment operator
//

   AliRsnCut::operator=(copy);
   if (this == &copy)
      return *this;
   fSpecies = copy.fSpecies;
   fDetector = copy.fDetector;
   fMyPID = copy.fMyPID;
   fRanges = copy.fRanges;
   f1 = copy.f1;
   f1 = copy.f2;

   
   return (*this);
}

//__________________________________________________________________________________________________
void AliRsnCutPIDkink::InitMyPID(Bool_t isMC, Bool_t isESD)
{
//
// Initialize manual PID object
//

   if (isESD)
      fMyPID = new AliESDpid(isMC);
   else
      fMyPID = new AliAODpidUtil(isMC);

}

//_________________________________________________________________________________________________
Bool_t AliRsnCutPIDkink::IsSelected(TObject *object)
{
//
// Cut checker.
// As usual, there are 'kFALSE' exit points whenever one of the conditions is not passed,
// and at the end, it returns kTRUE since it bypassed all possible exit points.
//

   // coherence check
   if (!TargetOK(object)) return kFALSE;

   // check initialization of PID object
   // if manual PID is used, use that, otherwise get from source event
   AliPIDResponse *pid = 0x0;
   if (fMyPID)
      pid = fMyPID;
   else
      pid = fEvent->GetPIDResponse();
   if (!pid) {
      AliFatal("NULL PID response");
      return kFALSE;
   }

   // convert input object into AliVTrack
   // if this fails, the cut cannot be checked
   AliVTrack *vtrack = fDaughter->Ref2Vtrack();
   if (!vtrack) {
      AliDebugClass(2, "Referenced daughter is not a track");
      return kFALSE;
   }

   // get reference momentum
   fTrackMom = (fDetector == kTPC) ? vtrack->GetTPCmomentum() : vtrack->P();

   Int_t indexKinkPos=vtrack->GetKinkIndex(0);
   if(indexKinkPos<0) {
     AliESDEvent *esd=fEvent->GetRefESD();

     AliESDkink *kink=esd->GetKink(TMath::Abs(indexKinkPos)-1);

     const TVector3 vposKink(kink->GetPosition()); //reco position
     
     Double_t etracK= TMath::Sqrt(vtrack->P()*vtrack->P() + 0.493677 *0.493677); //assuming Kaon
     Double_t rapiditK = 0.5 * (TMath::Log((etracK + vtrack->Pz()) / (etracK - vtrack->Pz())));
     if((TMath::Abs(rapiditK)) > 0.5) return kFALSE;

     Float_t kinkAngle=TMath::RadToDeg()*kink->GetAngle(2);
     if( (kinkAngle<2.)) return kFALSE;
     Int_t tpcNClHigh = -31.67+ (11./12.)*(kink->GetR());
     Int_t tpcNClMin  = -85.5 + (65./95.)*(kink->GetR()) ;
     if (vtrack->GetTPCNcls() > tpcNClHigh) return kFALSE;
     if (vtrack->GetTPCNcls() < tpcNClMin) return kFALSE;

     Double_t maxDecAngKmu=f1->Eval(vtrack->P(),0.,0.,0.);
     Double_t maxDecAngpimu=f2->Eval(vtrack->P(),0.,0.,0.);
     if ((kinkAngle<maxDecAngpimu*1.2)) return kFALSE;
     if ((kinkAngle>maxDecAngKmu*.98) && (vtrack->P()>1.2)) return kFALSE;
     if(kink->GetR()<120.) return kFALSE;
     if(kink->GetR()>210.) return kFALSE;

     Float_t qT=kink->GetQt();
     if ((qT<0.12)||(qT>0.30)) return kFALSE;

     if (TMath::Abs(vposKink[2]) > 225.) return kFALSE ;
     if (TMath::Abs(vposKink[2]) < 0.5) return kFALSE ;

   
   // get number of sigmas
   switch (fDetector) {
      case kITS:
         fTrackNSigma = TMath::Abs(pid->NumberOfSigmasITS(vtrack, fSpecies));
         break;
      case kTPC:
         fTrackNSigma = TMath::Abs(pid->NumberOfSigmasTPC(vtrack, fSpecies));
         break;
      case kTOF:
         fTrackNSigma = TMath::Abs(pid->NumberOfSigmasTOF(vtrack, fSpecies));
         break;
      default:
         AliError("Bad detector chosen. Rejecting track");
         return kFALSE;
   }

   // loop on all ranges, and use the one which contains this momentum
   // if none is found, the cut is not passed
   Bool_t accept = kFALSE;
   Int_t  i, goodRange = -1, nRanges = fRanges.GetEntriesFast();
   for (i = 0; i < nRanges; i++) {
      AliRsnPIDRange *range = (AliRsnPIDRange *)fRanges[i];
      if (!range) continue;
      if (!range->IsInRange(fTrackMom)) continue;
      else {
	goodRange = i;
         accept = range->CutPass(fTrackNSigma);
         AliDebugClass(2, Form("[%s] NSigma = %.3f, max = %.3f, track %s", GetName(), fTrackNSigma, range->NSigmaCut(), (accept ? "accepted" : "rejected")));
         break;
      }
   }
   if (goodRange < 0) {
     AliDebugClass(2, Form("[%s] No good range found. Rejecting track", GetName()));
      return kFALSE;
   } else {
      AliDebugClass(2, Form("[%s] Mom = %.3f, good range found (#%d), track was %s", GetName(), fTrackMom, goodRange, (accept ? "accepted" : "rejected")));
      return accept;
   }
   }

   return kFALSE;
}

//_________________________________________________________________________________________________
void AliRsnCutPIDkink::Print(const Option_t *) const
{
//
// Print information on this cut
//

   Char_t mom[200], det[100], match[200];

   switch (fDetector) {
      case kITS: snprintf(det, 3, "ITS"); break;
      case kTPC: snprintf(det, 3, "TPC"); break;
      case kTOF: snprintf(det, 3, "TOF"); break;
      default  : snprintf(det, 3, "undefined");
   }

   AliInfo(Form("Cut name          : %s", GetName()));
   AliInfo(Form("--> PID detector  : %s", det));
   AliInfo(Form("--> match criteria: %s", match));
   AliInfo(Form("--> momentum range: %s", mom));
}

//_________________________________________________________________________________________________
void AliRsnCutPIDkink::AddPIDRange(Double_t nsigma, Double_t pmin, Double_t pmax)
{
//
// Add a new slot for checking PID
//

   Int_t n = fRanges.GetEntries();

   new (fRanges[n]) AliRsnPIDRange(nsigma, pmin, pmax);
}

//_________________________________________________________________________________________________
void AliRsnCutPIDkink::SinglePIDRange(Double_t nsigma)
{
//
// Clear all slots and sets a unique one
//

   fRanges.Delete();

   new (fRanges[0]) AliRsnPIDRange(nsigma, 0.0, 1E20);
}
