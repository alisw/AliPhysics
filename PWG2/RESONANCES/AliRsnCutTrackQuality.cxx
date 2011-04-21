//
// Class AliRsnCutTrackQuality
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
#include <TFormula.h>
#include <TBits.h>

#include "AliLog.h"
#include "AliESDtrackCuts.h"

#include "AliRsnEvent.h"
#include "AliRsnDaughter.h"
#include "AliRsnCutTrackQuality.h"

ClassImp(AliRsnCutTrackQuality)

//_________________________________________________________________________________________________
AliRsnCutTrackQuality::AliRsnCutTrackQuality(const char *name) :
   AliRsnCut(name, AliRsnTarget::kDaughter, 0.0, 0.0),
   fFlagsOn(0x0),
   fFlagsOff(0x0),
   fRejectKinkDaughters(kTRUE),
   fDCARfixed(kTRUE),
   fDCARptFormula(""),
   fDCARmax(1E20),
   fDCAZfixed(kTRUE),
   fDCAZptFormula(""),
   fDCAZmax(1E20),
   fSPDminNClusters(0),
   fITSminNClusters(0),
   fITSmaxChi2(1E20),
   fTPCminNClusters(0),
   fTPCmaxChi2(1E20),
   fAODTestFilterBit(-1)
{
//
// Default constructor.
// Initializes all cuts in such a way that all of them are disabled.
//

   SetPtRange(0.0, 1E20);
   SetEtaRange(-1E20, 1E20);
}

//_________________________________________________________________________________________________
AliRsnCutTrackQuality::AliRsnCutTrackQuality(const AliRsnCutTrackQuality &copy) :
   AliRsnCut(copy),
   fFlagsOn(copy.fFlagsOn),
   fFlagsOff(copy.fFlagsOff),
   fRejectKinkDaughters(copy.fRejectKinkDaughters),
   fDCARfixed(copy.fDCARfixed),
   fDCARptFormula(copy.fDCARptFormula),
   fDCARmax(copy.fDCARmax),
   fDCAZfixed(copy.fDCAZfixed),
   fDCAZptFormula(copy.fDCAZptFormula),
   fDCAZmax(copy.fDCAZmax),
   fSPDminNClusters(copy.fSPDminNClusters),
   fITSminNClusters(copy.fITSminNClusters),
   fITSmaxChi2(copy.fITSmaxChi2),
   fTPCminNClusters(copy.fTPCminNClusters),
   fTPCmaxChi2(copy.fTPCmaxChi2),
   fAODTestFilterBit(copy.fAODTestFilterBit)
{
//
// Copy constructor.
// Just copy all data member values.
//

   SetPtRange(copy.fPt[0], copy.fPt[1]);
   SetEtaRange(copy.fEta[0], copy.fEta[1]);
}

//_________________________________________________________________________________________________
AliRsnCutTrackQuality& AliRsnCutTrackQuality::operator=(const AliRsnCutTrackQuality &copy)
{
//
// Assignment operator.
// Just copy all data member values.
//


   fFlagsOn = copy.fFlagsOn;
   fFlagsOff = copy.fFlagsOff;
   fRejectKinkDaughters = copy.fRejectKinkDaughters;
   fDCARfixed = copy.fDCARfixed;
   fDCARptFormula = copy.fDCARptFormula;
   fDCARmax = copy.fDCARmax;
   fDCAZfixed = copy.fDCAZfixed;
   fDCAZptFormula = copy.fDCAZptFormula;
   fDCAZmax = copy.fDCAZmax;
   fSPDminNClusters = copy.fSPDminNClusters;
   fITSminNClusters = copy.fITSminNClusters;
   fITSmaxChi2 = copy.fITSmaxChi2;
   fTPCminNClusters = copy.fTPCminNClusters;
   fTPCmaxChi2 = copy.fTPCmaxChi2;
   fAODTestFilterBit = copy.fAODTestFilterBit;

   SetPtRange(copy.fPt[0], copy.fPt[1]);
   SetEtaRange(copy.fEta[0], copy.fEta[1]);

   return (*this);
}

//_________________________________________________________________________________________________
void AliRsnCutTrackQuality::DisableAll()
{
//
// Disable all cuts
//

   fFlagsOn = 0x0;
   fFlagsOff = 0x0;
   fRejectKinkDaughters = kFALSE;
   fDCARfixed = kTRUE;
   fDCARptFormula = "";
   fDCARmax = 1E20;
   fDCAZfixed = kTRUE;
   fDCAZptFormula = "";
   fDCAZmax = 1E20;
   fSPDminNClusters = 0;
   fITSminNClusters = 0;
   fITSmaxChi2 = 1E20;
   fTPCminNClusters = 0;
   fTPCmaxChi2 = 1E20;
   fAODTestFilterBit = -1;

   SetPtRange(0.0, 1E20);
   SetEtaRange(-1E20, 1E20);
}

//_________________________________________________________________________________________________
Bool_t AliRsnCutTrackQuality::IsSelected(TObject *object)
{
//
// Cut checker.
// Checks the type of object being evaluated
// and then calls the appropriate sub-function (for ESD or AOD)
//

   // coherence check
   if (!TargetOK(object)) return kFALSE;

   // status is checked in the same way for all tracks, using AliVTrack
   // as a convention, if a the collection of 'on' flags is '0x0', it
   // is assumed that no flags are required, and this check is skipped;
   // for the collection of 'off' flags this is not needed
   AliVTrack *vtrack = fDaughter->Ref2Vtrack();
   if (!vtrack) {
      AliDebug(AliLog::kDebug + 2, "This object is not either an ESD nor AOD track");
      return kFALSE;
   }
   ULong_t status   = (ULong_t)vtrack->GetStatus();
   ULong_t checkOn  = status & fFlagsOn;
   ULong_t checkOff = status & fFlagsOff;
   if (fFlagsOn != 0x0 && checkOn != fFlagsOn) {
      AliDebug(AliLog::kDebug + 2, Form("Failed flag check: required  %s", Binary(fFlagsOn)));
      AliDebug(AliLog::kDebug + 2, Form("                   track has %s", Binary(status  )));
      return kFALSE;
   }
   if (checkOff != 0) {
      AliDebug(AliLog::kDebug + 2, Form("Failed flag check: forbidden %s", Binary(fFlagsOff)));
      AliDebug(AliLog::kDebug + 2, Form("                   track has %s", Binary(status  )));
      return kFALSE;
   }
   AliDebug(AliLog::kDebug + 3, Form("Flag check OK: required  %s", Binary(fFlagsOn)));
   AliDebug(AliLog::kDebug + 3, Form("               forbidden %s", Binary(fFlagsOff)));
   AliDebug(AliLog::kDebug + 3, Form("               track has %s", Binary(status  )));

   // retrieve real object type
   AliESDtrack *esdTrack = fDaughter->Ref2ESDtrack();
   AliAODTrack *aodTrack = fDaughter->Ref2AODtrack();
   if (esdTrack) {
      AliDebug(AliLog::kDebug + 2, "Checking an ESD track");
      return CheckESD(esdTrack);
   } else if (aodTrack) {
      AliDebug(AliLog::kDebug + 2, "Checking an AOD track");
      return CheckAOD(aodTrack);
   } else {
      AliDebug(AliLog::kDebug + 2, Form("This object is not either an ESD nor AOD track, it is an %s", fDaughter->GetRef()->ClassName()));
      return kFALSE;
   }
}

//_________________________________________________________________________________________________
Bool_t AliRsnCutTrackQuality::CheckESD(AliESDtrack *track)
{
//
// Check an ESD track.
// This is done using the default track checker for ESD.
// It is declared static, not to recreate it every time.
//

   static AliESDtrackCuts cuts;

   // general acceptance/pt cuts
   cuts.SetPtRange(fPt[0], fPt[1]);
   cuts.SetEtaRange(fEta[0], fEta[1]);

   // transverse DCA cuts
   if (fDCARfixed)
      cuts.SetMaxDCAToVertexXY(fDCARmax);
   else
      cuts.SetMaxDCAToVertexXYPtDep(fDCARptFormula.Data());

   // longitudinal DCA cuts
   if (fDCAZfixed)
      cuts.SetMaxDCAToVertexZ(fDCAZmax);
   else
      cuts.SetMaxDCAToVertexZPtDep(fDCAZptFormula.Data());

   // these options are always disabled in current version
   cuts.SetDCAToVertex2D(kFALSE);
   cuts.SetRequireSigmaToVertex(kFALSE);

   // TPC related cuts for TPC+ITS tracks
   cuts.SetMinNClustersTPC(fTPCminNClusters);
   cuts.SetMaxChi2PerClusterTPC(fTPCmaxChi2);
   cuts.SetAcceptKinkDaughters(!fRejectKinkDaughters);

   // ITS related cuts for TPC+ITS tracks
   if (fSPDminNClusters > 0)
      cuts.SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);

   // now that all is initialized, do the check
   return cuts.IsSelected(track);
}

//_________________________________________________________________________________________________
Bool_t AliRsnCutTrackQuality::CheckAOD(AliAODTrack *track)
{
//
// Check an AOD track.
// This is done doing directly all checks, since there is not
// an equivalend checker for AOD tracks
//

   // if a test bit is used, check it and skip the following
   if (fAODTestFilterBit >= 0) {
      UInt_t bit = (UInt_t)fAODTestFilterBit;
      AliDebugClass(2, Form("Required a test filter bit for AOD check: %u (result: %s)", bit, (track->TestFilterBit(bit) ? "accept" : "reject")));
      if (!track->TestFilterBit(bit)) 
         return kFALSE;
      else {
         if (track->Pt() < fPt[0] || track->Pt() > fPt[1]) return kFALSE;
         if (track->Eta() < fEta[0] || track->Eta() > fEta[1]) return kFALSE;
         return kTRUE;
      }
   }

   // try to retrieve the reference AOD event
   AliAODEvent *aodEvent = 0x0;
   if (fEvent) aodEvent = fEvent->GetRefAOD();
   if (!aodEvent) {
      AliError("AOD reference event is not initialized!");
      return kFALSE;
   }

   // step #0: check SPD and ITS clusters
   Int_t nSPD = 0;
   nSPD  = TESTBIT(track->GetITSClusterMap(), 0);
   nSPD += TESTBIT(track->GetITSClusterMap(), 1);
   if (nSPD < fSPDminNClusters) {
      AliDebug(AliLog::kDebug + 2, "Not enough SPD clusters in this track. Rejected");
      return kFALSE;
   }

   // step #1: check number of clusters in TPC
   if (track->GetTPCNcls() < fTPCminNClusters) {
      AliDebug(AliLog::kDebug + 2, "Too few TPC clusters. Rejected");
      return kFALSE;
   }
   if (track->GetITSNcls() < fITSminNClusters) {
      AliDebug(AliLog::kDebug + 2, "Too few ITS clusters. Rejected");
      return kFALSE;
   }

   // step #2: check chi square
   if (track->Chi2perNDF() > fTPCmaxChi2) {
      AliDebug(AliLog::kDebug + 2, "Bad chi2. Rejected");
      return kFALSE;
   }
   if (track->Chi2perNDF() > fITSmaxChi2) {
      AliDebug(AliLog::kDebug + 2, "Bad chi2. Rejected");
      return kFALSE;
   }

   // step #3: reject kink daughters
   AliAODVertex *vertex = track->GetProdVertex();
   if (vertex && fRejectKinkDaughters) {
      if (vertex->GetType() == AliAODVertex::kKink) {
         AliDebug(AliLog::kDebug + 2, "Kink daughter. Rejected");
         return kFALSE;
      }
   }

   // step #4: DCA cut (transverse)
   // --> reject all tracks not ITS refitted
   Double_t b[2], cov[3];
   vertex = aodEvent->GetPrimaryVertex();
   if (!vertex) {
      AliDebug(AliLog::kDebug + 2, "NULL vertex");
      return kFALSE;
   }
   if ((track->GetStatus() & AliESDtrack::kITSrefit) == 0) {
      AliDebug(AliLog::kDebug + 2, "Not ITS refitted");
      return kFALSE;
   }
   if (!track->PropagateToDCA(vertex, aodEvent->GetMagneticField(), kVeryBig, b, cov)) {
      AliDebug(AliLog::kDebug + 2, "Failed propagation to vertex");
      return kFALSE;
   }
   // if the DCA cut is not fixed, compute current value
   if (!fDCARfixed) {
      static TString str(fDCARptFormula);
      str.ReplaceAll("pt", "x");
      static const TFormula dcaXY(Form("%s_dcaXY", GetName()), str.Data());
      fDCARmax = dcaXY.Eval(track->Pt());
   }
   // check the cut
   if (TMath::Abs(b[0]) > fDCARmax) {
      AliDebug(AliLog::kDebug + 2, "Too large transverse DCA");
      return kFALSE;
   }

   // step #5: DCA cut (longitudinal)
   // the DCA has already been computed above
   // if the DCA cut is not fixed, compute current value
   if (!fDCAZfixed) {
      static TString str(fDCAZptFormula);
      str.ReplaceAll("pt", "x");
      static const TFormula dcaZ(Form("%s_dcaXY", GetName()), str.Data());
      fDCAZmax = dcaZ.Eval(track->Pt());
   }
   // check the cut
   if (TMath::Abs(b[1]) > fDCAZmax) {
      AliDebug(AliLog::kDebug + 2, "Too large longitudinal DCA");
      return kFALSE;
   }

   // step #6: check eta/pt range
   if (track->Eta() < fEta[0] || track->Eta() > fEta[1]) {
      AliDebug(AliLog::kDebug + 2, "Outside ETA acceptance");
      return kFALSE;
   }
   if (track->Pt() < fPt[0] || track->Pt() > fPt[1]) {
      AliDebug(AliLog::kDebug + 2, "Outside PT acceptance");
      return kFALSE;
   }

   // if we are here, all cuts were passed and no exit point was got
   AliDebug(AliLog::kDebug + 2, "============================= ACCEPTED TRACK =====================================================");
   return kTRUE;
}

//_________________________________________________________________________________________________
void AliRsnCutTrackQuality::Print(const Option_t *) const
{
//
// Print information on this cut
//

   AliInfo(Form("Cut name                : %s", GetName()));
   AliInfo(Form("Required flags (off, on): %lx %lx", fFlagsOn, fFlagsOff));
   AliInfo(Form("Ranges in eta, pt       : %.2f - %.2f, %.2f - %.2f", fEta[0], fEta[1], fPt[0], fPt[1]));
   AliInfo(Form("Kink daughters are      : %s", (fRejectKinkDaughters ? "rejected" : "accepted")));
   AliInfo(Form("TPC requirements        : min. cluster = %d, max chi2 = %f", fTPCminNClusters, fTPCmaxChi2));
   AliInfo(Form("ITS requirements        : min. cluster = %d (all), %d (SPD), max chi2 = %f", fITSminNClusters, fSPDminNClusters, fITSmaxChi2));

   if (fDCARfixed) {
      AliInfo(Form("DCA r cut               : fixed to %f cm", fDCARmax));
   } else {
      AliInfo(Form("DCA r cut formula       : %s", fDCARptFormula.Data()));
   }

   if (fDCAZfixed) {
      AliInfo(Form("DCA z cut               : fixed to %f cm", fDCAZmax));
   } else {
      AliInfo(Form("DCA z cut formula       : %s", fDCAZptFormula.Data()));
   }
}
