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
   fDCARmaxfixed(kTRUE),
   fDCARminfixed(kTRUE),
   fDCARptFormula(""),
   fDCARptFormulaMin(""),
   fDCARmax(1E20),
   fDCARmin(0),
   fDCAZfixed(kTRUE),
   fDCAZptFormula(""),
   fDCAZmax(1E20),
   fDCA2D(kFALSE),
   fSPDminNClusters(0),
   fITSminNClusters(0),
   fITSmaxChi2(1E20),
   fTPCminNClusters(0),
   fTPCmaxChi2(1E20),
   fCutMaxChi2TPCConstrainedVsGlobal(1E20),
   fTrackMaxChi2(1E20),
   fIsUseCrossedRowsCut(kFALSE),
   fTPCminNCrossedRows(0),
   fTPCminCrossedRowsOverFindableCls(0),
   fIsUseLengthActiveVolumeTPCCut(kFALSE),
   fCutMinLengthActiveVolumeTPC(0),
   fAODTestFilterBit(-1),
   fCheckOnlyFilterBit(kTRUE),
   fESDtrackCuts(0x0)
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
   fDCARmaxfixed(copy.fDCARmaxfixed),
   fDCARminfixed(copy.fDCARminfixed),
   fDCARptFormula(copy.fDCARptFormula),
   fDCARptFormulaMin(copy.fDCARptFormulaMin),
   fDCARmax(copy.fDCARmax),
   fDCARmin(copy.fDCARmin),
   fDCAZfixed(copy.fDCAZfixed),
   fDCAZptFormula(copy.fDCAZptFormula),
   fDCAZmax(copy.fDCAZmax),
   fDCA2D(copy.fDCA2D),
   fSPDminNClusters(copy.fSPDminNClusters),
   fITSminNClusters(copy.fITSminNClusters),
   fITSmaxChi2(copy.fITSmaxChi2),
   fTPCminNClusters(copy.fTPCminNClusters),
   fTPCmaxChi2(copy.fTPCmaxChi2),
   fCutMaxChi2TPCConstrainedVsGlobal(copy.fCutMaxChi2TPCConstrainedVsGlobal),
   fTrackMaxChi2(copy.fTrackMaxChi2),
   fIsUseCrossedRowsCut(copy.fIsUseCrossedRowsCut),
   fTPCminNCrossedRows(copy.fTPCminNCrossedRows),
   fTPCminCrossedRowsOverFindableCls(copy.fTPCminCrossedRowsOverFindableCls),
   fIsUseLengthActiveVolumeTPCCut(copy.fIsUseLengthActiveVolumeTPCCut),
   fCutMinLengthActiveVolumeTPC(copy.fCutMinLengthActiveVolumeTPC),
   fAODTestFilterBit(copy.fAODTestFilterBit),
   fCheckOnlyFilterBit(copy.fCheckOnlyFilterBit),
   fESDtrackCuts(copy.fESDtrackCuts)
{
//
// Copy constructor.
// Just copy all data member values.
//

   SetPtRange(copy.fPt[0], copy.fPt[1]);
   SetEtaRange(copy.fEta[0], copy.fEta[1]);
}

//_________________________________________________________________________________________________
AliRsnCutTrackQuality &AliRsnCutTrackQuality::operator=(const AliRsnCutTrackQuality &copy)
{
//
// Assignment operator.
// Just copy all data member values.
//

   if (this == &copy)
      return *this;

   fFlagsOn = copy.fFlagsOn;
   fFlagsOff = copy.fFlagsOff;
   fRejectKinkDaughters = copy.fRejectKinkDaughters;
   fDCARmaxfixed = copy.fDCARmaxfixed;
   fDCARminfixed = copy.fDCARminfixed;
   fDCARptFormula = copy.fDCARptFormula;
   fDCARptFormulaMin = copy.fDCARptFormulaMin;
   fDCARmax = copy.fDCARmax;
   fDCARmin = copy.fDCARmin;
   fDCAZfixed = copy.fDCAZfixed;
   fDCAZptFormula = copy.fDCAZptFormula;
   fDCAZmax = copy.fDCAZmax;
   fDCA2D = copy.fDCA2D;
   fSPDminNClusters = copy.fSPDminNClusters;
   fITSminNClusters = copy.fITSminNClusters;
   fITSmaxChi2 = copy.fITSmaxChi2;
   fTPCminNClusters = copy.fTPCminNClusters;
   fTPCmaxChi2 = copy.fTPCmaxChi2;
   fCutMaxChi2TPCConstrainedVsGlobal = copy.fCutMaxChi2TPCConstrainedVsGlobal;
   fTrackMaxChi2 = copy.fTrackMaxChi2;
   fIsUseCrossedRowsCut=copy.fIsUseCrossedRowsCut;
   fTPCminNCrossedRows = copy.fTPCminNCrossedRows;
   fTPCminCrossedRowsOverFindableCls = copy.fTPCminCrossedRowsOverFindableCls;
   fIsUseLengthActiveVolumeTPCCut=copy.fIsUseLengthActiveVolumeTPCCut;
   fCutMinLengthActiveVolumeTPC = copy.fCutMinLengthActiveVolumeTPC;
   
   fAODTestFilterBit = copy.fAODTestFilterBit;
   fCheckOnlyFilterBit = copy.fCheckOnlyFilterBit;
   fESDtrackCuts = copy.fESDtrackCuts;
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
   fDCARmaxfixed = kTRUE;
   fDCARminfixed = kTRUE;
   fDCARptFormula = "";
   fDCARptFormulaMin = "";
   fDCARmax = 1E20;
   fDCARmin = 0;
   fDCAZfixed = kTRUE;
   fDCAZptFormula = "";
   fDCAZmax = 1E20;
   fDCA2D = kFALSE;
   fSPDminNClusters = 0;
   fITSminNClusters = 0;
   fITSmaxChi2 = 1E20;
   fTPCminNClusters = 0;
   fTPCmaxChi2 = 1E20;
   fAODTestFilterBit = -1;
   fCutMaxChi2TPCConstrainedVsGlobal = 1E20;
   fTrackMaxChi2 = 1E20;
   fIsUseCrossedRowsCut = 0;
   fTPCminNCrossedRows = 0;
   fTPCminCrossedRowsOverFindableCls = 0;
   fIsUseLengthActiveVolumeTPCCut = 0;
   fCutMinLengthActiveVolumeTPC = 0.0;
 
   if (fESDtrackCuts) {
      const char *cutsName = fESDtrackCuts->GetName();
      const char *cutsTitle = fESDtrackCuts->GetTitle();
      delete fESDtrackCuts;
      fESDtrackCuts = new AliESDtrackCuts(cutsName,cutsTitle);
   }
   SetPtRange(0.0, 1E20);
   SetEtaRange(-1E20, 1E20);
}

//_________________________________________________________________________________________________
void AliRsnCutTrackQuality::SetPtRange(Double_t a, Double_t b)
{
  //Set Pt range cut
  fPt[0] = TMath::Min(a, b); 
  fPt[1] = TMath::Max(a, b);
  if (fESDtrackCuts) fESDtrackCuts->SetPtRange(fPt[0], fPt[1]);
  return;
}

//_________________________________________________________________________________________________
void AliRsnCutTrackQuality::SetEtaRange(Double_t a, Double_t b)   
{
  //Set Pt range cut
  fEta[0] = TMath::Min(a, b);
  fEta[1] = TMath::Max(a, b);
  if (fESDtrackCuts) fESDtrackCuts->SetEtaRange(fEta[0], fEta[1]);
  return;
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
      if (fESDtrackCuts)
         return fESDtrackCuts->IsSelected(esdTrack);
      else
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
  //static AliESDtrackCuts cuts;
  AliESDtrackCuts cuts;

   // general acceptance/pt cuts
   cuts.SetPtRange(fPt[0], fPt[1]);
   cuts.SetEtaRange(fEta[0], fEta[1]);

   // transverse DCA cuts
   if (fDCARmaxfixed)
      cuts.SetMaxDCAToVertexXY(fDCARmax);
   else
      cuts.SetMaxDCAToVertexXYPtDep(fDCARptFormula.Data());
      
   if (fDCARminfixed)
      cuts.SetMinDCAToVertexXY(fDCARmin);
   else
      cuts.SetMinDCAToVertexXYPtDep(fDCARptFormulaMin.Data());

   // longitudinal DCA cuts
   if (fDCAZfixed)
      cuts.SetMaxDCAToVertexZ(fDCAZmax);
   else
      cuts.SetMaxDCAToVertexZPtDep(fDCAZptFormula.Data());

   // 2D DCA
   cuts.SetDCAToVertex2D(fDCA2D);

   // theis option is always disabled in current version
   cuts.SetRequireSigmaToVertex(kFALSE);

   // TPC related cuts for TPC+ITS tracks
   if (fIsUseCrossedRowsCut) {
     cuts.SetMinNCrossedRowsTPC(fTPCminNCrossedRows);
     cuts.SetMinRatioCrossedRowsOverFindableClustersTPC(fTPCminCrossedRowsOverFindableCls);
   } else {
     cuts.SetMinNClustersTPC(fTPCminNClusters);
   }
   cuts.SetMaxChi2PerClusterTPC(fTPCmaxChi2);
   cuts.SetAcceptKinkDaughters(!fRejectKinkDaughters);
   cuts.SetMaxChi2TPCConstrainedGlobal(fCutMaxChi2TPCConstrainedVsGlobal);

   if (fIsUseLengthActiveVolumeTPCCut)
     cuts.SetMinLengthActiveVolumeTPC(fCutMinLengthActiveVolumeTPC);

   // ITS related cuts for TPC+ITS tracks
   if (fSPDminNClusters > 0)
      cuts.SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
   cuts.SetMaxChi2PerClusterITS(fITSmaxChi2);

   // now that all is initialized, do the check
   if (!track) {
     AliError("Invalid track object. Rejected.");
     return kFALSE;
   }
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
      UInt_t bit = 1 << fAODTestFilterBit;
      AliDebugClass(2, Form("Required a test filter bit for AOD check: %u (result: %s)", bit, (track->TestFilterBit(bit) ? "accept" : "reject")));
      if (!track->TestFilterBit(bit))
         return kFALSE;
      else {
         if (track->Pt() < fPt[0] || track->Pt() > fPt[1]) return kFALSE;
         if (track->Eta() < fEta[0] || track->Eta() > fEta[1]) return kFALSE;
         if (fCheckOnlyFilterBit) return kTRUE;
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


   //step #1: check number of clusters 
   if ((!fIsUseCrossedRowsCut) && (track->GetTPCNcls() < fTPCminNClusters)) {
      AliDebug(AliLog::kDebug + 2, "Too few TPC clusters. Rejected");
      return kFALSE;
   }
 
   if (track->GetITSNcls() < fITSminNClusters) {
      AliDebug(AliLog::kDebug + 2, "Too few ITS clusters. Rejected");
      return kFALSE;
   }

   //check track chi square
   if (track->Chi2perNDF() > fTrackMaxChi2) {
      AliDebug(AliLog::kDebug + 2, "Bad chi2. Rejected");
      return kFALSE;
   }

   //step #2a: check number of crossed rows in TPC
   if (fIsUseCrossedRowsCut) {
     Float_t nCrossedRowsTPC = track->GetTPCNCrossedRows();
     if (nCrossedRowsTPC < fTPCminNCrossedRows) {
       AliDebug(AliLog::kDebug + 2, "Too few TPC crossed rows. Rejected");
       return kFALSE;
     }
     if (track->GetTPCNclsF()>0) {
       Float_t ratioCrossedRowsOverFindableClustersTPC = nCrossedRowsTPC / track->GetTPCNclsF();
       if (ratioCrossedRowsOverFindableClustersTPC < fTPCminCrossedRowsOverFindableCls){
	 AliDebug(AliLog::kDebug + 2, "Too few TPC crossed rows/findable clusters. Rejected");
	 return kFALSE;
       }
     } else {
       AliDebug(AliLog::kDebug + 2, "Negative value for TPC crossed rows/findable clusters. Rejected");
       return kFALSE;
     }
   }
   //step #2b: check on track length in active volume of TPC implemented only for ESD tracks
   //if (fIsUseLengthActiveVolumeTPCCut) { // not yet implemented in AODs}
 
   //step #3: reject kink daughters
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
   if (!fDCARmaxfixed) {
      TString str(fDCARptFormula);
      str.ReplaceAll("pt", "x");
      TFormula dcaXY(Form("%s_dcaXY", GetName()), str.Data());
      fDCARmax = dcaXY.Eval(track->Pt());
   }
   if (!fDCARminfixed) {   
      TString str2(fDCARptFormulaMin);
      str2.ReplaceAll("pt", "x");
      TFormula dcaXY_2(Form("%s_dcaXY_2", GetName()), str2.Data());
      fDCARmin = dcaXY_2.Eval(track->Pt());
   }
   // check the cut
   if (TMath::Abs(b[0]) > fDCARmax) {
      AliDebug(AliLog::kDebug + 2, "Too large transverse DCA");
      return kFALSE;
   }
   
   if (TMath::Abs(b[0]) < fDCARmin) {
      AliDebug(AliLog::kDebug + 2, "Too short transverse DCA");
      return kFALSE;
   }

   // step #5: DCA cut (longitudinal)
   // the DCA has already been computed above
   // if the DCA cut is not fixed, compute current value
   if (!fDCAZfixed) {
      TString str(fDCAZptFormula);
      str.ReplaceAll("pt", "x");
      TFormula dcaZ(Form("%s_dcaXY", GetName()), str.Data());
      fDCAZmax = dcaZ.Eval(track->Pt());
   }
   // check the cut
   if (TMath::Abs(b[1]) > fDCAZmax) {
      AliDebug(AliLog::kDebug + 2, "Too large longitudinal DCA");
      return kFALSE;
   }

   //NOTE: 2D DCA cut not implemented for AODs

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

   if(!fESDtrackCuts){
      AliInfo(Form("Ranges in eta, pt       : %.2f - %.2f, %.2f - %.2f", fEta[0], fEta[1], fPt[0], fPt[1]));
      AliInfo(Form("Kink daughters are      : %s", (fRejectKinkDaughters ? "rejected" : "accepted")));
      AliInfo(Form("TPC requirements (clusters)       : min. cluster = %i, max chi2 = %f", fTPCminNClusters, fTPCmaxChi2));
      AliInfo(Form("TPC requirements (crossed rows)   : min. crossed rows = %f, min. crossed rows/findable clusters = %f", fTPCminNCrossedRows, fTPCminCrossedRowsOverFindableCls));
      AliInfo(Form("TPC requirements (track length)   : min. track length in active volume TPC = %f", fCutMinLengthActiveVolumeTPC));

      AliInfo(Form("ITS requirements        : min. cluster = %d (all), %d (SPD), max chi2 = %f", fITSminNClusters, fSPDminNClusters, fITSmaxChi2));

      if (fDCARmaxfixed) {
	 AliInfo(Form("Max DCA r cut               : fixed to %f cm", fDCARmax));
      } else {
	 AliInfo(Form("Max DCA r cut formula       : %s", fDCARptFormula.Data()));
      }
   
      if (fDCARminfixed) {
         AliInfo(Form("Min DCA r cut               : fixed to %f cm", fDCARmin));
      } else {
         AliInfo(Form("Min DCA r cut formula       : %s", fDCARptFormulaMin.Data()));
      }

      if (fDCAZfixed) {
	 AliInfo(Form("DCA z cut               : fixed to %f cm", fDCAZmax));
      } else {
	 AliInfo(Form("DCA z cut formula       : %s", fDCAZptFormula.Data()));
      }

   }else{
      Float_t eta1,eta2,pt1,pt2;
      fESDtrackCuts->GetEtaRange(eta1,eta2);
      fESDtrackCuts->GetPtRange(pt1,pt2);
      AliInfo(Form("Ranges in eta, pt       : %.2f - %.2f, %.2f - %.2f", eta1, eta2, pt1, pt2));
      AliInfo(Form("Kink daughters are      : %s", (!(fESDtrackCuts->GetAcceptKinkDaughters()) ? "rejected" : "accepted")));
      AliInfo(Form("TPC requirements (clusters)       : min. cluster = %i, max chi2 = %f", fESDtrackCuts->GetMinNClusterTPC(), fESDtrackCuts->GetMaxChi2PerClusterTPC()));
      AliInfo(Form("TPC requirements (crossed rows)   : min. crossed rows = %f, min. crossed rows/findable clusters = %f", fESDtrackCuts->GetMinNCrossedRowsTPC(), fESDtrackCuts->GetMinRatioCrossedRowsOverFindableClustersTPC()));
      AliInfo(Form("TPC requirements (track length)   : min. track length in active volume TPC = %f", fESDtrackCuts->GetMinLengthActiveVolumeTPC()));

      Int_t i=fESDtrackCuts->GetClusterRequirementITS(AliESDtrackCuts::kSPD);
      Int_t n=-1;
      char s[500];
      if(i==AliESDtrackCuts::kOff) n=0;
      else if(i==AliESDtrackCuts::kAny) n=1;
      else if(i==AliESDtrackCuts::kBoth) n=2;
      else if(i==AliESDtrackCuts::kNone) sprintf(s,"no clusters");
      else if(i==AliESDtrackCuts::kFirst) sprintf(s,"cluster in first layer");
      else if(i==AliESDtrackCuts::kOnlyFirst) sprintf(s,"cluster only in first layer");
      else if(i==AliESDtrackCuts::kSecond) sprintf(s,"cluster in second layer");
      else if(i==AliESDtrackCuts::kOnlySecond) sprintf(s,"cluster only in second layer");
      if(n>-1) AliInfo(Form("ITS requirements        : min. cluster = %d (all), %d (SPD), max chi2 = %f", fESDtrackCuts->GetMinNClustersITS(), n, fESDtrackCuts->GetMaxChi2PerClusterITS()));
      else AliInfo(Form("ITS requirements        : min. cluster = %d (all), SPD: %s, max chi2 = %f", fESDtrackCuts->GetMinNClustersITS(), s, fESDtrackCuts->GetMaxChi2PerClusterITS()));

      sprintf(s,"%s",fESDtrackCuts->GetMaxDCAToVertexXYPtDep());
      if (!strcmp(s,"")) {
         AliInfo(Form("Max DCA r cut               : fixed to %f cm", fESDtrackCuts->GetMaxDCAToVertexXY()));
      } else {
         AliInfo(Form("Max DCA r cut formula       : %s", s));
      }
   
      sprintf(s,"%s",fESDtrackCuts->GetMinDCAToVertexXYPtDep());
      if (!strcmp(s,"")) {
         AliInfo(Form("Min DCA r cut               : fixed to %f cm", fESDtrackCuts->GetMinDCAToVertexXY()));
      } else {
         AliInfo(Form("Min DCA r cut formula       : %s", s));
      }

      sprintf(s,"%s",fESDtrackCuts->GetMaxDCAToVertexZPtDep());
      if (!strcmp(s,"")) {
	 AliInfo(Form("DCA z cut               : fixed to %f cm", fESDtrackCuts->GetMaxDCAToVertexZ()));
      } else {
         AliInfo(Form("DCA z cut formula       : %s", s));
      }

   }

   AliInfo(Form("fAODTestFilterBit       : filter bit %i",fAODTestFilterBit));
   AliInfo(Form("fCheckOnlyFilterBit     : %i",((int) fCheckOnlyFilterBit)));

   return;
}
//__________________________________________________________________________________________________
void AliRsnCutTrackQuality::SetDefaults2010(Bool_t useTPCCrossedRows, Bool_t useDefaultKinematicCuts)
{
//
// Default settings for cuts used in 2010
//
  
  fIsUseCrossedRowsCut=useTPCCrossedRows;
  fESDtrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE, fIsUseCrossedRowsCut);
  if (useDefaultKinematicCuts) {
    SetPtRange(0.15, 1E+20);
    SetEtaRange(-0.8, 0.8);
  } 
  SetAODTestFilterBit(5);
  return;
}

//__________________________________________________________________________________________________
void AliRsnCutTrackQuality::SetDefaultsHighPt2011(Bool_t useTPCCrossedRows, Bool_t useDefaultKinematicCuts)
{
//
// Default settings for cuts used in 2011 (for high-pT)
//
  fIsUseCrossedRowsCut=useTPCCrossedRows;
  fESDtrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kTRUE, fIsUseCrossedRowsCut);
  fESDtrackCuts->SetMinNCrossedRowsTPC(120); //default is min 70 crossed rows -> use 120 to go to higher pt
  fESDtrackCuts->SetMaxFractionSharedTPCClusters(0.4);//default is not set --> use to go to higher pt
  if (useDefaultKinematicCuts) {
    SetPtRange(0.15, 1E+20);
    SetEtaRange(-0.8, 0.8);
  } 
  SetAODTestFilterBit(10);
  return;
}

//__________________________________________________________________________________________________
void AliRsnCutTrackQuality::SetDefaults2011(Bool_t useTPCCrossedRows, Bool_t useDefaultKinematicCuts)
{
//
// Default std cuts 2011 with crossed rows (=70)
//
  fIsUseCrossedRowsCut=useTPCCrossedRows;
  fESDtrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kTRUE,fIsUseCrossedRowsCut);
  if (useDefaultKinematicCuts) {
    SetPtRange(0.15, 1E+20);
    SetEtaRange(-0.8, 0.8);
  } 
  SetAODTestFilterBit(5);
  return;
}

//__________________________________________________________________________________________________
void AliRsnCutTrackQuality::SetDefaultsTPCOnly(Bool_t useDefaultKinematicCuts)
{
//
// Default std cuts TPC-only
//
  fESDtrackCuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
  if (useDefaultKinematicCuts) {
    SetPtRange(0.15, 1E+20);
    SetEtaRange(-0.8, 0.8);
  } 
  SetAODTestFilterBit(-1);
  return;
}
//__________________________________________________________________________________________________
const char *AliRsnCutTrackQuality::Binary(UInt_t number)
{
//
// Convert an integer in binary
//

   static char b[50];
   b[0] = '\0';

   UInt_t z;
   for (z = 512; z > 0; z >>= 1)
      strncat(b, ((number & z) == z) ? "1" : "0", 1);

   return b;
}
