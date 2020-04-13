/**************************************************************************
 * Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

/////////////////////////////////////////////////////////////
//
// Class for cuts on AOD reconstructed D+->K0S+pi
//
// Author: J.Hamon, julien.hamon@cern.ch (IPHC)
//
/////////////////////////////////////////////////////////////

#include <Riostream.h>
#include <TDatabasePDG.h>
#include <TMath.h>

#include "AliAnalysisManager.h"
#include "AliAODPidHF.h"
#include "AliAODRecoCascadeHF.h"
#include "AliAODTrack.h"
#include "AliAODv0.h"
#include "AliAODVertex.h"
#include "AliESDtrack.h"
#include "AliESDv0.h"
#include "AliESDVertex.h"
#include "AliInputEventHandler.h"
#include "AliRDHFCutsDplustoK0spi.h"



/// \cond CLASSIMP
ClassImp(AliRDHFCutsDplustoK0spi);
/// \endcond




//--------------------------------------------------------------------------
AliRDHFCutsDplustoK0spi::AliRDHFCutsDplustoK0spi(const char* name) :
AliRDHFCuts(name)
   , fExcludedCut(-1)
   , fV0Type(0)
   , fV0daughtersCuts(0)
{
   //
   /// Default Constructor
   //

   const Int_t nvars = 16;
   SetNVars(nvars);
   TString varNames[nvars] = { "Inv. Mass D+ -> K0s+pi [GeV/c2]",                //  0
                               "Inv. Mass K0s [GeV/c2]",                         //  1 -------- K0S reconstruction --------
                               "K0s reco: min pT V0-daughters [GeV/c]",          //  2
                               "K0s reco: max DCA V0 (prong-to-prong) [sigma]",  //  3
                               "K0s reco: min d0(xy) V0-daughters [cm]",         //  4
                               "K0s reco: min cosPointing (3D) V0",              //  5
                               "K0s reco: min fiducial radius (xy) [cm]",        //  6
                               "min pT bachelor track [GeV/c]",                  //  7 ------ Cascade reconstruction ------
                               "min pT V0 track [GeV/c]",                        //  8
                               "max d0(xy) bachelor wrt PV [cm]",                //  9
                               "max d0(xy) V0 wrt PV [cm]",                      // 10
                               "max cosThetaStar D+ daughters",                  // 11
                               "min cosPointing (3D) cascade daughters",         // 12
                               "min |cosPointing| (xy) cascade daughters",       // 13
                               "min decay length (xy) cascade [cm]",             // 14
                               "V0 type"                                         // 15
                             };


   Bool_t isUpperCut[nvars] = { kTRUE,     //  0
                                kTRUE,     //  1 -------- K0S reconstruction --------
                                kFALSE,    //  2
                                kTRUE,     //  3
                                kFALSE,    //  4
                                kFALSE,    //  5
                                kFALSE,    //  6
                                kFALSE,    //  7 ------ Cascade reconstruction ------
                                kFALSE,    //  8
                                kTRUE,     //  9
                                kTRUE,     // 10
                                kTRUE,     // 11
                                kFALSE,    // 12
                                kFALSE,    // 13
                                kFALSE,    // 14
                                kFALSE     // 15
   };
   SetVarNames(nvars, varNames, isUpperCut);


   Float_t limits[2] = {0, 999999999.};
   SetPtBins(2, limits);


   // PID settings
   Double_t nsigma[5] = {3., 3., 3., 3., 0.}; // 0-2 for TPC, 3 for TOF, 4 for ITS

   if (fPidHF) delete fPidHF;
   fPidHF = new AliAODPidHF();

   fPidHF -> SetMatch(1);          // switch to combine the info from more detectors: 1 = || , 2 = &, 3 = p region
   fPidHF -> SetAsym(kFALSE);      // asimmetric PID required (different sigmas for different p bins)
   fPidHF -> SetSigma(nsigma);     // sigma for the raw signal PID: 0-2 for TPC, 3 for TOF, 4 for ITS
   fPidHF -> SetCompat(kTRUE);     // compatibility region : useful only if fMatch=1
   fPidHF -> SetTPC(1);
   fPidHF -> SetTOF(1);
   fPidHF -> SetITS(0);
   fPidHF -> SetTRD(0);
}




//--------------------------------------------------------------------------
AliRDHFCutsDplustoK0spi::AliRDHFCutsDplustoK0spi(const AliRDHFCutsDplustoK0spi& source) :
AliRDHFCuts(source)
   , fExcludedCut(source.fExcludedCut)
   , fV0Type(source.fV0Type)
   , fV0daughtersCuts(0)
{
   //
   /// Standard constructor
   //

   if (source.fV0daughtersCuts) AddTrackCutsV0daughters(source.fV0daughtersCuts);
   else fV0daughtersCuts = new AliESDtrackCuts();
}




//--------------------------------------------------------------------------
AliRDHFCutsDplustoK0spi &AliRDHFCutsDplustoK0spi::operator=(const AliRDHFCutsDplustoK0spi &source)
{
   //
   /// Assignment operator
   //

   if (this != &source) {

      AliRDHFCuts::operator = (source);
      fV0Type      = source.fV0Type;
      fExcludedCut = source.fExcludedCut;

      delete fV0daughtersCuts;
      fV0daughtersCuts = new AliESDtrackCuts(*(source.fV0daughtersCuts));

   }

   return *this;
}




//--------------------------------------------------------------------------
AliRDHFCutsDplustoK0spi::~AliRDHFCutsDplustoK0spi()
{
   //
   ///  Default destructor
   //

   if (fV0daughtersCuts) { delete fV0daughtersCuts; fV0daughtersCuts=0; }
}




//---------------------------------------------------------------------------
void AliRDHFCutsDplustoK0spi::GetCutVarsForOpt(AliAODRecoDecayHF* obj, Float_t* vars, Int_t nvars, Int_t *pdgdaughters, AliAODEvent *aod)
{
   //
   /// Fills in the array 'vars' the selection cut values
   //

   if(nvars!=fnVarsForOpt) {
      AliDebug(2, "wrong number of variables");
      return;
   }


   Double_t mDplusPDG = TDatabasePDG::Instance()->GetParticle(411)->Mass();
   Double_t mk0sPDG   = TDatabasePDG::Instance()->GetParticle(310)->Mass();


   AliAODRecoCascadeHF *d = (AliAODRecoCascadeHF*) obj;
   if (!d) {
      AliDebug(2, "AliAODRecoCascadeHF null");
      return;
   }

   // Get the v0 and all daughter tracks
   AliAODTrack *bachelorTrack = (AliAODTrack*) d->GetBachelor();
   AliAODv0    *v0            = (AliAODv0*)    d->Getv0();
   if (!v0 || !bachelorTrack) {
      AliDebug(2, "Missing V0 or missing bachelor for the current cascade");
      return;
   }

   AliAODTrack *v0positiveTrack = (AliAODTrack*) d->Getv0PositiveTrack();
   AliAODTrack *v0negativeTrack = (AliAODTrack*) d->Getv0NegativeTrack();
   if (!v0positiveTrack || !v0negativeTrack) {
      AliDebug(2, "No V0 daughters' objects");
      return;
   }


   //recalculate vertex w/o daughters
   Bool_t cleanvtx = kFALSE;
   AliAODVertex *origownvtx = 0x0;
   if (fRemoveDaughtersFromPrimary) {
      if (d->GetOwnPrimaryVtx()) origownvtx = new AliAODVertex(*d->GetOwnPrimaryVtx());
      cleanvtx = kTRUE;
      if (!RecalcOwnPrimaryVtx(d, aod)) {
         CleanOwnPrimaryVtx(d, aod, origownvtx);
         cleanvtx = kFALSE;
      }
   }


   Int_t iter=-1;

   if (fVarsForOpt[0]) {
      iter++;
      if (TMath::Abs(pdgdaughters[0]==411)) vars[iter] = TMath::Abs(d->InvMassDplustoK0spi()-mDplusPDG);
   }


   //___ K0s reconstruction  _________________________
   //_________________________________________________

   // - Cut on the K0s invariant mass
   if (fVarsForOpt[1]) {
      iter++;
      vars[iter] = TMath::Abs(v0->MassK0Short()-mk0sPDG);
   }

   // - Cut on V0-daughters pT
   if (fVarsForOpt[2]) {
      iter++;
      vars[iter] = v0positiveTrack->Pt();
   }

   // - Cut on V0-daughters DCA (prong-to-prong)
   if (fVarsForOpt[3]) {
      iter++;
      vars[iter] = TMath::Abs(v0->GetDCA());
   }

   // - Cut on V0-daughters transverse impact parameter d0(xy)
   if (fVarsForOpt[4]) {
      iter++;
      vars[iter] = TMath::Abs(v0->DcaPosToPrimVertex());
   }

   // - Cut on V0 cosine of pointing angle (3D) wrt PV
   if (fVarsForOpt[5]) {
      iter++;
      vars[iter] = d->CosV0PointingAngle();
   }

   // - Cut on V0 fiducial transverse radius
   if (fVarsForOpt[6]) {
      iter++;
      vars[iter] = TMath::Abs(d->DecayLengthXYV0());
   }


   //___ Cascade reconstruction  _____________________
   //_________________________________________________

   // - Cut on bachelor pT
   if (fVarsForOpt[7]) {
      iter++;
      vars[iter] = bachelorTrack->Pt();
   }

   // - Cut on V0 pT
   if (fVarsForOpt[8]) {
      iter++;
      vars[iter] = v0->Pt();
   }

   // - Cut on bachelor transverse impact parameter d0(xy)
   if (fVarsForOpt[9]) {
      iter++;
      vars[iter] = TMath::Abs(d->Getd0Prong(0));
   }

   // - Cut on V0 transverse impact parameter d0(xy)
   if (fVarsForOpt[10]) {
      iter++;
      vars[iter] = TMath::Abs(d->Getd0Prong(1));
   }

   // - Cut on cascade-daughters cosThetaStar
   if (fVarsForOpt[11]) {
      iter++;
      vars[iter] = TMath::Abs(d->CosThetaStar(0, 411, 211, 310));
   }

   // - Cut on cascade cosPointingAngle (3D) wrt PV
   if (fVarsForOpt[12]) {
      iter++;
      vars[iter] = d->CosPointingAngle();
   }

   // - Cut on cascade cosPointingAngle (xy) wrt PV
   if (fVarsForOpt[13]) {
      iter++;
      vars[iter] = TMath::Abs(d->CosPointingAngleXY());
   }

   // - Cut on cascade decay length (xy)
   if (fVarsForOpt[14]) {
      iter++;
      vars[iter] = TMath::Abs(d->DecayLengthXY());
   }


   if (cleanvtx) CleanOwnPrimaryVtx(d, aod, origownvtx);

   return;
}


//--------------------------------------------------------------------------
Bool_t AliRDHFCutsDplustoK0spi::PreSelect(TObject* obj, AliAODv0 *v0, AliVTrack *bachelorTrack){
  //
  // Apply pre-selections, used in the AOD filtering
  //
  
  if (!fCutsRD) {
    AliFatal("Cut matrix not inizialized. Exit...");
    return 0;
  }

  AliAODRecoCascadeHF* d = (AliAODRecoCascadeHF*)obj;
  if (!d) {
    AliDebug(2,"AliAODRecoCascadeHF null");
    return 0;
  }
  
  Double_t pt = d->Pt();
  Int_t ptbin = PtBin(pt);
  if (ptbin<0) return 0;
  
  if ( v0 && ((v0->GetOnFlyStatus() == kTRUE  && GetV0Type() == AliRDHFCuts::kOnlyOfflineV0s) ||
	      (v0->GetOnFlyStatus() == kFALSE && GetV0Type() == AliRDHFCuts::kOnlyOnTheFlyV0s)) ) return 0;

  // cut on V0 pT min
  if (v0->Pt() < fCutsRD[GetGlobalIndex(8,ptbin)]) return 0;
  
  // cut on the minimum pt of the bachelor
  if (bachelorTrack->Pt() < fCutsRD[GetGlobalIndex(7,ptbin)]) return 0;


  // Cut on the K0s invariant mass
  Double_t  mk0s   = v0->MassK0Short();
  Double_t mk0sPDG   = TDatabasePDG::Instance()->GetParticle(310)->Mass();
  if (TMath::Abs(mk0s-mk0sPDG) > fCutsRD[GetGlobalIndex(1,ptbin)]) return 0;
    
  // cut on the D+ invariant mass
  Double_t  mDplus = d->InvMassDplustoK0spi();
  Double_t mDplusPDG = TDatabasePDG::Instance()->GetParticle(411)->Mass();
  if (TMath::Abs(mDplus - mDplusPDG) > fCutsRD[GetGlobalIndex(0,ptbin)]) return 0;

  // - Cut on V0-daughters DCA (prong-to-prong)
  if (TMath::Abs(v0->GetDCA()) > fCutsRD[GetGlobalIndex(3,ptbin)]) return 0;
  
  return kTRUE;
}
//--------------------------------------------------------------------------
Int_t AliRDHFCutsDplustoK0spi::IsSelected(TObject* obj, Int_t selectionLevel, AliAODEvent* aod)
{
   //
   /// Apply selection cuts
   /// Return value:  0: wrong candidates
   ///                1: D+->K0s+pi
   //

   if (!fCutsRD) {
      AliFatal("Cut matrice not inizialized. Exit...");
      return 0;
   }


   AliAODRecoCascadeHF *d = (AliAODRecoCascadeHF*) obj;
   if (!d) {
      AliDebug(2, "AliAODRecoCascadeHF null");
      return 0;
   }

   if (!d->GetSecondaryVtx()) {
      AliDebug(2, "No secondary vertex for cascade");
      return 0;
   }

   if (d->GetNDaughters() != 2) {
      AliDebug(2, Form("No 2 daughters for current cascade (nDaughters=%d)", d->GetNDaughters()));
      return 0;
   }



   //___ 1./ Check daughter tracks  __________________
   //_________________________________________________
   AliAODv0    *v0            = (AliAODv0*)    d->Getv0();
   AliAODTrack *bachelorTrack = (AliAODTrack*) d->GetBachelor();


   if (!v0 || !bachelorTrack) {
      AliDebug(2, "Missing V0 or missing bachelor for the current cascade");
      return 0;
   }

   if ((v0->GetOnFlyStatus()==kTRUE  && GetV0Type()==AliRDHFCuts::kOnlyOfflineV0s) ||
       (v0->GetOnFlyStatus()==kFALSE && GetV0Type()==AliRDHFCuts::kOnlyOnTheFlyV0s)) {
      AliDebug(2, "Wrong V0 status (offline or on-the-fly)");
      return 0;
   }

   if (bachelorTrack->GetID()<0) {
      AliDebug(2, Form("Bachelor has negative ID %d", bachelorTrack->GetID()));
      return 0;
   }

   if (!v0->GetSecondaryVtx()) {
      AliDebug(2, "No secondary vertex for V0 by cascade");
      return 0;
   }

   if (v0->GetNDaughters()!=2) {
      AliDebug(2, Form("More than 2 daughters for V0 of current cascade (onTheFly=%d, nDaughters=%d)", v0->GetOnFlyStatus(), v0->GetNDaughters()));
      return 0;
   }


   // - Get the V0 daughter tracks
   AliAODTrack *v0positiveTrack = (AliAODTrack*) d->Getv0PositiveTrack();
   AliAODTrack *v0negativeTrack = (AliAODTrack*) d->Getv0NegativeTrack();

   if (!v0positiveTrack || !v0negativeTrack) {
      AliDebug(2, "No V0 daughters' objects");
      return 0;
   }

   if (v0positiveTrack->GetID()<0 || v0negativeTrack->GetID()<0) {
      AliDebug(2, Form("At least one of V0 daughters has negative ID %d %d", v0positiveTrack->GetID(), v0negativeTrack->GetID()));
      return 0;
   }

   if (fUseTrackSelectionWithFilterBits && d->HasBadDaughters()) {
      AliDebug(2, "Check on the bachelor FilterBit: no BIT(4). Candidate rejected.");
      return 0;
   }


   // - Selection on daughter tracks
   if ( selectionLevel == AliRDHFCuts::kAll ||
        selectionLevel == AliRDHFCuts::kTracks )
   { if (!AreDtoK0sDaughtersSelected(d)) return 0; }




   //___ 2./ Check for cascade D+  ___________________
   //_________________________________________________

   // - Selection on candidate
   if (selectionLevel == AliRDHFCuts::kAll ||
       selectionLevel == AliRDHFCuts::kCandidate) {


      Double_t mDplusPDG = TDatabasePDG::Instance()->GetParticle(411)->Mass();
      Double_t mk0sPDG   = TDatabasePDG::Instance()->GetParticle(310)->Mass();

      Double_t  mDplus = d->InvMassDplustoK0spi();
      Double_t  mk0s   = v0->MassK0Short();


      Double_t pt = d->Pt();
      Int_t ptbin = PtBin(pt);
      if (ptbin==-1) {
         return 0;
      }


      // Recalculate vertex w/o daughters
      AliAODVertex *origownvtx = 0x0;
      if (fRemoveDaughtersFromPrimary && !fUseMCVertex) {
         if (d->GetOwnPrimaryVtx()) origownvtx = new AliAODVertex(*d->GetOwnPrimaryVtx());
         if (!RecalcOwnPrimaryVtx(d, aod)) {
            CleanOwnPrimaryVtx(d, aod, origownvtx);
            return 0;
         }
      }


      if (fUseMCVertex) {
         if (d->GetOwnPrimaryVtx()) origownvtx = new AliAODVertex(*d->GetOwnPrimaryVtx());
         if (!SetMCPrimaryVtx(d, aod)) {
            CleanOwnPrimaryVtx(d, aod, origownvtx);
            return 0;
         }
      }



      //___ Invariant mass compatibility  _______________
      //_________________________________________________

      // - Check invariant mass of cascade D+
      if (TMath::Abs(mDplus - mDplusPDG) > fCutsRD[GetGlobalIndex(0,ptbin)] && fExcludedCut!=0) {
         AliDebug(4, "cascade not compatible with D+ invariant mass");
         return 0;
      }



      //___ K0s reconstruction  _________________________
      //_________________________________________________

      // - Cut on the K0s invariant mass
      if (TMath::Abs(mk0s-mk0sPDG) > fCutsRD[GetGlobalIndex(1,ptbin)] && fExcludedCut!=1) {
         AliDebug(4, Form(" V0 mass is %2.2e and does not correspond to K0S cut", mk0s));
         CleanOwnPrimaryVtx(d, aod, origownvtx);
         return 0;
      }


      // - Cut on V0-daughters pT
      if (v0positiveTrack->Pt() < fCutsRD[GetGlobalIndex(2,ptbin)] && fExcludedCut!=2) {
         AliDebug(4, Form(" v0 positive track Pt=%2.2e < %2.2e", v0positiveTrack->Pt(), fCutsRD[GetGlobalIndex(2,ptbin)]));
         CleanOwnPrimaryVtx(d, aod, origownvtx);
         return 0;
      }
      if (v0negativeTrack->Pt() < fCutsRD[GetGlobalIndex(2,ptbin)] && fExcludedCut!=2) {
         AliDebug(4, Form(" v0 positive track Pt=%2.2e < %2.2e", v0negativeTrack->Pt(), fCutsRD[GetGlobalIndex(2,ptbin)]));
         CleanOwnPrimaryVtx(d, aod, origownvtx);
         return 0;
      }


      // - Cut on V0-daughters DCA (prong-to-prong)
      if (TMath::Abs(v0->GetDCA()) > fCutsRD[GetGlobalIndex(3,ptbin)] && fExcludedCut!=3) {
         AliDebug(4, " V0-daugters DCA doesn't pass the cut");
         CleanOwnPrimaryVtx(d, aod, origownvtx);
         return 0;
      }


      // - Cut on V0-daughters transverse impact parameter d0(xy)
      if (TMath::Abs(v0->DcaPosToPrimVertex()) < fCutsRD[GetGlobalIndex(4,ptbin)] && fExcludedCut!=4) {
         AliDebug(4, " V0-daugters impact parameter (xy) don't pass the cut");
         CleanOwnPrimaryVtx(d, aod, origownvtx);
         return 0;
      }
      if (TMath::Abs(v0->DcaNegToPrimVertex()) < fCutsRD[GetGlobalIndex(4,ptbin)] && fExcludedCut!=4) {
         AliDebug(4, " V0-daugters impact parameter (xy) don't pass the cut");
         CleanOwnPrimaryVtx(d, aod, origownvtx);
         return 0;
      }


      // - Cut on V0 cosine of pointing angle (3D) wrt PV
      if (d->CosV0PointingAngle() < fCutsRD[GetGlobalIndex(5,ptbin)] && fExcludedCut!=5) {
         AliDebug(4, " V0 cosine of pointing angle (3D) doesn't pass the cut");
         CleanOwnPrimaryVtx(d, aod, origownvtx);
         return 0;
      }


      // - Cut on V0 fiducial transverse radius
      if (TMath::Abs(d->DecayLengthXYV0()) < fCutsRD[GetGlobalIndex(6,ptbin)] && fExcludedCut!=6) {
         AliDebug(4, " V0 fiducial transverse radius min doesn't pass the cut");
         CleanOwnPrimaryVtx(d, aod, origownvtx);
         return 0;
      }



      //___ Cascade reconstruction  _____________________
      //_________________________________________________

      // - Cut on bachelor pT
      if (bachelorTrack->Pt() < fCutsRD[GetGlobalIndex(7,ptbin)] && fExcludedCut!=7) {
         AliDebug(4, Form(" bachelor track Pt=%2.2e < %2.2e", bachelorTrack->Pt(), fCutsRD[GetGlobalIndex(7,ptbin)]));
         CleanOwnPrimaryVtx(d, aod, origownvtx);
         return 0;
      }


      // - Cut on V0 pT
      if (v0->Pt() < fCutsRD[GetGlobalIndex(8,ptbin)] && fExcludedCut!=8) {
         AliDebug(4, Form(" V0 track Pt=%2.2e < %2.2e", v0->Pt(), fCutsRD[GetGlobalIndex(8,ptbin)]));
         CleanOwnPrimaryVtx(d, aod, origownvtx);
         return 0;
      }


      // - Cut on bachelor transverse impact parameter d0(xy)
      if (TMath::Abs(d->Getd0Prong(0)) > fCutsRD[GetGlobalIndex(9,ptbin)] && fExcludedCut!=9) {
         AliDebug(4, " bachelor transverse impact parameter doesn't pass the cut");
         CleanOwnPrimaryVtx(d, aod, origownvtx);
         return 0;
      }


      // - Cut on V0 transverse impact parameter d0(xy)
      if (TMath::Abs(d->Getd0Prong(1)) > fCutsRD[GetGlobalIndex(10,ptbin)] && fExcludedCut!=10) {
         AliDebug(4, " V0 transverse impact parameter doesn't pass the cut");
         CleanOwnPrimaryVtx(d, aod, origownvtx);
         return 0;
      }


      // - Cut on D+-daughters cosThetaStar
      if (TMath::Abs(d->CosThetaStar(0, 411, 211, 310)) > fCutsRD[GetGlobalIndex(11,ptbin)] && fExcludedCut!=11) {
         AliDebug(4, " D+ daughter cosThetaStar doesn't pass the cut");
         CleanOwnPrimaryVtx(d, aod, origownvtx);
         return 0;
      }


      // - Cut on cascade cosPointingAngle (3D) wrt PV
      if (d->CosPointingAngle() < fCutsRD[GetGlobalIndex(12,ptbin)] && fExcludedCut!=12) {
         AliDebug(4, " Cascade cosPointingAngle (3D) doesn't pass the cut");
         CleanOwnPrimaryVtx(d, aod, origownvtx);
         return 0;
      }

      // - Cut on cascade cosPointingAngle (xy) wrt PV
      if (TMath::Abs(d->CosPointingAngleXY()) < fCutsRD[GetGlobalIndex(13,ptbin)] && fExcludedCut!=13) {
         AliDebug(4, " Cascade cosPointingAngle (xy) doesn't pass the cut");
         CleanOwnPrimaryVtx(d, aod, origownvtx);
         return 0;
      }


      // - Cut on cascade decay length (xy)
      if (TMath::Abs(d->DecayLengthXY()) < fCutsRD[GetGlobalIndex(14,ptbin)] && fExcludedCut!=14) {
         AliDebug(4, " D+ normalised decay length (XY) doesn't pass the cut");
         CleanOwnPrimaryVtx(d, aod, origownvtx);
         return 0;
      }


      // Unset recalculated primary vertex when not needed any more
      CleanOwnPrimaryVtx(d, aod, origownvtx);

   } // end if (kAll || kCandidate)




   //___ 3./ PID selection  __________________________
   //_________________________________________________
   Int_t returnvalueTopo = 1;
   Int_t returnvaluePID  = 1;

   if (selectionLevel == AliRDHFCuts::kAll       ||
       selectionLevel == AliRDHFCuts::kCandidate ||
       selectionLevel == AliRDHFCuts::kPID)
   { returnvaluePID = IsSelectedPID(d); }


   Int_t returnvalueTot = returnvalueTopo & returnvaluePID;

   return returnvalueTot;
}




//---------------------------------------------------------------------------
Int_t AliRDHFCutsDplustoK0spi::IsSelectedPID(AliAODRecoDecayHF* obj)
{
   //
   /// Apply PID selections on the bachelor track (D+->K0s+pi)
   /// Return value:  0: not compatible with Pion
   ///                1: compatible with Pion
   //

   if (!fUsePID) {
      AliDebug(2, "PID selection inactive. Candidate accepted.");
      return 1;
   }

   AliAODRecoCascadeHF *d = (AliAODRecoCascadeHF*) obj;
   if (!d) {
      AliDebug(2, "AliAODRecoCascadeHF null");
      return 0;
   }

   AliAODTrack *bachelorTrack = (AliAODTrack*) d->GetBachelor();
   if (!bachelorTrack) {
      AliDebug(2, "Missing bachelor for the current cascade");
      return 0;
   }


   // - Combined TPC/TOF PID
   Int_t isPion = fPidHF -> MakeRawPid(bachelorTrack, AliPID::kPion);

   if (isPion>=0) return 1;
   else           return 0;
}




//---------------------------------------------------------------------------
Bool_t AliRDHFCutsDplustoK0spi::IsInFiducialAcceptance(Double_t pt, Double_t y) const
{
   //
   /// Checking if the cascade is in the fiducial acceptance region
   //

   if (fMaxRapidityCand > -998.) {
      if (TMath::Abs(y) > fMaxRapidityCand) return kFALSE;
      else return kTRUE;
   }

   if (pt > 5.) {
      // applying cut for pt > 5 GeV
      AliDebug(2,Form("pt of cascade = %f (> 5), cutting at |y| < 0.8", pt));
      if (TMath::Abs(y) > 0.8) return kFALSE;
   } else {
      // appliying smooth cut for pt < 5 GeV
      Double_t maxFiducialY = -0.2/15*pt*pt + 1.9/15*pt + 0.5;
      Double_t minFiducialY =  0.2/15*pt*pt - 1.9/15*pt - 0.5;
      AliDebug(2,Form("pt of cascade = %f (< 5), cutting  according to the fiducial zone [%f, %f]\n", pt, minFiducialY, maxFiducialY));
      if (y < minFiducialY || y > maxFiducialY) return kFALSE;
   }

   return kTRUE;
}




//---------------------------------------------------------------------------
Bool_t AliRDHFCutsDplustoK0spi::AreDtoK0sDaughtersSelected(AliAODRecoDecayHF *rd) const
{
   //
   /// Daughter track selections.
   //

   AliAODRecoCascadeHF* d = (AliAODRecoCascadeHF*) rd;
   if (!d) {
     AliDebug(2, "AliAODRecoCascadeHF null");
     return kFALSE;
   }

   if (!fTrackCuts) {
     AliFatal("Cut object is not defined for bachelor. Candidate accepted.");
     return kTRUE;
   }


   //___ 1./ Check bachelor tracks  __________________
   //_________________________________________________
   AliAODTrack *bachelorTrack = (AliAODTrack*) d->GetBachelor();
   if (!bachelorTrack) return kFALSE;


   if (fIsCandTrackSPDFirst && d->Pt()<fMaxPtCandTrackSPDFirst) {
      if(!bachelorTrack->HasPointOnITSLayer(0)) return kFALSE;
   }


   if (fKinkReject != (!(fTrackCuts->GetAcceptKinkDaughters())) ) {
      AliError(Form("Not compatible setting: fKinkReject=%1d - fTrackCuts->GetAcceptKinkDaughters()=%1d",fKinkReject, fTrackCuts->GetAcceptKinkDaughters()));
      return kFALSE;
   }


   AliAODVertex *vAOD = (AliAODVertex*) d->GetPrimaryVtx();
   Double_t pos[3]; vAOD->GetXYZ(pos);
   Double_t cov[6]; vAOD->GetCovarianceMatrix(cov);
   const AliESDVertex vESD(pos,cov,100.,100);

   if (!IsDaughterSelected(bachelorTrack, &vESD, fTrackCuts)) return kFALSE;

   if (!fV0daughtersCuts) {
     AliFatal("Cut object is not defined for V0daughters. Candidate accepted.");
     return kTRUE;
   }



   //___ 2./ Check V0  _______________________________
   //_________________________________________________
   AliAODv0 *v0 = (AliAODv0*) d->Getv0();
   if (!v0) return kFALSE;

   AliAODTrack *v0positiveTrack = (AliAODTrack*) d->Getv0PositiveTrack();
   AliAODTrack *v0negativeTrack = (AliAODTrack*) d->Getv0NegativeTrack();
   if (!v0positiveTrack || !v0negativeTrack) return kFALSE;


   Float_t etaMin=0, etaMax=0; fV0daughtersCuts->GetEtaRange(etaMin,etaMax);
   if ( (v0positiveTrack->Eta()<=etaMin || v0positiveTrack->Eta()>=etaMax) ||
        (v0negativeTrack->Eta()<=etaMin || v0negativeTrack->Eta()>=etaMax) ) return kFALSE;
   Float_t ptMin=0, ptMax=0; fV0daughtersCuts->GetPtRange(ptMin,ptMax);
   if ( (v0positiveTrack->Pt()<=ptMin || v0positiveTrack->Pt()>=ptMax) ||
        (v0negativeTrack->Pt()<=ptMin || v0negativeTrack->Pt()>=ptMax) ) return kFALSE;


   // Condition on nTPCclusters
   if (fV0daughtersCuts->GetMinNClusterTPC()>0) {
      if ( ((v0positiveTrack->GetTPCClusterInfo(2,1)) < fV0daughtersCuts->GetMinNClusterTPC()) ||
           ((v0negativeTrack->GetTPCClusterInfo(2,1)) < fV0daughtersCuts->GetMinNClusterTPC()) ) return kFALSE;
   }


   // kTPCrefit status
   if (v0->GetOnFlyStatus()==kFALSE) { // only for offline V0s
      if (fV0daughtersCuts->GetRequireTPCRefit()) {
         if( !(v0positiveTrack->GetStatus() & AliESDtrack::kTPCrefit)) return kFALSE;
         if( !(v0negativeTrack->GetStatus() & AliESDtrack::kTPCrefit)) return kFALSE;
      }
   }


   // kink condition
   if (!fV0daughtersCuts->GetAcceptKinkDaughters()) {
      AliAODVertex *maybeKinkPos = (AliAODVertex*)v0positiveTrack->GetProdVertex();
      AliAODVertex *maybeKinkNeg = (AliAODVertex*)v0negativeTrack->GetProdVertex();
      if (maybeKinkPos->GetType()==AliAODVertex::kKink ||
      maybeKinkNeg->GetType()==AliAODVertex::kKink) return kFALSE;
   }


  // Findable clusters > 0 condition - from V0 analysis
  // if( v0positiveTrack->GetTPCNclsF()<=0 || v0negativeTrack->GetTPCNclsF()<=0 ) return kFALSE;
  /*
      Float_t lPosTrackCrossedRows = v0positiveTrack->GetTPCClusterInfo(2,1);
      Float_t lNegTrackCrossedRows = v0positiveTrack->GetTPCClusterInfo(2,1);
      fTreeVariableLeastNbrCrossedRows = (Int_t) lPosTrackCrossedRows;
      if( lNegTrackCrossedRows < fTreeVariableLeastNbrCrossedRows )
      fTreeVariableLeastNbrCrossedRows = (Int_t) lNegTrackCrossedRows;
      //Compute ratio Crossed Rows / Findable clusters
      //Note: above test avoids division by zero!
      Float_t lPosTrackCrossedRowsOverFindable = lPosTrackCrossedRows / ((double)(pTrack->GetTPCNclsF()));
      Float_t lNegTrackCrossedRowsOverFindable = lNegTrackCrossedRows / ((double)(nTrack->GetTPCNclsF()));
      fTreeVariableLeastRatioCrossedRowsOverFindable = lPosTrackCrossedRowsOverFindable;
      if( lNegTrackCrossedRowsOverFindable < fTreeVariableLeastRatioCrossedRowsOverFindable )
      fTreeVariableLeastRatioCrossedRowsOverFindable = lNegTrackCrossedRowsOverFindable;
      //Lowest Cut Level for Ratio Crossed Rows / Findable = 0.8, set here
      if ( fTreeVariableLeastRatioCrossedRowsOverFindable < 0.8 ) return kFALSE;
   */

  return kTRUE;
}




//---------------------------------------------------------------------------
Int_t AliRDHFCutsDplustoK0spi::GetV0Type()
{
   //
   /// Get the V0 type: offline or on-the-fly
   //

   const Int_t nvars = this->GetNVars() ;
   fV0Type = (this->GetCuts())[nvars-1];
   TString *sVarNames = GetVarNames();

   if (sVarNames[nvars-1].Contains("V0 type")) {
      return (Int_t)fV0Type;
   } else {
      AliInfo("AliRDHFCutsDplustoK0spi Last variable is not the V0 type!!!");
      return -999;
   }
}
