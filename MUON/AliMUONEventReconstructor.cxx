/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

/*
$Log$
Revision 1.3  2000/06/16 07:27:08  gosset
To remove problem in running RuleChecker, like in MUON-dev

Revision 1.1.2.5  2000/06/16 07:00:26  gosset
To remove problem in running RuleChecker

Revision 1.1.2.4  2000/06/12 08:00:07  morsch
Dummy streamer to solve CINT compilation problem (to be investigated !)

Revision 1.1.2.3  2000/06/09 20:59:57  morsch
Make includes consistent with new file structure.

Revision 1.1.2.2  2000/06/09 12:58:05  gosset
Removed comment beginnings in Log sections of .cxx files
Suppressed most violations of coding rules

Revision 1.1.2.1  2000/06/07 14:44:53  gosset
Addition of files for track reconstruction in C++
*/

//__________________________________________________________________________
//
// MUON event reconstructor in ALICE
//
// This class contains as data:
// * the parameters for the event reconstruction
// * a pointer to the array of hits to be reconstructed (the event)
// * a pointer to the array of segments made with these hits inside each station
// * a pointer to the array of reconstructed tracks
//
// It contains as methods, among others:
// * MakeEventToBeReconstructed to build the array of hits to be reconstructed
// * MakeSegments to build the segments
// * MakeTracks to build the tracks
//__________________________________________________________________________

#include <iostream.h>

#include <TRandom.h>
#include <TFile.h>

#include "AliCallf77.h" 
#include "AliMUONEventReconstructor.h"
#include "AliMUON.h"
#include "AliMUONHitForRec.h"
#include "AliMUONSegment.h"
#include "AliMUONHit.h"
#include "AliMUONRawCluster.h"
#include "AliMUONTrack.h"
#include "AliMUONChamber.h"
#include "AliMUONTrackHit.h"
#include "AliRun.h"

#ifndef WIN32 
# define initfield initfield_
# define reco_gufld reco_gufld_
#else 
# define initfield INITFIELD
# define reco_gufld RECO_GUFLD
#endif 

extern "C"
{
void type_of_call initfield();
void type_of_call reco_gufld(Double_t *Coor, Double_t *Field);
}

//************* Defaults parameters for reconstruction
static const Double_t kDefaultMinBendingMomentum = 3.0;
static const Double_t kDefaultMaxSigma2Distance = 16.0;
static const Double_t kDefaultBendingResolution = 0.01;
static const Double_t kDefaultNonBendingResolution = 0.144;
static const Double_t kDefaultChamberThicknessInX0 = 0.03;
// Simple magnetic field:
// Value taken from macro MUONtracking.C: 0.7 T, hence 7 kG
// Length and Position from reco_muon.F, with opposite sign:
// Length = ZMAGEND-ZCOIL
// Position = (ZMAGEND+ZCOIL)/2
// to be ajusted differently from real magnetic field ????
static const Double_t kDefaultSimpleBValue = 7.0;
static const Double_t kDefaultSimpleBLength = 428.0;
static const Double_t kDefaultSimpleBPosition = 1019.0;
static const Int_t kDefaultRecGeantHits = 0;
static const Double_t kDefaultEfficiency = 0.95;

static const Int_t kDefaultPrintLevel = 0;

ClassImp(AliMUONEventReconstructor) // Class implementation in ROOT context

  //__________________________________________________________________________
AliMUONEventReconstructor::AliMUONEventReconstructor(void)
{
  // Constructor for class AliMUONEventReconstructor
  SetReconstructionParametersToDefaults();
  // Memory allocation for the TClonesArray of hits for reconstruction
  // Is 10000 the right size ????
  fHitsForRecPtr = new TClonesArray("AliMUONHitForRec", 10000);
  fNHitsForRec = 0; // really needed or GetEntriesFast sufficient ????
  // Memory allocation for the TClonesArray's of segments in stations
  // Is 2000 the right size ????
  for (Int_t st = 0; st < kMaxMuonTrackingStations; st++) {
    fSegmentsPtr[st] = new TClonesArray("AliMUONSegment", 2000);
    fNSegments[st] = 0; // really needed or GetEntriesFast sufficient ????
  }
  // Memory allocation for the TClonesArray of reconstructed tracks
  // Is 10 the right size ????
  fRecTracksPtr = new TClonesArray("AliMUONTrack", 10);
  fNRecTracks = 0; // really needed or GetEntriesFast sufficient ????

  // Initialize magnetic field
  // using Fortran subroutine INITFIELD in "reco_muon.F".
  // Should rather use AliMagF ???? and remove prototyping ...
  initfield();
  // Impression de quelques valeurs
  Double_t coor[3], field[3];
  coor[0] = 50.0;
  coor[1] = 50.0;
  coor[2] = 950.0;
  reco_gufld(coor, field);
  cout << "coor: " << coor[0] << ", " << coor[1] << ", " << coor[2] << endl;
  cout << "field: " << field[0] << ", " << field[1] << ", " << field[2] << endl;
  coor[2] = -950.0;
  reco_gufld(coor, field);
  cout << "coor: " << coor[0] << ", " << coor[1] << ", " << coor[2] << endl;
  cout << "field: " << field[0] << ", " << field[1] << ", " << field[2] << endl;
  coor[2] = -950.0;

  if (fPrintLevel >= 0) {
    cout << "AliMUONEventReconstructor constructed with defaults" << endl; Dump();}
  return;
}

AliMUONEventReconstructor::AliMUONEventReconstructor (const AliMUONEventReconstructor& Reconstructor)
{
  // Dummy copy constructor
}

AliMUONEventReconstructor & AliMUONEventReconstructor::operator=(const AliMUONEventReconstructor& Reconstructor)
{
  // Dummy assignment operator
    return *this;
}

  //__________________________________________________________________________
AliMUONEventReconstructor::~AliMUONEventReconstructor(void)
{
  // Destructor for class AliMUONEventReconstructor
  delete fHitsForRecPtr; // Correct destruction of everything ???? or delete [] ????
  for (Int_t st = 0; st < kMaxMuonTrackingStations; st++)
    delete fSegmentsPtr[st]; // Correct destruction of everything ????
  return;
}

  //__________________________________________________________________________
void AliMUONEventReconstructor::SetReconstructionParametersToDefaults(void)
{
  // Set reconstruction parameters to default values
  // Would be much more convenient with a structure (or class) ????
  fMinBendingMomentum = kDefaultMinBendingMomentum;
  fMaxSigma2Distance = kDefaultMaxSigma2Distance;

  AliMUON *pMUON = (AliMUON*) gAlice->GetModule("MUON");
  // ******** Parameters for making HitsForRec
  // minimum radius,
  // like in TRACKF_STAT:
  // 2 degrees for stations 1 and 2, or ch(0...) from 0 to 3;
  // 30 cm for stations 3 to 5, or ch(0...) from 4 to 9
  for (Int_t ch = 0; ch < kMaxMuonTrackingChambers; ch++) {
    if (ch < 4) fRMin[ch] = TMath::Abs((&(pMUON->Chamber(ch)))->Z()) *
		  2.0 * TMath::Pi() / 180.0;
    else fRMin[ch] = 30.0;
  }
  // maximum radius
  // like in TRACKF_STAT (10 degrees ????)
  fRMax[0] = fRMax[1] = 91.5;
  fRMax[2] = fRMax[3] = 122.5;
  fRMax[4] = fRMax[5] = 158.3;
  fRMax[6] = fRMax[7] = 260.0;
  fRMax[8] = fRMax[9] = 260.0;

  // ******** Parameters for making segments
  // should be parametrized ????
  // according to interval between chambers in a station ????
  // Maximum distance in non bending plane
  // 5 * 0.22 just to remember the way it was made in TRACKF_STAT
  // SIGCUT*DYMAX(IZ)
  for (Int_t st = 0; st < kMaxMuonTrackingStations; st++)
    fSegmentMaxDistNonBending[st] = 5. * 0.22;
  // Maximum distance in bending plane
  // values from TRACKF_STAT corresponding to (J psi 20cm)
  fSegmentMaxDistBending[0] = 1.5;
  fSegmentMaxDistBending[1] = 1.5;
  fSegmentMaxDistBending[2] = 3.0;
  fSegmentMaxDistBending[3] = 6.0;
  fSegmentMaxDistBending[4] = 6.0;
  
  fBendingResolution = kDefaultBendingResolution;
  fNonBendingResolution = kDefaultNonBendingResolution;
  fChamberThicknessInX0 = kDefaultChamberThicknessInX0;
  fSimpleBValue = kDefaultSimpleBValue;
  fSimpleBLength = kDefaultSimpleBLength;
  fSimpleBPosition = kDefaultSimpleBPosition;
  fRecGeantHits = kDefaultRecGeantHits;
  fEfficiency = kDefaultEfficiency;
  fPrintLevel = kDefaultPrintLevel;
  return;
}

//__________________________________________________________________________
Double_t AliMUONEventReconstructor::GetImpactParamFromBendingMomentum(Double_t BendingMomentum)
{
  // Returns impact parameter at vertex in bending plane (cm),
  // from the signed bending momentum "BendingMomentum" in bending plane (GeV/c),
  // using simple values for dipole magnetic field.
  // The sign is the sign of the charge.
  return (-0.0003 * fSimpleBValue * fSimpleBLength * fSimpleBPosition /
	  BendingMomentum);
}

//__________________________________________________________________________
Double_t AliMUONEventReconstructor::GetBendingMomentumFromImpactParam(Double_t ImpactParam)
{
  // Returns signed bending momentum in bending plane (GeV/c),
  // from the impact parameter "ImpactParam" at vertex in bending plane (cm),
  // using simple values for dipole magnetic field.
  // The sign is the sign of the charge.
  return (-0.0003 * fSimpleBValue * fSimpleBLength * fSimpleBPosition /
	  ImpactParam);
}

//__________________________________________________________________________
void AliMUONEventReconstructor::SetBkgGeantFile(Text_t *BkgGeantFileName)
{
  // Set background file ... for GEANT hits
  // Must be called after having loaded the firts signal event
  if (fPrintLevel >= 0) {
    cout << "Enter SetBkgGeantFile with BkgGeantFileName ``"
	 << BkgGeantFileName << "''" << endl;}
  if (strlen(BkgGeantFileName)) {
    // BkgGeantFileName not empty: try to open the file
    if (fPrintLevel >= 2) {cout << "Before File(Bkg)" << endl; gDirectory->Dump();}
    fBkgGeantFile = new TFile(BkgGeantFileName);
    if (fPrintLevel >= 2) {cout << "After File(Bkg)" << endl; gDirectory->Dump();}
    if (fBkgGeantFile-> IsOpen()) {
      if (fPrintLevel >= 0) {
	cout << "Background for GEANT hits in file: ``" << BkgGeantFileName
	     << "'' successfully opened" << endl;}
    }
    else {
      cout << "Background for GEANT hits in file: " << BkgGeantFileName << endl;
      cout << "NOT FOUND: EXIT" << endl;
      exit(0); // right instruction for exit ????
    }
    // Arrays for "particles" and "hits"
    fBkgGeantParticles = new TClonesArray("TParticle", 200);
    fBkgGeantHits = new TClonesArray("AliMUONHit", 2000);
    // Event number to -1 for initialization
    fBkgGeantEventNumber = -1;
    // Back to the signal file:
    // first signal event must have been loaded previously,
    // otherwise, Segmentation violation at the next instruction
    // How is it possible to do smething better ????
    ((gAlice->TreeK())->GetCurrentFile())->cd();
    if (fPrintLevel >= 2) {cout << "After cd(gAlice)" << endl; gDirectory->Dump();}
  }
  return;
}

//__________________________________________________________________________
void AliMUONEventReconstructor::NextBkgGeantEvent(void)
{
  // Get next event in background file for GEANT hits
  // Goes back to event number 0 when end of file is reached
  char treeName[20];
  TBranch *branch;
  if (fPrintLevel >= 0) {
    cout << "Enter NextBkgGeantEvent" << endl;}
  // Clean previous event
  if(fBkgGeantTK) delete fBkgGeantTK;
  fBkgGeantTK = NULL;
  if(fBkgGeantParticles) fBkgGeantParticles->Clear();
  if(fBkgGeantTH) delete fBkgGeantTH;
  fBkgGeantTH = NULL;
  if(fBkgGeantHits) fBkgGeantHits->Clear();
  // Increment event number
  fBkgGeantEventNumber++;
  // Get access to Particles and Hits for event from background file
  if (fPrintLevel >= 2) {cout << "Before cd(Bkg)" << endl; gDirectory->Dump();}
  fBkgGeantFile->cd();
  if (fPrintLevel >= 2) {cout << "After cd(Bkg)" << endl; gDirectory->Dump();}
  // Particles: TreeK for event and branch "Particles"
  sprintf(treeName, "TreeK%d", fBkgGeantEventNumber);
  fBkgGeantTK = (TTree*)gDirectory->Get(treeName);
  if (!fBkgGeantTK) {
    if (fPrintLevel >= 0) {
      cout << "Cannot find Kine Tree for background event: " <<
	fBkgGeantEventNumber << endl;
      cout << "Goes back to event 0" << endl;
    }
    fBkgGeantEventNumber = 0;
    sprintf(treeName, "TreeK%d", fBkgGeantEventNumber);
    fBkgGeantTK = (TTree*)gDirectory->Get(treeName);
    if (!fBkgGeantTK) {
      cout << "ERROR: cannot find Kine Tree for background event: " <<
	fBkgGeantEventNumber << endl;
      exit(0);
    }
  }
  if (fBkgGeantTK) 
    fBkgGeantTK->SetBranchAddress("Particles", &fBkgGeantParticles);
  fBkgGeantTK->GetEvent(0); // why event 0 ???? necessary ????
  // Hits: TreeH for event and branch "MUON"
  sprintf(treeName, "TreeH%d", fBkgGeantEventNumber);
  fBkgGeantTH = (TTree*)gDirectory->Get(treeName);
  if (!fBkgGeantTH) {
    cout << "ERROR: cannot find Hits Tree for background event: " <<
      fBkgGeantEventNumber << endl;
      exit(0);
  }
  if (fBkgGeantTH && fBkgGeantHits) {
    branch = fBkgGeantTH->GetBranch("MUON");
    if (branch) branch->SetAddress(&fBkgGeantHits);
  }
  fBkgGeantTH->GetEntries(); // necessary ????
  // Back to the signal file
  ((gAlice->TreeK())->GetCurrentFile())->cd();
  if (fPrintLevel >= 2) {cout << "After cd(gAlice)" << endl; gDirectory->Dump();}
  return;
}

//__________________________________________________________________________
void AliMUONEventReconstructor::EventReconstruct(void)
{
  // To reconstruct one event
  if (fPrintLevel >= 1) cout << "enter EventReconstruct" << endl;
  MakeEventToBeReconstructed();
  MakeSegments();
  MakeTracks();
  return;
}

  //__________________________________________________________________________
void AliMUONEventReconstructor::ResetHitsForRec(void)
{
  // To reset the array and the number of HitsForRec,
  // and also the number of HitsForRec
  // and the index of the first HitForRec per chamber
  if (fHitsForRecPtr) fHitsForRecPtr->Clear();
  fNHitsForRec = 0;
  for (Int_t ch = 0; ch < kMaxMuonTrackingChambers; ch++)
    fNHitsForRecPerChamber[ch] = fIndexOfFirstHitForRecPerChamber[ch] = 0;
  return;
}

  //__________________________________________________________________________
void AliMUONEventReconstructor::ResetSegments(void)
{
  // To reset the TClonesArray of segments and the number of Segments
  // for all stations
  for (Int_t st = 0; st < kMaxMuonTrackingStations; st++) {
    if (fSegmentsPtr[st]) fSegmentsPtr[st]->Clear();
    fNSegments[st] = 0;
  }
  return;
}

  //__________________________________________________________________________
void AliMUONEventReconstructor::ResetTracks(void)
{
  // To reset the TClonesArray of reconstructed tracks
  if (fRecTracksPtr) fRecTracksPtr->Clear();
  fNRecTracks = 0;
  return;
}

  //__________________________________________________________________________
void AliMUONEventReconstructor::MakeEventToBeReconstructed(void)
{
  // To make the list of hits to be reconstructed,
  // either from the GEANT hits or from the raw clusters
  // according to the parameter set for the reconstructor
  if (fPrintLevel >= 1) cout << "enter MakeEventToBeReconstructed" << endl;
  ResetHitsForRec();
  if (fRecGeantHits == 1) {
    // Reconstruction from GEANT hits
    // Back to the signal file
    ((gAlice->TreeK())->GetCurrentFile())->cd();
    // Signal hits
    // AliMUON *MUON  = (AliMUON*) gAlice->GetModule("MUON"); // necessary ????
    // Security on MUON ????
    AddHitsForRecFromGEANT(gAlice->TreeH());
    // Background hits
    AddHitsForRecFromBkgGEANT(fBkgGeantTH, fBkgGeantHits);
    // Sort HitsForRec in increasing order with respect to chamber number
    SortHitsForRecWithIncreasingChamber();
  }
  else {
    // Reconstruction from raw clusters
    // AliMUON *MUON  = (AliMUON*) gAlice->GetModule("MUON"); // necessary ????
    // Security on MUON ????
    // TreeR assumed to be be "prepared" in calling function
    // by "MUON->GetTreeR(nev)" ????
    TTree *treeR = gAlice->TreeR();
    AddHitsForRecFromRawClusters(treeR);
    // No sorting: it is done automatically in the previous function
  }
  if (fPrintLevel >= 10) {
    cout << "end of MakeEventToBeReconstructed" << endl;
    cout << "NHitsForRec: " << fNHitsForRec << endl;
    for (Int_t ch = 0; ch < kMaxMuonTrackingChambers; ch++) {
      cout << "chamber(0...): " << ch
	   << "  NHitsForRec: " << fNHitsForRecPerChamber[ch]
	   << "  index(first HitForRec): " << fIndexOfFirstHitForRecPerChamber[ch]
	   << endl;
      for (Int_t hit = fIndexOfFirstHitForRecPerChamber[ch];
	   hit < fIndexOfFirstHitForRecPerChamber[ch] + fNHitsForRecPerChamber[ch];
	   hit++) {
	cout << "HitForRec index(0...): " << hit << endl;
	((*fHitsForRecPtr)[hit])->Dump();
      }
    }
  }
  return;
}

  //__________________________________________________________________________
void AliMUONEventReconstructor::AddHitsForRecFromGEANT(TTree *TH)
{
  // To add to the list of hits for reconstruction
  // the GEANT signal hits from a hit tree TH.
  if (fPrintLevel >= 2)
    cout << "enter AddHitsForRecFromGEANT with TH: " << TH << endl;
  if (TH == NULL) return;
  AliMUON *pMUON  = (AliMUON*) gAlice->GetModule("MUON"); // necessary ????
  // Security on MUON ????
  // See whether it could be the same for signal and background ????
  // Loop over tracks in tree
  Int_t ntracks = (Int_t) TH->GetEntries();
  if (fPrintLevel >= 2)
    cout << "ntracks: " << ntracks << endl;
  for (Int_t track = 0; track < ntracks; track++) {
    gAlice->ResetHits();
    TH->GetEvent(track);
    // Loop over hits
    Int_t hit = 0;
    for (AliMUONHit* mHit = (AliMUONHit*) pMUON->FirstHit(-1); 
	 mHit;
	 mHit = (AliMUONHit*) pMUON->NextHit(), hit++) {
      NewHitForRecFromGEANT(mHit,track, hit, 1);
    } // end of hit loop
  } // end of track loop
  return;
}

  //__________________________________________________________________________
void AliMUONEventReconstructor::AddHitsForRecFromBkgGEANT(TTree *TH, TClonesArray *Hits)
{
  // To add to the list of hits for reconstruction
  // the GEANT background hits from a hit tree TH and a pointer Hits to a hit list.
  // How to have only one function "AddHitsForRecFromGEANT" ????
  if (fPrintLevel >= 2)
    cout << "enter AddHitsForRecFromBkgGEANT with TH: " << TH << endl;
  if (TH == NULL) return;
  // Loop over tracks in tree
  Int_t ntracks = (Int_t) TH->GetEntries();
  if (fPrintLevel >= 2)
    cout << "ntracks: " << ntracks << endl;
  for (Int_t track = 0; track < ntracks; track++) {
    if (Hits) Hits->Clear();
    TH->GetEvent(track);
    // Loop over hits
    for (Int_t hit = 0; hit < Hits->GetEntriesFast(); hit++) {
      NewHitForRecFromGEANT((AliMUONHit*) (*Hits)[hit], track, hit, 0);
    } // end of hit loop
  } // end of track loop
  return;
}

  //__________________________________________________________________________
AliMUONHitForRec* AliMUONEventReconstructor::NewHitForRecFromGEANT(AliMUONHit* Hit, Int_t TrackNumber, Int_t HitNumber, Int_t Signal)
{
  // To make a new hit for reconstruction from a GEANT hit pointed to by "Hit",
  // with hit number "HitNumber" in the track numbered "TrackNumber",
  // either from signal ("Signal" = 1) or background ("Signal" = 0) event.
  // Selects hits in tracking (not trigger) chambers.
  // Takes into account the efficiency (fEfficiency)
  // and the smearing from resolution (fBendingResolution and fNonBendingResolution).
  // Adds a condition on the radius between RMin and RMax
  // to better simulate the real chambers.
  // Returns the pointer to the new hit for reconstruction,
  // or NULL in case of inefficiency or non tracking chamber or bad radius.
  // No condition on at most 20 cm from a muon from a resonance
  // like in Fortran TRACKF_STAT.
  AliMUONHitForRec* hitForRec;
  Double_t bendCoor, nonBendCoor, radius;
  Int_t chamber = Hit->fChamber - 1; // chamber(0...)
  // only in tracking chambers (fChamber starts at 1)
  if (chamber >= kMaxMuonTrackingChambers) return NULL;
  // only if hit is efficient (keep track for checking ????)
  if (gRandom->Rndm() > fEfficiency) return NULL;
  // only if radius between RMin and RMax
  bendCoor = Hit->fY;
  nonBendCoor = Hit->fX;
  radius = TMath::Sqrt((bendCoor * bendCoor) + (nonBendCoor * nonBendCoor));
  if ((radius < fRMin[chamber]) || (radius > fRMax[chamber])) return NULL;
  // new AliMUONHitForRec from GEANT hit and increment number of AliMUONHitForRec's
  hitForRec = new ((*fHitsForRecPtr)[fNHitsForRec]) AliMUONHitForRec(Hit);
  fNHitsForRec++;
  // add smearing from resolution
  hitForRec->SetBendingCoor(bendCoor + gRandom->Gaus(0., fBendingResolution));
  hitForRec->SetNonBendingCoor(nonBendCoor
			       + gRandom->Gaus(0., fNonBendingResolution));
  // more information into HitForRec
  //  resolution: angular effect to be added here ????
  hitForRec->SetBendingReso2(fBendingResolution * fBendingResolution);
  hitForRec->SetNonBendingReso2(fNonBendingResolution * fNonBendingResolution);
  //  GEANT track info
  hitForRec->SetHitNumber(HitNumber);
  hitForRec->SetTHTrack(TrackNumber);
  hitForRec->SetGeantSignal(Signal);
  if (fPrintLevel >= 10) {
    cout << "track: " << TrackNumber << " hit: " << HitNumber << endl;
    Hit->Dump();
    cout << "AliMUONHitForRec number (1...): " << fNHitsForRec << endl;
    hitForRec->Dump();}
  return hitForRec;
}

  //__________________________________________________________________________
void AliMUONEventReconstructor::SortHitsForRecWithIncreasingChamber()
{
  // Sort HitsForRec's in increasing order with respect to chamber number.
  // Uses the function "Compare".
  // Update the information for HitsForRec per chamber too.
  Int_t ch, nhits, prevch;
  fHitsForRecPtr->Sort();
  for (ch = 0; ch < kMaxMuonTrackingChambers; ch++) {
    fNHitsForRecPerChamber[ch] = 0;
    fIndexOfFirstHitForRecPerChamber[ch] = 0;
  }
  prevch = 0; // previous chamber
  nhits = 0; // number of hits in current chamber
  // Loop over HitsForRec
  for (Int_t hit = 0; hit < fNHitsForRec; hit++) {
    // chamber number (0...)
    ch = ((AliMUONHitForRec*)  ((*fHitsForRecPtr)[hit]))->GetChamberNumber();
    // increment number of hits in current chamber
    (fNHitsForRecPerChamber[ch])++;
    // update index of first HitForRec in current chamber
    // if chamber number different from previous one
    if (ch != prevch) {
      fIndexOfFirstHitForRecPerChamber[ch] = hit;
      prevch = ch;
    }
  }
  return;
}

//   //__________________________________________________________________________
// void AliMUONEventReconstructor::AddHitsForRecFromCathodeCorrelations(TTree* TC)
// {
//   // OLD VERSION WHEN ONE ONE WAS USING SO CALLED CATHODE CORRELATIONS
//   // To add to the list of hits for reconstruction
//   // the (cathode correlated) raw clusters
//   // No condition added, like in Fortran TRACKF_STAT,
//   // on the radius between RMin and RMax.
//   AliMUONHitForRec *hitForRec;
//   if (fPrintLevel >= 1) cout << "enter AddHitsForRecFromCathodeCorrelations" << endl;
//   AliMUON *MUON  = (AliMUON*) gAlice->GetModule("MUON"); // necessary ????
//   // Security on MUON ????
//   // Loop over tracking chambers
//   for (Int_t ch = 0; ch < kMaxMuonTrackingChambers; ch++) {
//     // number of HitsForRec to 0 for the chamber
//     fNHitsForRecPerChamber[ch] = 0;
//     // index of first HitForRec for the chamber
//     if (ch == 0) fIndexOfFirstHitForRecPerChamber[ch] = 0;
//     else fIndexOfFirstHitForRecPerChamber[ch] = fNHitsForRec;
//     TClonesArray *reconst_hits  = MUON->ReconstHitsAddress(ch);
//     MUON->ResetReconstHits();
//     TC->GetEvent();
//     Int_t ncor = (Int_t)reconst_hits->GetEntries();
//     // Loop over (cathode correlated) raw clusters
//     for (Int_t cor = 0; cor < ncor; cor++) {
//       AliMUONReconstHit * mCor = 
// 	(AliMUONReconstHit*) reconst_hits->UncheckedAt(cor);
//       // new AliMUONHitForRec from (cathode correlated) raw cluster
//       // and increment number of AliMUONHitForRec's (total and in chamber)
//       hitForRec = new ((*fHitsForRecPtr)[fNHitsForRec]) AliMUONHitForRec(mCor);
//       fNHitsForRec++;
//       (fNHitsForRecPerChamber[ch])++;
//       // more information into HitForRec
//       hitForRec->SetChamberNumber(ch);
//       hitForRec->SetHitNumber(cor);
//       // Z coordinate of the chamber (cm) with sign opposite to GEANT sign
//       // could (should) be more exact from chamber geometry ???? 
//       hitForRec->SetZ(-(&(MUON->Chamber(ch)))->Z());
//       if (fPrintLevel >= 10) {
// 	cout << "chamber (0...): " << ch <<
// 	  " cathcorrel (0...): " << cor << endl;
// 	mCor->Dump();
// 	cout << "AliMUONHitForRec number (1...): " << fNHitsForRec << endl;
// 	hitForRec->Dump();}
//     } // end of cluster loop
//   } // end of chamber loop
//   return;
// }

  //__________________________________________________________________________
void AliMUONEventReconstructor::AddHitsForRecFromRawClusters(TTree* TR)
{
  // To add to the list of hits for reconstruction all the raw clusters
  // No condition added, like in Fortran TRACKF_STAT,
  // on the radius between RMin and RMax.
  AliMUONHitForRec *hitForRec;
  AliMUONRawCluster *clus;
  Int_t iclus, nclus;
  TClonesArray *rawclusters;
  if (fPrintLevel >= 1) cout << "enter AddHitsForRecFromRawClusters" << endl;
  AliMUON *pMUON  = (AliMUON*) gAlice->GetModule("MUON"); // necessary ????
  // Security on MUON ????
  // Loop over tracking chambers
  for (Int_t ch = 0; ch < kMaxMuonTrackingChambers; ch++) {
    // number of HitsForRec to 0 for the chamber
    fNHitsForRecPerChamber[ch] = 0;
    // index of first HitForRec for the chamber
    if (ch == 0) fIndexOfFirstHitForRecPerChamber[ch] = 0;
    else fIndexOfFirstHitForRecPerChamber[ch] = fNHitsForRec;
    rawclusters = pMUON->RawClustAddress(ch);
    pMUON->ResetRawClusters();
    TR->GetEvent((Int_t) (TR->GetEntries()) - 1); // to be checked ????
    nclus = (Int_t) (rawclusters->GetEntries());
    // Loop over (cathode correlated) raw clusters
    for (iclus = 0; iclus < nclus; iclus++) {
      clus = (AliMUONRawCluster*) rawclusters->UncheckedAt(iclus);
      // new AliMUONHitForRec from raw cluster
      // and increment number of AliMUONHitForRec's (total and in chamber)
      hitForRec = new ((*fHitsForRecPtr)[fNHitsForRec]) AliMUONHitForRec(clus);
      fNHitsForRec++;
      (fNHitsForRecPerChamber[ch])++;
      // more information into HitForRec
      //  resolution: info should be already in raw cluster and taken from it ????
      hitForRec->SetBendingReso2(fBendingResolution * fBendingResolution);
      hitForRec->SetNonBendingReso2(fNonBendingResolution * fNonBendingResolution);
      //  original raw cluster
      hitForRec->SetChamberNumber(ch);
      hitForRec->SetHitNumber(iclus);
      // Z coordinate of the chamber (cm) with sign opposite to GEANT sign
      // could (should) be more exact from chamber geometry ???? 
      hitForRec->SetZ(-(&(pMUON->Chamber(ch)))->Z());
      if (fPrintLevel >= 10) {
	cout << "chamber (0...): " << ch <<
	  " raw cluster (0...): " << iclus << endl;
	clus->Dump();
	cout << "AliMUONHitForRec number (1...): " << fNHitsForRec << endl;
	hitForRec->Dump();}
    } // end of cluster loop
  } // end of chamber loop
  return;
}

  //__________________________________________________________________________
void AliMUONEventReconstructor::MakeSegments(void)
{
  // To make the list of segments in all stations,
  // from the list of hits to be reconstructed
  if (fPrintLevel >= 1) cout << "enter MakeSegments" << endl;
  ResetSegments();
  // Loop over stations
  for (Int_t st = 0; st < kMaxMuonTrackingStations; st++)
    MakeSegmentsPerStation(st); 
  if (fPrintLevel >= 10) {
    cout << "end of MakeSegments" << endl;
    for (Int_t st = 0; st < kMaxMuonTrackingStations; st++) {
      cout << "station(0...): " << st
	   << "  Segments: " << fNSegments[st]
	   << endl;
      for (Int_t seg = 0;
	   seg < fNSegments[st];
	   seg++) {
	cout << "Segment index(0...): " << seg << endl;
	((*fSegmentsPtr[st])[seg])->Dump();
      }
    }
  }
  return;
}

  //__________________________________________________________________________
void AliMUONEventReconstructor::MakeSegmentsPerStation(Int_t Station)
{
  // To make the list of segments in station number "Station" (0...)
  // from the list of hits to be reconstructed.
  // Updates "fNSegments"[Station].
  // Segments in stations 4 and 5 are sorted
  // according to increasing absolute value of "impact parameter"
  AliMUONHitForRec *hit1Ptr, *hit2Ptr;
  AliMUONSegment *segment;
  Bool_t last2st;
  Double_t bendingSlope, distBend, distNonBend, extBendCoor, extNonBendCoor,
    impactParam, maxImpactParam;
  AliMUON *pMUON  = (AliMUON*) gAlice->GetModule("MUON"); // necessary ????
  if (fPrintLevel >= 1)
    cout << "enter MakeSegmentsPerStation (0...) " << Station << endl;
  // first and second chambers (0...) in the station
  Int_t ch1 = 2 * Station;
  Int_t ch2 = ch1 + 1;
  // variable true for stations downstream of the dipole:
  // Station(0..4) equal to 3 or 4
  if ((Station == 3) || (Station == 4)) {
    last2st = kTRUE;
    // maximum impact parameter (cm) according to fMinBendingMomentum...
    maxImpactParam =
      TMath::Abs(GetImpactParamFromBendingMomentum(fMinBendingMomentum));
  }
  else last2st = kFALSE;
  // extrapolation factor from Z of first chamber to Z of second chamber
  // dZ to be changed to take into account fine structure of chambers ????
  Double_t extrapFact =
    (&(pMUON->Chamber(ch2)))->Z() / (&(pMUON->Chamber(ch1)))->Z();
  // index for current segment
  Int_t segmentIndex = 0;
  // Loop over HitsForRec in the first chamber of the station
  for (Int_t hit1 = fIndexOfFirstHitForRecPerChamber[ch1];
       hit1 < fIndexOfFirstHitForRecPerChamber[ch1] + fNHitsForRecPerChamber[ch1];
       hit1++) {
    // pointer to the HitForRec
    hit1Ptr = (AliMUONHitForRec*) ((*fHitsForRecPtr)[hit1]);
    // extrapolation,
    // on the straight line joining the HitForRec to the vertex (0,0,0),
    // to the Z of the second chamber of the station
    extBendCoor = extrapFact * hit1Ptr->GetBendingCoor();
    extNonBendCoor = extrapFact * hit1Ptr->GetNonBendingCoor();
    // Loop over HitsForRec in the second chamber of the station
    for (Int_t hit2 = fIndexOfFirstHitForRecPerChamber[ch2];
	 hit2 < fIndexOfFirstHitForRecPerChamber[ch2] + fNHitsForRecPerChamber[ch2];
	 hit2++) {
      // pointer to the HitForRec
      hit2Ptr = (AliMUONHitForRec*) ((*fHitsForRecPtr)[hit2]);
      // absolute values of distances, in bending and non bending planes,
      // between the HitForRec in the second chamber
      // and the previous extrapolation
      distBend = TMath::Abs(hit2Ptr->GetBendingCoor() - extBendCoor);
      distNonBend = TMath::Abs(hit2Ptr->GetNonBendingCoor() - extNonBendCoor);
      if (last2st) {
	// bending slope
	bendingSlope = (hit1Ptr->GetBendingCoor() - hit2Ptr->GetBendingCoor()) /
	  (hit1Ptr->GetZ() - hit2Ptr->GetZ());
	// absolute value of impact parameter
	impactParam =
	  TMath::Abs(hit1Ptr->GetBendingCoor() - hit2Ptr->GetZ() * bendingSlope);
      }
      // check for distances not too large,
      // and impact parameter not too big if stations downstream of the dipole.
      // Conditions "distBend" and "impactParam" correlated for these stations ????
      if ((distBend < fSegmentMaxDistBending[Station]) &&
	  (distNonBend < fSegmentMaxDistNonBending[Station]) &&
	  (!last2st || (impactParam < maxImpactParam))) {
	// make new segment
	segment = new ((*fSegmentsPtr[Station])[segmentIndex])
	  AliMUONSegment(hit1Ptr, hit2Ptr);
	// update "link" to this segment from the hit in the first chamber
	if (hit1Ptr->GetNSegments() == 0)
	  hit1Ptr->SetIndexOfFirstSegment(segmentIndex);
	hit1Ptr->SetNSegments(hit1Ptr->GetNSegments() + 1);
	if (fPrintLevel >= 10) {
	  cout << "segmentIndex(0...): " << segmentIndex
	       << "  distBend: " << distBend
	       << "  distNonBend: " << distNonBend
	       << endl;
	  segment->Dump();
	  cout << "HitForRec in first chamber" << endl;
	  hit1Ptr->Dump();
	  cout << "HitForRec in second chamber" << endl;
	  hit2Ptr->Dump();
	};
	// increment index for current segment
	segmentIndex++;
      }
    } //for (Int_t hit2
  } // for (Int_t hit1...
  fNSegments[Station] = segmentIndex;
  // Sorting according to "impact parameter" if station(1..5) 4 or 5,
  // i.e. Station(0..4) 3 or 4, using the function "Compare".
  // After this sorting, it is impossible to use
  // the "fNSegments" and "fIndexOfFirstSegment"
  // of the HitForRec in the first chamber to explore all segments formed with it.
  // Is this sorting really needed ????
  if ((Station == 3) || (Station == 4)) (fSegmentsPtr[Station])->Sort();
  if (fPrintLevel >= 1) cout << "Station: " << Station << "  NSegments: "
			     << fNSegments[Station] << endl;
  return;
}

  //__________________________________________________________________________
void AliMUONEventReconstructor::MakeTracks(void)
{
  // To make the tracks,
  // from the list of segments and points in all stations
  if (fPrintLevel >= 1) cout << "enter MakeTracks" << endl;
  ResetTracks();
  // Look for candidates from at least 3 aligned points in stations(1..) 4 and 5
  MakeTrackCandidates();
  // Follow tracks in stations(1..) 3, 2 and 1
  FollowTracks();
  return;
}

  //__________________________________________________________________________
Int_t AliMUONEventReconstructor::MakeTrackCandidatesWithTwoSegments(AliMUONSegment *BegSegment)
{
  // To make track candidates with two segments in stations(1..) 4 and 5,
  // the first segment being pointed to by "BegSegment".
  // Returns the number of such track candidates.
  Int_t endStation, iEndSegment, nbCan2Seg;
  AliMUONSegment *endSegment, *extrapSegment;
  AliMUONTrack *recTrack;
  Double_t mcsFactor;
  if (fPrintLevel >= 1) cout << "enter MakeTrackCandidatesWithTwoSegments" << endl;
  // Station for the end segment
  endStation = 7 - (BegSegment->GetHitForRec1())->GetChamberNumber() / 2;
  // multiple scattering factor corresponding to one chamber
  mcsFactor = 0.0136 /
    GetBendingMomentumFromImpactParam(BegSegment->GetBendingImpact());
  mcsFactor	= fChamberThicknessInX0 * mcsFactor * mcsFactor;
  // linear extrapolation to end station
  extrapSegment =
    BegSegment->CreateSegmentFromLinearExtrapToStation(endStation, mcsFactor);
  // number of candidates with 2 segments to 0
  nbCan2Seg = 0;
  // Loop over segments in the end station
  for (iEndSegment = 0; iEndSegment < fNSegments[endStation]; iEndSegment++) {
    // pointer to segment
    endSegment = (AliMUONSegment*) ((*fSegmentsPtr[endStation])[iEndSegment]);
    // test compatibility between current segment and "extrapSegment"
    if ((endSegment->
	 NormalizedChi2WithSegment(extrapSegment,
				   fMaxSigma2Distance)) <= 4.0) {
      // both segments compatible:
      // make track candidate from "begSegment" and "endSegment"
      if (fPrintLevel >= 2)
	cout << "TrackCandidate with Segment " << iEndSegment <<
	  " in Station(0..) " << endStation << endl;
      // flag for both segments in one track:
      // to be done in track constructor ????
      BegSegment->SetInTrack(kTRUE);
      endSegment->SetInTrack(kTRUE);
      recTrack = new ((*fRecTracksPtr)[fNRecTracks])
	AliMUONTrack(BegSegment, endSegment, this);
      fNRecTracks++;
      if (fPrintLevel >= 10) recTrack->RecursiveDump();
      // increment number of track candidates with 2 segments
      nbCan2Seg++;
    }
  } // for (iEndSegment = 0;...
  delete extrapSegment; // should not delete HitForRec's it points to !!!!
  return nbCan2Seg;
}

  //__________________________________________________________________________
Int_t AliMUONEventReconstructor::MakeTrackCandidatesWithOneSegmentAndOnePoint(AliMUONSegment *BegSegment)
{
  // To make track candidates with one segment and one point
  // in stations(1..) 4 and 5,
  // the segment being pointed to by "BegSegment".
  Int_t ch, ch1, ch2, endStation, iHit, iHitMax, iHitMin, nbCan1Seg1Hit;
  AliMUONHitForRec *extrapHitForRec, *hit;
  AliMUONTrack *recTrack;
  Double_t mcsFactor;
  if (fPrintLevel >= 1)
    cout << "enter MakeTrackCandidatesWithOneSegmentAndOnePoint" << endl;
  // station for the end point
  endStation = 7 - (BegSegment->GetHitForRec1())->GetChamberNumber() / 2;
  // multiple scattering factor corresponding to one chamber
  mcsFactor = 0.0136 /
    GetBendingMomentumFromImpactParam(BegSegment->GetBendingImpact());
  mcsFactor	= fChamberThicknessInX0 * mcsFactor * mcsFactor;
  // first and second chambers(0..) in the end station
  ch1 = 2 * endStation;
  ch2 = ch1 + 1;
  // number of candidates to 0
  nbCan1Seg1Hit = 0;
  // Loop over chambers of the end station
  for (ch = ch2; ch >= ch1; ch--) {
    // linear extrapolation to chamber
    extrapHitForRec =
      BegSegment->CreateHitForRecFromLinearExtrapToChamber(ch, mcsFactor);
    // limits for the hit index in the loop
    iHitMin = fIndexOfFirstHitForRecPerChamber[ch];
    iHitMax = iHitMin + fNHitsForRecPerChamber[ch];
    // Loop over HitForRec's in the chamber
    for (iHit = iHitMin; iHit < iHitMax; iHit++) {
      // pointer to HitForRec
      hit = (AliMUONHitForRec*) ((*fHitsForRecPtr)[iHit]);
      // test compatibility between current HitForRec and "extrapHitForRec"
      if ((hit->
	   NormalizedChi2WithHitForRec(extrapHitForRec,
				       fMaxSigma2Distance)) <= 2.0) {
	// both HitForRec's compatible:
	// make track candidate from begSegment and current HitForRec
	if (fPrintLevel >= 2)
	  cout << "TrackCandidate with HitForRec " << iHit <<
	    " in Chamber(0..) " << ch << endl;
	// flag for beginning segments in one track:
	// to be done in track constructor ????
	BegSegment->SetInTrack(kTRUE);
	recTrack = new ((*fRecTracksPtr)[fNRecTracks])
	  AliMUONTrack(BegSegment, hit, this);
	// the right place to eliminate "double counting" ???? how ????
	fNRecTracks++;
	if (fPrintLevel >= 10) recTrack->RecursiveDump();
	// increment number of track candidates
	nbCan1Seg1Hit++;
      }
    } // for (iHit = iHitMin;...
    delete extrapHitForRec;
  } // for (ch = ch2;...
  return nbCan1Seg1Hit;
}

  //__________________________________________________________________________
void AliMUONEventReconstructor::MakeTrackCandidates(void)
{
  // To make track candidates
  // with at least 3 aligned points in stations(1..) 4 and 5
  // (two Segment's or one Segment and one HitForRec)
  Int_t begStation, iBegSegment, nbCan1Seg1Hit, nbCan2Seg;
  AliMUONSegment *begSegment;
  if (fPrintLevel >= 1) cout << "enter MakeTrackCandidates" << endl;
  // Loop over stations(1..) 5 and 4 for the beginning segment
  for (begStation = 4; begStation > 2; begStation--) {
    // Loop over segments in the beginning station
    for (iBegSegment = 0; iBegSegment < fNSegments[begStation]; iBegSegment++) {
      // pointer to segment
      begSegment = (AliMUONSegment*) ((*fSegmentsPtr[begStation])[iBegSegment]);
      if (fPrintLevel >= 2)
	cout << "look for TrackCandidate's with Segment " << iBegSegment <<
	  " in Station(0..) " << begStation << endl;
      // Look for track candidates with two segments,
      // "begSegment" and all compatible segments in other station.
      // Only for beginning station(1..) 5
      // because candidates with 2 segments have to looked for only once.
      if (begStation == 4)
	nbCan2Seg = MakeTrackCandidatesWithTwoSegments(begSegment);
      // Look for track candidates with one segments and one point,
      // "begSegment" and all compatible HitForRec's in other station.
      // Only if "begSegment" does not belong already to a track candidate.
      // Is that a too strong condition ????
      if (!(begSegment->GetInTrack()))
	nbCan1Seg1Hit = MakeTrackCandidatesWithOneSegmentAndOnePoint(begSegment);
    } // for (iBegSegment = 0;...
  } // for (begStation = 4;...
  return;
}

  //__________________________________________________________________________
void AliMUONEventReconstructor::FollowTracks(void)
{
  // Follow tracks in stations(1..) 3, 2 and 1
  AliMUONHitForRec *bestHit, *extrapHit, *hit;
  AliMUONSegment *bestSegment, *extrapSegment, *segment;
  AliMUONTrack *track;
  AliMUONTrackParam *trackParam1, trackParam[2];
  Int_t ch, chBestHit, iHit, iSegment, station, trackIndex;
  Double_t bestChi2, chi2, dZ1, dZ2, maxSigma2Distance, mcsFactor;
  AliMUON *pMUON = (AliMUON*) gAlice->GetModule("MUON"); // necessary ????
  maxSigma2Distance = 4. * fMaxSigma2Distance; // sigma2cut increased by 4 !!!!
  if (fPrintLevel >= 2)
    cout << "enter FollowTracks" << endl;
  // Loop over track candidates
  for (trackIndex = 0; trackIndex < fNRecTracks; trackIndex++) {
    if (fPrintLevel >= 2)
      cout << "FollowTracks: track candidate(0..): " << trackIndex << endl;
    // function for each track candidate ????
    track = (AliMUONTrack*) ((*fRecTracksPtr)[trackIndex]);
    // Fit track candidate from vertex at X = Y = 0 
    track->Fit(track->GetTrackParamAtVertex(), 3);
    if (fPrintLevel >= 10) {
      cout << "FollowTracks: track candidate(0..): " << trackIndex
	   << " after fit in stations(1..) 4 and 5" << endl;
      track->RecursiveDump();
    }
    // Loop over stations(1..) 3, 2 and 1
    // something SPECIAL for stations 2 and 1 for majority coincidence ????
    for (station = 2; station >= 0; station--) {
      // Track parameters at first track hit (smallest Z)
      trackParam1 = ((AliMUONTrackHit*)
		     (track->GetTrackHitsPtr()->First()))->GetTrackParam();
      // extrapolation to station
      trackParam1->ExtrapToStation(station, trackParam);
      extrapSegment = new AliMUONSegment(); //  empty segment
      // multiple scattering factor corresponding to one chamber
      // and momentum in bending plane (not total)
      mcsFactor = 0.0136 * trackParam1->GetInverseBendingMomentum();
      mcsFactor	= fChamberThicknessInX0 * mcsFactor * mcsFactor;
      // Z difference from previous station
      dZ1 = (&(pMUON->Chamber(2 * station)))->Z() -
	(&(pMUON->Chamber(2 * station + 2)))->Z();
      // Z difference between the two previous stations
      dZ2 = (&(pMUON->Chamber(2 * station + 2)))->Z() -
	(&(pMUON->Chamber(2 * station + 4)))->Z();
      extrapSegment->SetBendingCoorReso2(fBendingResolution);
      extrapSegment->SetNonBendingCoorReso2(fNonBendingResolution);
      extrapSegment->UpdateFromStationTrackParam(trackParam, mcsFactor, dZ1, dZ2);
      bestChi2 = 5.0;
      bestSegment = NULL;
      if (fPrintLevel >= 10) {
	cout << "FollowTracks: track candidate(0..): " << trackIndex
	     << " Look for segment in station(0..): " << station << endl;
      }
      // Loop over segments in station
      for (iSegment = 0; iSegment < fNSegments[station]; iSegment++) {
	// Look for best compatible Segment in station
	// should consider all possibilities ????
	// multiple scattering ????
	// separation in 2 functions: Segment and HitForRec ????
	segment = (AliMUONSegment*) ((*fSegmentsPtr[station])[iSegment]);
	chi2 = segment->NormalizedChi2WithSegment(extrapSegment, maxSigma2Distance);
	if (chi2 < bestChi2) {
	  // update best Chi2 and Segment if better found
	  bestSegment = segment;
	  bestChi2 = chi2;
	}
      }
      if (bestSegment) {
	// best segment found: add it to track candidate
	track->AddSegment(bestSegment);
	// set track parameters at these two TrakHit's
	track->SetTrackParamAtHit(track->GetNTrackHits() - 2, &(trackParam[0]));
	track->SetTrackParamAtHit(track->GetNTrackHits() - 1, &(trackParam[1]));
	if (fPrintLevel >= 10) {
	  cout << "FollowTracks: track candidate(0..): " << trackIndex
	       << " Added segment in station(0..): " << station << endl;
	  track->RecursiveDump();
	}
      }
      else {
	// No best segment found:
	// Look for best compatible HitForRec in station:
	// should consider all possibilities ????
	// multiple scattering ???? do about like for extrapSegment !!!!
	extrapHit = new AliMUONHitForRec(); //  empty hit
	bestChi2 = 3.0;
	bestHit = NULL;
	if (fPrintLevel >= 10) {
	  cout << "FollowTracks: track candidate(0..): " << trackIndex
	       << " Segment not found, look for hit in station(0..): " << station
	       << endl;
	}
	// Loop over chambers of the station
	for (ch = 0; ch < 2; ch++) {
	  // coordinates of extrapolated hit
	  extrapHit->SetBendingCoor((&(trackParam[ch]))->GetBendingCoor());
	  extrapHit->SetNonBendingCoor((&(trackParam[ch]))->GetNonBendingCoor());
	  // resolutions from "extrapSegment"
	  extrapHit->SetBendingReso2(extrapSegment->GetBendingCoorReso2());
	  extrapHit->SetNonBendingReso2(extrapSegment->GetNonBendingCoorReso2());
	  // Loop over hits in the chamber
	  for (iHit = fIndexOfFirstHitForRecPerChamber[ch];
	       iHit < fIndexOfFirstHitForRecPerChamber[ch] +
		 fNHitsForRecPerChamber[ch];
	       iHit++) {
	    hit = (AliMUONHitForRec*) ((*fHitsForRecPtr)[iHit]);
	    // condition for hit not already in segment ????
	    chi2 = hit->NormalizedChi2WithHitForRec(extrapHit, maxSigma2Distance);
	    if (chi2 < bestChi2) {
	      // update best Chi2 and HitForRec if better found
	      bestHit = hit;
	      bestChi2 = chi2;
	      chBestHit = ch;
	    }
	  }
	}
	if (bestHit) {
	  // best hit found: add it to track candidate
	  track->AddHitForRec(bestHit);
	  // set track parameters at these two TrakHit's
	  track->SetTrackParamAtHit(track->GetNTrackHits() - 1,
				    &(trackParam[chBestHit]));
	  if (fPrintLevel >= 10) {
	    cout << "FollowTracks: track candidate(0..): " << trackIndex
		 << " Added hit in station(0..): " << station << endl;
	    track->RecursiveDump();
	  }
	}
	else {
	  fRecTracksPtr->RemoveAt(trackIndex); // Remove candidate
	  break; // stop the search for this candidate:
	  // exit from the loop over station
	}
      }
      // Sort track hits according to increasing Z
      track->GetTrackHitsPtr()->Sort();
      // Update track parameters at first track hit (smallest Z)
      trackParam1 = ((AliMUONTrackHit*)
		     (track->GetTrackHitsPtr()->First()))->GetTrackParam();
      // Track fit from first track hit varying X and Y
      track->Fit(trackParam1, 5);
      if (fPrintLevel >= 10) {
	cout << "FollowTracks: track candidate(0..): " << trackIndex
	     << " after fit from station(0..): " << station << " to 4" << endl;
	track->RecursiveDump();
      }
      delete extrapSegment;
    } // for (station = 2;...
  } // for (trackIndex = 0;...
  return;
}

void AliMUONEventReconstructor::Streamer(TBuffer &R__b)
{
    ;
}
