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

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD trigger class                                                        //
//                                                                           //
//  Author:                                                                  //
//    Bogdan Vulpescu                                                        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TTree.h>
#include <TBranch.h>
#include <TMatrixD.h>
#include <TClonesArray.h>
#include <TObjArray.h>

#include "AliLog.h"
#include "AliRun.h"
#include "AliLoader.h"

#include "AliTRDdigitsManager.h"
#include "AliTRDgeometry.h"
#include "AliTRDdataArrayI.h"
#include "AliTRDcalibDB.h"
#include "AliTRDCommonParam.h"
#include "AliTRDrawData.h"
#include "AliTRDtrigger.h"
#include "AliTRDmodule.h"
#include "AliTRDmcmTracklet.h"
#include "AliTRDgtuTrack.h"
#include "AliTRDtrigParam.h"
#include "AliTRDmcm.h"
#include "AliTRDzmaps.h"
#include "AliTRDCalibra.h"
#include "Cal/AliTRDCalPIDLQ.h"

ClassImp(AliTRDtrigger)

//_____________________________________________________________________________
AliTRDtrigger::AliTRDtrigger()
  :TNamed()
  ,fField(0)
  ,fGeo(NULL)
  ,fCalib(NULL)
  ,fCParam(NULL)
  ,fTrigParam(NULL)
  ,fRunLoader(NULL)
  ,fDigitsManager(NULL)
  ,fTrackletTree(NULL)
  ,fTracklets(NULL)
  ,fNROB(0)
  ,fMCM(NULL)
  ,fTrk(NULL)
  ,fTrkTest(NULL)
  ,fModule(NULL)
  ,fGTUtrk(NULL)
  ,fNtracklets(0)
  ,fDigits(NULL)
  ,fTrack0(NULL)
  ,fTrack1(NULL)
  ,fTrack2(NULL)
  ,fNPrimary(0)
  ,fTracks(NULL)
{
  //
  // AliTRDtrigger default constructor
  //

}

//_____________________________________________________________________________
AliTRDtrigger::AliTRDtrigger(const Text_t *name, const Text_t *title)
  :TNamed(name,title)
  ,fField(0)
  ,fGeo(NULL)
  ,fCalib(NULL)
  ,fCParam(NULL)
  ,fTrigParam(NULL)
  ,fRunLoader(NULL)
  ,fDigitsManager(new AliTRDdigitsManager())
  ,fTrackletTree(NULL)
  ,fTracklets(new TObjArray(400))
  ,fNROB(0)
  ,fMCM(NULL)
  ,fTrk(NULL)
  ,fTrkTest(NULL)
  ,fModule(NULL)
  ,fGTUtrk(NULL)
  ,fNtracklets(0)
  ,fDigits(NULL)
  ,fTrack0(NULL)
  ,fTrack1(NULL)
  ,fTrack2(NULL)
  ,fNPrimary(0)
  ,fTracks(new TClonesArray("AliTRDgtuTrack",1000))
{
  //
  // AliTRDtrigger constructor
  //

}

//_____________________________________________________________________________
AliTRDtrigger::AliTRDtrigger(const AliTRDtrigger &p)
  :TNamed(p)
  ,fField(p.fField)
  ,fGeo(NULL)
  ,fCalib(NULL)
  ,fCParam(NULL)
  ,fTrigParam(NULL)
  ,fRunLoader(NULL)
  ,fDigitsManager(NULL)
  ,fTrackletTree(NULL)
  ,fTracklets(NULL)
  ,fNROB(p.fNROB)
  ,fMCM(NULL)
  ,fTrk(NULL)
  ,fTrkTest(NULL)
  ,fModule(NULL)
  ,fGTUtrk(NULL)
  ,fNtracklets(p.fNtracklets)
  ,fDigits(NULL)
  ,fTrack0(NULL)
  ,fTrack1(NULL)
  ,fTrack2(NULL)
  ,fNPrimary(p.fNPrimary)
  ,fTracks(NULL)
{
  //
  // AliTRDtrigger copy constructor
  //

}

///_____________________________________________________________________________
AliTRDtrigger::~AliTRDtrigger()
{
  //
  // AliTRDtrigger destructor
  //

  if (fTracklets) {
    fTracklets->Delete();
    delete fTracklets;
  }

  if (fTracks) {
    fTracks->Delete();
    delete fTracks;
  }

}

//_____________________________________________________________________________
AliTRDtrigger &AliTRDtrigger::operator=(const AliTRDtrigger &p)
{
  //
  // Assignment operator
  //

  if (this != &p) ((AliTRDtrigger &) p).Copy(*this);
  return *this;

}

//_____________________________________________________________________________
void AliTRDtrigger::Copy(TObject &) const
{
  //
  // Copy function
  //

  AliFatal("Not implemented");

}

//_____________________________________________________________________________
void AliTRDtrigger::Init()
{

  fModule = new AliTRDmodule(fTrigParam); 
  fTracks->Clear();

  fField  = fTrigParam->GetField();
  fGeo    = (AliTRDgeometry*)AliTRDgeometry::GetGeometry(fRunLoader);

  fCalib  = AliTRDcalibDB::Instance();
  if (!fCalib) {
    AliError("No instance of AliTRDcalibDB.");
    return;  
  }

  fCParam = AliTRDCommonParam::Instance();
  if (!fCParam) {
    AliError("No common parameters.");
    return;
  }

}

//_____________________________________________________________________________
Bool_t AliTRDtrigger::Open(const Char_t *name, Int_t nEvent)
{
  //
  // Opens the AliROOT file.
  //

  TString evfoldname = AliConfig::GetDefaultEventFolderName();
  fRunLoader         = AliRunLoader::GetRunLoader(evfoldname);

  if (!fRunLoader) {
    fRunLoader = AliRunLoader::Open(name);
  }
  if (!fRunLoader) {
    AliError(Form("Can not open session for file %s.",name));
    return kFALSE;
  }

  // Import the Trees for the event nEvent in the file
  fRunLoader->GetEvent(nEvent);

  // Open output
  TObjArray *ioArray = 0;
  AliLoader* loader  = fRunLoader->GetLoader("TRDLoader");
  loader->MakeTree("T");
  fTrackletTree = loader->TreeT();
  fTrackletTree->Branch("TRDmcmTracklet","TObjArray",&ioArray,32000,0);
  Init();

  return kTRUE;

}

//_____________________________________________________________________________
Bool_t AliTRDtrigger::ReadDigits() 
{
  //
  // Reads the digits arrays from the input aliroot file
  //

  if (!fRunLoader) {
    AliError("Can not find the Run Loader");
    return kFALSE;
  }

  AliLoader* loader = fRunLoader->GetLoader("TRDLoader");
  if (!loader->TreeD()) {
    loader->LoadDigits();
  }
  if (!loader->TreeD()) {
     return kFALSE;
  }

  return (fDigitsManager->ReadDigits(loader->TreeD()));

}

//_____________________________________________________________________________
Bool_t AliTRDtrigger::ReadDigits(AliRawReader* rawReader)
{
  //
  // Reads the digits arrays from the ddl file
  //

  AliTRDrawData *raw = new AliTRDrawData();
  fDigitsManager     = raw->Raw2Digits(rawReader);

  return kTRUE;

}

//_____________________________________________________________________________
Bool_t AliTRDtrigger::ReadTracklets(AliRunLoader *rl) 
{
  //
  // Reads the tracklets find the tracks
  //

  Int_t idet;

  AliLoader *loader = rl->GetLoader("TRDLoader");
  loader->LoadTracks();
  fTrackletTree     = loader->TreeT();

  TBranch *branch   = fTrackletTree->GetBranch("TRDmcmTracklet");
  if (!branch) {
    AliError("Can't get the branch !");
    return kFALSE;
  }
  TObjArray *tracklets = new TObjArray(400);
  branch->SetAddress(&tracklets);

  Int_t nEntries   = (Int_t) fTrackletTree->GetEntries();
  Int_t iEntry;
  Int_t itrk;
  Int_t iStack;
  Int_t iStackPrev = -1;
  
  for (iEntry = 0; iEntry < nEntries; iEntry++) {    

    fTrackletTree->GetEvent(iEntry);
    
    for (itrk = 0; itrk < tracklets->GetEntriesFast(); itrk++) {

      fTrk   = (AliTRDmcmTracklet*)tracklets->UncheckedAt(itrk);
      idet   = fTrk->GetDetector();
      iStack = idet / (AliTRDgeometry::Nplan());

      if (iStackPrev != iStack) {
	if (iStackPrev == -1) {
	  iStackPrev = iStack;
	} 
        else {
	  MakeTracks(idet - AliTRDgeometry::Nplan());
	  ResetTracklets();
	  iStackPrev = iStack;
	}
      }
      
      Tracklets()->Add(fTrk);

      if ((iEntry == (nEntries-1)) && 
          (itrk   == (tracklets->GetEntriesFast() - 1))) {
	idet++;
	MakeTracks(idet-AliTRDgeometry::Nplan());
	ResetTracklets();
      }

    }

  }

  loader->UnloadTracks();

  return kTRUE;

}

//_____________________________________________________________________________
Bool_t AliTRDtrigger::MakeTracklets(Bool_t makeTracks)
{
  //
  // Create tracklets from digits
  //

  Int_t chamBeg = 0;
  Int_t chamEnd = AliTRDgeometry::Ncham();
  Int_t planBeg = 0;
  Int_t planEnd = AliTRDgeometry::Nplan();
  Int_t sectBeg = 0;
  Int_t sectEnd = AliTRDgeometry::Nsect();

  fTrkTest = new AliTRDmcmTracklet(0,0,0);
  fMCM     = new AliTRDmcm(fTrigParam,0);

  Int_t   time;
  Int_t   col;
  Int_t   row;
  Int_t   col1;
  Int_t   col2;
  Int_t   idet       = -1;
  Int_t   iStack     = -1;
  Int_t   iStackPrev = -1;
  Float_t amp;

  for (Int_t isect = sectBeg; isect < sectEnd; isect++) {

    for (Int_t icham = chamBeg; icham < chamEnd; icham++) {

      // Number of ROBs in the chamber
      if(icham == 2) {
	fNROB = 6;
      } 
      else {
	fNROB = 8;
      }

      for (Int_t iplan = planBeg; iplan < planEnd; iplan++) {

        idet = fGeo->GetDetector(iplan,icham,isect);
	ResetTracklets();
	
	if (makeTracks) {
	  iStack = idet / (AliTRDgeometry::Nplan());
	  if (iStackPrev != iStack) {
	    if (iStackPrev == -1) {
	      iStackPrev = iStack;
	    } 
            else {
	      MakeTracks(idet-AliTRDgeometry::Nplan());
	      ResetTracklets();
	      iStackPrev = iStack;
	    }
	  }
	}

        Int_t nRowMax    = fCParam->GetRowMax(iplan,icham,isect);
	Int_t nColMax    = fCParam->GetColMax(iplan);
        Int_t nTimeTotal = fCalib->GetNumberOfTimeBins();

        // Get the digits
        fDigits = fDigitsManager->GetDigits(idet);
	if (!fDigits) return kFALSE;
	// This is to take care of switched off super modules
        if (fDigits->GetNtime() == 0) {
          continue;
	}
        fDigits->Expand();
        fTrack0 = fDigitsManager->GetDictionary(idet,0);
	if (!fTrack0) return kFALSE;
        fTrack0->Expand();
        fTrack1 = fDigitsManager->GetDictionary(idet,1);
	if (!fTrack1) return kFALSE;
        fTrack1->Expand();
        fTrack2 = fDigitsManager->GetDictionary(idet,2); 
	if (!fTrack2) return kFALSE;
        fTrack2->Expand();

	for (Int_t iRob = 0; iRob < fNROB; iRob++) {

	  for (Int_t iMcm = 0; iMcm < kNMCM; iMcm++) {

	    fMCM->Reset();
	    fMCM->SetRobId(iRob);
	    fMCM->SetChaId(idet);

	    SetMCMcoordinates(iMcm);

	    row = fMCM->GetRow();

	    if ((row < 0) || (row >= nRowMax)) {
	      AliError("MCM row number out of range.");
	      continue;
	    }

	    fMCM->GetColRange(col1,col2);
	    
            for (time = 0; time < nTimeTotal; time++) {
	      for (col = col1; col < col2; col++) {
		if ((col >= 0) && (col < nColMax)) {
		  amp = TMath::Abs(fDigits->GetDataUnchecked(row,col,time));
		} 
                else {
		  amp = 0.0;
		}
		fMCM->SetADC(col-col1,time,amp);
	      }
	    }

	    if (fTrigParam->GetTailCancelation()) {
	      fMCM->Filter(fTrigParam->GetNexponential(),fTrigParam->GetFilterType());
	    }
	    
	    if (fMCM->Run()) {

	      for (Int_t iSeed = 0; iSeed < kMaxTrackletsPerMCM; iSeed++) {
		
		if (fMCM->GetSeedCol()[iSeed] < 0) {
                  continue;
		}

		if (fTrigParam->GetDebugLevel()   > 1) { 
		  AliInfo(Form("Add tracklet %d in col %02d \n",fNtracklets,fMCM->GetSeedCol()[iSeed]));
		}

		if (TestTracklet(idet,row,iSeed,0)) {
		  AddTracklet(idet,row,iSeed,fNtracklets++);
		}

	      }

	    }

	  }

      
	}

	// Compress the arrays
        fDigits->Compress(1,0);
        fTrack0->Compress(1,0);
        fTrack1->Compress(1,0);
        fTrack2->Compress(1,0);

	WriteTracklets(idet);

     }
    }
  }

  if (makeTracks) {
    idet++;
    MakeTracks(idet - AliTRDgeometry::Nplan());
    ResetTracklets();
  }

  return kTRUE;

}

//_____________________________________________________________________________
void AliTRDtrigger::SetMCMcoordinates(Int_t imcm)
{
  //
  // Configure MCM position in the pad plane
  //

  Int_t robid = fMCM->GetRobId();

  // setting the Row and Col range

  const Int_t kNcolRob = 2;  // number of ROBs per chamber in column direction
  const Int_t kNmcmRob = 4;  // number of MCMs per ROB in column/row direction

  Int_t mcmid = imcm%(kNmcmRob*kNmcmRob);

  if (robid%kNcolRob == 0) {

    if (mcmid%kNmcmRob == 0) {
      fMCM->SetColRange(18*0-1,18*1-1+2+1);
    }
    if (mcmid%kNmcmRob == 1) {
      fMCM->SetColRange(18*1-1,18*2-1+2+1);
    }
    if (mcmid%kNmcmRob == 2) {
      fMCM->SetColRange(18*2-1,18*3-1+2+1);
    }
    if (mcmid%kNmcmRob == 3) {
      fMCM->SetColRange(18*3-1,18*4-1+2+1);
    }

  } 
  else {

    if (mcmid%kNmcmRob == 0) {
      fMCM->SetColRange(18*4-1,18*5-1+2+1);
    }
    if (mcmid%kNmcmRob == 1) {
      fMCM->SetColRange(18*5-1,18*6-1+2+1);
    }
    if (mcmid%kNmcmRob == 2) {
      fMCM->SetColRange(18*6-1,18*7-1+2+1);
    }
    if (mcmid%kNmcmRob == 3) {
      fMCM->SetColRange(18*7-1,18*8-1+2+1);
    }

  } 

  fMCM->SetRow(kNmcmRob*(robid/kNcolRob)+mcmid/kNmcmRob);

}

//_____________________________________________________________________________
Bool_t AliTRDtrigger::TestTracklet(Int_t det, Int_t row, Int_t seed, Int_t n)
{
  //
  // Check first the tracklet pt
  //

  Int_t nTimeTotal  = fCalib->GetNumberOfTimeBins();

  // Calibration fill 2D
  AliTRDCalibra *calibra = AliTRDCalibra::Instance();
  if (!calibra) {
    AliInfo("Could not get Calibra instance\n");
  }

  fTrkTest->Reset();

  fTrkTest->SetDetector(det);
  fTrkTest->SetRow(row);
  fTrkTest->SetN(n);

  Int_t iCol, iCol1, iCol2, track[3];
  iCol = fMCM->GetSeedCol()[seed];  // 0....20 (MCM)
  fMCM->GetColRange(iCol1,iCol2);   // range in the pad plane
	    
  Float_t amp[3];
  for (Int_t iTime = 0; iTime < nTimeTotal; iTime++) {

    amp[0] = fMCM->GetADC(iCol-1,iTime);
    amp[1] = fMCM->GetADC(iCol  ,iTime);
    amp[2] = fMCM->GetADC(iCol+1,iTime);

    // extract track contribution only from the central pad
    track[0] = fTrack0->GetDataUnchecked(row,iCol+iCol1,iTime);
    track[1] = fTrack1->GetDataUnchecked(row,iCol+iCol1,iTime);
    track[2] = fTrack2->GetDataUnchecked(row,iCol+iCol1,iTime);

    if      (fMCM->IsCluster(iCol,iTime)) {

      fTrkTest->AddCluster(iCol+iCol1,iTime,amp,track);

    } 
    else if ((iCol+1+1) < kMcmCol) {

      amp[0] = fMCM->GetADC(iCol-1+1,iTime);
      amp[1] = fMCM->GetADC(iCol  +1,iTime);
      amp[2] = fMCM->GetADC(iCol+1+1,iTime);

      if (fMCM->IsCluster(iCol+1,iTime)) {

	// extract track contribution only from the central pad
	track[0] = fTrack0->GetDataUnchecked(row,iCol+1+iCol1,iTime);
	track[1] = fTrack1->GetDataUnchecked(row,iCol+1+iCol1,iTime);
	track[2] = fTrack2->GetDataUnchecked(row,iCol+1+iCol1,iTime);

	fTrkTest->AddCluster(iCol+1+iCol1,iTime,amp,track);

      }

    } 

  }

  fTrkTest->CookLabel(0.8);  
  /*
  if (fTrkTest->GetLabel() >= fNPrimary) {
    Info("AddTracklet","Only primaries are stored!");
    return;
  }
  */
  // LTU Pt cut
  fTrkTest->MakeTrackletGraph(fGeo,fField);

  // TRD Online calibration
  if (calibra->Getmcmtracking()) {
    calibra->UpdateHistogramcm(fTrkTest);
  }

  fTrkTest->MakeClusAmpGraph();

  if (TMath::Abs(fTrkTest->GetPt()) < fTrigParam->GetLtuPtCut()) {
    return kFALSE;
  }
  
  return kTRUE;  

}

//_____________________________________________________________________________
void AliTRDtrigger::AddTracklet(Int_t det, Int_t row, Int_t seed, Int_t n)
{
  //
  // Add a found tracklet
  //

  Int_t nTimeTotal  = fCalib->GetNumberOfTimeBins();

  fTrk = new AliTRDmcmTracklet(det,row,n);

  Int_t iCol, iCol1, iCol2, track[3];
  iCol = fMCM->GetSeedCol()[seed];  // 0....20 (MCM)
  fMCM->GetColRange(iCol1,iCol2);   // range in the pad plane
	    
  Float_t amp[3];
  for (Int_t iTime = 0; iTime < nTimeTotal; iTime++) {

    amp[0] = fMCM->GetADC(iCol-1,iTime);
    amp[1] = fMCM->GetADC(iCol  ,iTime);
    amp[2] = fMCM->GetADC(iCol+1,iTime);

    // extract track contribution only from the central pad
    track[0] = fTrack0->GetDataUnchecked(row,iCol+iCol1,iTime);
    track[1] = fTrack1->GetDataUnchecked(row,iCol+iCol1,iTime);
    track[2] = fTrack2->GetDataUnchecked(row,iCol+iCol1,iTime);

    if      (fMCM->IsCluster(iCol,iTime)) {

      fTrk->AddCluster(iCol+iCol1,iTime,amp,track);

    } 
    else if ((iCol+1+1) < kMcmCol) {

      amp[0] = fMCM->GetADC(iCol-1+1,iTime);
      amp[1] = fMCM->GetADC(iCol  +1,iTime);
      amp[2] = fMCM->GetADC(iCol+1+1,iTime);

      if (fMCM->IsCluster(iCol+1,iTime)) {

	// extract track contribution only from the central pad
	track[0] = fTrack0->GetDataUnchecked(row,iCol+1+iCol1,iTime);
	track[1] = fTrack1->GetDataUnchecked(row,iCol+1+iCol1,iTime);
	track[2] = fTrack2->GetDataUnchecked(row,iCol+1+iCol1,iTime);

	fTrk->AddCluster(iCol+1+iCol1,iTime,amp,track);

      }

    }

  }

  fTrk->CookLabel(0.8);  
  /*
  if (fTrk->GetLabel() >= fNPrimary) {
    Info("AddTracklet","Only primaries are stored!");
    return;
  }
  */
  // LTU Pt cut
  fTrk->MakeTrackletGraph(fGeo,fField);
  fTrk->MakeClusAmpGraph();
  if (TMath::Abs(fTrk->GetPt()) < fTrigParam->GetLtuPtCut()) {
    return;
  }
      
  Tracklets()->Add(fTrk);

}

//_____________________________________________________________________________
Bool_t AliTRDtrigger::WriteTracklets(Int_t det) 
{
  //
  // Fills TRDmcmTracklet branch in the tree with the Tracklets 
  // found in detector = det. For det=-1 writes the tree. 
  //

  if ((det < -1) || (det >= AliTRDgeometry::Ndet())) {
    AliError(Form("Unexpected detector index %d.",det));
    return kFALSE;
  }

  TBranch *branch = fTrackletTree->GetBranch("TRDmcmTracklet");
  if (!branch) {
    TObjArray *ioArray = 0;
    branch = fTrackletTree->Branch("TRDmcmTracklet","TObjArray",&ioArray,32000,0);
  }

  if ((det >= 0) && (det < AliTRDgeometry::Ndet())) {

    Int_t nTracklets = Tracklets()->GetEntriesFast();
    TObjArray *detTracklets = new TObjArray(400);

    for (Int_t i = 0; i < nTracklets; i++) {

      AliTRDmcmTracklet *trk = (AliTRDmcmTracklet *) Tracklets()->UncheckedAt(i);
      
      if (det == trk->GetDetector()) {
        detTracklets->AddLast(trk);
      }

    }

    branch->SetAddress(&detTracklets);
    fTrackletTree->Fill();

    delete detTracklets;

    return kTRUE;

  }

  if (det == -1) {

    AliInfo(Form("Writing the Tracklet tree %s for event %d."
	        ,fTrackletTree->GetName(),fRunLoader->GetEventNumber()));

    AliLoader* loader = fRunLoader->GetLoader("TRDLoader");
    loader->WriteTracks("OVERWRITE");
    
    return kTRUE;  

  }

  return kFALSE;

}

//_____________________________________________________________________________
void AliTRDtrigger::MakeTracks(Int_t det)
{
  //
  // Create GTU tracks per module (stack of 6 chambers)
  //
  
  fModule->Reset();

  Int_t nRowMax, iplan, icham, isect, row;

  if ((det < 0) || (det >= AliTRDgeometry::Ndet())) {
    AliError(Form("Unexpected detector index %d.",det));
    return;
  }
  
  Int_t nTracklets = Tracklets()->GetEntriesFast();
  
  AliTRDmcmTracklet *trk;
  for (Int_t i = 0; i < nTracklets; i++) {
    
    trk = (AliTRDmcmTracklet *) Tracklets()->UncheckedAt(i);
    
    iplan = fGeo->GetPlane(trk->GetDetector());
    icham = fGeo->GetChamber(trk->GetDetector());
    isect = fGeo->GetSector(trk->GetDetector());

    nRowMax = fCParam->GetRowMax(iplan,icham,isect);
    row = trk->GetRow();

    fModule->AddTracklet(trk->GetDetector(),
			 row,
			 trk->GetRowz(),
			 trk->GetSlope(),
			 trk->GetOffset(),
			 trk->GetTime0(),
			 trk->GetNclusters(),
			 trk->GetLabel(),
			 trk->GetdQdl());
    
  }

  fModule->SortTracklets();
  fModule->RemoveMultipleTracklets();
  fModule->SortZ((Int_t)fGeo->GetChamber(det));
  fModule->FindTracks();
  fModule->SortTracks();
  fModule->RemoveMultipleTracks();

  Int_t nModTracks = fModule->GetNtracks();
  AliTRDgtuTrack *gtutrk;
  for (Int_t i = 0; i < nModTracks; i++) {
    gtutrk = (AliTRDgtuTrack*)fModule->GetTrack(i);
    if (TMath::Abs(gtutrk->GetPt()) < fTrigParam->GetGtuPtCut()) continue;
    gtutrk->CookLabel();
    gtutrk->MakePID();
    AddTrack(gtutrk,det);
  }
  
}

//_____________________________________________________________________________
void AliTRDtrigger::AddTrack(const AliTRDgtuTrack *t, Int_t det)  
{ 
  //
  // Add a track to the list
  //

  AliTRDgtuTrack *track = new(fTracks->operator[](fTracks->GetEntriesFast())) 
                          AliTRDgtuTrack(*t);
  track->SetDetector(det); 

}

//_____________________________________________________________________________
TObjArray* AliTRDtrigger::Tracklets()
{ 
  //
  // Returns list of tracklets
  //

  if (!fTracklets) {
    fTracklets = new TObjArray(400); 
  }
  return fTracklets;       

}

//_____________________________________________________________________________
void AliTRDtrigger::ResetTracklets()                              
{
  //
  // Resets the list of tracklets
  //
 
  if (fTracklets) {
    fTracklets->Delete();             
  }

}

//_____________________________________________________________________________
Int_t AliTRDtrigger::GetNumberOfTracks() const                     
{ 
  //
  // Returns number of tracks
  //

  return fTracks->GetEntriesFast(); 
                
}

//_____________________________________________________________________________
AliTRDgtuTrack* AliTRDtrigger::GetTrack(Int_t i) const               
{ 
  //
  // Returns a given track from the list
  //

  return (AliTRDgtuTrack *) fTracks->UncheckedAt(i); 

}
