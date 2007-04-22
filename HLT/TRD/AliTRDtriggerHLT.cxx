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
//#include "AliLoader.h"

#include "AliTRDdigitsManager.h"
#include "AliTRDgeometry.h"
#include "AliTRDdataArrayI.h"
#include "AliTRDcalibDB.h"
#include "AliTRDCommonParam.h"
#include "AliTRDrawData.h"
#include "AliTRDmodule.h"
#include "AliTRDmcmTracklet.h"
#include "AliTRDgtuTrack.h"
#include "AliTRDtrigParam.h"
#include "AliTRDmcm.h"
#include "AliTRDzmaps.h"
#include "AliHLTTRDCalibra.h"
#include "Cal/AliTRDCalPIDLQ.h"

//#include "AliTRDtrigger.h"
#include "AliTRDtriggerHLT.h"

ClassImp(AliTRDtriggerHLT)

//_____________________________________________________________________________
AliTRDtriggerHLT::AliTRDtriggerHLT()
  :AliTRDtrigger(),
  fTreeCreatedHere(kFALSE)
{
  //
  // AliTRDtriggerHLT default constructor
  //  
}

//_____________________________________________________________________________
AliTRDtriggerHLT::AliTRDtriggerHLT(const Text_t *name, const Text_t *title)
  :AliTRDtrigger(name,title),
  fTreeCreatedHere(kFALSE)   
{
  //
  // AliTRDtriggerHLT constructor
  //
}

//_____________________________________________________________________________
AliTRDtriggerHLT::AliTRDtriggerHLT(const AliTRDtriggerHLT &p)
  :AliTRDtrigger(p),
  fTreeCreatedHere(kFALSE)   
{
  //
  // AliTRDtriggerHLT copy constructor
  //
}

///_____________________________________________________________________________
AliTRDtriggerHLT::~AliTRDtriggerHLT()
{
  //
  // AliTRDtrigger destructor
  //

  //   if (fTracklets) {
  //     fTracklets->Delete();
  //     delete fTracklets;
  //   }
  
  //   if (fTracks) {
  //     fTracks->Delete();
  //     delete fTracks;
  //   }

  if (fGeo)
    delete fGeo;
}

//_____________________________________________________________________________
AliTRDtriggerHLT &AliTRDtriggerHLT::operator=(const AliTRDtriggerHLT &p)
{
  //
  // Assignment operator
  //

  if (this != &p) ((AliTRDtriggerHLT &) p).Copy(*this);
  return *this;

}

//_____________________________________________________________________________
void AliTRDtriggerHLT::Copy(TObject &) const
{
  //
  // Copy function
  //

  AliFatal("Not implemented");

}

//_____________________________________________________________________________
void AliTRDtriggerHLT::Init()
{
  //
  // Init: module, field, geometry, calibDB
  //

  fModule = new AliTRDmodule(fTrigParam); 
  fTracks->Clear();
  fField  = fTrigParam->GetField();

  //AliInfo("Recreating default geometry!");
  AliDebug(1, Form("Recreating default geometry!"));
  //Modified by MP:
  //fGeo    = (AliTRDgeometry*)AliTRDgeometry::GetGeometry(fRunLoader);
  fGeo    = new AliTRDgeometry();

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
Bool_t AliTRDtriggerHLT::ReadDigits(AliRawReaderMemory* rawReader)
{
  //
  // Reads the digits arrays from the ddl file
  //

  AliTRDrawData raw;
  fDigitsManager = raw.Raw2Digits((AliRawReader*)rawReader);
  //AliInfo(Form("Digits manager at 0x%x", fDigitsManager));
  AliDebug(1, Form("Digits manager at 0x%x", fDigitsManager));
  if (fDigitsManager)
    return kTRUE;
  else
    return kFALSE;
}

//_____________________________________________________________________________
Bool_t AliTRDtriggerHLT::MakeTracklets(Bool_t makeTracks)
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

		if (fTrigParam->GetDebugLevel()   > 1) 
		  { 
		    //AliInfo(Form("Add tracklet %d in col %02d \n",fNtracklets,fMCM->GetSeedCol()[iSeed]));
		    AliDebug(5, Form("Add tracklet %d in col %02d \n",fNtracklets,fMCM->GetSeedCol()[iSeed]));
		  }

		if (TestTracklet(idet,row,iSeed,0)) {
		  AddTracklet(idet,row,iSeed,fNtracklets++);
		  //AliInfo(Form("TESTED OK tracklet %d in col %02d \n",fNtracklets,fMCM->GetSeedCol()[iSeed]));
		  AliDebug(5, Form("TESTED OK tracklet %d in col %02d \n",fNtracklets,fMCM->GetSeedCol()[iSeed]));
		}
		else
		  {
		    //AliInfo(Form("REJECTED tracklet %d in col %02d \n",fNtracklets,fMCM->GetSeedCol()[iSeed]));
		    AliDebug(5, Form("REJECTED tracklet %d in col %02d \n",fNtracklets,fMCM->GetSeedCol()[iSeed]));
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

	//WriteTracklets(idet);
	TreeTracklets(idet);

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
Bool_t AliTRDtriggerHLT::ResetTree()
{
  //
  // Recreate the cluster tree
  //


  //   if (fTrackletTree != 0)
  //     fTrackletTree->Reset();
  // well we'd better delete the whole tree and branches
  delete fTrackletTree;
  fTrackletTree = NULL;
  fTrackletTree = new TTree("TRDclusterTree", "TRDclusterTree");
  if (fTrackletTree)
    {
      fTreeCreatedHere = kTRUE;
      //AliInfo("Tree Reset Successful");
      AliDebug(1, Form("Tree Reset Successful"));
    }
  else
    {
      fTreeCreatedHere = kFALSE;
      AliError("Reset Tree Error!\n");
    }
  
  return fTreeCreatedHere;
}

//_____________________________________________________________________________
Bool_t AliTRDtriggerHLT::TreeTracklets(Int_t idet)
{
  //
  // Failsafe ::WriteTracklets - recreate the cluster tree if necessary
  //

  if (fTrackletTree == 0)
    {
      fTrackletTree = new TTree("TRDtrackletTree", "TRDtrackletTree");
      fTreeCreatedHere = kTRUE;
    }
  
  if (idet < 0)
    {
      AliError("Attempt to write tracklet tree on HLT - bad idea!");
      AliError("HLT does not write any data here!\n");    
      return kFALSE;
    }
  else
    return WriteTracklets(idet);
}

//_____________________________________________________________________________
Bool_t AliTRDtriggerHLT::TestTracklet(Int_t det, Int_t row, Int_t seed, Int_t n)
{
  //
  // Check first the tracklet pt
  //

  Int_t nTimeTotal  = fCalib->GetNumberOfTimeBins();

  // Calibration fill 2D
  AliHLTTRDCalibra *calibra = AliHLTTRDCalibra::Instance();
  if (!calibra) 
    {
      AliError("Could not get Calibra instance\n");
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
  //AliInfo(Form("Geom at 0x%x Field = %f", fGeo, fField));
  fTrkTest->MakeTrackletGraph(fGeo,fField);

  // TRD Online calibration
  //if (calibra->GetMcmTracking()) {
  //AliInfo(Form("Updating tracklet histos : %d", calibra->GetMcmTracking()));
  //calibra->UpdateHistogramcm(fTrkTest);
  //}

  //AliInfo("Updating tracklet histos");
  AliDebug(2, Form("Updating tracklet histos"));
  calibra->UpdateHistogramcm(fTrkTest);

  fTrkTest->MakeClusAmpGraph();
  //AliInfo("MakeClusAmpGraph Done.");

  if (TMath::Abs(fTrkTest->GetPt()) < fTrigParam->GetLtuPtCut()) 
    {
      //AliInfo(Form("Cut tracklet - pT is %f", TMath::Abs(fTrkTest->GetPt())));
      AliDebug(5, Form("Cut tracklet - pT is %f", TMath::Abs(fTrkTest->GetPt())));
      return kFALSE;
    }
  else
    {
      //AliInfo(Form("Tracklet passed - pT is %f", TMath::Abs(fTrkTest->GetPt())));
      AliDebug(5, Form("Tracklet passed - pT is %f", TMath::Abs(fTrkTest->GetPt())));
    }
  return kTRUE;  

}
