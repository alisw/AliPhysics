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

//Class to calculate the intrinsic efficiency of the detection elements of the
//MUON tracking chambers in function of the position in the detection element.
//WOrk on ESD only
//Author:  Nicolas LE BRIS - SUBATECH Nantes
// Modified by Matthieu LENHARDT - SUBATECH Nantes


//ROOT includes
#include <TChain.h>
#include <TFile.h>
#include <TList.h>
#include <TH2F.h>

//ANALYSIS includes
#include "AliAnalysisManager.h"

//STEER includes
#include "AliESDEvent.h"
#include "AliESDMuonTrack.h"
#include "AliESDInputHandler.h"
#include "AliMagF.h"
#include "AliTracker.h"
#include "AliAnalysisManager.h"
#include "AliGeomManager.h"
#include "AliCDBManager.h"


//PWG3/muondep includes
#include "AliAnalysisTaskMuonTrackingEff.h"
#include "AliCheckMuonDetEltResponse.h"

//MUON includes
#include "AliMUONTrackExtrap.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMUONConstants.h"

ClassImp(AliAnalysisTaskMuonTrackingEff)

const Int_t AliAnalysisTaskMuonTrackingEff::fgkTotNbrOfDetectionElt  = 156;
const Int_t AliAnalysisTaskMuonTrackingEff::fgkTotNbrOfChamber  = AliMUONConstants::NTrackingCh();

//________________________________________________________________________
AliAnalysisTaskMuonTrackingEff::AliAnalysisTaskMuonTrackingEff()
  :
  AliAnalysisTask(),
  fIsInit(kFALSE),
  fIsLoaded(kFALSE),
  fOCDBpath(0),
  fUsableTracks(0),
  fESD(0x0),
  fTransformer(0x0),
  fDetEltTDHistList(0),
  fDetEltTTHistList(0),
  fChamberTDHistList(0),
  fChamberTTHistList(0),
  fChamberEff(0x0)
{
/// Default constructor
}
//________________________________________________________________________
AliAnalysisTaskMuonTrackingEff::AliAnalysisTaskMuonTrackingEff(const AliAnalysisTaskMuonTrackingEff& src)
  :
  AliAnalysisTask(src),
  fIsInit(kFALSE),
  fIsLoaded(kFALSE),
  fOCDBpath(0),
  fUsableTracks(0),
  fESD(0x0),
  fTransformer(0x0),
  fDetEltTDHistList(0),
  fDetEltTTHistList(0),
  fChamberTDHistList(0),
  fChamberTTHistList(0),
  fChamberEff(0x0)
{
  /// copy ctor
  src.Copy(*this);
}
//________________________________________________________________________
AliAnalysisTaskMuonTrackingEff& AliAnalysisTaskMuonTrackingEff::operator=(const AliAnalysisTaskMuonTrackingEff& src)
{
  /// assignement operator
  if ( this != &src ) 
  {
    src.Copy(*this);
  }
  return *this;
}

//________________________________________________________________________
AliAnalysisTaskMuonTrackingEff::AliAnalysisTaskMuonTrackingEff(TString name, TString path)
  :
  AliAnalysisTask(name, "AnalysisTaskESD"),
  fIsInit(kFALSE),
  fIsLoaded(kFALSE),
  fOCDBpath(path),
  fUsableTracks(0),
  fESD(0x0),
  fTransformer(0x0),
  fDetEltTDHistList(0),
  fDetEltTTHistList(0),
  fChamberTDHistList(0),
  fChamberTTHistList(0),
  fChamberEff(0x0)
{
//Constructor
//-----------
//Define input & output
//---------------------
// -Input slot 0 works with a TChain:
    DefineInput(0, TChain::Class());

// -Output slots 0 to 5 writes into a TClonesArray:
    DefineOutput(0, TList::Class());
    DefineOutput(1, TList::Class());
    DefineOutput(2, TList::Class());
    DefineOutput(3, TList::Class());
}



//______________________________________________________________________________
AliAnalysisTaskMuonTrackingEff::~AliAnalysisTaskMuonTrackingEff()
{
// Destructor.
  delete fTransformer;
  delete fDetEltTDHistList;
  delete fDetEltTTHistList;
  delete fChamberTDHistList;
  delete fChamberTTHistList;
  delete fChamberEff;
}



//________________________________________________________________________
void AliAnalysisTaskMuonTrackingEff::CreateOutputObjects()
{
  if (!fIsInit)
    Init();

  OpenFile(0);
}



//________________________________________________________________________
void AliAnalysisTaskMuonTrackingEff::ConnectInputData(Option_t */*option*/)

{
  //Connect input
    AliESDInputHandler* esdHandler = (AliESDInputHandler*) ((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler());
    if (!esdHandler)
      Printf("ERROR: Could not get ESDInputHandler");
     
    else fESD = esdHandler->GetEvent();
}


//________________________________________________________________________
void AliAnalysisTaskMuonTrackingEff::Exec(Option_t */*option*/)
{
//Execute analysis for current event

  // Initialize the object first
  if (!fIsInit)
    Init();

  // Then load the OCDB and Geometry
  if (!fIsLoaded)
    LoadOCDBandGeometry();

  // Fill the histograms
  if (fChamberEff == 0x0)
    fChamberEff = new AliCheckMuonDetEltResponse(fTransformer, fESD, fDetEltTDHistList, fDetEltTTHistList, fChamberTDHistList, fChamberTTHistList);
  fChamberEff->CheckDetEltResponse();
  fUsableTracks += fChamberEff -> GetNbrUsableTracks();

//Post the output data:
    PostData(0, fDetEltTDHistList);  
    PostData(1, fDetEltTTHistList);  
    PostData(2, fChamberTDHistList);
    PostData(3, fChamberTTHistList);
}



//________________________________________________________________________
void AliAnalysisTaskMuonTrackingEff::Terminate(Option_t */*option*/)
{
//Terminate analysis
  fDetEltTDHistList = (TList*) GetOutputData(0);
  fDetEltTTHistList = (TList*) GetOutputData(1);
  fChamberTDHistList = (TList*) GetOutputData(2);
  fChamberTTHistList = (TList*) GetOutputData(3);
}



//________________________________________________________________________
void AliAnalysisTaskMuonTrackingEff::Init()
{
  //Define detection element efficiency histograms
  //----------------------------------------------
  fUsableTracks = 0;


  fDetEltTDHistList  = new TList();
  fDetEltTTHistList  = new TList();

  fChamberTDHistList = new TList();
  fChamberTTHistList = new TList();
  
  TString histName; 
  TString histTitle;

  
  for (Int_t ii = 0; ii < fgkTotNbrOfDetectionElt; ++ii)
    {
      Int_t iDetElt = 0;
      if      (ii<16) iDetElt = 100*(ii/4+1)       + ii -  int(ii/4)*4;
      else if (ii<52) iDetElt = 100*((ii-16)/18+5) + ii - (int((ii-16)/18)*18 + 16);
      if      (ii>51) iDetElt = 100*((ii-52)/26+7) + ii - (int((ii-52)/26)*26 + 52); 

      
      histTitle.Form("detEltNbr %d",iDetElt);
      histName.Form("TD_detEltNbr%d",iDetElt);
      if(ii<16)
	{// Chambers 1 -> 4
	  TH2F *histo = new TH2F(histName, histTitle, 12, -10.0 , 110.0, 12, -10.0, 110.0);
	  fDetEltTDHistList->AddAt(histo, ii);
	}

      else 
	{// Chambers 5 -> 10
	  TH2F *histo = new TH2F(histName, histTitle, 28, -140.0, 140.0, 8, -40.0, 40.0);
	  fDetEltTDHistList->AddAt(histo, ii);	  
	}


      histName.Form("TT_detEltNbr%d",iDetElt);
      if(ii<16)
	{// Chambers 1 -> 4
	  TH2F *histo = new TH2F(histName, histTitle, 12, -10.0 , 110.0, 12, -10.0, 110.0);
	  fDetEltTTHistList->AddAt(histo, ii);
	}

      else 
	{// Chambers 5 -> 10
	  TH2F *histo = new TH2F(histName, histTitle, 28, -140.0, 140.0, 8, -40.0, 40.0);
	  fDetEltTTHistList->AddAt(histo, ii);	  
	}
    }


  for (Int_t jj = 0; jj < fgkTotNbrOfChamber+1; jj++)
    {
      histName.Form("TD_ChamberNbr%d", jj+1);
      histTitle.Form("ChamberNbr %d", jj+1);
      if (jj<4)
	{//Chambers 1 -> 4
	  TH1F *histo = new TH1F(histName, histTitle, 4, 0.0, 4.0);
	  fChamberTDHistList->AddAt(histo, jj);
	}
      else 
	{
	  if (jj<6) 
	    {//Chambers 5 -> 6
	      TH1F *histo = new TH1F(histName, histTitle, 18, 0.0, 18.0);
	      fChamberTDHistList->AddAt(histo, jj);
	    }

	  if (jj>=6 && jj < 10) 
	    {//Chambers 7 -> 10
	      TH1F *histo = new TH1F(histName, histTitle, 26, 0.0, 26.0);
	      fChamberTDHistList->AddAt(histo, jj);
	    }
	}
      
      histName.Form("TT_ChamberNbr%d",jj+1);
      if (jj<4)
	{//Chambers 1 -> 4
	  TH1F *histo = new TH1F(histName, histTitle, 4, 0.0, 4.0);
	  fChamberTTHistList->AddAt(histo, jj);
	}
      else 
	{
	  if (jj<6) 
	    {//Chambers 5 -> 6
	      TH1F *histo = new TH1F(histName, histTitle, 18, 0.0, 18.0);
	      fChamberTTHistList->AddAt(histo, jj);
	    }

	  if (jj>=6 && jj < 10) 
	    {//Chambers 7 -> 10
	      TH1F *histo = new TH1F(histName, histTitle, 26, 0.0, 26.0);
	      fChamberTTHistList->AddAt(histo, jj);
	    }
	}
      
      if (jj == 10)
	{
	  histName.Form("TD_Chambers", jj+1);
	  histTitle.Form("Chambers", jj+1);
	  TH1F *histoTDchambers = new TH1F(histName, histTitle, 10, 0.5, 10.5);
	  fChamberTDHistList->AddAt(histoTDchambers, jj);
	  
	  histName.Form("TT_Chambers", jj+1);
	  TH1F *histoTTchambers = new TH1F(histName, histTitle, 10, 0.5, 10.5);
	  fChamberTTHistList->AddAt(histoTTchambers, jj);
	}

    }
  
  fIsInit = kTRUE;
}


//________________________________________________________________________
void AliAnalysisTaskMuonTrackingEff::LoadOCDBandGeometry()
{
  // Load the OCDB and the Geometry

  // OCDB
  AliCDBManager* man = AliCDBManager::Instance();
  TString ocdbPath= "alien://folder=/alice/data/2010/OCDB";
  man->SetDefaultStorage(ocdbPath.Data());
  man->SetRun(fESD->GetRunNumber());

  //Geometry
  if (!AliGeomManager::GetGeometry())
    AliGeomManager::LoadGeometry();
  fTransformer = new AliMUONGeometryTransformer();
  fTransformer->LoadGeometryData();

  //Set Field Map for track extrapolation
  AliMUONTrackExtrap::SetField();

  fIsLoaded = kTRUE;
}

