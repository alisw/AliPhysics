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

/* $Id$ */
 
//------------------------------------------------------------------------
// Jet reader base class
// manages the reading of input for jet algorithms
// Authors: jgcn@mda.cinvestav.mx
//          magali.estienne@subatech.in2p3.fr
//          alexandre.shabetai@cern.ch
//
// **February 2011
// implemented  standard geometry (AliEMCALGeometry) (was AliJetDummyGeo implented separately in ESDReader and AODReader
// local2master matrices are now get from $ALICE_ROOT/OADB/PWG4/JetReconstruction/EMCALlocal2master.root
// you can choose the geometry (EMCAL_COMPLETE, EMCAL_FIRSTYEARv1, etc) via SetEMCALgeo2bLoad('Name_of_Geometry') in the Readerheader
// different options for survey(ed) matrice are provided too
// ** August 2011
// OADB path changed from  '../OADB/PWG4/JetReconstruction/'  to   '../OADB/EMCAL/'
// marco.bregant@subatech.in2p3.fr
// ** 2011
// - AliJetESD/AODReader classes removed from JETAN. Reader now ESD/AOD independent. It uses VEvent in the AliJetFill* classes.
// - EMCal utilities added for bad cells id, calibration, etc.
//------------------------------------------------------------------------

// root
#include <TFile.h>

//AliRoot
#include "AliVEvent.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliJetReader.h"
#include "AliJetReaderHeader.h"
#include "AliAnalysisManager.h"
#include "AliJetFillCalTrkTrack.h" 
#include "AliJetFillCalTrkTrackMC.h"
#include "AliJetCalTrk.h"

using std::cout;
using std::endl;
ClassImp(AliJetReader)

////////////////////////////////////////////////////////////////////////


AliJetReader::AliJetReader():
  fCalTrkEvent(0x0),
  fFillEvent(0x0),
  fReaderHeader(0),
  fFillEventwTrks(0x0), 
  fDebug(0),
  fVEvent(0x0),
  fMCEvent(0x0),
  fOpt(0)
{
  // Default constructor
}

//-----------------------------------------------------------------------
AliJetReader::~AliJetReader()
{
  // Destructor
  if (fCalTrkEvent) {
    fCalTrkEvent->Delete();
    delete fCalTrkEvent;
  }

  if (fFillEventwTrks) {
    delete fFillEventwTrks;
  }

}

//-----------------------------------------------------------------------
void AliJetReader::InitTasks()
{
  // Initialization
  fOpt = GetReaderHeader()->GetDetector();

  TString datatype = fReaderHeader->GetDataType();
  datatype.ToUpper();
  Bool_t kIsKine = kFALSE;
  if((!datatype.Contains("AOD") && datatype.Contains("MC")) ||
     (!datatype.Contains("AOD") && datatype.Contains("MC2")) ){
    kIsKine = kTRUE;
  } 
  Bool_t kIsHighMult = GetReaderHeader()->GetIsHighMult();
  fCalTrkEvent = new AliJetCalTrkEvent(fOpt,kIsKine,kIsHighMult);

  // Initialize jet analysis
  CreateTasks();

}

//-----------------------------------------------------------------------
Bool_t AliJetReader::ProcessEvent()
{
  // Process one event
  // Charged only or charged+neutral jets

  Bool_t ok = ExecTasks();

  if(!ok) {return kFALSE;}

  return kTRUE;

}

//-----------------------------------------------------------------------
void AliJetReader::SetInputEvent(const TObject* esd, const TObject* aod, const AliMCEvent* mc)
{
  // set input event pointers
  if( fReaderHeader->GetDataType().Contains("AOD") && aod) {fVEvent = (AliAODEvent*) aod;}
  else if( fReaderHeader->GetDataType().Contains("ESD") && esd) {fVEvent = (AliESDEvent*) esd;}
  else if ( fReaderHeader->GetDataType().Contains("MC") || fReaderHeader->GetDataType().Contains("MC2")) { fMCEvent = (AliMCEvent*) mc;}
  else {printf("No input event ! ");}

}

//-----------------------------------------------------------------------
Bool_t AliJetReader::CreateTasks()
{
  // For reader task initialization

  fDebug = fReaderHeader->GetDebug();

  fFillEvent = new AliJetFillCalTrkEvent();
  if (fOpt>0) {
    // Tracks
    if(fOpt%2==!0 && fOpt!=0){
      fFillEventwTrks = new AliJetFillCalTrkTrack();
      fFillEventwTrks->SetReaderHeader(fReaderHeader);
    }
  }
  else { // MC/Kine cases
    fFillEventwTrks = new AliJetFillCalTrkTrackMC();
    fFillEventwTrks->SetReaderHeader(fReaderHeader);
  }

  if(fDebug>1) cout << "Tasks instantiated at that stage ! " << endl;

  return kTRUE;

}

//-----------------------------------------------------------------------
Bool_t AliJetReader::ExecTasks()
{
  // Main function
  // Fill the reader part
  
  fDebug = fReaderHeader->GetDebug();

  if(!fVEvent && !fMCEvent) {
    return kFALSE;
  }

  // TPC only or Digits+TPC or Clusters+TPC
  if(fOpt%2==!0 || fOpt==0){
    fFillEventwTrks->SetVEvent(fVEvent);
    fFillEventwTrks->SetMCEvent(fMCEvent);
    fFillEventwTrks->SetCalTrkEvent(fCalTrkEvent);    
    fFillEventwTrks->Exec("tpc");
  }

  return kTRUE;

}

//-----------------------------------------------------------------------
void AliJetReader::WriteRHeaderToFile() const
{
  // write reader header
  AliJetReaderHeader *rh = GetReaderHeader();
  rh->Write();

}

//-----------------------------------------------------------------------
void AliJetReader::WriteReaderHeader()
{
  // Write the Headers
  TFile* f = new TFile("jets_local.root", "recreate");
  WriteRHeaderToFile();
  f->Close();

}

