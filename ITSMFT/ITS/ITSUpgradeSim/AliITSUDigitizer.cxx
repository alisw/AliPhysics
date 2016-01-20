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
 
/* $Id: AliITSUDigitizer.cxx 52261 2011-10-23 15:46:57Z hristov $ */
///////////////////////////////////////////////////////////////////////////
//Piotr.Skowronski@cern.ch :                                             //
//Corrections applied in order to compile (only)                         // 
//   with new I/O and folder structure                                   //
//To be implemented correctly by responsible                             //
//                                                                       //
//  Class used to steer                                                  //
//  the digitization for ITS                                             //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include <stdlib.h>
#include <TClonesArray.h>
#include <TTree.h>
#include <TBranch.h>

#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliLoader.h"
#include "AliLog.h"
#include "AliDigitizationInput.h"
#include "AliITSUDigitizer.h"
#include "AliITSUGeomTGeo.h"
#include "AliITSMFTSimulation.h"
#include "AliITSMFTSDigit.h"

ClassImp(AliITSUDigitizer)

//______________________________________________________________________
AliITSUDigitizer::AliITSUDigitizer() 
:  fITS(0)
  ,fModActive(0)
  ,fInit(kFALSE)
  ,fRoif(-1)
  ,fRoiifile(0)
  ,fFlagFirstEv(kTRUE)
{
  // Default constructor.
}

//______________________________________________________________________
AliITSUDigitizer::AliITSUDigitizer(AliDigitizationInput* digInp) 
  :AliDigitizer(digInp)
  ,fITS(0)
  ,fModActive(0)
  ,fInit(kFALSE)
  ,fRoif(-1)
  ,fRoiifile(0)
  ,fFlagFirstEv(kTRUE)
{
  // Standard constructor.
}

//______________________________________________________________________
AliITSUDigitizer::~AliITSUDigitizer()
{
  // destructor. 
  fITS = 0; // don't delete fITS. Done else where.
  delete[] fModActive;
}

//______________________________________________________________________
Bool_t AliITSUDigitizer::Init()
{
  // Initialization. Set up region of interest, if switched on, and loads ITS and ITSgeom.
  //
  if (fInit) return kTRUE;
  //
  fInit = kTRUE; // Assume for now init will work.
  //
  if(!gAlice) {
    fITS      = 0;
    fRoiifile = 0;
    fInit     = kFALSE;
    AliFatal("gAlice not found");
  } // end if
  //
  fITS = (AliITSU *)(gAlice->GetDetector("ITS"));
  if(!fITS){
    fRoiifile = 0;
    fInit     = kFALSE;
    AliFatal("ITS not found");
  } 
  if (!fITS->IsSimInitDone()) fITS->InitSimulation();
  int nm = fITS->GetITSGeomTGeo()->GetNChips();
  fModActive = new Bool_t[nm];
  for (Int_t i=nm;i--;) fModActive[i] = kTRUE;

  return fInit;
}

//______________________________________________________________________
void AliITSUDigitizer::Digitize(Option_t* /*opt*/)
{
  // Main digitization function. 
  //
  if (!fInit) AliFatal("Init not successful, aborting.");
  //
  Int_t nfiles = GetDigInput()->GetNinputs();
  Int_t event  = GetDigInput()->GetOutputEventNr();
  //
  TString loadname = Form("%sLoader",fITS->GetName());
  //
  AliITSUGeomTGeo* geom = fITS->GetITSGeomTGeo();
  Int_t nChips  = geom->GetNChips();
  Bool_t lmod;
  Int_t *fl = new Int_t[nfiles];
  fl[0] = fRoiifile;
  int mask = 1;
  for (int id=0;id<nfiles;id++) if(id!=fRoiifile) fl[mask++] = id;
  //
  TClonesArray * sdig = new TClonesArray( "AliITSMFTSDigit",1000 );
  //
  AliRunLoader *inRL = 0x0, *outRL = 0x0;
  AliLoader *ingime = 0x0, *outgime = 0x0;    
  //
  outRL = AliRunLoader::GetRunLoader(fDigInput->GetOutputFolderName());    
  if (!outRL) AliFatal("Can not get Output Run Loader");
  //
  outRL->GetEvent(event);
  outgime = outRL->GetLoader(loadname);
  if ( outgime == 0x0) AliFatal("Can not get Output ITS Loader");
  //
  outgime->LoadDigits("update");
  if (outgime->TreeD() == 0x0) outgime->MakeTree("D");
  //
  // Digitize
  fITS->MakeBranchInTreeD(outgime->TreeD());
  if(fRoif!=0) AliDebug(1,"Region of Interest digitization selected");
  else         AliDebug(1,"No Region of Interest selected. Digitizing everything");
  //
  for (int ifiles=0; ifiles<nfiles; ifiles++ ) {
    inRL =  AliRunLoader::GetRunLoader(fDigInput->GetInputFolderName(fl[ifiles]));
    ingime = inRL->GetLoader(loadname);
    if (ingime->TreeS() == 0x0) ingime->LoadSDigits();
  }
  //
  for (int chip=0; chip<nChips; chip++ ) {
    //
    if (!fRoif && !fModActive[chip]) continue;
    int lr = geom->GetLayer(chip);
    AliITSMFTSimulation *sim = fITS->GetSimulationModel(lr);
    if (!sim) AliFatal(Form("The simulation model for layer %d is not available",lr));
    //
    // Fill the chip with the sum of SDigits
    sim->InitSimulationChip((AliITSMFTChip*)(fITS->GetChip(chip)), event, fITS->GetSegmentation(lr), fITS->GetResponseParam(lr));
    //
    for (int ifiles=0; ifiles<nfiles; ifiles++ ) {
      //
      if (!fRoif && !fModActive[chip]) continue;
      inRL =  AliRunLoader::GetRunLoader(fDigInput->GetInputFolderName(fl[ifiles]));
      ingime = inRL->GetLoader(loadname);
      //
      TTree *treeS = ingime->TreeS();
      fITS->SetTreeAddress();
      //
      if( !treeS  ) continue; 
      TBranch *brchSDigits = treeS->GetBranch( fITS->GetName() );
      if( brchSDigits )	brchSDigits->SetAddress( &sdig ); 
      else {
	AliError(Form("branch ITS not found in TreeS, input file %d ", ifiles));
	delete [] fl;
	return;
      } 
      //
      sdig->Clear();
      mask = GetDigInput()->GetMask(ifiles);
      brchSDigits->GetEvent( chip );
      lmod = sim->AddSDigitsToChip(sdig,mask);
      if(GetRegionOfInterest() && !ifiles) fModActive[chip] = lmod;
      //
    } 
    // Digitize current chip sum(SDigits)->Digits
    sim->FinishSDigitiseChip(fITS->GetDigits());
    //
    outgime->TreeD()->Fill();       // fills all branches - wasted disk space
    fITS->ResetDigits();
  } // end for chip
  //
  //  fITS->WriteFOSignals(); 
  outgime->TreeD()->AutoSave();
  outgime->WriteDigits("OVERWRITE");
  outgime->UnloadDigits();
  for(int ifiles=0; ifiles<nfiles; ifiles++ ) {
    inRL =  AliRunLoader::GetRunLoader(fDigInput->GetInputFolderName(fl[ifiles]));
    ingime = inRL->GetLoader(loadname);
    ingime->UnloadSDigits();
  }
  //
  delete[] fl;
  sdig->Clear();
  delete sdig;
  for (Int_t i=nChips;i--;) fModActive[i] = kTRUE;
  //
  return;
}

//______________________________________________________________________
void AliITSUDigitizer::SetByRegionOfInterest(TTree *ts)
{
  // Scans through the ITS branch of the SDigits tree, ts, for chips
  // which have SDigits in them. For these chips, a flag is set to
  // digitize only these chips. The value of fRoif determines how many
  // neighboring chips will also be turned on. fRoif=0 will turn on only
  // those chips with SDigits in them. fRoif=1 will turn on, in addition,
  // those chips that are +-1 chip from the one with the SDigits. And
  // So on. This last feature is not supported yet.
  // Inputs:
  //      TTree *ts  The tree in which the existing SDigits will define the
  //                 region of interest.
  if (fRoif==0) return;
  if (ts==0)    return;
  TBranch *brchSDigits = ts->GetBranch(fITS->GetName());
  TClonesArray * sdig = new TClonesArray( "AliITSMFTSDigit",1000 );
  //
  if( brchSDigits ) brchSDigits->SetAddress( &sdig );
  else  {AliError("Branch ITS not found in TreeS"); return;}
  //
  int nm = fITS->GetITSGeomTGeo()->GetNChips();
  for (int m=0;m<nm;m++) {
    fModActive[m] = kFALSE; // Not active by default
    sdig->Clear();
    brchSDigits->GetEvent(m);
    int ndig = sdig->GetEntries();
    for(int i=0;i<ndig;i++) {
      // activate the necessary chips
      if ( ((AliITSMFTSDigit*)sdig->At(m))->GetSumSignal()>0.0 ) { // Must have non zero signal.
	fModActive[m] = kTRUE;
	break;
      } // end if
    } // end if. end for i.
  } // end for m
  AliDebug(1,"Digitization by Region of Interest selected");
  sdig->Clear();
  delete sdig;
  return;
}
