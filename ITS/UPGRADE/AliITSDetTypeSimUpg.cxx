/***************************************************************************
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
 $Id: AliITSDetTypeSimUpg.cxx 56993 2012-06-08 08:24:17Z fca $
*/


/////////////////////////////////////////////////////////////////////
// Base simulation functions for ITS upgrade                       //
//                                                                 //
//                                                                 //
/////////////////////////////////////////////////////////////////////          
#include "TBranch.h"
#include "TClonesArray.h"
#include "TObjArray.h"
#include "TTree.h"

#include "AliRun.h"

#include "AliCDBManager.h"
#include "AliCDBId.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"
#include "AliCDBMetaData.h"
#include "AliITSdigit.h"
#include "AliITSdigitPixUpg.h"
#include "AliITSgeom.h"
#include "AliITSDetTypeSimUpg.h"
#include "AliITSpListItem.h"
#include "AliITSCalibration.h"
#include "AliITSsimulation.h"
#include "AliITSTriggerConditions.h"
#include "AliITSsegmentation.h"
#include "AliITSsegmentationPixUpg.h"
#include "AliITSsimulationPixUpg.h"
#include "AliBaseLoader.h"

using std::endl;
using std::cout;
ClassImp(AliITSDetTypeSimUpg)

const char* AliITSDetTypeSimUpg::fgkDetTypeName[AliITSDetTypeSimUpg::kNDetTypes] = {"PixUpg"};


//----------------------------------------------------------------------
AliITSDetTypeSimUpg::AliITSDetTypeSimUpg() 
{ 
  // Default Constructor
  // Note: the base class AliITSDetTypeSim is badly designed, since its 
  // default c-tor dynamically creates lot of stuff. 
  // To be derived directly from TObject later...
  // Inputs:
  //    none.
  // Outputs:
  //    none.
  // Return:
  //    A properly zero-ed AliITSDetTypeSimUpg class.
  if (kNDetTypes>3) AliFatal("The hack failed: kNDetTypes>fgkNdettypes");
  //
  for (int i=kNDetTypes;i--;) fDetNModules[i] = 0;
}

//----------------------------------------------------------------------
AliITSDetTypeSimUpg::~AliITSDetTypeSimUpg()
{
  // Destructor
  // RS: delegate everything to base class
}

//_______________________________________________________________________
void AliITSDetTypeSimUpg::SetDefaultSegmentation(Int_t idet)
{
  // Set default segmentation model objects
  AliITSsegmentation *seg = 0;
  //  
  if (GetSegmentationModel(idet)) delete GetSegmentation()->RemoveAt(idet);
  //
  switch (idet) 
    {
    case kDetPixUpg: seg = new AliITSsegmentationPixUpg(); break;
    default        : AliFatal(Form("Uknown detector type %d",idet));
  };
  SetSegmentationModel(idet,seg);
  //
}

//_______________________________________________________________________
AliITSCalibration* AliITSDetTypeSimUpg::GetCalibrationModel(Int_t iMod) const 
{
  //Get response model for module number iMod 
  //
  TObjArray* calarr = GetCalibrationArray();
  if (!calarr) {
    AliError("Calibration array is not initialized");
    return 0;
  }
  return (AliITSCalibration*)calarr->At(iMod);
}

//_______________________________________________________________________
void AliITSDetTypeSimUpg::SetDefaults()
{
  //Set defaults for segmentation and response
  
  if(GetITSgeom()==0){
    Warning("SetDefaults","GetITSgeom() is 0!");
    return;
  } // end if
  if (GetCalibrationArray()==0) {
    CreateCalibrationArray();
  } // end if
  //
  ResetSegmentation();
  if(!GetCalibration()){AliFatal("Exit"); exit(0);}
  
  SetDigitClassName(0,"AliITSdigitPixUpg");
  //  
  for(Int_t idet=0;idet<kNDetTypes;idet++){
    if(!GetSegmentationModel(idet)) SetDefaultSegmentation(idet);
  }
  //
}

//______________________________________________________________________
Bool_t AliITSDetTypeSimUpg::GetCalibration() 
{
  // Get Default calibration if a storage is not defined.
  /*
  //
  if(!fFirstcall){
    AliITSCalibration* cal = GetCalibrationModel(0);
    if(cal)return kTRUE;
  }
  else {
    fFirstcall = kFALSE;
  }
  
  SetRunNumber((Int_t)AliCDBManager::Instance()->GetRun());
  Int_t run=GetRunNumber();
  
  Bool_t isCacheActive = kTRUE;
  if (run) {
    isCacheActive=kFALSE;
    fCalibration->SetOwner(kTRUE);
  }
  else{
    fCalibration->SetOwner(kFALSE);
  }
  
  AliCDBManager::Instance()->SetCacheFlag(isCacheActive);
  //
  // TO DO, RS
  */
  return kTRUE;
}

//_______________________________________________________________________
void AliITSDetTypeSimUpg::SetDefaultSimulation()
{
  //Set default simulation for detector type
  //
  if (GetITSgeom()==0) {
    Warning("SetDefaultSimulation","GetITSgeom() is 0!");
    return;
  }
  if (GetCalibrationArray()==0) {
    Warning("SetDefaultSimulation","fCalibration is 0!");
    return;
  }
  if (GetSegmentation()==0) {
    Warning("SetDefaultSimulation","fSegmentation is 0!");
    for (Int_t i=0;i<kNDetTypes;i++) SetDefaultSegmentation(i);
  } else for(Int_t i=0;i<kNDetTypes;i++) if(!GetSegmentationModel(i)){
	Warning("SetDefaultSimulation",
		"Segmentation not defined for det %d - Default taken\n!",i);
	SetDefaultSegmentation(i);
      }
  //
  AliITSsimulation* sim;
  //
  for (Int_t idet=0;idet<kNDetTypes;idet++) {
    //PixUpg
    if (idet==kDetPixUpg) {
      sim = GetSimulationModel(idet); 
      if (!sim) {
	sim = new AliITSsimulationPixUpg(this);
	SetSimulationModel(idet,sim);
      }
    }
  }
}


//__________________________________________________________
void AliITSDetTypeSimUpg::AddSimDigit(Int_t branch, const AliITSdigit* d)
{  
  // Add a simulated digit.
  TClonesArray* ldigits = DigitsAddress(branch);
  if (!ldigits) AliFatal(Form("Digits array for detector %d is not initialized",branch));
  int* ndig = GetNDigitArray();
  //
  switch(branch){
  case kDetPixUpg:
    new( (*ldigits)[ndig[branch]++]) AliITSdigitPixUpg(*((AliITSdigitPixUpg*)d));
    break;
  default:
    AliFatal(Form("Digit for unknown detector type %d",branch));
  }  
}

//______________________________________________________________________
void AliITSDetTypeSimUpg::AddSimDigit(Int_t branch,Float_t phys,Int_t *digits,
				      Int_t *tracks,Int_t *hits,Float_t *charges, 
				      Int_t sigexpanded)
{
  //   Add a simulated digit to the list.

  TClonesArray *ldigits = DigitsAddress(branch);
  if (!ldigits) AliFatal(Form("Digits array for detector %d is not initialized",branch));
  int* ndig = GetNDigitArray();
  //
  switch(branch){
  case kDetPixUpg:
    new((*ldigits)[ndig[branch]++]) AliITSdigitPixUpg(digits,tracks,hits);
    break;
  default:
    AliFatal(Form("Digit for unknown detector type %d",branch));
  } 
}
