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

//.
// ITSupgrade base class to reconstruct an event
//.
//.
//.
#include "TObjArray.h"
#include "TTree.h"
#include "AliITSRecPoint.h"
#include "AliITSReconstructor.h"
#include "AliITSupgrade.h"
#include "AliITSUpgradeReconstructor.h" //class header
#include "AliITSDetTypeRec.h"
#include "AliITS.h"              //Reconstruct() 
#include "AliESDEvent.h"           //FillEsd()
#include "AliRawReader.h"          //Reconstruct() for raw digits
#include "AliRun.h"
#include "AliLog.h"                //
#include "AliITSRawStream.h"     //ConvertDigits()
#include "AliITSsegmentationUpgrade.h"
#include "AliITSUpgradeClusterFinder.h"
#include "AliITStrackerUpgrade.h"
#include "AliITStrackerU.h"
ClassImp(AliITSUpgradeReconstructor)

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  AliITSUpgradeReconstructor::AliITSUpgradeReconstructor():
    AliReconstructor(), 
    fDigits(0x0),
    fNlayers(0)
{
  //
  //ctor
  //

  AliITSsegmentationUpgrade *s = new AliITSsegmentationUpgrade();
  fNlayers = s->GetNLayers();
  delete s;
  
  fDigits = new TObjArray(fNlayers);
  fDigits->SetOwner(kTRUE);
  for(Int_t iLay =0; iLay<fNlayers; iLay++) fDigits->AddAt(new TClonesArray("AliITSDigitUpgrade",5000),iLay);
  AliInfo("    ************* Using the upgrade reconstructor! ****** ");



}//AliITSReconstructor
//-----------------------------------------------------------------------
AliITSUpgradeReconstructor::~AliITSUpgradeReconstructor(){
  //Destructor
  if(fDigits){
    fDigits->Delete();
    delete fDigits;
    fDigits=0;
  }
}

//_________________________________________________________________
void AliITSUpgradeReconstructor::Init() {
  // Initalize this constructor bet getting/creating the objects
  // nesseary for a proper ITS reconstruction.
  // Inputs:
  //   none.
  // Output:
  //   none.
  // Return:
  //   none.
  return;
}
//_____________________________________________________________________
void AliITSUpgradeReconstructor::ResetDigits(){
  // Reset number of digits and the digits array for the ITS detector.

  if(!fDigits) return;
  for(Int_t i=0;i<fNlayers;i++){
    ResetDigits(i);
  }
}
//____________________________________________________________________
void AliITSUpgradeReconstructor::ResetDigits(Int_t branch){
  // Reset number of digits and the digits array for this branch.

  if(fDigits->At(branch)) ((TClonesArray*)fDigits->At(branch))->Clear();

}
//___________________________________________________________________
void AliITSUpgradeReconstructor::Reconstruct(TTree* digitsTree, TTree* clustersTree) const {
  //
  // Clustering
  //
  for(Int_t iLay=0;iLay<fNlayers;iLay++) {
    digitsTree->SetBranchAddress(Form("Layer%d",iLay),&(*fDigits)[iLay]);
  }
 
  digitsTree->GetEntry(0);
  AliITSUpgradeClusterFinder *clf = new AliITSUpgradeClusterFinder();
  clf->MakeRecPointBranch(clustersTree);
  clf->SetRecPointTreeAddress(clustersTree);
  clf->DigitsToRecPoints(fDigits);     //cluster finder
  clustersTree->Fill();                //fill tree for current event
  delete clf;
  for(Int_t iLay=0;iLay<fNlayers;iLay++){
    fDigits->At(iLay)->Clear();
  }

}
//_______________________________________________________________________________________________________________
AliTracker* AliITSUpgradeReconstructor::CreateTracker() const
{
  //
  // Create the ITSUpgrade tracker switching betqween the general case and the StandAlone case. 
  // AliITStrackerUpgrade is used only for SA tracks.
  // 

  if(GetRecoParam()->GetTrackerSAOnly()){
    AliITStrackerUpgrade *trackUp = new AliITStrackerUpgrade(fNlayers);
    if(GetRecoParam()->GetTrackerSAOnly()) trackUp->SetSAFlag(kTRUE);
    if(trackUp->GetSAFlag())AliDebug(1,"Tracking Performed in ITS only\n");
    if(GetRecoParam()->GetInwardFindingSA()){
      trackUp->SetInwardFinding();
      trackUp->SetInnerStartLayer(GetRecoParam()->GetInnerStartLayerSA());
    }else{
      trackUp->SetOutwardFinding();
      trackUp->SetOuterStartLayer(GetRecoParam()->GetOuterStartLayerSA());
    }
    trackUp->SetMinNPoints(GetRecoParam()->GetMinNPointsSA());
    return trackUp;
  } else {
    AliITStrackerU *t = new AliITStrackerU();
    return t;
  }
}
//_______________________________________________________________________
