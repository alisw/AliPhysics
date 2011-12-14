/**************************************************************************
 * Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
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

///////////////////////////////////////////////////////////////////
//                                                               //
// Implementation of the class to apply SDD map corrections      //
// Origin: F.Prino, Torino, prino@to.infn.it                     //
//                                                               //
///////////////////////////////////////////////////////////////////

#include "TObjArray.h"
#include "TString.h"
#include "TFile.h"
#include "AliLog.h"
#include "AliCDBEntry.h"
#include "AliITSCorrMapSDD.h"
#include "AliITSCorrectSDDPoints.h"
#include "AliITSsegmentationSDD.h"

ClassImp(AliITSCorrectSDDPoints)

//______________________________________________________________________
AliITSCorrectSDDPoints::AliITSCorrectSDDPoints():
  TObject(),
  fArrayOfMaps(0),
  fSegmentationSDD(0)
{
  // default constructor
  TFile* fil=new TFile("$ALICE_ROOT/OCDB/ITS/Calib/MapsTimeSDD/Run0_9999999_v0_s0.root");
  AliCDBEntry* e=(AliCDBEntry*)fil->Get("AliCDBEntry");
  fArrayOfMaps=(TObjArray*)e->GetObject();
  e->SetOwner(kTRUE);
  fil->Close();
  AliInfo(Form("%d AliITSCorrMapSDD objects in file %s",fArrayOfMaps->GetEntries(),fil->GetName()));
  fSegmentationSDD=new AliITSsegmentationSDD();
}

//______________________________________________________________________
AliITSCorrectSDDPoints::AliITSCorrectSDDPoints(TObjArray* maps):
  TObject(),
  fArrayOfMaps(maps),
  fSegmentationSDD(new AliITSsegmentationSDD())
{
  // constructor from external array
}

//______________________________________________________________________
AliITSCorrectSDDPoints::AliITSCorrectSDDPoints(TString filname):
  TObject(),
  fArrayOfMaps(0),
  fSegmentationSDD(0)
{
  // standard constructor
  TFile* fil=new TFile(filname.Data());
  AliCDBEntry* e=(AliCDBEntry*)fil->Get("AliCDBEntry");
  fArrayOfMaps=(TObjArray*)e->GetObject();
  e->SetOwner(kTRUE);
  fil->Close();
  AliInfo(Form("%d AliITSCorrMapSDD objects in file %s",fArrayOfMaps->GetEntries(),fil->GetName()));
  fSegmentationSDD=new AliITSsegmentationSDD();
}

//______________________________________________________________________
AliITSCorrectSDDPoints::~AliITSCorrectSDDPoints(){
  //
  if(fArrayOfMaps) delete fArrayOfMaps;
}

//______________________________________________________________________
void AliITSCorrectSDDPoints::SetCorrectionMaps(const TObjArray *arr)
{
  // replace the maps
  delete fArrayOfMaps;
  fArrayOfMaps = (TObjArray*)arr;
}

//______________________________________________________________________
Float_t AliITSCorrectSDDPoints::GetCorrection(Int_t modId, Float_t zloc, Float_t xloc) const{
  // returns correction to SDD drift corrdinate in cm
  Int_t nSide=fSegmentationSDD->GetSideFromLocalX(xloc);
  Int_t iSide=2*(modId-240)+nSide;
  if(iSide<0 || iSide >= 520){ 
    AliError(Form("Side out of range %d",iSide));
    return 0.;
  }
  AliITSCorrMapSDD* m=(AliITSCorrMapSDD*)fArrayOfMaps->At(iSide);
  return m->GetCorrection(zloc,xloc,fSegmentationSDD);
}
