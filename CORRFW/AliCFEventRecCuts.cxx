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
// Cut on the Event at reconstructed level: for the moment 
// the requirements on the number of charged tracks and on 
// the vertex position and resolution are implemented
// The argument of IsSelected member function (passed object) is cast into 
// an AliESDEvent. In the future may be modified to use AliVEvent interface
// and include more cut variables.
// The class derives from AliCFCutBase
// Author:S.Arcelli Silvia.Arcelli@cern.ch
//
//
#include "TBits.h"
#include "AliLog.h"
#include "AliESDEvent.h"
#include "AliESDVertex.h"
#include "AliCFEventRecCuts.h"
ClassImp(AliCFEventRecCuts) 
//____________________________________________________________________
AliCFEventRecCuts::AliCFEventRecCuts() : 
  AliCFCutBase(),
  fNTracksMin(-1),
  fNTracksMax(1000000),
  fRequireVtxCuts(kFALSE),
  fVtxXMax(1.e99),
  fVtxYMax(1.e99),
  fVtxZMax(1.e99),
  fVtxXMin(-1.e99),
  fVtxYMin(-1.e99),
  fVtxZMin(-1.e99),
  fVtxXResMax(1.e99),
  fVtxYResMax(1.e99),
  fVtxZResMax(1.e99),
  fBitMap(0x0)
{
  //
  //ctor
  //
  fBitMap=new TBits(0);
}
//____________________________________________________________________
AliCFEventRecCuts::AliCFEventRecCuts(Char_t* name, Char_t* title) : 
  AliCFCutBase(name,title),
  fNTracksMin(-1),
  fNTracksMax(1000000),
  fRequireVtxCuts(kFALSE),
  fVtxXMax(1.e99),
  fVtxYMax(1.e99),
  fVtxZMax(1.e99),
  fVtxXMin(-1.e99),
  fVtxYMin(-1.e99),
  fVtxZMin(-1.e99),
  fVtxXResMax(1.e99),
  fVtxYResMax(1.e99),
  fVtxZResMax(1.e99),
  fBitMap(0x0)
 {
  //
  //ctor
  //
  fBitMap=new TBits(0);
 }
//____________________________________________________________________
AliCFEventRecCuts::AliCFEventRecCuts(const AliCFEventRecCuts& c) : 
  AliCFCutBase(c),
  fNTracksMin(c.fNTracksMin),
  fNTracksMax(c.fNTracksMax),
  fRequireVtxCuts(c.fRequireVtxCuts),
  fVtxXMax(c.fVtxXMax),
  fVtxYMax(c.fVtxYMax),
  fVtxZMax(c.fVtxZMax),
  fVtxXMin(c.fVtxXMin),
  fVtxYMin(c.fVtxYMin),
  fVtxZMin(c.fVtxZMin),
  fVtxXResMax(c.fVtxXResMax),
  fVtxYResMax(c.fVtxYResMax),
  fVtxZResMax(c.fVtxZResMax),
  fBitMap(c.fBitMap)
 
{
  //
  //copy constructor
  //
}
//____________________________________________________________________
AliCFEventRecCuts::~AliCFEventRecCuts() {
  //
  //dtor
  //

  if(fBitMap)delete fBitMap;
}
//____________________________________________________________________
AliCFEventRecCuts& AliCFEventRecCuts::operator=(const AliCFEventRecCuts& c)
{
  //
  // Assignment operator
  //
  if (this != &c) {
    AliCFCutBase::operator=(c) ;
    fNTracksMin=c.fNTracksMin;
    fNTracksMax=c.fNTracksMax;
    fRequireVtxCuts=c.fRequireVtxCuts;
    fVtxXMax=c.fVtxXMax;
    fVtxYMax=c.fVtxYMax;
    fVtxZMax=c.fVtxZMax;
    fVtxXMin=c.fVtxXMin;
    fVtxYMin=c.fVtxYMin;
    fVtxZMin=c.fVtxZMin;
    fVtxXResMax=c.fVtxXResMax;
    fVtxYResMax=c.fVtxYResMax;
    fVtxZResMax=c.fVtxZResMax;
    fBitMap=c.fBitMap;
  }
  return *this ;
}
//____________________________________________________________________
Bool_t AliCFEventRecCuts::IsSelected(TObject* obj) {
  //
  //Check if the requested cuts are passed
  //

  TBits *bitmap = SelectionBitMap(obj);

  Bool_t isSelected = kTRUE;

  for (UInt_t icut=0; icut<bitmap->GetNbits();icut++)
	if(!bitmap->TestBitNumber(icut)) isSelected = kFALSE;

  return isSelected;

}

//____________________________________________________________________
TBits *AliCFEventRecCuts::SelectionBitMap(TObject* obj) {
  //
  //cut on the number of charged tracks and on the event vertex.
  //so far specific to AliESDEvents
  //

  //Check if the requested cuts are passed and return a bitmap
  for(Int_t j=0;j<kNCuts;j++)fBitMap->SetBitNumber(j,kFALSE);
  AliESDEvent* esd = dynamic_cast<AliESDEvent *>(obj);
  if ( !esd ) return fBitMap ;

  //now start checking the cuts,
  //first assume the event will be accepted: 
  for(Int_t j=0;j<kNCuts;j++)fBitMap->SetBitNumber(j,kTRUE);

  //Number of charged tracks:
  Int_t nTracks = esd->GetNumberOfTracks();
  if(nTracks<fNTracksMin || nTracks>fNTracksMax)
    fBitMap->SetBitNumber(0,kFALSE); 
  
  if(fRequireVtxCuts){
    const AliESDVertex* vtxESD = esd->GetVertex();
    if(!vtxESD){
      for(Int_t j=1;j<kNCuts;j++)fBitMap->SetBitNumber(j,kFALSE); 
      return fBitMap;
    }
    // Require the vertex to have been reconstructed successfully
    if (strcmp(vtxESD->GetName(), "default")==0){
      AliWarning(Form(" No reconstructed vertex found, skip event"));    
      for(Int_t j=1;j<kNCuts;j++)fBitMap->SetBitNumber(j,kFALSE); 
      return fBitMap;
    }    
    // Pick up the position and uncertainties
    
    Double_t vtxPos[3];
    vtxPos[0] = vtxESD->GetXv();
    vtxPos[1] = vtxESD->GetYv();
    vtxPos[2] = vtxESD->GetZv();
    
    Double_t vtxRes[3];
    vtxRes[0] = vtxESD->GetXRes();
    vtxRes[1] = vtxESD->GetYRes();
    vtxRes[2] = vtxESD->GetZRes();
 
    // Apply the cut
    
    if (vtxPos[0]>fVtxXMax || vtxPos[0]<fVtxXMin)
      fBitMap->SetBitNumber(1,kFALSE); 
    if (vtxPos[1]>fVtxYMax || vtxPos[1]<fVtxYMin)
      fBitMap->SetBitNumber(2,kFALSE); 
    if (vtxPos[2]>fVtxZMax || vtxPos[2]<fVtxZMin)
      fBitMap->SetBitNumber(3,kFALSE); 
    if (vtxRes[0]==0 || vtxRes[0]>fVtxXResMax)
      fBitMap->SetBitNumber(4,kFALSE); 
    if (vtxRes[1]==0 || vtxRes[1]>fVtxYResMax)
      fBitMap->SetBitNumber(5,kFALSE); 
    if (vtxRes[2]==0 || vtxRes[2]>fVtxZResMax)
      fBitMap->SetBitNumber(6,kFALSE); 
  }  
  return fBitMap;
}
//__________________________________________________________________________________
void AliCFEventRecCuts::GetBitMap(TObject* obj, TBits *bitmap) {
  //
  // retrieve the pointer to the bitmap
  //
  bitmap = SelectionBitMap(obj);

}
