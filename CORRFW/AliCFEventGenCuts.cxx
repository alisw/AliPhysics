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
// Cut on the Event at generator level: for the moment just 
// the requirements on the MB process type, number of tracks and on 
// the 3-D vertex position are implemented
// The argument of IsSelected member function (passed object) is cast into 
// an AliMCEvent. In the future may be modified to use AliVEvent interface
// and include more cut variables.
// The class derives from AliCFCutBase
// Author:S.Arcelli Silvia.Arcelli@cern.ch

#include "TBits.h"
#include "TList.h"
#include "AliLog.h"
#include "AliMCEvent.h"
#include <AliGenEventHeader.h>
#include <AliGenPythiaEventHeader.h>
#include <AliGenCocktailEventHeader.h>
#include "AliCFEventGenCuts.h"

ClassImp(AliCFEventGenCuts) 
//____________________________________________________________________
AliCFEventGenCuts::AliCFEventGenCuts() : 
  AliCFCutBase(),
  fMBProcType(-1),
  fNTracksMin(-1),
  fNTracksMax(100000),
  fRequireVtxCuts(kFALSE),
  fVtxXMax(1.e99),
  fVtxYMax(1.e99),
  fVtxZMax(1.e99),
  fVtxXMin(-1.e99),
  fVtxYMin(-1.e99),
  fVtxZMin(-1.e99),
  fBitMap(0x0)
{
  //
  //ctor
  //
  fBitMap=new TBits(0);
}
//____________________________________________________________________
AliCFEventGenCuts::AliCFEventGenCuts(Char_t* name, Char_t* title) : 
  AliCFCutBase(name,title),
  fMBProcType(-1),
  fNTracksMin(-1),
  fNTracksMax(100000),
  fRequireVtxCuts(kFALSE),
  fVtxXMax(1.e99),
  fVtxYMax(1.e99),
  fVtxZMax(1.e99),
  fVtxXMin(-1.e99),
  fVtxYMin(-1.e99),
  fVtxZMin(-1.e99),
  fBitMap(0x0)
 {
  //
  //ctor
  //
  fBitMap=new TBits(0);
 }
//____________________________________________________________________
AliCFEventGenCuts::AliCFEventGenCuts(const AliCFEventGenCuts& c) : 
  AliCFCutBase(c),
  fMBProcType(c.fMBProcType),
  fNTracksMin(c.fNTracksMin),
  fNTracksMax(c.fNTracksMax),
  fRequireVtxCuts(c.fRequireVtxCuts),
  fVtxXMax(c.fVtxXMax),
  fVtxYMax(c.fVtxYMax),
  fVtxZMax(c.fVtxZMax),
  fVtxXMin(c.fVtxXMin),
  fVtxYMin(c.fVtxYMin),
  fVtxZMin(c.fVtxZMin),
  fBitMap(c.fBitMap)
 
{
  //
  //copy constructor
  //
}
//____________________________________________________________________
AliCFEventGenCuts::~AliCFEventGenCuts() {
  //
  //dtor
  //

  if(fBitMap)delete fBitMap;
}
//____________________________________________________________________
AliCFEventGenCuts& AliCFEventGenCuts::operator=(const AliCFEventGenCuts& c)
{
  //
  // Assignment operator
  //
  if (this != &c) {
    AliCFCutBase::operator=(c) ;
    fMBProcType=c.fMBProcType;
    fNTracksMin=c.fNTracksMin;
    fNTracksMax=c.fNTracksMax;
    fRequireVtxCuts=c.fRequireVtxCuts;
    fVtxXMax=c.fVtxXMax;
    fVtxYMax=c.fVtxYMax;
    fVtxZMax=c.fVtxZMax;
    fVtxXMin=c.fVtxXMin;
    fVtxYMin=c.fVtxYMin;
    fVtxZMin=c.fVtxZMin;
    fBitMap=c.fBitMap;
  }
  return *this ;
}
//____________________________________________________________________
Bool_t AliCFEventGenCuts::IsSelected(TObject* obj) {
  //
  //Check if the requested cuts are passed
  //

  SelectionBitMap(obj);

  Bool_t isSelected = kTRUE;

  for (UInt_t icut=0; icut<fBitMap->GetNbits();icut++)
	if(!fBitMap->TestBitNumber(icut)) isSelected = kFALSE;

  return isSelected;

}

//____________________________________________________________________
void AliCFEventGenCuts::SelectionBitMap(TObject* obj){
  //
  //cut on the MB process type, the number of charged and neutral 
  //tracks and on the event vertex. So far specific to AliMCEvents
  //

  //Check if the requested cuts are passed and return a bitmap
  for(Int_t j=0;j<kNCuts;j++)fBitMap->SetBitNumber(j,kFALSE);
  AliMCEvent* ev = dynamic_cast<AliMCEvent *>(obj);
  if ( !ev ) return;
  AliGenEventHeader*genHeader = ev->GenEventHeader();  


  fBitMap->SetBitNumber(0,kTRUE);
  if(fMBProcType>-1){
    Int_t process=ProcType(genHeader);
    if(process==-1){
      AliInfo(Form(" not a pythia event, not checking on the process type"));
    }else{

      switch (fMBProcType)  {
      case kND:
	{
	  if (!( process!=92 && process!=93 && process!=94))
	    fBitMap->SetBitNumber(0,kFALSE);
	  break;
	}
      case kSD:
	{
	  if (!( process==92 || process==93))
	    fBitMap->SetBitNumber(0,kFALSE);
	  break;
	}
      case kDD:
	{
	  if (!( process==94))
	    fBitMap->SetBitNumber(0,kFALSE);
	  break;
	}
      }
    }
  }


  //Number of charged+neutral tracks:
  Int_t nTracks = ev->GetNumberOfTracks();
  fBitMap->SetBitNumber(1,kTRUE); //assume it is ok...
  if(nTracks<fNTracksMin || nTracks>fNTracksMax)
    fBitMap->SetBitNumber(1,kFALSE); 

  //now check the vertex cuts
  for(Int_t j=2;j<kNCuts;j++)fBitMap->SetBitNumber(j,kTRUE);

  TArrayF vtxPos(3);
  genHeader->PrimaryVertex(vtxPos);    

  if(fRequireVtxCuts){
    // Apply the cut
    if (vtxPos[0]>fVtxXMax || vtxPos[0]<fVtxXMin)
      fBitMap->SetBitNumber(2,kFALSE); 
    if (vtxPos[1]>fVtxYMax || vtxPos[1]<fVtxYMin)
      fBitMap->SetBitNumber(3,kFALSE); 
    if (vtxPos[2]>fVtxZMax || vtxPos[2]<fVtxZMin)
      fBitMap->SetBitNumber(4,kFALSE); 
  }  
  return;
}

 //______________________________________________________________________
Bool_t AliCFEventGenCuts::IsMBProcType(AliMCEvent *ev, PrType iproc){
  //
  //returns the type of MB process (if pythia events)
  //

  if ( !ev ) return kFALSE ;

  AliGenEventHeader*genHeader = ev->GenEventHeader();  

  Int_t process=ProcType(genHeader);

  switch (iproc)  {
  case kND: //Non Diffractive: Actually what is checked is ALL - SD - DD
    {
      if ( process!=92 && process!=93 && process!=94)
	return kTRUE;
    }
    break;
  case kSD: //Single Diffractive
    {
      if ( process==92 || process==93)
	return kTRUE;
    }
    break;
  case kDD: //Double Diffractive
    {
      if ( process==94)
	return kTRUE;
    }
    break;
  default: return kFALSE; break;
  }
  
  return kFALSE;
}
 //____________________________________________________________________________
Int_t AliCFEventGenCuts::ProcType(AliGenEventHeader *genHeader) {

  //get the Pythia process type: if we are not dealing with pythia stuff, 
  //return -1 and we do not apply the cut
  //

  // can only read pythia headers, either directly or from cocktalil header
  AliGenPythiaEventHeader* pythiaGenHeader = dynamic_cast<AliGenPythiaEventHeader*>(genHeader);

  if (!pythiaGenHeader) {

    AliGenCocktailEventHeader* genCocktailHeader = dynamic_cast<AliGenCocktailEventHeader*>(genHeader);
    if (!genCocktailHeader) {
      return -1;
    }

    TList* headerList = genCocktailHeader->GetHeaders();
    if (!headerList) {
      return -1;
    }

    for (Int_t i=0; i<headerList->GetEntries(); i++) {
      pythiaGenHeader = dynamic_cast<AliGenPythiaEventHeader*>(headerList->At(i));
      if (pythiaGenHeader)
        break;
    }

    if (!pythiaGenHeader) {
      return -1;
    }
  }
  
  Int_t process=pythiaGenHeader->ProcessType();
  return process;
}

