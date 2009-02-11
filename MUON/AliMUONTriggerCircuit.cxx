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

//-----------------------------------------------------------------------------
/// \class AliMUONTriggerCircuit
/// Contains as data members the Y positions of the X declusturized strips and 
/// the X positions of the (doubled or not) Y strips.
/// This is used to associate the global positions to the fired strips of the 
/// local trigger output (see AliMUONTrackReconstructor::MakeTriggerTrack)
///
/// \author Philippe Crochet (LPCCFd)
//-----------------------------------------------------------------------------

#include "AliMUONTriggerCircuit.h"
#include "AliMUONConstants.h"
#include "AliMUONGeometryTransformer.h"

#include "AliMpTrigger.h"
#include "AliMpSlat.h"
#include "AliMpPCB.h"
#include "AliMpSegmentation.h"
#include "AliMpVSegmentation.h"
#include "AliMpCathodType.h"
#include "AliMpDDLStore.h"
#include "AliMpLocalBoard.h"
#include "AliMpConstants.h"
#include "AliMpPad.h"

#include "AliRun.h"
#include "AliLog.h"

#include <TMath.h>
#include <Riostream.h>

/// \cond CLASSIMP
ClassImp(AliMUONTriggerCircuit)
/// \endcond

//----------------------------------------------------------------------
AliMUONTriggerCircuit::AliMUONTriggerCircuit(const AliMUONGeometryTransformer* transformer)
    : TObject(),
      fkTransformer(transformer),
      fkCurrentSeg(0x0),
      fCurrentDetElem(0x0),
      fCurrentLocalBoard(0x0)
{
/// Constructor
  
    for (Int_t i = 1; i < AliMpConstants::NofLocalBoards()+1; ++i) {
      fXpos11[i].Set(16); 
      fYpos11[i].Set(31);
      fYpos21[i].Set(63);
    }

    for (Int_t i = 1; i < AliMpConstants::NofLocalBoards()+1; ++i) { // board begins at 1
  
    AliMpLocalBoard* localBoard = AliMpDDLStore::Instance()->GetLocalBoard(i);

    if (!localBoard)
    {
      AliError(Form("Did not get localboard %d",i));
      continue;
    }

    LoadXPos(localBoard);
    LoadYPos(localBoard);

  }

}

//----------------------------------------------------------------------
AliMUONTriggerCircuit::~AliMUONTriggerCircuit()
{
/// Destructor
   for (Int_t i = 1; i < AliMpConstants::NofLocalBoards()+1; ++i) {
     fXpos11[i].Reset();
     fYpos11[i].Reset();
     fYpos21[i].Reset();
    }

} 

//----------------------------------------------------------------------
AliMUONTriggerCircuit::AliMUONTriggerCircuit(const AliMUONTriggerCircuit& circuit)
    :  TObject(circuit),
       fkTransformer(circuit.fkTransformer), // do not copy, just pointed to
       fkCurrentSeg(circuit.fkCurrentSeg),
       fCurrentDetElem(circuit.fCurrentDetElem),
       fCurrentLocalBoard(circuit.fCurrentLocalBoard)
{
/// Copy constructor

    for (Int_t i = 1; i < AliMpConstants::NofLocalBoards()+1; ++i) {
      fXpos11[i] = circuit.fXpos11[i];
      fYpos11[i] = circuit.fYpos11[i];
      fYpos21[i] = circuit.fYpos21[i];
    }

}
//----------------------------------------------------------------------
AliMUONTriggerCircuit& AliMUONTriggerCircuit::operator=(const AliMUONTriggerCircuit& circuit) 
{
/// Assignment operator

  if (this == &circuit) return *this;

  fkTransformer      = circuit.fkTransformer;
  fkCurrentSeg       = circuit.fkCurrentSeg;
  fCurrentDetElem    = circuit.fCurrentDetElem;
  fCurrentLocalBoard = circuit.fCurrentLocalBoard;

  for (Int_t i = 1; i < AliMpConstants::NofLocalBoards()+1; ++i) {
    fXpos11[i] = circuit.fXpos11[i];
    fYpos11[i] = circuit.fYpos11[i];
    fYpos21[i] = circuit.fYpos21[i];
  }

  return *this;

}

//---------------------------------------------------------------------
void AliMUONTriggerCircuit::LoadYPos(AliMpLocalBoard* const localBoard)
{
/// fill fYpos11 and fYpos21 -> y position of X declusterized strips
  
  fCurrentLocalBoard = localBoard->GetId();
  Int_t ichamber = 0;
  Int_t icathode = 0;    
    
  Int_t zeroDown = localBoard->GetSwitch(AliMpLocalBoard::kZeroDown);
  Int_t zeroUp   = localBoard->GetSwitch(AliMpLocalBoard::kZeroUp);
 
  Int_t iline = localBoard->GetPosition().GetFirst();
  Int_t icol  = localBoard->GetPosition().GetSecond();
  if ( iline == 5 ) --icol;

  //--- first plane 
  ichamber = 10;
  fCurrentDetElem = AliMpDDLStore::Instance()->GetDEfromLocalBoard(fCurrentLocalBoard, ichamber);

  fkCurrentSeg = AliMpSegmentation::Instance()
        ->GetMpSegmentation(fCurrentDetElem, AliMp::GetCathodType(icathode));  

  Int_t iFirstStrip = FirstStrip(localBoard);
  Int_t iLastStrip  = iFirstStrip + 16;    
  Int_t iStripCircuit = 0;

  FillXstrips(icol, iFirstStrip, iLastStrip, 
	      iStripCircuit, fYpos11[fCurrentLocalBoard]);
  
  //--- second plane 
  ichamber = 12;
  fCurrentDetElem = AliMpDDLStore::Instance()->GetDEfromLocalBoard(fCurrentLocalBoard, ichamber);

  fkCurrentSeg = AliMpSegmentation::Instance()
        ->GetMpSegmentation(fCurrentDetElem, AliMp::GetCathodType(icathode));  

  // second plane middle part
  Int_t iFirstStripMiddle = FirstStrip(localBoard);
  Int_t iLastStripMiddle  = iFirstStrip + 16;
  iStripCircuit = 8;

  FillXstrips(icol, iFirstStripMiddle, iLastStripMiddle,
	      iStripCircuit, fYpos21[fCurrentLocalBoard]);
  
  // second plane upper part
  if (zeroUp == 0) { // something up
    Int_t iFirstStripUp;
    Int_t iLastStripUp;
    Int_t icolUp = icol;

    // check if we need to move to another detElemId
    AliMpPad pad = fkCurrentSeg->PadByIndices(AliMpIntPair(icol-1,iLastStripMiddle+1),kFALSE);

    if (pad.IsValid()) { // upper strips within same detElemId
      iFirstStripUp = iLastStripMiddle;
      iLastStripUp  = iFirstStripUp + 8;

    } else {             // upper strips in another detElemId
      fCurrentDetElem = AliMpDDLStore::Instance()->
	           GetNextDEfromLocalBoard(fCurrentLocalBoard, ichamber);

      fkCurrentSeg = AliMpSegmentation::Instance()
            ->GetMpSegmentation(fCurrentDetElem, AliMp::GetCathodType(icathode));  

      iFirstStripUp = 0;
      iLastStripUp  = iFirstStripUp + 8;
      if (iline == 4) icolUp = icol - 1; // special case
      if (iline == 5) icolUp = icol + 1; // special case
    }
    
    iStripCircuit = 24;
    FillXstrips(icolUp, iFirstStripUp, iLastStripUp,
		iStripCircuit, fYpos21[fCurrentLocalBoard]);
    
    // fill strip between middle and upper part
    fYpos21[fCurrentLocalBoard][47] = (fYpos21[fCurrentLocalBoard][46] + 
				       fYpos21[fCurrentLocalBoard][48])/2.;
  } // end of something up
  
  // restore current detElemId & segmentation
  fCurrentDetElem = AliMpDDLStore::Instance()->GetDEfromLocalBoard(fCurrentLocalBoard, ichamber);
  fkCurrentSeg = AliMpSegmentation::Instance()
      ->GetMpSegmentation(fCurrentDetElem, AliMp::GetCathodType(icathode));  

  // second plane lower part
  if (zeroDown == 0) { // something down
    Int_t iFirstStripDo;
    Int_t iLastStripDo;
    Int_t icolDo = icol;
    
    // check if we need to move to another detElemId      
    AliMpPad pad = fkCurrentSeg->PadByIndices(AliMpIntPair(icol-1,iFirstStripMiddle-1),kFALSE);
    if (pad.IsValid()) { // lower strips within same detElemId
      iFirstStripDo = iFirstStripMiddle - 8;
      iLastStripDo  = iFirstStripDo + 8;	      

    } else {             // lower strips in another detElemId 
      fCurrentDetElem = AliMpDDLStore::Instance()
	  ->GetPreviousDEfromLocalBoard(fCurrentLocalBoard, ichamber);

      fkCurrentSeg = AliMpSegmentation::Instance()
	  ->GetMpSegmentation(fCurrentDetElem, AliMp::GetCathodType(icathode));  

      // get iFirstStrip in this module 
      const AliMpTrigger* t = AliMpSegmentation::Instance()->GetTrigger(fkCurrentSeg);
      const AliMpSlat* slat = t->GetLayer(0);

      if (iline == 5) icolDo = icol + 1; // special case
      if (iline == 6) icolDo = icol - 1; // special case
	    
      const AliMpPCB* pcb = slat->GetPCB(icolDo-1);
      iFirstStripDo = (pcb->Iymax() + 1) - 8;
      iLastStripDo =  iFirstStripDo + 8;
    }  
    
    iStripCircuit = 0;
    FillXstrips(icolDo, iFirstStripDo, iLastStripDo,
		iStripCircuit, fYpos21[fCurrentLocalBoard]);
    
    // fill strip between middle and upper part
    fYpos21[fCurrentLocalBoard][15] = (fYpos21[fCurrentLocalBoard][14] + 
				       fYpos21[fCurrentLocalBoard][16])/2.;
  } // end of something down
  
}

//----------------------------------------------------------------------
void AliMUONTriggerCircuit::FillXstrips(const Int_t icol, 
					const Int_t iFirstStrip, const Int_t iLastStrip, 
					Int_t liStripCircuit, TArrayF& ypos)
{    
/// fill 
  Double_t xyGlobal[2] = {0.};
  for (Int_t istrip = iFirstStrip; istrip < iLastStrip; ++istrip) {

    AliMpPad pad = fkCurrentSeg->PadByIndices(AliMpIntPair(icol-1,istrip),kTRUE);
    if ( !pad.IsValid() ) {
      StdoutToAliError(cout << "Pad not found in seg " << endl;
                       fkCurrentSeg->Print();
                       cout << " ix,iy=" << icol-1 << "," << istrip << endl;
                       );
    }
    Float_t yDim = pad.Dimensions().Y(); // half size! 

    XYGlobal(pad,xyGlobal);
    
    ypos[2*liStripCircuit] = xyGlobal[1];
    if (istrip != (iLastStrip - 1)) ypos[2*liStripCircuit+1] = xyGlobal[1] + yDim;
    liStripCircuit++;
  }    
}


//----------------------------------------------------------------------
void AliMUONTriggerCircuit::LoadXPos(AliMpLocalBoard* const localBoard)
{
/// fill fXpos11 -> x position of Y strips for the first plane only
/// fXpos11 contains the x position of Y strip for the current circuit
/// taking into account whether or nor not part(s) of the circuit
/// (middle, up or down) has(have) 16 strips (handdled by means of switchs)

  fCurrentLocalBoard = localBoard->GetId();

  Int_t ichamber = 10;
  Int_t icathode = 1;
  
  Int_t x2u = localBoard->GetSwitch(AliMpLocalBoard::kX2u);
  Int_t x2m = localBoard->GetSwitch(AliMpLocalBoard::kX2m);
  Int_t x2d = localBoard->GetSwitch(AliMpLocalBoard::kX2d);
  Int_t zeroAllYLSB = localBoard->GetSwitch(AliMpLocalBoard::kZeroAllYLSB);

  Int_t  iStripCircuit = 0;
  Int_t  iFirstStrip   = 0;
  Int_t  iLastStrip    = 0;
  Bool_t doubling      = kFALSE;
  
  Int_t iline  = localBoard->GetPosition().GetFirst();
  Int_t icol   = localBoard->GetPosition().GetSecond();
  if ( iline == 5 ) --icol;

  fCurrentDetElem = AliMpDDLStore::Instance()->GetDEfromLocalBoard(fCurrentLocalBoard, ichamber);

  fkCurrentSeg = AliMpSegmentation::Instance()
        ->GetMpSegmentation(fCurrentDetElem, AliMp::GetCathodType(icathode));  

  // check if one needs a strip doubling or not
  if ( (x2u || x2m || x2d ) && x2m ) doubling = kTRUE;
  

  // check if one starts at strip = 0 or 8 (boards 26-29 and 143-146)
  if (zeroAllYLSB) iStripCircuit = 8;
  
  // get iFirstStrip in this module 
  const AliMpTrigger* t = AliMpSegmentation::Instance()->GetTrigger(fkCurrentSeg);
  const AliMpSlat* slat = t->GetLayer(0);
  
  const AliMpPCB* pcb = slat->GetPCB(icol-1);
  iFirstStrip = pcb->Ixmin();
  

  if (doubling || zeroAllYLSB == 1) iLastStrip = iFirstStrip + 8;
  else iLastStrip = iFirstStrip + 16;
  
  FillYstrips(iFirstStrip, iLastStrip, iStripCircuit, doubling);  
}

//----------------------------------------------------------------------
void AliMUONTriggerCircuit::FillYstrips(const Int_t iFirstStrip, const Int_t iLastStrip, 
					Int_t liStripCircuit,
					const Bool_t doubling)
{    
/// fill
  Double_t xyGlobal[2] = {0.};

  for (Int_t istrip = iFirstStrip; istrip < iLastStrip; ++istrip) {

    AliMpPad pad = fkCurrentSeg->PadByIndices(AliMpIntPair(istrip,0),kTRUE);

    if ( !pad.IsValid() )
    {
	StdoutToAliError(cout << "Pad not found in seg " << endl;
			 fkCurrentSeg->Print();
			 cout << " ix,iy=" << istrip << "," << 0 << endl;
			 );
    }
    Float_t xDim = pad.Dimensions().X(); // half size!

    XYGlobal(pad,xyGlobal);
    
    if (!doubling) {	
      fXpos11[fCurrentLocalBoard].AddAt(xyGlobal[0], liStripCircuit);
    } else if (doubling) {

      fXpos11[fCurrentLocalBoard].AddAt(TMath::Sign(1.,xyGlobal[0]) * 
				  (TMath::Abs(xyGlobal[0]) - xDim/2.), 2*liStripCircuit);

      fXpos11[fCurrentLocalBoard].AddAt(TMath::Sign(1.,xyGlobal[0]) *
				  (TMath::Abs(fXpos11[fCurrentLocalBoard][2*liStripCircuit]) + xDim),
					2*liStripCircuit + 1); 
    }

    liStripCircuit++;
  }    
}

//----------------------------------------------------------------------
void AliMUONTriggerCircuit::XYGlobal(const AliMpPad& pad,
				     Double_t* xyGlobal)
{
/// returns pad x & y positions and x & y pad dimensions in global coordinates
/// note: no need for transformation for pad dimensions
  
  // get the pad position and dimensions
  Double_t xl1 = pad.Position().X();
  Double_t yl1 = pad.Position().Y();
  Double_t zg1 = 0;
  
  // positions from local to global 
  fkTransformer->Local2Global(fCurrentDetElem, xl1, yl1, 0, 
                                 xyGlobal[0], xyGlobal[1], zg1);
}


//----------------------------------------------------------------------
//--- methods which return member data related info
//----------------------------------------------------------------------
Float_t AliMUONTriggerCircuit::GetY11Pos(Int_t localBoardId, Int_t istrip) const 
{
/// returns Y position of X strip istrip in MC11
  return fYpos11[localBoardId][istrip];
}
//----------------------------------------------------------------------
Float_t AliMUONTriggerCircuit::GetY21Pos(Int_t localBoardId, Int_t istrip) const 
{
/// returns Y position of X strip istrip in MC21
  return fYpos21[localBoardId][istrip];
}
//----------------------------------------------------------------------
Float_t AliMUONTriggerCircuit::GetX11Pos(Int_t localBoardId, Int_t istrip) const 
{
/// returns X position of Y strip istrip in MC11
  return fXpos11[localBoardId][istrip];
}

//----------------------------------------------------------------------
Int_t AliMUONTriggerCircuit::FirstStrip(AliMpLocalBoard* localBoard)
{
/// returns the first strip from mapping for board boardName
/// take care of special case for boards RC1L6B12 & LC1L6B12
  Int_t iFirstStrip = -1;
  Int_t boardNumber = atoi(localBoard->GetName()+6);

  Int_t iline = localBoard->GetPosition().GetFirst();
  Int_t icol  = localBoard->GetPosition().GetSecond();
  if ( iline == 5 ) --icol;

  switch (boardNumber)
  {
    case 12:
      iFirstStrip = 0;
      break;		
    case 34:
      iFirstStrip = 16;
      break;		
    case 56:
      iFirstStrip = 32;
      break;		
    case 78:
      iFirstStrip = 48;
      break;		
  }
  if (icol == 1 && iline == 6) iFirstStrip += 16; // special case
  return iFirstStrip;
}

//----------------------------------------------------------------------
Float_t AliMUONTriggerCircuit::PtCal(Int_t localBoardId, Int_t istripX, Int_t idev, Int_t istripY) const{
/// returns calculated pt for circuit/istripX/idev/istripY according 
/// to the formula of the TRD. Note : idev (input) is in [0+30]

  //  Int_t jdev = idev - 15;        // jdev in [-15+15]
  Int_t istripX2=istripX+idev+1; // find istripX2 using istripX and idev

  Float_t yPosX1=fYpos11[localBoardId][istripX];
  Float_t yPosX2=fYpos21[localBoardId][istripX2];
  Float_t xPosY1=fXpos11[localBoardId][istripY];

// Z distance between IP and center of dipole
  Float_t zf= TMath::Abs(0.5 *(AliMUONConstants::CoilZ() + AliMUONConstants::YokeZ()));
  Float_t z1=AliMUONConstants::DefaultChamberZ(10);
  Float_t z2=AliMUONConstants::DefaultChamberZ(12);
  Float_t thetaDev=(1./zf)*(yPosX1*z2-yPosX2*z1)/(z2-z1);
  Float_t xf=xPosY1*zf/z1; 
  Float_t yf=yPosX2-((yPosX2-yPosX1)*(z2-zf))/(z2-z1);
  return (3.*0.3/TMath::Abs(thetaDev)) * TMath::Sqrt(xf*xf+yf*yf)/zf;
}
