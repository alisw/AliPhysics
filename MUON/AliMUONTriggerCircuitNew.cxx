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

// --------------------
// Class AliMUONTriggerCircuitNew
// --------------------
// Contains as data members the Y positions of the X declusturized strips and 
// the X positions of the (doubled or not) Y strips.
// This is used to associate the global positions to the fired strips of the 
// local trigger output (see AliMUONTrackReconstructor::MakeTriggerTrack)

#include <TMath.h>
#include "AliMUONTriggerCircuitNew.h"
#include "AliRun.h"
#include "AliMUON.h"
#include "AliMUONTriggerConstants.h"
#include "AliMUONConstants.h"
#include "AliLog.h"
#include "AliMUONLocalTriggerBoard.h"
#include "AliMUONTriggerCrateStore.h"
#include "AliMUONTriggerCrate.h"
#include "AliMpTriggerSegmentation.h"
#include "AliMpTrigger.h"
#include "AliMpSlat.h"
#include "AliMpPCB.h"
#include "AliMpVSegmentation.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMpSegFactory.h"

/// \cond CLASSIMP
ClassImp(AliMUONTriggerCircuitNew)
/// \endcond

//----------------------------------------------------------------------
AliMUONTriggerCircuitNew::AliMUONTriggerCircuitNew()
: TObject(),
  fILocalBoard(0),
  fSegFactory(0x0),
  fTransformer(0x0)
{
/// Constructor
  
  Int_t i;  
  for (i=0; i<16; i++) { fXpos11[i]=0.; }
  for (i=0; i<31; i++) { fYpos11[i]=0.; }
  for (i=0; i<63; i++) { fYpos21[i]=0.; }
}

//----------------------------------------------------------------------
AliMUONTriggerCircuitNew::~AliMUONTriggerCircuitNew()
{
/// Destructor
} 

//----------------------------------------------------------------------
AliMUONTriggerCircuitNew::AliMUONTriggerCircuitNew(const AliMUONTriggerCircuitNew& circuit)
   :  TObject(circuit),
      fILocalBoard(circuit.fILocalBoard)
{
/// Copy constructor

  for (Int_t i = 0; i < 16; ++i)
    fXpos11[i] = circuit.fXpos11[i];

  for (Int_t i = 0; i < 31; ++i)
    fYpos11[i] = circuit.fYpos11[i];

  for (Int_t i = 0; i < 63; ++i)
    fYpos21[i] = circuit.fYpos21[i];

}
//----------------------------------------------------------------------
AliMUONTriggerCircuitNew& AliMUONTriggerCircuitNew::operator=(const AliMUONTriggerCircuitNew& circuit) 
{
/// Assignment operator

  if (this == &circuit) return *this;

  fILocalBoard = circuit.fILocalBoard;

  for (Int_t i = 0; i < 16; ++i)
    fXpos11[i] = circuit.fXpos11[i];

  for (Int_t i = 0; i < 31; ++i)
    fYpos11[i] = circuit.fYpos11[i];

  for (Int_t i = 0; i < 63; ++i)
    fYpos21[i] = circuit.fYpos21[i];

  return *this;

}
//----------------------------------------------------------------------
void AliMUONTriggerCircuitNew::Init(Int_t iCircuit, const AliMUONTriggerCrateStore& crates) 
{
/// initialize circuit characteristics
  fILocalBoard=iCircuit+1;//AliMUONTriggerConstants::CircuitId(iCircuit);
 
  LoadXPos(crates);
  LoadYPos(crates);
  
}

//---------------------------------------------------------------------
void AliMUONTriggerCircuitNew::LoadYPos(const AliMUONTriggerCrateStore& crates)
{
/// fill fYpos11 and fYpos21 -> y position of X declusterized strips
  
  const AliMUONLocalTriggerBoard* localBoard = crates.LocalBoard(fILocalBoard);
  
  if (!localBoard)
  {
    AliError(Form("Did not get localboard %d",fILocalBoard));
    return;
  }
  StdoutToAliDebug(1,localBoard->Print("CONF"););
  
  Int_t ichamber = 0;
  Int_t icathode = 0;    
  
  const AliMpVSegmentation* seg;
  
  Int_t zeroDown = localBoard->GetSwitch(AliMUONLocalTriggerBoard::kZeroDown);
  Int_t zeroUp = localBoard->GetSwitch(AliMUONLocalTriggerBoard::kZeroUp);
  
  //--- first plane 
  ichamber = 10;
  
  char side;
  Int_t iline, icol;
  DecodeBoardName(localBoard->GetName(),side,iline,icol);
  
  Int_t detElemId = DetElemId(ichamber,side,iline);
  seg = fSegFactory->CreateMpSegmentation(detElemId, icathode);  

  Int_t iFirstStrip = FirstStrip(localBoard->GetName());
  Int_t iLastStrip = iFirstStrip + 16;    
  Int_t iStripCircuit = 0;
  FillXstrips(seg,detElemId,icol, 
              iFirstStrip,iLastStrip,iStripCircuit,fYpos11);
  
  //--- second plane 
  ichamber = 12;
  
  detElemId = DetElemId(ichamber,side,iline);
  seg = fSegFactory->CreateMpSegmentation(detElemId, icathode);  

  // second plane middle part
  Int_t iFirstStripMiddle = FirstStrip(localBoard->GetName());
  Int_t iLastStripMiddle = iFirstStrip + 16;
  iStripCircuit=8;
  FillXstrips(seg,detElemId,icol, 
              iFirstStripMiddle,iLastStripMiddle,iStripCircuit,fYpos21);
  
  // second plane upper part
  if (zeroUp == 0) { // something up
    Int_t iFirstStripUp;
    Int_t iLastStripUp;
    Int_t icolUp=icol;
    // check if we need to move to another detElemId
    AliMpPad pad = seg->PadByIndices(AliMpIntPair(icol-1,iLastStripMiddle+1),kFALSE);
    if (pad.IsValid()) { // upper strips within same detElemId
      iFirstStripUp = iLastStripMiddle;
      iLastStripUp = iFirstStripUp + 8;
      //	    icolUp = icol;
    } else {             // upper strips in another detElemId
      detElemId = DetElemId(ichamber,side,iline+1); // get detElemId
      seg = fSegFactory->CreateMpSegmentation(detElemId, icathode);  

      iFirstStripUp = 0;
      iLastStripUp = iFirstStripUp + 8;
      if (iline == 4) icolUp = icol - 1; // special case
      if (iline == 5) icolUp = icol + 1; // special case
    }
    
    iStripCircuit=24;
    FillXstrips(seg,detElemId,icolUp, 
                iFirstStripUp,iLastStripUp,iStripCircuit,fYpos21);
    
    // fill strip between middle and upper part
    fYpos21[47]=(fYpos21[46]+fYpos21[48])/2.;
  } // end of something up
  
  // restore current detElemId & segmentation
  detElemId = DetElemId(ichamber,side,iline); 
  seg = fSegFactory->CreateMpSegmentation(detElemId, icathode);  

  // second plane lower part
  if (zeroDown == 0) { // something down
    Int_t iFirstStripDo;
    Int_t iLastStripDo;
    Int_t icolDo=icol;
    
    // check if we need to move to another detElemId      
    AliMpPad pad = seg->PadByIndices(AliMpIntPair(icol-1,iFirstStripMiddle-1),kFALSE);
    if (pad.IsValid()) { // lower strips within same detElemId
      iFirstStripDo = iFirstStripMiddle - 8;
      iLastStripDo = iFirstStripDo + 8;	      
      //	    icolDo = icol;
    } else {             // lower strips in another detElemId 
      detElemId = DetElemId(ichamber,side,iline-1); // get detElemId
      seg = fSegFactory->CreateMpSegmentation(detElemId, icathode);  

      // get iFirstStrip in this module 
      const AliMpTriggerSegmentation* trig = (AliMpTriggerSegmentation*)(seg);
      const AliMpTrigger* t = trig->Slat();
      const AliMpSlat* slat = t->GetLayer(0);
      if (iline == 5) icolDo = icol + 1; // special case
      if (iline == 6) icolDo = icol - 1; // special case	    
      const AliMpPCB* pcb = slat->GetPCB(icolDo-1);
      iFirstStripDo = (pcb->Iymax() + 1) - 8;
      iLastStripDo =  iFirstStripDo + 8;
    }  
    
    iStripCircuit=0;
    FillXstrips(seg,detElemId,icolDo, 
                iFirstStripDo,iLastStripDo,iStripCircuit,fYpos21);
    
    // fill strip between middle and upper part
    fYpos21[15]=(fYpos21[14]+fYpos21[16])/2.;
  } // end of something down
  
}

//----------------------------------------------------------------------
void AliMUONTriggerCircuitNew::LoadXPos(const AliMUONTriggerCrateStore& crates)
{
/// fill fXpos11 -> x position of Y strips for the first plane only
/// fXpos11 contains the x position of Y strip for the current circuit
/// taking into account whether or nor not part(s) of the circuit
/// (middle, up or down) has(have) 16 strips (handdled by means of switchs)
  
  const AliMUONLocalTriggerBoard* localBoard = crates.LocalBoard(fILocalBoard);
  
  if (!localBoard)
  {
    AliError(Form("Did not get localboard %d",fILocalBoard));
    return;
  }
  StdoutToAliDebug(1,localBoard->Print("CONF"););
  
  Int_t ichamber = 10;
  Int_t icathode = 1;
  
  Int_t x2u = localBoard->GetSwitch(AliMUONLocalTriggerBoard::kX2u);
  Int_t x2m = localBoard->GetSwitch(AliMUONLocalTriggerBoard::kX2m);
  Int_t x2d = localBoard->GetSwitch(AliMUONLocalTriggerBoard::kX2d);
  Int_t zeroAllYLSB = localBoard->GetSwitch(AliMUONLocalTriggerBoard::kZeroAllYLSB);
  Int_t iStripCircuit = 0;
  Int_t iFirstStrip = 0;
  Int_t iLastStrip = 0;
  Bool_t doubling = kFALSE;
  
  const AliMpVSegmentation* seg;

  char side;
  Int_t iline, icol;
  
  DecodeBoardName(localBoard->GetName(),side,iline,icol);
  
  Int_t detElemId=DetElemId(ichamber,side,iline); // get detElemId
  seg = fSegFactory->CreateMpSegmentation(detElemId, icathode);  

  // check if one needs a strip doubling or not
  if ( (x2u == 1 || x2m == 1 || x2d == 1) && x2m == 1) doubling = kTRUE;
  
  // check if one starts at strip = 0 or 8 (boards 26-29 and 143-146)
  if (zeroAllYLSB == 1) iStripCircuit = 8;
  
  // get iFirstStrip in this module 
  const AliMpTriggerSegmentation* trig = (AliMpTriggerSegmentation*)(seg);
  const AliMpTrigger* t = trig->Slat();
  const AliMpSlat* slat = t->GetLayer(0);
  
  const AliMpPCB* pcb = slat->GetPCB(icol-1);
  iFirstStrip = pcb->Ixmin();
  
  if (doubling) iLastStrip = iFirstStrip + 8;
  else iLastStrip = iFirstStrip + 16;
  
  FillYstrips(seg,detElemId, 
              iFirstStrip,iLastStrip,iStripCircuit,doubling);  
}

//----------------------------------------------------------------------
void AliMUONTriggerCircuitNew::FillYstrips(
                                           const AliMpVSegmentation* seg,
                                           const Int_t detElemId, 
                                           const Int_t iFirstStrip, const Int_t iLastStrip, Int_t liStripCircuit,
                                           const Bool_t doubling)
{    
/// fill
  Double_t xyGlobal[4]={0.,0.,0.,0.};
  for (Int_t istrip=iFirstStrip; istrip<iLastStrip; istrip++) {
    AliMpPad pad = seg->PadByIndices(AliMpIntPair(istrip,0),kTRUE);
    if ( !pad.IsValid() )
    {
	StdoutToAliError(cout << "Pad not found in seg " << endl;
			 seg->Print();
			 cout << " ix,iy=" << istrip << "," << 0 << endl;
			 );
    }
    XYGlobal(detElemId,pad,xyGlobal);
    
    if (!doubling) {	    
      fXpos11[liStripCircuit]=xyGlobal[0];
    } else if (doubling) {	    
      fXpos11[2*liStripCircuit]=TMath::Sign(1.,xyGlobal[0]) * 
	(TMath::Abs(xyGlobal[0]) - xyGlobal[2]/2.);
      fXpos11[2*liStripCircuit+1]=TMath::Sign(1.,xyGlobal[0]) *
        (TMath::Abs(fXpos11[2*liStripCircuit]) + xyGlobal[2]); 
    }	
    liStripCircuit++;
  }    
}

//----------------------------------------------------------------------
void AliMUONTriggerCircuitNew::FillXstrips(
                                           const AliMpVSegmentation* seg,
                                           const Int_t detElemId, const Int_t icol, 
                                           const Int_t iFirstStrip, const Int_t iLastStrip, 
                                           Int_t liStripCircuit, Float_t *tab)
{    
/// fill 
  Double_t xyGlobal[4]={0.,0.,0.,0.};
  for (Int_t istrip=iFirstStrip; istrip<iLastStrip; istrip++) {
    AliMpPad pad = seg->PadByIndices(AliMpIntPair(icol-1,istrip),kTRUE);
    if ( !pad.IsValid() )
    {
      StdoutToAliError(cout << "Pad not found in seg " << endl;
                       seg->Print();
                       cout << " ix,iy=" << icol-1 << "," << istrip << endl;
                       );
    }
    
    XYGlobal(detElemId,pad,xyGlobal);
    
    tab[2*liStripCircuit]=xyGlobal[1];
    if (istrip!=(iLastStrip-1)) tab[2*liStripCircuit+1]=xyGlobal[1]+xyGlobal[3];
    liStripCircuit++;
  }    
}

//----------------------------------------------------------------------
//--- methods which return member data related info
//----------------------------------------------------------------------
Float_t AliMUONTriggerCircuitNew::GetY11Pos(Int_t istrip) const {
/// returns Y position of X strip istrip in MC11
  return fYpos11[istrip];
}
//----------------------------------------------------------------------
Float_t AliMUONTriggerCircuitNew::GetY21Pos(Int_t istrip) const {
/// returns Y position of X strip istrip in MC21
  return fYpos21[istrip];
}
//----------------------------------------------------------------------
Float_t AliMUONTriggerCircuitNew::GetX11Pos(Int_t istrip) const {
/// returns X position of Y strip istrip in MC11
  return fXpos11[istrip];
}
//----------------------------------------------------------------------
//--- end of methods which return member data related info
//----------------------------------------------------------------------
/* removed tmp
void AliMUONTriggerCircuitNew::dump(const char* what, const Float_t* array, Int_t size)
{
  cout << what << " " << endl;
  for ( Int_t i = 0; i < size; ++i )
  {
    cout << array[i] << " , ";
  }
  cout << endl;
}

void AliMUONTriggerCircuitNew::dump(const char* what, const Int_t* array, Int_t size)
{
  cout << what << " " << endl;
  for ( Int_t i = 0; i < size; ++i )
  {
    cout << array[i] << " , ";
  }
  cout << endl;
}

//_____________________________________________________________________________
void AliMUONTriggerCircuitNew::Print(Option_t* ) const
{
  cout << "IdCircuit " << fILocalBoard << " X2m,X2ud=" << fX2m << ","
  << fX2ud;
  for ( Int_t i = 0; i < 2; ++i )
  {
    cout << " OrMud[" << i << "]=" << fOrMud[i];
  }
  cout << endl;
  dump("Xpos11",fXpos11,16);
  dump("Ypos11",fYpos11,31);
  dump("Ypos21",fYpos21,63);
  for ( Int_t i = 0; i < 4; ++i )
  {
    char s[80];
    sprintf(s,"Xcode[%d]",i);
    dump(s,fXcode[i],32);
    sprintf(s,"Ycode[%d]",i);
    dump(s,fYcode[i],32);
  }
  // Int_t fILocalBoard;            // circuit Id number
  //  Int_t fX2m;                  // internal info needed by TriggerDecision
  //  Int_t fX2ud;                 // internal info needed by TriggerDecision
  //  Int_t fOrMud[2];             // internal info needed by TriggerDecision
  //  Int_t fXcode[4][32];         // code of X strips
  //  Int_t fYcode[4][32];         // code of Y strips 
  //  Float_t fXpos11[16];         // X position of Y strips in MC11
  //  Float_t fYpos11[31];         // Y position of X strips in MC11
  //  Float_t fYpos21[63];         // Y position of X strips in MC21
  
}
removed tmp*/

//----------------------------------------------------------------------
Int_t AliMUONTriggerCircuitNew::DetElemId(Int_t ichamber, char side, Int_t iline)
{
/// returns detection element Id for chamber iChamber, side side and line iline
  Int_t itmp=0;    
  if ( side == 'R' ) {	       // right side 
    switch (iline) // (from 1 to 9, from bottom to top)
    {
      case 1:
        itmp = 14;
        break;	
      case 2:
        itmp = 15;
        break;		
      case 3:
        itmp = 16;
        break;
      case 4:
        itmp = 17;
        break;	    
      case 5:
        itmp = 0;
        break;		
      case 6:
        itmp = 1;
        break;		
      case 7:
        itmp = 2;
        break;		
      case 8:
        itmp = 3;
        break;		
      case 9:
        itmp = 4;
        break;
    }	
  } else if ( side == 'L' ) { // left side	    
    switch (iline) // (from 1 to 9, from bottom to top)
    {
      case 1:
        itmp = 13;
        break;		
      case 2:
        itmp = 12;
        break;		
      case 3:
        itmp = 11;
        break;		
      case 4:
        itmp = 10;
        break;		
      case 5:
        itmp = 9;
        break;		
      case 6:
        itmp = 8;
        break;		
      case 7:
        itmp = 7;
        break;		
      case 8:
        itmp = 6;
        break;
      case 9:
        itmp = 5;
        break;		
    }
  }    
  return ((ichamber+1)*100)+itmp;      
}

//----------------------------------------------------------------------
Int_t
AliMUONTriggerCircuitNew::DetElemId(Int_t iChamber, const char* boardName)
{
/// returns detection element Id for chamber iChamber and board boardName
  char side;
  Int_t iline;
  Int_t icol;

  DecodeBoardName(boardName, side, iline, icol);

  return DetElemId(iChamber,side,iline);
}

//----------------------------------------------------------------------
void
AliMUONTriggerCircuitNew::DecodeBoardName(const char* boardName,
                                          char& side,
                                          Int_t& iLine,
                                          Int_t& iCol)
{
/// get side, line and col from board boardName
/// note: icol = icol -1 for iline = 5 w.r.t other ilines
  side = boardName[0];
  iLine = boardName[4] - '0';
  iCol = boardName[2] - '0';
  if ( iLine == 5 ) --iCol;
}

//----------------------------------------------------------------------
Int_t
AliMUONTriggerCircuitNew::FirstStrip(const char* boardName)
{
/// returns the first strip from mapping for board boardName
/// take care of special case for boards RC1L6B12 & LC1L6B12
  Int_t iFirstStrip = -1;
  Int_t boardNumber = atoi(boardName+6);
  char side;
  Int_t iline, icol;
  DecodeBoardName(boardName,side,iline,icol);
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
  if (icol == 1 && iline == 6) iFirstStrip = iFirstStrip + 16; // special case
  return iFirstStrip;
}

//----------------------------------------------------------------------
void AliMUONTriggerCircuitNew::XYGlobal(
                                        Int_t detElemId, const AliMpPad& pad,
                                        Double_t xyGlobal[4])
{
/// returns pad x & y positions and x & y pad dimensions in global coordinates
/// note: no need for transformation for pad dimensions
  
  // get the pad position and dimensions
  Double_t xl1 = pad.Position().X();
  Double_t yl1 = pad.Position().Y();
  xyGlobal[2] = pad.Dimensions().X(); // half size!
  xyGlobal[3] = pad.Dimensions().Y(); // half size! 
  Double_t zg1=0;
  
  // positions from local to global 
  fTransformer->Local2Global(detElemId, xl1, yl1, 0, 
                                 xyGlobal[0], xyGlobal[1], zg1);
}
