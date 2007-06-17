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

//*-- Author: Rachid Guernane (LPCCFd)
//*   Manager class for muon trigger electronics
//*   Client of trigger board classes
//*
//*

#include "AliLoader.h"
#include "AliLog.h"
#include "AliMUON.h" 
#include "AliMUONCalibrationData.h"
#include "AliMUONVDigit.h"
#include "AliMUONVDigitStore.h"
#include "AliMUONGlobalTrigger.h"
#include "AliMUONGlobalTriggerBoard.h"
#include "AliMUONLocalTrigger.h"
#include "AliMUONLocalTriggerBoard.h"
#include "AliMUONRegionalTrigger.h"
#include "AliMUONRegionalTriggerBoard.h"
#include "AliMUONTriggerCrate.h"
#include "AliMUONTriggerCrateStore.h"
#include "AliMUONTriggerElectronics.h"
#include "AliMUONVTriggerStore.h"
#include "AliMUONVCalibParam.h"
#include "AliMpCathodType.h"
#include "AliMpDEManager.h"
#include "AliMpSegmentation.h"
#include "AliMpVSegmentation.h"
#include "AliMpCathodType.h"

#include "AliMUONTriggerGUIboard.h"

#include "AliLog.h"
#include "AliLoader.h"
#include "AliRun.h"
#include <TBits.h>
#include <TSystem.h>


/// \cond CLASSIMP
ClassImp(AliMUONTriggerElectronics)
/// \endcond

//___________________________________________
AliMUONTriggerElectronics::AliMUONTriggerElectronics(AliMUONCalibrationData* calibData) 
: TObject(),
  fSourceFileName(),
  fCrates(new AliMUONTriggerCrateStore),
  fGlobalTriggerBoard(new AliMUONGlobalTriggerBoard)
{
 /// CONSTRUCTOR
///
  SetDataSource();
  Factory(calibData);
  LoadMasks(calibData);
}

//___________________________________________
AliMUONTriggerElectronics::~AliMUONTriggerElectronics()
{
/// DESTRUCTOR
///
  delete fGlobalTriggerBoard;
  delete fCrates;
}

//___________________________________________
void AliMUONTriggerElectronics::Factory(AliMUONCalibrationData* calibData)
{  
 /// BUILD ALL ELECTRONICS
 ///

// get coinc44 from AliMUON (added 12/09/06)
  AliMUON *pMUON  = (AliMUON*)gAlice->GetModule("MUON");
  Int_t coinc44 = pMUON->GetTriggerCoinc44();
  if (coinc44 != 0 && coinc44 != 1) {
      AliFatal("Coinc 44 should be equal to 0 or 1");
      return;
  }

  fCrates->ReadFromFile(gSystem->ExpandPathName(fSourceFileName.Data()));
  
  if ( !calibData ) return;
  
  AliMUONTriggerLut* lut = calibData->TriggerLut();
  
  if (!lut) return;
  
  AliMUONLocalTriggerBoard* localBoard;
  
  fCrates->FirstLocalBoard();
  
  while ( (localBoard=fCrates->NextLocalBoard()) )
  {
    localBoard->SetLUT(lut);
    localBoard->SetCoinc44(coinc44);
  }
}

//___________________________________________
void AliMUONTriggerElectronics::Feed(const AliMUONVDigitStore& digitStore)
{
  /// FILL INPUTS
  ///
  
  TIter next(digitStore.CreateTriggerIterator());
  AliMUONVDigit* mdig;
  
  while ( ( mdig = static_cast<AliMUONVDigit*>(next()) ) )
  {      
    //       CHECKME ! The TrackCharge is not ok with new digitizerV3 !
    //			for (Int_t ichg=0; ichg<10; ichg++) schg += mdig->TrackCharge(ichg);
    Int_t ichamber = AliMpDEManager::GetChamberId(mdig->DetElemId());
    Int_t schg = (Int_t)(mdig->Charge() + 0.5);
    
    //       APPLY CONDITION ON SOFT BACKGROUND	
    Int_t tchg = schg - (Int_t(schg/10))*10;	
    
    if (schg<=10 || tchg>0) 
    {
      Int_t detElemId  = mdig->DetElemId();
      Int_t cathode    = mdig->Cathode();
      
      const AliMpVSegmentation *seg = 
        AliMpSegmentation::Instance()
        ->GetMpSegmentation(detElemId,AliMp::GetCathodType(cathode));
      
      Int_t ix = mdig->PadX(), iy = mdig->PadY();
      
      AliDebug(3,Form("cathode %d ix %d iy %d ",cathode,ix,iy));
      
      AliMpPad pad = seg->PadByIndices(AliMpIntPair(ix,iy),kTRUE);
      
      for (Int_t i=0; i<pad.GetNofLocations(); i++) 
      {
        AliMpIntPair location = pad.GetLocation(i);
        
        Int_t nboard = location.GetFirst();
        
        Int_t ibitxy = location.GetSecond();
        
        AliMUONLocalTriggerBoard *b = fCrates->LocalBoard(nboard);
        
        if (b) 
        {
          if (cathode && b->GetSwitch(6)) ibitxy += 8;
          
          b->SetbitM(ibitxy,cathode,ichamber-10);
        }
        else
        {
          AliError(Form("Could not get local board number %d",b->GetNumber()));
        }
      }
    }		
  }
  
  // Particular case of the columns with 22 local boards (2R(L) 3R(L))   
  AliMUONTriggerCrate *crate = 0x0; TObjArray *bs = 0x0;
  
  char *scratess[4] = {  "2R",   "2L",   "3L",   "3R"}; 
  char *scratesd[4] = {"2-3R", "2-3L", "2-3L", "2-3R"}; 
  Int_t    slotf[4] = {     2,      2,     10,     10}; 
  Int_t    slotd[4] = {     1,      1,      9,      9}; 
  
  for (Int_t i = 0; i < 4; i++)
  {
    crate = fCrates->Crate(scratess[i]); 
    bs = crate->Boards();
    AliMUONLocalTriggerBoard *desybb = (AliMUONLocalTriggerBoard*)bs->At(14);
    AliMUONLocalTriggerBoard *fromcb = (AliMUONLocalTriggerBoard*)bs->At(15);
    AliMUONLocalTriggerBoard *desxbb = (AliMUONLocalTriggerBoard*)bs->At(16);
    
    crate = fCrates->Crate(scratesd[i]); 
    bs = crate->Boards();
    AliMUONLocalTriggerBoard *frombb = (AliMUONLocalTriggerBoard*)bs->At(slotf[i]);
    AliMUONLocalTriggerBoard *desycb = (AliMUONLocalTriggerBoard*)bs->At(slotd[i]);
    
    UShort_t cX[2];
    
    //    COPY X3-4 FROM BOARD  2 OF CRATE 2-3 TO BOARD 16 OF CRATE 2
    //    COPY X3-4 FROM BOARD 10 OF CRATE 2-3 TO BOARD 16 OF CRATE 3
    frombb->GetX34(cX); desxbb->SetX34(cX);
    
    //    COPY X3-4 FROM BOARD 15 OF CRATE 2 TO BOARD 1 OF CRATE 2-3
    //    COPY X3-4 FROM BOARD 15 OF CRATE 3 TO BOARD 9 OF CRATE 2-3
    fromcb->GetX34(cX); desycb->SetX34(cX);
    
    UShort_t cY[4];
    
    desybb->GetY(cY); frombb->SetY(cY);
    
    frombb->GetY(cY); desxbb->SetY(cY);
    fromcb->GetY(cY); desycb->SetY(cY);
  }
  
  // FILL UP/DOWN OF CURRENT BOARD (DONE VIA J3 BUS IN REAL LIFE)
  AliMUONTriggerCrate* cr;
  
  fCrates->FirstCrate();
  
  while ( ( cr = fCrates->NextCrate() ) )
  {            
    TObjArray *boards = cr->Boards();
    
    for (Int_t j = 1; j < boards->GetEntries()-1; j++)
    {
      TObject *o = boards->At(j);
			
      if (!o) break;
			
      AliMUONLocalTriggerBoard *currboard = (AliMUONLocalTriggerBoard*)o;
			
      AliMUONLocalTriggerBoard *neighbour = (AliMUONLocalTriggerBoard*)boards->At(j+1);
			
      UShort_t cXY[2][4];
			
      if (j==1) {neighbour->GetXY(cXY); currboard->SetXYU(cXY);}
			
      //       LAST BOARD IN THE CRATE HAS NO UP EXCEPT FOR CRATES 2 & 3
      if (j < boards->GetEntries()-2)  
      {
	      AliMUONLocalTriggerBoard *nextboard = (AliMUONLocalTriggerBoard*)boards->At(j+2);
				
	      currboard->GetXY(cXY); neighbour->SetXYD(cXY);
	      nextboard->GetXY(cXY); neighbour->SetXYU(cXY);
				
	      if (j==boards->GetEntries()-3) {neighbour->GetXY(cXY); nextboard->SetXYD(cXY);}
      }
    }
  }
}


//___________________________________________
void AliMUONTriggerElectronics::Feed(UShort_t pattern[2][4])
{
  /// FILL INPUTS
  ///
  AliMUONTriggerCrate* cr;
   
   fCrates->FirstCrate();
   
   while ( ( cr = fCrates->NextCrate() ) )
   {                 
     TObjArray *boards = cr->Boards();
     
     for (Int_t j = 1; j < boards->GetEntries(); j++)
     {
       TObject *o = boards->At(j);
       
       if (!o) break;
       
       AliMUONLocalTriggerBoard *board = (AliMUONLocalTriggerBoard*)o;
       
       board->SetXY(pattern);
     }
   }
}

//___________________________________________
void AliMUONTriggerElectronics::DumpOS()
{
/// DUMP IN THE OLD WAY
///
   for (Int_t i= 0; i < 234;i++)
   {
      AliMUONLocalTriggerBoard *board = fCrates->LocalBoard(i);

      if (board) board->Scan("ALL");
   }
}

//___________________________________________
void AliMUONTriggerElectronics::Scan(Option_t *option)
{
  /// SCAN
  ///

  AliMUONTriggerCrate* cr;
  
  fCrates->FirstCrate();
  
  while ( ( cr = fCrates->NextCrate() ) )
  {                
    TObjArray *boards = cr->Boards();
    
    for (Int_t j = 0; j < boards->GetEntries(); j++)
    {
      TObject *o = boards->At(j);
      
      TString op = option;
      
      Bool_t cdtion = kFALSE;
      
      if (op.Contains("LOCAL"))    cdtion = o->IsA() == AliMUONLocalTriggerBoard::Class();
      if (op.Contains("REGIONAL")) cdtion = o->IsA() == AliMUONRegionalTriggerBoard::Class();
      if (op.Contains("GLOBAL"))   cdtion = o->IsA() == AliMUONGlobalTriggerBoard::Class();
      
      if (!o || !cdtion) continue;
      
      AliMUONLocalTriggerBoard *board = (AliMUONLocalTriggerBoard*)o;
      
      board->Scan();
    }
  }
}

//___________________________________________
void AliMUONTriggerElectronics::Reset()
{
  /// RESET
  ///
  
   AliMUONTriggerCrate* cr;
   
   fCrates->FirstCrate();
   
   while ( ( cr = fCrates->NextCrate() ) )
   {            
      TObjArray *boards = cr->Boards();
            
      for (Int_t j=0; j<boards->GetEntries(); j++)
      {     
         AliMUONTriggerBoard *b = (AliMUONTriggerBoard*)boards->At(j);

         if (b) b->Reset();
      }
   }
}

//_______________________________________________________________________
void AliMUONTriggerElectronics::LoadMasks(AliMUONCalibrationData* calibData)
{
  /// LOAD MASKS FROM CDB
  

  // SET MASKS
  
  AliMUONTriggerCrate* cr;
  
  fCrates->FirstCrate();
  
  Int_t irb(0);
  
  while ( ( cr = fCrates->NextCrate() ) )
  {            
    TObjArray *boards = cr->Boards();
    
    AliMUONRegionalTriggerBoard *regb =
      (AliMUONRegionalTriggerBoard*)boards->At(0);

    AliMUONVCalibParam* regionalBoardMasks = calibData->RegionalTriggerBoardMasks(irb);
    
    for ( Int_t i = 0; i < regionalBoardMasks->Size(); ++i )
    {
      UShort_t rmask = static_cast<UShort_t>(regionalBoardMasks->ValueAsInt(i) & 0x3F);
      regb->Mask(i,rmask);
    }
    
    for (Int_t j = 1; j < boards->GetEntries(); j++)
    {
      AliMUONLocalTriggerBoard *b = (AliMUONLocalTriggerBoard*)boards->At(j);
      
      Int_t cardNumber = b->GetNumber();
      
      if (cardNumber) // interface board are not interested
      {
        AliMUONVCalibParam* localBoardMasks = calibData->LocalTriggerBoardMasks(cardNumber);
        for ( Int_t i = 0; i < localBoardMasks->Size(); ++i )
        {
          UShort_t lmask = static_cast<UShort_t>(localBoardMasks->ValueAsInt(i) & 0xFFFF);
          b->Mask(i,lmask);
        }
      }
    }
    ++irb;
  }
  
  AliMUONVCalibParam* globalBoardMasks = calibData->GlobalTriggerBoardMasks();
  for ( Int_t i = 0; i < globalBoardMasks->Size(); ++i )
  {
    UShort_t gmask = static_cast<UShort_t>(globalBoardMasks->ValueAsInt(i) & 0xFFF);
    fGlobalTriggerBoard->Mask(i,gmask);
  }
}


//___________________________________________
void AliMUONTriggerElectronics::LocalResponse()
{
/// \todo add comment
	
  AliMUONTriggerCrate* cr;
  
  fCrates->FirstCrate();
  
  while ( ( cr = fCrates->NextCrate() ) )
  {            
    
    TObjArray *boards = cr->Boards();
    
    AliMUONRegionalTriggerBoard *regb = (AliMUONRegionalTriggerBoard*)boards->At(0);
    
    UShort_t thisl[16]; for (Int_t j=0; j<16; j++) thisl[j] = 0;
  
    for (Int_t j = 1; j < boards->GetEntries(); j++)
    {     
	TObject *o = boards->At(j);
      
	if (!o) break;
      
	AliMUONLocalTriggerBoard *board = (AliMUONLocalTriggerBoard*)o;
      
	if (board) // check if empty slot
	{
	  board->Response();
				
	  UShort_t response = board->GetResponse();            
        
	  // CRATE CONTAINING INTERFACE BOARD
	  if (board->GetNumber() == 0) // copy boards
	  {
	    if ( response != 0 ) 
	      AliWarning(Form("Interface board %s in slot %d of crate %s has a non zero response",
					  board->GetName(),j,cr->GetName()));
	    AliDebug(1, Form("local slot %d, number %d in crate %s\n", j, board->GetNumber(), cr->GetName()));

	  }
        
	  thisl[j-1] = response;
	}
      }
    
    regb->SetLocalResponse(thisl);
  }
}

//___________________________________________
void AliMUONTriggerElectronics::RegionalResponse()
{
  /// Compute the response for all regional cards.
  AliMUONTriggerCrate* cr;
  
  fCrates->FirstCrate();
  
  while ( ( cr = fCrates->NextCrate() ) )
  {            
      TObjArray *boards = cr->Boards();

      AliMUONRegionalTriggerBoard *regb = (AliMUONRegionalTriggerBoard*)boards->At(0);
      
      if (regb) 
      {
         regb->Response();
      }  
   }
}

//___________________________________________
void AliMUONTriggerElectronics::GlobalResponse()
{
  /// Compute the global response

  UShort_t regional[16];
  
  AliMUONTriggerCrate* cr;
  
  fCrates->FirstCrate();
  Int_t irb(0);
  
  if ( !fCrates->NumberOfCrates() >= 16 ) 
  {
    AliFatal(Form("Something is wrong : too many crates %d",
                  fCrates->NumberOfCrates()));
  }
  
  while ( ( cr = fCrates->NextCrate() ) )
  {            
    AliMUONTriggerBoard* rb = 
      static_cast<AliMUONTriggerBoard*>(cr->Boards()->At(0));
    regional[irb] = rb->GetResponse();
    ++irb;
  }
  
  fGlobalTriggerBoard->SetRegionalResponse(regional);
  fGlobalTriggerBoard->Response();
}

//_______________________________________________________________________
void AliMUONTriggerElectronics::Digits2Trigger(const AliMUONVDigitStore& digitStore,
                                               AliMUONVTriggerStore& triggerStore)
{
  /// Main method to go from digits to trigger decision
  AliMUONRegionalTrigger pRegTrig;
  
  triggerStore.Clear();
  
  // RUN THE FULL BEE CHAIN
  Feed(digitStore);
  LocalResponse();
  RegionalResponse();      
  GlobalResponse();
  //    DumpOS();
	
  AliMUONTriggerCrate* cr;
  AliMUONLocalTrigger localTrigger;
  
  // stored in right order
  // do not used iterator order
  
  for (Int_t iSide = 0; iSide < 2; iSide++) // right & left side
  {            
    for (Int_t iReg = 0; iReg < 8; iReg++) // 8 crates/regional boards for each side.
    {
      cr = fCrates->Crate(iSide, iReg);     
      TObjArray *boards = cr->Boards();
      
      UInt_t regInpLpt = 0;
      UInt_t regInpHpt = 0;
      
      AliMUONRegionalTriggerBoard *regBoard = (AliMUONRegionalTriggerBoard*)boards->At(0);
      
      for (Int_t j = 1; j < boards->GetEntries(); j++)
      {     
        TObject *o = boards->At(j);
        
        if (!o) break;
        
        AliMUONLocalTriggerBoard *board = (AliMUONLocalTriggerBoard*)o;
        
        if (board) 
        {
          //          L0 TRIGGER
          // pcrochet 181206: MOOD needs ALL boards
          //	  if (board->Triggered())
          //	  {
          
          Int_t icirc = board->GetNumber();
          if (icirc != 0) { // pcrochet 181206: MOOD needs ALL boards
            
            localTrigger.SetLoCircuit(icirc);
            localTrigger.SetLoStripX(board->GetStripX11());
            localTrigger.SetLoDev(board->GetDev());
            localTrigger.SetLoSdev(board->GetSdev());
            localTrigger.SetLoTrigY(board->GetTrigY());
            localTrigger.SetLoStripY(board->GetStripY11());
            
            //             SAVE LUT OUTPUT 
            UShort_t response = board->GetResponse();
            localTrigger.SetLoHpt((response & 12) >> 2);
            localTrigger.SetLoLpt(response &  3);
            
            // calculates regional inputs from local for the moment
            UInt_t hPt = (response >> 2) & 0x3;
            UInt_t lPt =  response       & 0x3;
            
            regInpHpt |= hPt << (30 - (j-1)*2);
            regInpLpt |= lPt << (30 - (j-1)*2);
            
            TBits rrr;
            rrr.Set(6,&response);	  
            
            //             SAVE BIT PATTERN
            localTrigger.SetX1Pattern(board->GetXY(0,0));
            localTrigger.SetX2Pattern(board->GetXY(0,1));
            localTrigger.SetX3Pattern(board->GetXY(0,2));
            localTrigger.SetX4Pattern(board->GetXY(0,3));
            
            localTrigger.SetY1Pattern(board->GetXY(1,0));
            localTrigger.SetY2Pattern(board->GetXY(1,1));
            localTrigger.SetY3Pattern(board->GetXY(1,2));
            localTrigger.SetY4Pattern(board->GetXY(1,3));
            
            //             ADD A NEW LOCAL TRIGGER          
            triggerStore.Add(localTrigger);  
            
          }
          }
        }
      pRegTrig.SetId(iReg + 8*iSide);
      pRegTrig.SetLocalOutput(regInpLpt, 0);
      pRegTrig.SetLocalOutput(regInpHpt, 1);
      pRegTrig.SetOutput(regBoard->GetResponse());
      
      triggerStore.Add(pRegTrig);  
      }
    }
  
  // GLOBAL TRIGGER INFORMATION
  UShort_t global = fGlobalTriggerBoard->GetResponse();
  
  AliMUONGlobalTrigger globalTrigger;
  
  globalTrigger.SetFromGlobalResponse(global);
  // ADD A LOCAL TRIGGER IN THE LIST 
  triggerStore.SetGlobal(globalTrigger);
  
  // NOW RESET ELECTRONICS
  Reset();
}

//_______________________________________________________________________
void AliMUONTriggerElectronics::FeedBoardsGUI(TObjArray *guibs)
{
  /// feed digits from board objects from the TriggerGUI, with values
  /// read from a file or set interactively in the GUI
  /// 

  // adaptated from FeedM()

  AliMUONTriggerGUIboard* board;
  Int_t cathode, nstripX, nstripY, ix, iy, detElemId0, detElemId, charge;
  Int_t iX1, iY1, schg, tchg;
  Bool_t triggerBgn;

  for (Int_t ib = 0; ib < 234; ib++) {

    board = (AliMUONTriggerGUIboard*)guibs->At(ib);
    if (board == 0) continue;

    detElemId0 = board->GetDetElemId();

    nstripX = board->GetNStripX();
    nstripY = board->GetNStripY();

    for (Int_t ichamber = 11; ichamber <= 14; ichamber++) {
  
      detElemId = ichamber * 100 + detElemId0;

      // x strips
      cathode = 0;
      for (Int_t isx = 0; isx < nstripX; isx++) {

	charge = (Int_t)board->GetXDig(ichamber-11,isx);
	if (charge) {

	  triggerBgn = kFALSE;
	  schg = (Int_t)(charge + 0.5);
	  // APPLY CONDITION ON SOFT BACKGROUND	
	  tchg = schg - (Int_t(schg/10))*10;	
	  if (schg<=10 || tchg>0) {
	    triggerBgn = kFALSE;
	  } else {
	    triggerBgn = kTRUE;
	  }
	  if (triggerBgn) continue;

	  //printf("MT %2d SX %2d \n",ichamber,isx);

	  ix  = board->GetXSix();
	  iY1 = board->GetXSiy1();
	  iy  = isx + iY1;

	  //printf("X: CH %1d B %3d ID %4d ix %2d iy %2d \n",ichamber,ib,detElemId,ix,iy);
	  
	  const AliMpVSegmentation* seg = 
	    AliMpSegmentation::Instance()->GetMpSegmentation(detElemId,AliMp::GetCathodType(cathode));

	  AliMpPad pad = seg->PadByIndices(AliMpIntPair(ix,iy),kTRUE);

	  if (!pad.IsValid()) printf("Invalid pad! \n");

	  for (Int_t i=0; i<pad.GetNofLocations(); i++) {

	    AliMpIntPair location = pad.GetLocation(i);
	    Int_t nboard = location.GetFirst();
	    if (nboard != board->GetIdCircuit()) continue;
	    Int_t ibitxy = location.GetSecond();
	    
	    //printf("FeedGUI x (%2d): ix %d iy %d detElemId %d \n",ichamber,ix,iy,detElemId);

	    AliMUONLocalTriggerBoard *b = fCrates->LocalBoard(nboard);
	    
	    if (b) {
	      if (cathode && b->GetSwitch(6)) ibitxy += 8;
	      
	      //printf("Feed x-digits in board: %d cha %d s %d \n",nboard,ichamber,isx);
	      b->SetbitM(ibitxy,cathode,ichamber-11);
	      
	    } else {
	      AliError(Form("Could not get local board number %d",b->GetNumber()));
	    }	      
	    
	  }

	}

      }

      // y strips
      cathode = 1;
      for (Int_t isy = 0; isy < nstripY; isy++) {

	charge = board->GetYDig(ichamber-11,isy);
	if (charge) {

	  triggerBgn = kFALSE;
	  schg = (Int_t)(charge + 0.5);
	  // APPLY CONDITION ON SOFT BACKGROUND	
	  tchg = schg - (Int_t(schg/10))*10;	
	  if (schg<=10 || tchg>0) {
	    triggerBgn = kFALSE;
	  } else {
	    triggerBgn = kTRUE;
	  }
	  if (triggerBgn) continue;

	  //printf("MT %2d SY %2d \n",ichamber,isy);

	  iX1 = board->GetYSix1();
	  ix  = isy + iX1;
	  iy  = board->GetYSiy();

	  //printf("Y: CH %1d B %3d ID %4d ix %2d iy %2d \n",ichamber,ib,detElemId,ix,iy);
	  
	  const AliMpVSegmentation* seg = 
	    AliMpSegmentation::Instance()->GetMpSegmentation(detElemId,AliMp::GetCathodType(cathode));

	  AliMpPad pad = seg->PadByIndices(AliMpIntPair(ix,iy),kTRUE);

	  if (!pad.IsValid()) printf("Invalid pad! \n");;

	  for (Int_t i=0; i<pad.GetNofLocations(); i++) {

	    AliMpIntPair location = pad.GetLocation(i);
	    Int_t nboard = location.GetFirst();
	    if (nboard != board->GetIdCircuit()) continue;
	    Int_t ibitxy = location.GetSecond();
	    
	    //printf("FeedGUI y (%2d): ix %d iy %d detElemId %d \n",ichamber,ix,iy,detElemId);

	    AliMUONLocalTriggerBoard *b = fCrates->LocalBoard(nboard);
	    
	    if (b) {
	      if (cathode && b->GetSwitch(6)) ibitxy += 8;
	      
	      //printf("Feed y-digits in board: %d cha %d s %d \n",nboard,ichamber,isy);
	      b->SetbitM(ibitxy,cathode,ichamber-11);
	      
	    } else {
	      AliError(Form("Could not get local board number %d",b->GetNumber()));
	    }

	  }

	}

      }

    }
    
  }

  // ... the rest from FeedM()

  // Particular case of the columns with 22 local boards (2R(L) 3R(L))   
  AliMUONTriggerCrate *crate = 0x0; TObjArray *bs = 0x0;

  char *scratess[4] = {  "2R",   "2L",   "3L",   "3R"}; 
  char *scratesd[4] = {"2-3R", "2-3L", "2-3L", "2-3R"}; 
  Int_t    slotf[4] = {     2,      2,     10,     10}; 
  Int_t    slotd[4] = {     1,      1,      9,      9}; 

  for (Int_t i = 0; i < 4; i++)
  {
      crate = fCrates->Crate(scratess[i]); 
      bs = crate->Boards();
      AliMUONLocalTriggerBoard *desybb = (AliMUONLocalTriggerBoard*)bs->At(14);
      AliMUONLocalTriggerBoard *fromcb = (AliMUONLocalTriggerBoard*)bs->At(15);
      AliMUONLocalTriggerBoard *desxbb = (AliMUONLocalTriggerBoard*)bs->At(16);

      crate = fCrates->Crate(scratesd[i]); 
      bs = crate->Boards();
      AliMUONLocalTriggerBoard *frombb = (AliMUONLocalTriggerBoard*)bs->At(slotf[i]);
      AliMUONLocalTriggerBoard *desycb = (AliMUONLocalTriggerBoard*)bs->At(slotd[i]);

      UShort_t cX[2];

      //    COPY X3-4 FROM BOARD  2 OF CRATE 2-3 TO BOARD 16 OF CRATE 2
      //    COPY X3-4 FROM BOARD 10 OF CRATE 2-3 TO BOARD 16 OF CRATE 3
      frombb->GetX34(cX); desxbb->SetX34(cX);

      //    COPY X3-4 FROM BOARD 15 OF CRATE 2 TO BOARD 1 OF CRATE 2-3
      //    COPY X3-4 FROM BOARD 15 OF CRATE 3 TO BOARD 9 OF CRATE 2-3
      fromcb->GetX34(cX); desycb->SetX34(cX);

      UShort_t cY[4];

      desybb->GetY(cY); frombb->SetY(cY);

      frombb->GetY(cY); desxbb->SetY(cY);
      fromcb->GetY(cY); desycb->SetY(cY);
  }

  // FILL UP/DOWN OF CURRENT BOARD (DONE VIA J3 BUS IN REAL LIFE)
  AliMUONTriggerCrate* cr;
 
  fCrates->FirstCrate();
 
  while ( ( cr = fCrates->NextCrate() ) )
  {            
      TObjArray *boards = cr->Boards();
		
      for (Int_t j = 1; j < boards->GetEntries()-1; j++)
      {
	  TObject *o = boards->At(j);
			
	  if (!o) break;
			
	  AliMUONLocalTriggerBoard *currboard = (AliMUONLocalTriggerBoard*)o;
			
	  AliMUONLocalTriggerBoard *neighbour = (AliMUONLocalTriggerBoard*)boards->At(j+1);
			
	  UShort_t cXY[2][4];
			
	  if (j==1) {neighbour->GetXY(cXY); currboard->SetXYU(cXY);}
			
	  //       LAST BOARD IN THE CRATE HAS NO UP EXCEPT FOR CRATES 2 & 3
	  if (j < boards->GetEntries()-2)  
	  {
	      AliMUONLocalTriggerBoard *nextboard = (AliMUONLocalTriggerBoard*)boards->At(j+2);
				
	      currboard->GetXY(cXY); neighbour->SetXYD(cXY);
	      nextboard->GetXY(cXY); neighbour->SetXYU(cXY);
				
	      if (j==boards->GetEntries()-3) {neighbour->GetXY(cXY); nextboard->SetXYD(cXY);}
	  }
      }
  }
}

//_______________________________________________________________________
Int_t AliMUONTriggerElectronics::TriggerGUI(Int_t *trigInfo, Bool_t patt)
{
  /// trigger with digits from TriggerGUI and return local trigger information
  /// and optionally the strips pattern
  /// 

  Int_t nlo = 0;

  LocalResponse();
  RegionalResponse();      
  GlobalResponse();

  AliMUONTriggerCrate* cr;
 
  // stored in right order
  // do not used iterator order

  for (Int_t iSide = 0; iSide < 2; iSide++) // right & left side
  {            
    for (Int_t iReg = 0; iReg < 8; iReg++) // 8 crates/regional boards for each side.
    {
      cr = fCrates->Crate(iSide, iReg);     
      TObjArray *boards = cr->Boards();

      for (Int_t j = 1; j < boards->GetEntries(); j++)
      {     
	TObject *o = boards->At(j);
      
	if (!o) break;
      
	AliMUONLocalTriggerBoard *board = (AliMUONLocalTriggerBoard*)o;
      
	if (board) 
	{
	  //          L0 TRIGGER
	  if (board->Triggered())
	  {
          
	    if (patt) {
	      cout << "                                   " << endl;
	      cout << "Local trigger board P A T T E R N :" << endl;
	      cout << "-----------------------------------" << endl;
	      cout << "                                   " << endl;
	      board->Pattern();
	      board->Scan("RESPF");
	    }

	    Int_t icirc    = board->GetNumber();
	    Int_t loStripX = board->GetStripX11();
	    Int_t loStripY = board->GetStripY11();
	    Int_t loDev    = board->GetDev();

	    UShort_t response = board->GetResponse();
	    Int_t loHpt = (response & 12) >> 2;
	    Int_t loLpt = response &  3;
	    /*
	    cout << "TriggerGUI done!" << endl;

	    cout << "Circuit = "  << icirc    << endl;
	    cout << "LoStripX = " << loStripX << endl;
	    cout << "LoStripY = " << loStripY << endl;
	    cout << "LoDev = "    << loDev    << endl;
	    cout                              << endl;
	    */
	    trigInfo[6*nlo+0] = icirc;
	    trigInfo[6*nlo+1] = loStripX;
	    trigInfo[6*nlo+2] = loStripY;
	    trigInfo[6*nlo+3] = loDev;

	    trigInfo[6*nlo+4] = loLpt;
	    trigInfo[6*nlo+5] = loHpt;

	    nlo++;

	  }
	}
      }
    }
  }

  Reset();

  return nlo;
          
}

