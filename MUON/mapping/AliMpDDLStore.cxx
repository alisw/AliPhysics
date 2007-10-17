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

// $Id$
// $MpId: AliMpDDLStore.cxx,v 1.4 2006/05/24 13:58:34 ivana Exp $
// Category: management

//-----------------------------------------------------------------------------
// Class AliMpDDLStore
// --------------------
// The top container class for DDLs, det elements and bus patched
// It provides acces to DDL, det element and bus patches objects
// via various characteristics.
// Authors: Ivana Hrivnacova, IPN Orsay
//          Christian Finck, SUBATECH Nantes
//-----------------------------------------------------------------------------

#include "AliMpDDLStore.h"
#include "AliMpConstants.h"
#include "AliMpDEStore.h"
#include "AliMpDDL.h"
#include "AliMpFiles.h"
#include "AliMpHelper.h"
#include "AliMpDEManager.h"
#include "AliMpDetElement.h"
#include "AliMpBusPatch.h"
#include "AliMpTriggerCrate.h"
#include "AliMpLocalBoard.h"
#include "AliMpSegmentation.h"
#include "AliMpVSegmentation.h"
#include "AliMpStringObjMap.h"
#include "AliLog.h"

#include <Riostream.h>
#include <TList.h>
#include <TObjArray.h>
#include <TString.h>
#include <TObjString.h>
#include <TClass.h>

/// \cond CLASSIMP
ClassImp(AliMpDDLStore)
/// \endcond

AliMpDDLStore* AliMpDDLStore::fgInstance = 0;
const Int_t    AliMpDDLStore::fgkNofDDLs = 20;
const Int_t    AliMpDDLStore::fgkNofTriggerDDLs = 2;

//
// static methods
//

//______________________________________________________________________________
AliMpDDLStore* AliMpDDLStore::Instance(Bool_t warn)
{
/// Create the DDL store if it does not yet exist
/// and return its instance

  if ( ! fgInstance && warn  ) {
    AliWarningClass("DDL Store has not been loaded");
  }  
     
  return fgInstance;
}    

//______________________________________________________________________________
AliMpDDLStore* AliMpDDLStore::ReadData(Bool_t warn)
{
/// Load the DDL store from ASCII data files
/// and return its instance

  if ( fgInstance ) {
    if ( warn )
      AliWarningClass("DDL Store has been already loaded");
    return fgInstance;
  }  
  
  AliInfoClass("Reading DDL Store from ASCII files.");

  fgInstance = new AliMpDDLStore();
  return fgInstance;
}    

//
// ctors, dtor
//

//______________________________________________________________________________
AliMpDDLStore::AliMpDDLStore()
: TObject(),
  fDDLs(fgkNofDDLs+fgkNofTriggerDDLs), // FIXEME
  fBusPatches(true),
  fTriggerCrates(true),
  fLocalBoards(true),
  fManuList12(),
  fManuBridge2()
{  
/// Standard constructor

  AliDebug(1,"");
  fDDLs.SetOwner(true);
  fBusPatches.SetOwner(true);
  fBusPatches.SetSize(900);

  fTriggerCrates.SetOwner(true);
  fTriggerCrates.SetSize(16);

  fLocalBoards.SetOwner(true);
  fLocalBoards.SetSize(242); // included non-identied board


  // Load segmentation & DE store data
  if ( ! AliMpSegmentation::Instance(false) )
    AliMpSegmentation::ReadData(true);

  // Create all detection elements
  ReadDDLs();
  ReadTriggerDDLs();
  SetManus();
  SetPatchModules();
  SetBusPatchLength();
}

//______________________________________________________________________________
AliMpDDLStore::AliMpDDLStore(TRootIOCtor* /*ioCtor*/)
: TObject(),
  fDDLs(),
  fBusPatches(),
  fTriggerCrates(true),
  fLocalBoards(true)
{  
/// Constructor for IO

  AliDebug(1,"");

  fgInstance = this;
}


//______________________________________________________________________________
AliMpDDLStore::~AliMpDDLStore()
{
/// Destructor

  AliDebug(1,"");

  // DDL objects are deleted with fDDLs 
  // Bus patches objects are deleted with fBusPatches 
  
  fgInstance = 0;
}

//
// private methods
//

//______________________________________________________________________________
Int_t  AliMpDDLStore::GetManuListIndex(Int_t detElemId) const
{
/// Return the index of the manu list for given detElemId

  return AliMpDEManager::GetChamberId(detElemId)*4 + (detElemId % 100);
}  
  

//______________________________________________________________________________
Int_t AliMpDDLStore::GetBusPatchIndex(Int_t detElemId, Int_t manuId) const
{
/// Calculate the index of the buspatch from manuId

    Int_t pos = 0;
    AliMp::StationType stationType = AliMpDEManager::GetStationType(detElemId);
    static Int_t manuMask = AliMpConstants::ManuMask(AliMp::kNonBendingPlane) - 1;

    if( stationType == AliMp::kStation345) {
      pos = (manuId & manuMask)/100;
    } else {
      Int_t idx = GetManuListIndex(detElemId);

      // using array defined from DetElemIdToBusPatch.dat file
      for (pos = fManuList12[idx].GetSize()-1; pos >= 0; --pos)
	  if ( manuId >= fManuList12[idx].At(pos)) break;
    }

    return pos;
}


//______________________________________________________________________________
Bool_t AliMpDDLStore::ReadDDLs()
{
/// Read ddl <-> bus patch file
 
    
  TString infile = AliMpFiles::BusPatchFilePath();

  ifstream in(infile, ios::in);
  if (!in) {
    AliErrorStream() << "Data file " << infile << " not found.";
    return false;
  }  
      
  char line[255];

  while ( in.getline(line,255) ) {

     if ( line[0] == '#' ) continue;

     TString tmp(AliMpHelper::Normalize(line));

     TObjArray* stringList = tmp.Tokenize(TString(" "));

     TString sDE = ((TObjString*)stringList->At(0))->GetString();
     Int_t idDE  = atoi(sDE.Data());

     if ( ! AliMpDEManager::IsValidDetElemId(idDE, false) ) {
       AliErrorStream() << "DetElemId "<< idDE << " not valid." << endl;
       return false;
     }  

     TString busPatch = ((TObjString*)stringList->At(1))->GetString();


     TString sDDL = ((TObjString*)stringList->At(2))->GetString();
     Int_t  iDDL  = atoi(sDDL.Data());

     if ( iDDL < 0 || iDDL >= fgkNofDDLs ) {
       AliErrorStream() << "DDL id "<< iDDL << " outside limits." << endl;
       return false;
     }
  
     AliDebugStream(3)
	 << "idDE " << idDE << " buspatch " << busPatch.Data() << " iDDL " << iDDL 
	 << endl;

     // reading 1st manu Id for each bus patch (station 1 & 2)
     if(AliMpDEManager::GetStationType(idDE) != AliMp::kStation345) {
     
       TString sManu = ((TObjString*)stringList->At(3))->GetString();
       AliMpHelper::DecodeName(sManu,',',fManuList12[GetManuListIndex(idDE)]);
  
       if(AliMpDEManager::GetStationType(idDE) == AliMp::kStation2) {
	 TString sManuBridge = ((TObjString*)stringList->At(4))->GetString();
	 AliMpHelper::DecodeName(sManuBridge,',',fManuBridge2[GetManuListIndex(idDE)]);
       }

     }

     delete stringList;

     AliMpDDL* ddl = GetDDL(iDDL, false);
     if ( !ddl) {
       ddl = new AliMpDDL(iDDL);
       fDDLs.AddAt(ddl, iDDL);
     }  
     ddl->AddDE(idDE);
     
     TArrayI busPatchList;
     // decoding range of buspatch
     AliMpHelper::DecodeName(busPatch,';',busPatchList);
     
     // Update DE
     AliMpDetElement* de = AliMpDEManager::GetDetElement(idDE);
     de->SetDdlId(iDDL);
     // filling buspatch -> idDE
     for (Int_t i = 0; i < busPatchList.GetSize(); i++) {
       fBusPatches.Add(busPatchList[i], 
                       new AliMpBusPatch(busPatchList[i], idDE, iDDL));
       de->AddBusPatch(busPatchList[i]);
     }  
   }
   
   // Fill bus patch Ids array in DDLs now
   for ( Int_t i=0; i<fDDLs.GetEntriesFast(); i++ ) {
     AliMpDDL* ddl = (AliMpDDL*) fDDLs.At(i);
     ddl->FillBusPatchIds();
   }  
   
   in.close();
   return true;
}

//______________________________________________________________________________
Bool_t  AliMpDDLStore::ReadTriggerDDLs()
{
/// create trigger DDL object ddl<->Crate<->local board
  
   AliMpStringObjMap inputXfromMap;
   AliMpStringObjMap inputXtoMap;
   AliMpStringObjMap inputYfromMap;
   AliMpStringObjMap inputYtoMap;

    Int_t nonNotified = 0;

    Int_t iDDL = -1;

    TString infile = AliMpFiles::LocalTriggerBoardMapping();
  
    ifstream in(infile, ios::in);
    
    if (!in) {
      AliError(Form("Local Trigger Board Mapping File %s not found", infile.Data()));
      return kFALSE;
    }
  
    AliMpLocalBoard* board = 0x0;
    AliMpTriggerCrate* crate = 0x0;

    TString localNumber;
    Int_t localNumberId = 0;
    TArrayI list;

    char line[255];
    while (in.getline(line, 255))
    {

      // files contains lines with some blank caracters ???
      if (line[0] == ' ' && strlen(line) < 12) continue;

      TString tmp(AliMpHelper::Normalize(line));
      if (tmp.IsNull()) continue; // Ignore empty lines
 
      if (tmp.Contains("Board")) 
      {
	// extract id, name, slot & crate name
	TObjArray* stringList = tmp.Tokenize(TString(" "));
  
	TString localNumber = ((TObjString*)stringList->At(1))->GetString();

	if (localNumber.CompareTo("nn") != 0) // Notified cards
	    localNumberId = localNumber.Atoi();
	else 
	    localNumberId = AliMpConstants::NofLocalBoards()  + (++nonNotified);

	TString localName = ((TObjString*)stringList->At(4))->GetString();
	TString slotName  = ((TObjString*)stringList->At(10))->GetString();

	AliDebug(3, Form("%d %s %d\n", localNumberId, localName.Data(), slotName.Atoi()));

	board = new AliMpLocalBoard(localNumberId, localName.Data(), slotName.Atoi());

	// not notified board
	if (localNumber.CompareTo("nn") == 0)
	    board->SetNotified(false);

	// compose name with number and side (Right or Left)
	TString crateNumber = ((TObjString*)stringList->At(6))->GetString();
	TString crateSide = ((TObjString*)stringList->At(7))->GetString();
	TString crateName = crateNumber + crateSide;

	// set crate name
	board->SetCrate(crateName);

	// determine ddl number vs crate side
	if (crateName.Contains("R")) 
	    iDDL = fgkNofDDLs; // starts where tracker ends
	else 
	    iDDL = fgkNofDDLs + 1;
	
	AliMpDDL* ddl = GetDDL(iDDL, false);
	if ( !ddl) {
	  ddl = new AliMpDDL(iDDL);
	  fDDLs.AddAt(ddl, iDDL);
	}  
	
	Int_t crateId; 
	if (crateNumber.CompareTo("2-3") != 0)
	    crateId = crateNumber.Atoi();
	else
	    crateId = 23;

	// add trigger crate number for given ddl if not present
	if ( !ddl->HasTriggerCrateId(crateId) )
	    ddl->AddTriggerCrate(crateId);

	// if crate not existing then creating
	if (!GetTriggerCrate(crateName, false)) {
	    crate = new AliMpTriggerCrate(crateName.Data(), iDDL);
	    fTriggerCrates.Add(crateName.Data(), crate);
	}

	// create list of local board in this crate
	crate->AddLocalBoard(localNumberId);

	delete stringList;
      }

      if (tmp.Contains("DetElemId")) 
      {
	list.Reset();
	AliMpDDL* ddl = GetDDL(iDDL, true);

	// add  DE for local board and DDL
	TString detElem = &tmp[tmp.First(" ")+1];
	AliMpHelper::DecodeName(detElem,' ',list);

	for (Int_t i = 0; i < list.GetSize(); ++i) {

	  if (list[i]) { // skip copy card
	      AliMpDetElement* de = AliMpDEManager::GetDetElement(list[i]);
	      if (de->GetDdlId() == -1)
		  de->SetDdlId(iDDL);
	  }
	    board->AddDE(list[i]);
	    if(!ddl->HasDEId(list[i]))
		ddl->AddDE(list[i]);
	}
      }

      // map the different copies between cards in X and Y inputs
      // when copying Xinput, Yinput are copied as well
      if (tmp.Contains("X3input1")) 
      {
	TObjArray* stringList = tmp.Tokenize(TString(" "));
  
	TString copyX = ((TObjString*)stringList->At(3))->GetString();
	TString input = ((TObjString*)stringList->At(2))->GetString();;

	if (copyX.Contains("(1)")) { 
	  inputXfromMap.Add(input, board);
	  inputYfromMap.Add(input, board);

	}
	if (copyX.Contains("(2)")) {
	  inputXtoMap.Add(input, board);
	  inputYtoMap.Add(input, board);
	}
	delete stringList;
      }

      if (tmp.Contains("Y1input1")) 
      {
	TObjArray* stringList = tmp.Tokenize(TString(" "));

	TString copyY = ((TObjString*)stringList->At(3))->GetString();
	TString input = ((TObjString*)stringList->At(2))->GetString();;

	if (copyY.Contains("(1)")) 
	    inputYfromMap.Add(input, board);
	if (copyY.Contains("(2)")) 
	    inputYtoMap.Add(input, board);
	delete stringList;
      }


      if (tmp.Contains("transv")) 
      {
	// set transverse connector
	Int_t blankPos = tmp.Last(' ');

	TString transv(tmp(blankPos+1, tmp.Length()-blankPos));
	transv.ToUpper();
	if (transv.CompareTo("NONE") == 0)
	    board->SetTC(false);
      }

      if (tmp.Contains("Switch")) 
      {
	list.Reset();
	in.getline(line, 255);// skip one line

	in.getline(line, 255);
	TString tmp(AliMpHelper::Normalize(line));

	// decode switches
	AliMpHelper::DecodeName(tmp,' ',list);

	// store switches
	AliMpArrayI switches;

	for (Int_t i = 0; i < list.GetSize(); ++i)
	    board->AddSwitch(list[i]);

	// add local board into map
	fLocalBoards.Add(board->GetId(), board);

      }
    }
    
    // set copy card number to where the X-Y inputs are copied and
    // from where the X-Y inputs come.
    // deleting the first item (TString) done by AliMpStringObjMap itself
    // keep AliMpLocalBoard object undelete

    TString value;

    for (inputXfromMap.First(); !inputXfromMap.IsDone(); inputXfromMap.Next()) {
      
      value = inputXfromMap.CurrentKey();
      AliMpLocalBoard* boardFrom = (AliMpLocalBoard*)inputXfromMap.CurrentItem();
      AliMpLocalBoard* boardTo   = (AliMpLocalBoard*)inputXtoMap.Get(value);
      boardFrom->SetInputXto(boardTo->GetId());
      boardTo->SetInputXfrom(boardFrom->GetId());
      AliDebug(3, Form("copy xInputs from local id %d_%s_%d to %d_%s_%d\n", 
		       boardTo->GetInputXfrom(), boardFrom->GetCrate().Data(), boardFrom->GetSlot(),
		       boardFrom->GetInputXto(), boardTo->GetCrate().Data(), boardTo->GetSlot()));
    }

    for (inputYfromMap.First(); !inputYfromMap.IsDone(); inputYfromMap.Next()) {

      value = inputYfromMap.CurrentKey();
      AliMpLocalBoard* boardFrom = (AliMpLocalBoard*)inputYfromMap.CurrentItem();
      AliMpLocalBoard* boardTo   = (AliMpLocalBoard*)inputYtoMap.Get(value);
      boardFrom->SetInputYto(boardTo->GetId());
      boardTo->SetInputYfrom(boardFrom->GetId());
      AliDebug(3, Form("copy yInputs from local id %d_%s_%d to %d_%s_%d\n", 
		       boardTo->GetInputYfrom(), boardFrom->GetCrate().Data(), boardFrom->GetSlot(),
		       boardFrom->GetInputYto(), boardTo->GetCrate().Data(), boardTo->GetSlot()));
    }

    return kTRUE;
}

//______________________________________________________________________________
Bool_t AliMpDDLStore::SetManus()
{
/// Set manus for each bus patch

    Int_t manuMask = AliMpConstants::ManuMask(AliMp::kNonBendingPlane) - 1;

    // loop over DDL 
    for (Int_t iDDL = 0; iDDL < fgkNofDDLs; ++iDDL) {

      AliDebug(3, Form("DDL # %d\n", iDDL));

      AliMpDDL* ddl = GetDDL(iDDL);

      // loop over DE in the given DDL
      for (Int_t detElemIdx = 0; detElemIdx < ddl->GetNofDEs(); ++detElemIdx) {

	Int_t detElemId = ddl->GetDEId(detElemIdx);

	AliMpDetElement* detElement = GetDetElement(detElemId);

	AliMp::StationType stationType = AliMpDEManager::GetStationType(detElemId);


	// list of manu per DE on both cathode
	TList manuList;
	for ( Int_t cath = 0; cath < 2 ; ++cath ) {
	  const AliMpVSegmentation* seg 
	      = AliMpSegmentation::Instance()->GetMpSegmentation(detElemId,AliMp::GetCathodType(cath));
        
	  AliMp::PlaneType planeType = detElement->GetPlaneType(AliMp::GetCathodType(cath));

	  TArrayI manus;
	  seg->GetAllElectronicCardIDs(manus);
          
	  // filling TList manu
	  for ( Int_t im = 0; im < manus.GetSize(); ++im ) {
      
	    AliMpIntPair* manu = 0x0;
	    if( stationType == AliMp::kStation345) 
        manu = new AliMpIntPair((manus[im] & manuMask), planeType, kTRUE); //remove offset for NB
	    else 
        manu = new AliMpIntPair(manus[im], planeType, kTRUE); //keep offset for NB
	    
	    manuList.Add(manu);
      
      detElement->AddManu(manus[im]);
	  }        
	}// cathode
  
	manuList.Sort(); // sort regardless B or NB plane

	// filling manu to the corresponding buspatch
	for (Int_t iEntry = 0; iEntry < manuList.GetEntries(); ++iEntry) {

	  AliMpIntPair* manuPtr = (AliMpIntPair*)manuList.At(iEntry);

	  Int_t manuId = manuPtr->GetFirst();
	  Int_t pos    = GetBusPatchIndex(detElemId, manuId);

	  if (pos > detElement->GetNofBusPatches()) {
	    AliError(Form("pos greater %d than size %d manuId %d detElemId %d \n", 
			  pos, detElement->GetNofBusPatches(), manuId, detElemId));
	    return false;
	  }
	  
	  // get buspatch and fill manus
	  Int_t busPatchId = detElement->GetBusPatchId(pos);
	  AliMpBusPatch* busPatch = GetBusPatch(busPatchId);

	  if( stationType == AliMp::kStation345) {
      
	    if (manuPtr->GetSecond()) 
        busPatch->AddManu(manuId+manuMask+1); // add offset again after sorted
	    else 
        busPatch->AddManu(manuId); 
      
	  } else {
      
	    busPatch->AddManu(manuId); 
      
	  }
	}

	manuList.Delete();

	if (AliDebugLevel() == 3) {

	  // print out for checking
	  for(Int_t pos = 0; pos < detElement->GetNofBusPatches(); ++pos) {
	    Int_t busPatchId = detElement->GetBusPatchId(pos);
	    AliMpBusPatch* busPatch = GetBusPatch(busPatchId);
	    printf("BusPatch: %d\n", busPatch->GetId());
	    for (Int_t iEntry = 0; iEntry < busPatch->GetNofManus(); ++iEntry) 
		printf("manu Id: %d\n", busPatch->GetManuId(iEntry));
	  }
	}

      } // detection element loop
    }// DDL loop

    return true;
}

//______________________________________________________________________________
Bool_t AliMpDDLStore::SetPatchModules()
{
/// Compute the number of manu per PCB for each buspatch 

    AliMpDEIterator it;
    Bool_t result = true;

    for ( it.First(); !it.IsDone(); it.Next() ) {

      AliMpDetElement* detElement = it.CurrentDE();
      
      for (Int_t i = 0; i < detElement->GetNofBusPatches(); ++i) {
	AliMpBusPatch* busPatch = GetBusPatch(detElement->GetBusPatchId(i));
	Bool_t newResult = false;
	Int_t idDE = busPatch->GetDEId();

	if (AliMpDEManager::GetStationType(idDE) == AliMp::kStation2 ) 
	    newResult = busPatch->SetNofManusPerModule(fManuBridge2[GetManuListIndex(idDE)].At(i));
	else 
	    newResult = busPatch->SetNofManusPerModule();
      }
    }

    return result;
}

//______________________________________________________________________________
Bool_t AliMpDDLStore::SetBusPatchLength()
{
/// read the buspatch length file and set buspatch length 

    TString infile = AliMpFiles::BusPatchLengthFilePath();
    ifstream in(infile, ios::in);
    if (!in) {
      AliErrorStream() << "Data file " << infile << " not found.";
      return false;
    }  
    char line[255];

    for (Int_t iDDL = 0; iDDL < fgkNofDDLs; ++iDDL ) {
      AliMpDDL* ddl = GetDDL(iDDL);

      for (Int_t iBusPatch = 0; iBusPatch < ddl->GetNofBusPatches(); ++iBusPatch) {
  
	do {
	  if (!in.getline(line,255)) {
	    AliWarning(Form("Wrong size in bus patch length file; index %d DDL %d", 
			    iBusPatch, iDDL));
	    return false;
	  }
	} while(line[0] == '#');

	TString tmp(AliMpHelper::Normalize(line));

	TObjArray* stringList = tmp.Tokenize(TString(" "));
	
	TString sLocalBusId = ((TObjString*)stringList->At(0))->GetString();
	Int_t   localBusId  = sLocalBusId.Atoi();

	TString sLength = ((TObjString*)stringList->At(1))->GetString();
	Float_t length  = sLength.Atof();

	delete stringList;

	if (localBusId != iBusPatch + 1)
	    AliWarning(Form("Wrong local buspatch id %d instead of %d", iBusPatch+1, localBusId));

	Int_t busPatchId = ddl->GetBusPatchId(iBusPatch);
	AliMpBusPatch* busPatch = GetBusPatch(busPatchId);
	busPatch->SetCableLength(length);
      }
    }

    return true;
}


//________________________________________________________________
Int_t AliMpDDLStore::GetLocalBoardId(TString name) const
{
/// return the first board with a given side and line

   TExMapIter i = fLocalBoards.GetIterator();
    Long_t key, value;
    while ( i.Next(key, value) ) {
      AliMpLocalBoard* local = (AliMpLocalBoard*)value;

      TString tmp(&local->GetName()[4], 2);
      if (name.Contains(tmp))
	  if (name[0] == local->GetName()[0])
	      return local->GetId();
    }

    return 0;

}

//
// public methods
//


//______________________________________________________________________________
AliMpDDL* AliMpDDLStore::GetDDL(Int_t ddlId, Bool_t warn) const
{
/// Return DDL for given ddlId

  AliMpDDL* ddl 
    = (AliMpDDL*)fDDLs.At(ddlId);
    
  if ( ! ddl && warn ) {  
    AliErrorStream() 
        << "DDL with Id = " << ddlId << " not defined." << endl;
  }	

  return ddl;
}    

//______________________________________________________________________________
AliMpDetElement*  AliMpDDLStore::GetDetElement(Int_t detElemId, Bool_t warn) const
{
/// Return detection element with given detElemId

  if ( ! AliMpDEStore::Instance() ) {
    AliFatal("DE Store has not been loaded.");
    return 0; 
  }  

  return AliMpDEStore::Instance()->GetDetElement(detElemId, warn);
}  

//______________________________________________________________________________
AliMpBusPatch* AliMpDDLStore::GetBusPatch(Int_t busPatchId, Bool_t warn) const
{
/// Return bus patch with given Id

  AliMpBusPatch* busPatch
    = (AliMpBusPatch*) fBusPatches.GetValue(busPatchId);
    
  if ( ! busPatch && warn ) {  
    AliErrorStream() 
        << "Bus patch with Id = " << busPatchId << " not defined." << endl;
  }	

  return busPatch;
}    

//______________________________________________________________________________
AliMpLocalBoard* AliMpDDLStore::GetLocalBoard(Int_t localBoardId, Bool_t warn) const
{
/// Return bus patch with given Id

  AliMpLocalBoard* localBoard
    = (AliMpLocalBoard*) fLocalBoards.GetValue(localBoardId);
    
  if ( ! localBoard && warn ) {  
    AliErrorStream() 
        << "Local board with Id = " << localBoardId << " not defined." << endl;
  }	

  return localBoard;
}    

//______________________________________________________________________________
AliMpTriggerCrate* AliMpDDLStore::GetTriggerCrate(TString name, Bool_t warn) const
{
/// Return trigger crate with given name

  AliMpTriggerCrate* crate
     = (AliMpTriggerCrate*) fTriggerCrates.GetValue(name.Data());
    
  if ( ! crate && warn ) {  
    AliErrorStream() 
        << "Trigger crate with name = " << name.Data() << " not defined." << endl;
  }	

  return crate;
}  

//______________________________________________________________________________
AliMpTriggerCrate* AliMpDDLStore::GetTriggerCrate(Int_t ddlId, Int_t index, Bool_t warn) const
{
/// Return trigger crate with given ddl and index crate

  if (ddlId == 0 || ddlId == 1)
      ddlId += fgkNofDDLs;

  AliMpDDL* ddl = GetDDL(ddlId, warn);
  if ( ! ddl ) return 0; 
  
  if ( index >= ddl->GetNofTriggerCrates() ) {
      AliError(Form("crate id %d greater than array[%d]", index, ddl->GetNofTriggerCrates()));
      return 0;
  }    
   
  Int_t crateId = ddl->GetTriggerCrateId(index);
  TString name = AliMpTriggerCrate::GenerateName(crateId, ddlId, fgkNofDDLs);
  
  return GetTriggerCrate(name, warn);
}  

//______________________________________________________________________________
Int_t  AliMpDDLStore::GetDEfromBus(Int_t busPatchId) const
{
/// Return detection element Id for given busPatchId

  AliMpBusPatch* busPatch = GetBusPatch(busPatchId);
  
  if ( ! busPatch ) {
    AliErrorStream() 
        << "Bus patch with Id = " << busPatchId << " not defined." << endl;
    return 0;    
  }	

  return busPatch->GetDEId();
}  

//______________________________________________________________________________
Int_t  AliMpDDLStore::GetDEfromLocalBoard(Int_t localBoardId, Int_t chamberId) const
{
/// Return detElemId for local board Id and chamber id.

  AliMpLocalBoard* localBoard = GetLocalBoard(localBoardId);
  
  if ( ! localBoard ) {
    AliErrorStream() 
        << "Loacl board with Id = " << localBoardId << " not defined." << endl;
    return 0;    
  }	

   return localBoard->GetDEIdByChamber(chamberId);
}

//______________________________________________________________________________
Int_t  AliMpDDLStore::GetDDLfromBus(Int_t busPatchId) const
{
/// Return DDL Id for given busPatchId

  AliMpBusPatch* busPatch = GetBusPatch(busPatchId);
  
  if ( ! busPatch ) {
    AliErrorStream() 
        << "Bus patch with Id = " << busPatchId << " not defined." << endl;
    return 0;    
  }	

  return busPatch->GetDdlId();
}  

//______________________________________________________________________________
Int_t AliMpDDLStore::GetBusPatchId(Int_t detElemId, Int_t manuId) const
{
/// Return bus patch for a given manuId

  AliMpDetElement* detElement = GetDetElement(detElemId);
  Int_t pos = GetBusPatchIndex(detElemId, manuId);

  if ( pos > detElement->GetNofBusPatches() ) {
    AliErrorStream() 
      << "Pos = " << pos 
      << " greater than the size = " <<  detElement->GetNofBusPatches()
      << " for detElemId = " << detElemId 
      << " manuId = " << manuId << endl;
    return -1;
  }

  return detElement->GetBusPatchId(pos);
}    

//______________________________________________________________________________
AliMpIntPair  AliMpDDLStore::GetDetElemIdManu(Int_t manuSerial) const
{
/// Return the detElemId and manuId for given serial manu number

  if ( ! AliMpDEStore::Instance() ) {
    AliFatal("DE Store has not been loaded.");
    return AliMpIntPair::Invalid(); 
  }  

  return AliMpDEStore::Instance()->GetDetElemIdManu(manuSerial);
}  

//______________________________________________________________________________
void AliMpDDLStore::PrintAllManu() const
{
/// Print all manu Ids and their serial numbers sorted by detection element
/// and bus patch.                                                            \n
/// As serial manu numbers are filled in a different way than manu Ids this
/// printing allows to check that both ways are consistent

  // Loop over DE
  AliMpDEIterator it;
  for ( it.First(); ! it.IsDone(); it.Next() ) {
     AliMpDetElement* de = it.CurrentDE();
     cout << "DE: " << de->GetId() << endl; 
     
     // Loop over bus patches in this DE
     for ( Int_t i=0; i< de->GetNofBusPatches(); ++i ) {

       AliMpBusPatch* busPatch = GetBusPatch(de->GetBusPatchId(i));    
       cout << "  busPatch: " << busPatch->GetId() << endl; 

       cout << "    Manu       : ";
       for ( Int_t j=0; j<busPatch->GetNofManus(); ++j ) {
         cout << std::setw(6) << busPatch->GetManuId(j) << " ";
       }
       cout << endl;
       
       cout << "    Manu serial: ";
       for ( Int_t k=0; k<busPatch->GetNofManus(); ++k ) {
         cout << std::setw(6) << de->GetManuSerialFromId(busPatch->GetManuId(k)) << " ";
       }        
       cout << endl;
     }
  }
}  

//________________________________________________________________
Int_t  AliMpDDLStore::GetNextDEfromLocalBoard(Int_t localBoardId, Int_t chamberId ) const
{
/// return the next detection element in line

    AliMpLocalBoard* localBoard  =  GetLocalBoard(localBoardId);

    TString name(localBoard->GetName());

    Int_t line = localBoard->GetPosition().GetFirst();
    ++line;
    
    name.Replace(4,1,Form("%d", line));

   Int_t nextLocalId;
    if ((nextLocalId = GetLocalBoardId(name)))
	return GetDEfromLocalBoard(nextLocalId, chamberId);
    else
	return 0;

    return 0;
}

//________________________________________________________________
Int_t  AliMpDDLStore::GetPreviousDEfromLocalBoard(Int_t localBoardId, Int_t chamberId) const
{
/// return the previous detection element in line

    AliMpLocalBoard* localBoard  =  GetLocalBoard(localBoardId);

    TString name(localBoard->GetName());

    Int_t line = localBoard->GetPosition().GetFirst();
    --line;
    
    name.Replace(4,1,Form("%d", line));

    Int_t prevLocalId;
    if ((prevLocalId = GetLocalBoardId(name)))
	return GetDEfromLocalBoard(prevLocalId, chamberId);
    else
	return 0;

}

