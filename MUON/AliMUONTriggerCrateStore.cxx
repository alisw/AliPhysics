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

#include "AliMUONTriggerCrateStore.h"

#include "AliMUONTriggerCrate.h"
#include "AliMUONLocalTriggerBoard.h"
#include "AliMUONRegionalTriggerBoard.h"
#include "AliMpExMap.h"
#include "TString.h"
#include "TSystem.h"
#include "AliLog.h"
#include "Riostream.h"

/// 
/// \class AliMUONTriggerCrateStore
/// 
/// A container of trigger crate objects that offers iteration
/// over both the crates themselves and the local boards they contain
///
/// \author Laurent Aphecetche

/// \cond CLASSIMP
ClassImp(AliMUONTriggerCrateStore)
/// \endcond

//_____________________________________________________________________________
AliMUONTriggerCrateStore::AliMUONTriggerCrateStore()
: TObject(),
fCrates(0x0),
fLocalBoards(0x0),
fCrateIterator(0x0),
fLBIterator(0x0),
fCurrentCrate(0x0),
fCurrentLocalBoard(-1)
{
/// Default constructor
}

//_____________________________________________________________________________
AliMUONTriggerCrateStore::~AliMUONTriggerCrateStore()
{
/// Destructor
  delete fCrateIterator;
  delete fLBIterator;
  delete fCrates;
  delete fLocalBoards;
}

//_____________________________________________________________________________
void 
AliMUONTriggerCrateStore::AddCrate(const char *name)
{
  /// create and add a crate to our map
  if (!fCrates)
  {
    AliError("Object not properly initialized");
    return;
  }

  AliDebug(1,Form("Adding crate %s",name));
  TObject* there = fCrates->GetValue(name);
  if (there)
  {
    AliError(Form("Cannot add crate %s because it's already there !",name));
  }
  else
  {
    fCrates->Add(name,new AliMUONTriggerCrate(name,17));
  }
}

//_____________________________________________________________________________
AliMUONLocalTriggerBoard* 
AliMUONTriggerCrateStore::LocalBoard(Int_t boardNumber) const
{
  /// return a board by number
  
  if ( !fLocalBoards )
  {
    AliError("Object not properly initialized");
    return 0x0;
  }

  return static_cast<AliMUONLocalTriggerBoard*>(fLocalBoards->GetValue(boardNumber));
}

//_____________________________________________________________________________
AliMUONTriggerCrate* 
AliMUONTriggerCrateStore::Crate(const char *name) const
{
  /// return a crate by name
  if ( !fCrates )
  {
    AliError("Object not properly initialized");
    return 0x0;
  }
  return static_cast<AliMUONTriggerCrate*>(fCrates->GetValue(name));
}

//_____________________________________________________________________________
AliMUONTriggerCrate* 
AliMUONTriggerCrateStore::Crate(Int_t ddl, Int_t reg) const
{
  /// return a crate by name
  if ( !fCrates )
  {
    AliError("Object not properly initialized");
    return 0x0;
  }
  TString name = GetCrateName(ddl, reg);
  return static_cast<AliMUONTriggerCrate*>(fCrates->GetValue(name.Data()));
}
// //____________________________________________________________________
TString AliMUONTriggerCrateStore::GetCrateName(Int_t ddl, Int_t reg) const
{
  // set crate name from DDL & reg number

  Char_t name[10];
  switch(reg) {
      case 0:
      case 1:
	sprintf(name,"%d", reg+1);
	break;
      case 2:
	strcpy(name, "2-3");
	break;
      case 3:
      case 4:
      case 5:
      case 6:
      case 7:
	sprintf(name,"%d", reg);
	break;
  }

  // crate Right for first DDL
  if (ddl == 0)
    strcat(name, "R");
  else 
    strcat(name, "L"); 

  return TString(name);
}
//_____________________________________________________________________________
void
AliMUONTriggerCrateStore::FirstCrate()
{
  /// initialize iteration
  if ( !fCrates )
  {
    AliError("Object not properly initialized");
    return;
  }
  if (!fCrateIterator)
  {
    fCrateIterator = new TExMapIter(fCrates->GetIterator());
  }
  fCrateIterator->Reset();
}

//_____________________________________________________________________________
void
AliMUONTriggerCrateStore::FirstLocalBoard()
{
  /// Initialize iterator on local boards.
  /// Please note that we're not using directly the FirstCrate() and
  /// NextCrate() methods here to avoid mix and match between crate iterator
  /// and local board iterator
  fCurrentCrate = 0x0;
  fCurrentLocalBoard = 0;

  if ( !fLBIterator ) 
  {
    fLBIterator = new TExMapIter(fCrates->GetIterator());
  }
  fLBIterator->Reset();
  Long_t key, value;
  Bool_t ok = fLBIterator->Next(key,value);
  if ( ok )
  {
    fCurrentCrate = reinterpret_cast<AliMUONTriggerCrate*>(value);
    fCurrentLocalBoard = 1;
  }
}

//_____________________________________________________________________________
AliMUONTriggerCrate*
AliMUONTriggerCrateStore::NextCrate()
{
  /// Return the next crate in iteration, or 0 if iteration is ended.
  if (!fCrateIterator) return 0x0;
  
  Long_t key, value;
  Bool_t ok = fCrateIterator->Next(key,value);
  if (ok)
  {
    return reinterpret_cast<AliMUONTriggerCrate*>(value);
  }
  else
  {
    return 0x0;
  }
}

//_____________________________________________________________________________
AliMUONLocalTriggerBoard*
AliMUONTriggerCrateStore::NextLocalBoard()
{  
  /// Return the next local board in iteration, or 0 if iteration is ended.
  if ( !fLBIterator ) return 0x0;

  if ( fCurrentLocalBoard >= fCurrentCrate->Boards()->GetLast() +1)
//  if ( fCurrentLocalBoard >= fCurrentCrate->Boards()->GetLast() )
  {
    // try to go to next crate, if some are left
    Long_t key, value;
    Bool_t ok = fLBIterator->Next(key,value);
    if ( ok )
    {
      fCurrentCrate = reinterpret_cast<AliMUONTriggerCrate*>(value);
      fCurrentLocalBoard = 1;
    }
    else
    {
      fCurrentLocalBoard = 0;
      return 0x0;
    }
  }

  AliMUONLocalTriggerBoard* lb = static_cast<AliMUONLocalTriggerBoard*>
    (fCurrentCrate->Boards()->At(fCurrentLocalBoard));
  
  ++fCurrentLocalBoard;
  
  return lb;
}

//_____________________________________________________________________________
Int_t
AliMUONTriggerCrateStore::NumberOfCrates() const
{
  /// Number of crates we're holding
  if ( fCrates ) return fCrates->GetSize();
  return 0;
}

//_____________________________________________________________________________
Int_t
AliMUONTriggerCrateStore::NumberOfLocalBoards() const
{
  /// Number of local boards we're holding
  if ( fLocalBoards ) return fLocalBoards->GetSize();
  return 0;
}

//_____________________________________________________________________________
void
AliMUONTriggerCrateStore::ReadFromFile(const char* file)
{
    /// Read crate and local board information from file.
    fCrates = new AliMpExMap(kTRUE);
    fLocalBoards = new AliMpExMap(kFALSE);
  
    ifstream myInputFile(gSystem->ExpandPathName(file), ios::in);
  
    string sLine, sValue;
  
    if ( !myInputFile ) 
    {
      AliError("TRIGGER ELECTRONICS CONFIGURATION FILE COULD NOT BE OPENED");
    }
    else
    {
      while (getline(myInputFile,sLine))
      {
	if (sLine.empty()) continue; // Ignore empty lines
	else
	{
	  const Int_t kMaxfields = 15; char **fields = new char*[kMaxfields];
        
	  char s[100]; 
        
	  if (sLine.find("Board",0) != string::npos) 
	  {   
	    strcpy(s,sLine.c_str());
          
	    Int_t numlines = 0;
          
	    for (char *token = strtok(s, " ");
		 token != NULL;
		 token = strtok(NULL, " "))
	    {
	      fields[numlines] = new char[strlen(token)+1];
	      strcpy(fields[numlines++],token);
	    }
          
	    char str[10]; strcpy(str, fields[6]); strcat(str, fields[7]);
          
	    AliMUONTriggerCrate *crate = Crate(str); 
          
	    if (!crate) 
	    {
	      AddCrate(str); crate = Crate(str);
            
	      AliMUONRegionalTriggerBoard *rboard = new AliMUONRegionalTriggerBoard();
	      crate->AddBoard(rboard, 0);
	    }               
          
	    //             CONVENTION: SLOT 0 HOLDS THE REGIONAL BOARD
	    Int_t sl = atoi(fields[10]);
          
//          AliMUONTriggerLut* lut = calibData->TriggerLut();
          
	    AliMUONLocalTriggerBoard *board = 
		new AliMUONLocalTriggerBoard(fields[4], sl, 0x0);
          
	    if (strcmp(fields[1],"nn")) 
	    {
	      Int_t sboard = atoi(fields[1]);
            
	      board->SetNumber(sboard);
	      fLocalBoards->Add(sboard,board);
            
// 	      fCrateMap[sboard-1] = new char[strlen(str)+1]; strcpy(fCrateMap[sboard-1], str);
// 	      fBoardMap[sboard-1] = sl;
	    }
          
	    board->SetCrate(str);
          
	    crate->AddBoard(board, sl);
          
	    while (getline(myInputFile,sLine)) if (sLine.find("transv",0) != string::npos) break;
          
	    strcpy(s,sLine.c_str());
          
	    for (char *token = strtok(s, " ");
		 token != NULL;
		 token = strtok(NULL, " ")) if (!strcmp(token,"NONE")) board->SetTC(kFALSE);
          
	    while (getline(myInputFile,sLine)) if (sLine.find("Switch",0) != string::npos) break;
          
	    while (getline(myInputFile,sLine)) if (!sLine.empty()) break;   
          
	    strcpy(s,sLine.c_str());
          
	    Int_t lines = 0;
          
	    for (char *token = strtok(s, " ");
		 token != NULL;
		 token = strtok(NULL, " ")) board->SetSwitch(lines++, atoi(token));
          
	    for (Int_t i = 0; i<numlines; i++) 
		if (fields[i]) {delete [] fields[i]; fields[i] = 0;}
              
	    delete [] fields; fields = 0;
	  }
	}
      }
    }
}
