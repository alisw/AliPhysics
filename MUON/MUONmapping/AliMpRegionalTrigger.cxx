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
// $MpId: AliMpTrigger.cxx,v 1.4 2006/05/24 13:58:52 ivana Exp $

//-----------------------------------------------------------------------------
// Class AliMpRegionalTrigger
// --------------------
// The class defines the properties of regional trigger crate
// Author: Ch. Finck, Subatech Nantes
//-----------------------------------------------------------------------------

#include "AliMpRegionalTrigger.h"
#include "AliMpExMapIterator.h"
#include "AliMpTriggerCrate.h"
#include "AliMpLocalBoard.h"
#include "AliMpConstants.h"
#include "AliMpFiles.h"
#include "AliMpDataStreams.h"
#include "AliMpHelper.h"

#include "AliLog.h"

#include <TArrayI.h>
#include <Riostream.h>
#include <TClass.h>
#include <TSystem.h>


/// \cond CLASSIMP
ClassImp(AliMpRegionalTrigger)
/// \endcond


//______________________________________________________________________________
AliMpRegionalTrigger::AliMpRegionalTrigger()
  : TObject(),
    fTriggerCrates(),
    fLocalBoardMap(),
    fLocalBoardArray(AliMpConstants::TotalNofLocalBoards()+1) // included non-notified boards
{
      /// Standard constructor
  
    fTriggerCrates.SetOwner(true);
    fTriggerCrates.SetSize(AliMpConstants::LocalBoardNofChannels());
}

//______________________________________________________________________________
AliMpRegionalTrigger::AliMpRegionalTrigger(const AliMpRegionalTrigger& rhs)
  : TObject(rhs),
    fTriggerCrates(rhs.fTriggerCrates),
    fLocalBoardMap(rhs.fLocalBoardMap),
    fLocalBoardArray(rhs.fLocalBoardArray)
{
/// Copy constructor
}  

//______________________________________________________________________________
AliMpRegionalTrigger::AliMpRegionalTrigger(TRootIOCtor* ioCtor)
  : TObject(),
    fTriggerCrates(ioCtor),
    fLocalBoardMap(ioCtor),
    fLocalBoardArray()
{
/// Constructor for I0
}

//______________________________________________________________________________
AliMpRegionalTrigger& AliMpRegionalTrigger::operator=(const AliMpRegionalTrigger& rhs)
{
/// Assignment operator

  // check assignment to self
  if (this == &rhs) return *this;

  // base class assignment
  TObject::operator=(rhs);

  // assignment operator
  fTriggerCrates = rhs.fTriggerCrates;
  fLocalBoardArray = rhs.fLocalBoardArray;
  
  return *this;
}  

//______________________________________________________________________________
AliMpRegionalTrigger::~AliMpRegionalTrigger()
{
/// Destructor
}


//
// private methods
//

//______________________________________________________________________________
Bool_t AliMpRegionalTrigger::ReadData(istream& in)
{
/// Load the Regional trigger from ASCII data files
/// and fill objects. Return false if reading fails
  
  if ( !in.good() ) return kFALSE;
   
  Int_t localBoardId = 0;
  TArrayI listInt;
  UShort_t crateId;
  Int_t nofBoards;
  Int_t localBoardIndex(0);
  char line[80];
 
  // decode file and store in objects
  while (!in.eof())
  {
    in.getline(line,80);
    if (!strlen(line)) break;
    TString crateName(AliMpHelper::Normalize(line));
    
    in.getline(line,80);    
    sscanf(line,"%hx",&crateId);

    // skip data which are not stored in mapping object
    // (mode, coincidence, mask)
    in.getline(line,80);
    in.getline(line,80);
    in.getline(line,80);
    
    // read # local board
    in.getline(line,80);
    sscanf(line,"%d",&nofBoards);
    
    AliMpTriggerCrate* crate 
      = (AliMpTriggerCrate*)(fTriggerCrates.GetValue(crateName.Data()));
    if (!crate)  {
      crate = new AliMpTriggerCrate(crateName.Data(), crateId);
      fTriggerCrates.Add(crateName.Data(), crate);
    }

    Char_t localBoardName[20];
    Int_t slot;
    UInt_t switches;
    
    for ( Int_t i = 0; i < nofBoards; ++i ) 
    {
        in.getline(line,80);
        sscanf(line,"%02d %19s %03d %03x",&slot,localBoardName,&localBoardId,&switches);
        AliMpLocalBoard* board = new AliMpLocalBoard(localBoardId, localBoardName, slot); 
        board->SetSwitch(switches);
        board->SetCrate(crateName);
        
        if (localBoardId > AliMpConstants::NofLocalBoards())
          board->SetNotified(false); // copy cards
        
        crate->AddLocalBoard(localBoardId);
        
        // add  list of DEs for local board
        listInt.Reset();
        in.getline(line,80);
        TString tmp(AliMpHelper::Normalize(line));
        AliMpHelper::DecodeName(tmp,' ',listInt);
        for (Int_t ii = 0; ii < listInt.GetSize(); ++ii) { 
          if ( listInt[ii] ) board->AddDE(listInt[ii]);
        }  
         
        // set copy number and transverse connector
        in.getline(line,80);
        TString tmp1 = AliMpHelper::Normalize(line);
        AliMpHelper::DecodeName(tmp1,' ',listInt);
        
        board->SetInputXfrom(listInt[0]);
        board->SetInputXto(listInt[1]);
        
        board->SetInputYfrom(listInt[2]);
        board->SetInputYto(listInt[3]);
        
        board->SetTC(listInt[4]);
        
        // add local board into array
        fLocalBoardArray.AddAt(board,localBoardIndex);
        fLocalBoardMap.Add(board->GetId(),board);
      
      ++localBoardIndex;
    }
  }

  AliDebug(1,Form("%d trigger crate created",fTriggerCrates.GetSize()));
  AliDebug(1,Form("%d local board added to the map",fLocalBoardMap.GetSize()));
  AliDebug(1,Form("%d local board referenced from the array",fLocalBoardArray.GetLast()+1));
  
  return kTRUE;
}  

//
// public methods
//

//______________________________________________________________________________
Bool_t AliMpRegionalTrigger::ReadData(const TString& fileName)
{
/// Load the Regional trigger from ASCII data files
/// and return its instance
    
    AliDebugStream(2) << "Read data from file " << fileName.Data() << endl;
    
    TString inFileName(gSystem->ExpandPathName(fileName.Data()));
    ifstream inFile(inFileName.Data(), ios::in);
    if ( ! inFile.good() ) {
      AliErrorStream()
         << "Local Trigger Board Mapping File " << inFileName.Data() << " not found bordel de merde" << endl;
      return kFALSE;
    }
    
    return ReadData(inFile);  
}

//______________________________________________________________________________
Bool_t AliMpRegionalTrigger::ReadData(const AliMpDataStreams& dataStreams)
{
/// Load the Regional trigger from ASCII data files
/// and return its instance
    
    AliDebugStream(2) << "Read data from stream " << endl;
    istream& in
       = dataStreams.
           CreateDataStream(AliMpFiles::LocalTriggerBoardMapping());
           
    Bool_t result = ReadData(in);
    
    delete &in;
    return result;        
}

//______________________________________________________________________________
AliMpLocalBoard* AliMpRegionalTrigger::FindLocalBoard(Int_t localBoardId, 
                                                      Bool_t warn) const {
    /// Return local board with given Id

    AliMpLocalBoard* localBoard
      = static_cast<AliMpLocalBoard*>(fLocalBoardMap.GetValue(localBoardId));
    
    if ( ! localBoard && warn ) {
        AliErrorStream()
        << "Loacl board with localBoardId = " << localBoardId << " not found." << endl;
    }      

    return localBoard;
}

//______________________________________________________________________________
AliMpTriggerCrate* AliMpRegionalTrigger::FindTriggerCrate(TString name, 
                                                          Bool_t warn) const  {
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
Int_t AliMpRegionalTrigger::GetNofTriggerCrates() const 
{ 
    /// Return number of trigger crates

    return fTriggerCrates.GetSize(); 
}

//______________________________________________________________________________
Int_t AliMpRegionalTrigger::GetNofLocalBoards() const
{ 
    /// Return number of local boards
    
    return fLocalBoardArray.GetLast()+1; 
}

//______________________________________________________________________________
TIterator* 
AliMpRegionalTrigger::CreateCrateIterator() const
{
  /// Create iterator over crates

  return fTriggerCrates.CreateIterator();
}

//______________________________________________________________________________
TIterator* 
AliMpRegionalTrigger::CreateLocalBoardIterator() const
{
  /// Create iterator over local boards

  return fLocalBoardArray.MakeIterator();
}

//______________________________________________________________________________
Int_t 
AliMpRegionalTrigger::LocalBoardId(Int_t index) const
{
  /// Return local board Id for the local boards with a given index

  AliMpLocalBoard* lb = static_cast<AliMpLocalBoard*>(fLocalBoardArray.At(index));
  if (lb)
  {
    return lb->GetId();
  }
  AliError(Form("Could not get local board at index %d",index));
  return -1;
}

//______________________________________________________________________________
void AliMpRegionalTrigger::SetTriggerCratesOwner(Bool_t owner)
{
  /// Set ownership to trigger crates

  fTriggerCrates.SetOwner(owner);
}  
