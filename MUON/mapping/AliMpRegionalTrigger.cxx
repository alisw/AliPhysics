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
/// and return its instance
    

    AliMpLocalBoard* board = 0x0;
    AliMpTriggerCrate* crate = 0x0;


    Int_t localBoardId = 0;
    TArrayI list;
    UShort_t crateId, mask;
    
    char line[80];
   
    while (!in.eof())
    {
      in.getline(line,80);
      if (!strlen(line)) break;
      TString crateName(AliMpHelper::Normalize(line));
      
      in.getline(line,80);    
      sscanf(line,"%hx",&crateId);
  
      // skip mode
      in.getline(line,80);
      
      // skip coincidence
      in.getline(line,80);
      
      // skip mask
      in.getline(line,80);
      sscanf(line,"%hx",&mask);
      
      crate = (AliMpTriggerCrate*)(fTriggerCrates.GetValue(crateName.Data()));
      if (!crate) 
      {
        // cout << "Creating crate: " << crateName.Data() << endl;
        crate = new AliMpTriggerCrate(crateName.Data(), crateId);
        fTriggerCrates.Add(crateName.Data(), crate);
      }
      
      Char_t localBoardName[20];
      Int_t slot;
      UInt_t switches;
      
      for ( Int_t i = 0; i < AliMpConstants::LocalBoardNofChannels(); ++i ) 
      {
        if ( (mask >> i ) & 0x1 )
        {
          in.getline(line,80);
          sscanf(line,"%02d %s %03d %03x",&slot,localBoardName,&localBoardId,&switches);
          // cout << "  Creating local board: " << localBoardId << endl;
          board = new AliMpLocalBoard(localBoardId, localBoardName, slot); 
          board->SetSwitch(switches);
          board->SetCrate(crateName);
          
          if (localBoardId > AliMpConstants::NofLocalBoards())
            board->SetNotified(false); // copy cards
          
          crate->AddLocalBoard(localBoardId);

          // add  list of DEs for local board
          list.Reset();
          in.getline(line,80);
          TString tmp(AliMpHelper::Normalize(line));
          AliMpHelper::DecodeName(tmp,' ',list);
          for (Int_t ii = 0; ii < list.GetSize(); ++ii) { 
            if ( list[ii] ) board->AddDE(list[ii]);
          }  

          // set copy number and transverse connector
          in.getline(line,80);
          tmp = AliMpHelper::Normalize(line);
          AliMpHelper::DecodeName(tmp,' ',list);
    
          board->SetInputXfrom(list[0]);
          board->SetInputXto(list[1]);
          
          board->SetInputYfrom(list[2]);
          board->SetInputYto(list[3]);
          
          board->SetTC(list[4]);
          
          // add local board into array
          fLocalBoardArray.AddAt(board,board->GetId());
          fLocalBoardMap.Add(board->GetId(),board);
        }
      }
    }
    
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
    
    if ( fileName != "" ) {
      AliDebugStream(2) << "Read data from file " << fileName.Data() << endl;
    
      TString inFileName(fileName);
      inFileName = gSystem->ExpandPathName(inFileName.Data());
      ifstream inFile(inFileName.Data(), ios::in);
      if ( ! inFile.good() ) {
        AliErrorStream()
           << "Local Trigger Board Mapping File " << fileName.Data() << " not found" << endl;
        return kFALSE;
      }
      
      return ReadData(inFile);  
    } 
    else {
      AliDebugStream(2) << "Read data from stream " << fileName.Data() << endl;
      istream& in
         = AliMpDataStreams::Instance()
             ->CreateDataStream(AliMpFiles::LocalTriggerBoardMapping());
             
      Bool_t result = ReadData(in);
      
      delete &in;
      return result;        
    }
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
    
    return fLocalBoardArray.GetSize(); 
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


