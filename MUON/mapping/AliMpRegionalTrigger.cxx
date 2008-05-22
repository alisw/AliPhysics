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
#include "AliMpTriggerCrate.h"
#include "AliMpLocalBoard.h"
#include "AliMpConstants.h"
#include "AliMpFiles.h"
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
    fTriggerCrates(true),
    fLocalBoards(true)
{
/// Standard constructor
  
    fTriggerCrates.SetOwner(true);
    fTriggerCrates.SetSize(AliMpConstants::LocalBoardNofChannels());

    fLocalBoards.SetOwner(true);
    fLocalBoards.SetSize(AliMpConstants::TotalNofLocalBoards()); // included non-notified boards
}

//______________________________________________________________________________
AliMpRegionalTrigger::AliMpRegionalTrigger(const AliMpRegionalTrigger& rhs)
  : TObject(rhs),
    fTriggerCrates(rhs.fTriggerCrates),
    fLocalBoards(rhs.fLocalBoards)
{
/// Copy constructor
}  

//______________________________________________________________________________
AliMpRegionalTrigger::AliMpRegionalTrigger(TRootIOCtor* /*ioCtor*/)
  : TObject(),
    fTriggerCrates(),
    fLocalBoards()
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
  fLocalBoards = rhs.fLocalBoards;
  
  return *this;
}  

//______________________________________________________________________________
AliMpRegionalTrigger::~AliMpRegionalTrigger()
{
/// Destructor
}

//
// public methods
//

//______________________________________________________________________________
Bool_t AliMpRegionalTrigger::ReadData(const TString& fileName)
{
/// Load the Regional trigger from ASCII data files
/// and return its instance
    
    TString inFileName(fileName);
    if ( inFileName == "" )
      inFileName = AliMpFiles::LocalTriggerBoardMapping();
    
    inFileName = gSystem->ExpandPathName(inFileName.Data());

    ifstream in(inFileName.Data(), ios::in);

    if (!in) {
      AliErrorStream()
         << "Local Trigger Board Mapping File " << fileName.Data() << " not found" << endl;
      return kFALSE;
    }

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
          
          // add local board into map
          fLocalBoards.Add(board->GetId(), board);
        }
      }
    }
    return kTRUE;
}

//______________________________________________________________________________
AliMpLocalBoard* AliMpRegionalTrigger::FindLocalBoard(Int_t localBoardId, 
                                                      Bool_t warn) const {
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
AliMpTriggerCrate* AliMpRegionalTrigger::GetTriggerCrate(Int_t index) const
{ 
    /// Return the trigger crates with given index;

    return static_cast<AliMpTriggerCrate*>(fTriggerCrates.GetObject(index)); 
}

//______________________________________________________________________________
AliMpTriggerCrate* AliMpRegionalTrigger::GetTriggerCrateFast(Int_t index) const
{ 
    /// Return the trigger crates with given index;
    /// the index is not checked as we use the fast method in AliMpExMap.

    return static_cast<AliMpTriggerCrate*>(fTriggerCrates.GetObjectFast(index)); 
}

//______________________________________________________________________________
Int_t AliMpRegionalTrigger::GetNofLocalBoards() const
{ 
    /// Return number of local boards
    
    return fLocalBoards.GetSize(); 
}

//______________________________________________________________________________
AliMpLocalBoard* AliMpRegionalTrigger::GetLocalBoard(Int_t index) const
{ 
    /// Return local board with given index;
    
    return static_cast<AliMpLocalBoard*>(fLocalBoards.GetObject(index)); 
}

//______________________________________________________________________________
AliMpLocalBoard* AliMpRegionalTrigger::GetLocalBoardFast(Int_t index) const
{ 
    /// Return local board with given index;
    /// the index is not checked as we use the fast method in AliMpExMap.
    
    return static_cast<AliMpLocalBoard*>(fLocalBoards.GetObjectFast(index)); 
}




