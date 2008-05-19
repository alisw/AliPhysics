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
// Class AliMUONRegionalTriggerConfig
// --------------------
// The class defines the configuration of regional trigger crate
// Author: Ch. Finck, Subatech Nantes
//-----------------------------------------------------------------------------

#include "AliMUONRegionalTriggerConfig.h"
#include "AliMUONTriggerCrateConfig.h"
#include "AliMpConstants.h"
#include "AliMpFiles.h"
#include "AliMpHelper.h"

#include "AliLog.h"

#include <TArrayI.h>
#include <Riostream.h>
#include <TClass.h>
#include <TSystem.h>


/// \cond CLASSIMP
ClassImp(AliMUONRegionalTriggerConfig)
/// \endcond


//______________________________________________________________________________
AliMUONRegionalTriggerConfig::AliMUONRegionalTriggerConfig()
  : TObject(),
    fTriggerCrates(true)
{
/// Standard constructor
  
    fTriggerCrates.SetOwner(true);
    fTriggerCrates.SetSize(AliMpConstants::LocalBoardNofChannels());
}

//______________________________________________________________________________
AliMUONRegionalTriggerConfig::AliMUONRegionalTriggerConfig(const AliMUONRegionalTriggerConfig& rhs)
  : TObject(rhs),
    fTriggerCrates(rhs.fTriggerCrates)
{
/// Copy constructor
}  

//______________________________________________________________________________
AliMUONRegionalTriggerConfig& AliMUONRegionalTriggerConfig::operator=(const AliMUONRegionalTriggerConfig& rhs)
{
/// Assignment operator

  // check assignment to self
  if (this == &rhs) return *this;

  // base class assignment
  TObject::operator=(rhs);

  // assignment operator
  fTriggerCrates = rhs.fTriggerCrates;
  
  return *this;
}  

//______________________________________________________________________________
AliMUONRegionalTriggerConfig::~AliMUONRegionalTriggerConfig()
{
/// Destructor
}

//
// public methods
//

//______________________________________________________________________________
Int_t AliMUONRegionalTriggerConfig::ReadData(const TString& fileName)
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

    AliMUONTriggerCrateConfig* crate = 0x0;

    TArrayI list;
    UShort_t crateId, mask;
    Int_t mode, coincidence;
    Int_t localBoardId = 0;
    char line[80];
   
    while (!in.eof())
    {
      in.getline(line,80);
      if (!strlen(line)) break;
      TString crateName(AliMpHelper::Normalize(line));
      
      in.getline(line,80);    
      sscanf(line,"%hx",&crateId);
  
      in.getline(line,80);
      sscanf(line,"%d",&mode);
      
      in.getline(line,80);
      sscanf(line,"%d",&coincidence);
      
      in.getline(line,80);
      sscanf(line,"%hx",&mask);
      
      crate = (AliMUONTriggerCrateConfig*)(fTriggerCrates.GetValue(crateName.Data()));
      if (!crate) 
      {
        // cout << "Creating crate: " << crateName.Data() << endl;
        crate = new AliMUONTriggerCrateConfig(crateName.Data(), crateId, mask, mode, coincidence);
        fTriggerCrates.Add(crateName.Data(), crate);
      }
      
      Char_t localBoardName[20];
      Int_t slot;
      UInt_t switches;
      
      for ( Int_t i = 0; i < AliMpConstants::LocalBoardNofChannels(); ++i ) 
      {
        if ( (mask >> i ) & 0x1 )
        {
          // read local board
          in.getline(line,80);
          sscanf(line,"%02d %s %03d %03x",&slot,localBoardName,&localBoardId,&switches);
          crate->AddLocalBoard(localBoardId);

          // skip DEs for local board
          in.getline(line,80);

          // skip copy number and transverse connector
          in.getline(line,80);
         }
      }
    }
    return fTriggerCrates.GetSize();
}


//______________________________________________________________________________
AliMUONTriggerCrateConfig* AliMUONRegionalTriggerConfig::FindTriggerCrate(TString name, 
                                                          Bool_t warn) const  {
    /// Return trigger crate with given name

    AliMUONTriggerCrateConfig* crate
    = (AliMUONTriggerCrateConfig*) fTriggerCrates.GetValue(name.Data());

    if ( ! crate && warn ) {
        AliErrorStream()
        << "Trigger crate with name = " << name.Data() << " not defined." << endl;
    }

    return crate;
}

//______________________________________________________________________________
Int_t AliMUONRegionalTriggerConfig::GetNofTriggerCrates() const 
{ 
    /// Return number of trigger crates

    return fTriggerCrates.GetSize(); 
}

//______________________________________________________________________________
AliMUONTriggerCrateConfig* AliMUONRegionalTriggerConfig::GetTriggerCrate(Int_t index) const
{ 
    /// Return the trigger crates with given index;

    return static_cast<AliMUONTriggerCrateConfig*>(fTriggerCrates.GetObject(index)); 
}

//______________________________________________________________________________
AliMUONTriggerCrateConfig* AliMUONRegionalTriggerConfig::GetTriggerCrateFast(Int_t index) const
{ 
    /// Return the trigger crates with given index;
    /// the index is not checked as we use the fast method in AliMpExMap.

    return static_cast<AliMUONTriggerCrateConfig*>(fTriggerCrates.GetObjectFast(index)); 
}




