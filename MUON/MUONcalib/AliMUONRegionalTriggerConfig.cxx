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
#include "AliMpHelper.h"
#include "AliMpExMapIterator.h"
#include "AliMpRegionalTrigger.h"
#include "AliLog.h"

#include <TArrayI.h>
#include <Riostream.h>
#include <TClass.h>
#include <TSystem.h>
#include <TList.h>


using std::cout;
using std::endl;
using std::ifstream;
using std::ios;
/// \cond CLASSIMP
ClassImp(AliMUONRegionalTriggerConfig)
/// \endcond


//______________________________________________________________________________
AliMUONRegionalTriggerConfig::AliMUONRegionalTriggerConfig()
  : TObject(),
    fTriggerCrates()
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
    /// Load the Regional trigger from ASCII data file

    // Read first data contained in mapping object
    //
    AliMpRegionalTrigger mpRegionalTrigger;
    mpRegionalTrigger.SetTriggerCratesOwner(kFALSE); 
    if ( ! mpRegionalTrigger.ReadData(fileName) ) {
        AliErrorStream()
           << "Reading mapping regional trigger from file " << fileName.Data() << " failed." 
           << endl;
        return 0;
    }

    // Fill calibration object from mapping object
    //
    TIterator* it = mpRegionalTrigger.CreateCrateIterator();
    AliMpTriggerCrate* mpTriggerCrate;
    while ( ( mpTriggerCrate = (AliMpTriggerCrate*)it->Next() ) ) {
      fTriggerCrates.Add(
        mpTriggerCrate->GetName(), new AliMUONTriggerCrateConfig(mpTriggerCrate));
    }    
    delete it;     
        
    // 

    // Read remaining calibration data from file
    //
    ifstream in(gSystem->ExpandPathName(fileName.Data()), ios::in);
    if ( ! in.good() ) {
        AliErrorStream()
           << "Local Trigger Board Mapping File " << fileName.Data() << " not found" << endl;
        return 0;
    }

    UShort_t mask;
    Int_t mode, coincidence;
    Int_t nofBoards;
    char line[80];

    // decode file and store in objects
    while (!in.eof())
    {
      // Get name
      in.getline(line,80);
      if (!strlen(line)) break;
      TString crateName(AliMpHelper::Normalize(line));

      in.getline(line,80);    
  
      // read mode
      in.getline(line,80);
      sscanf(line,"%d",&mode);

      // read coincidence
      in.getline(line,80);
      sscanf(line,"%d",&coincidence);

      // read mask
      in.getline(line,80);
      sscanf(line,"%hx",&mask);
      
      // read # local board
      in.getline(line,80);
      sscanf(line,"%d",&nofBoards);

      AliMUONTriggerCrateConfig*  crateConfig 
        = (AliMUONTriggerCrateConfig*)(fTriggerCrates.GetValue(crateName.Data()));
        
      // This should never happen, but let's test it anyway  
      if ( ! crateConfig ) {
        AliErrorStream()
           << "Cannot find AliMUONTriggerCrateConfig " << crateName.Data() << endl;
        return 0;
      }
       
      crateConfig->SetMode(mode);
      crateConfig->SetCoinc(coincidence);
      crateConfig->SetMask(mask);
      
      // Skipp local board data
      for ( Int_t i = 0; i < 3*nofBoards; ++i ) 
          in.getline(line,80);
    }

    return fTriggerCrates.GetSize();
}

//______________________________________________________________________________
AliMUONTriggerCrateConfig* AliMUONRegionalTriggerConfig::FindTriggerCrate(TString name, 
                                                          Bool_t warn) const  
{
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
TIterator* 
AliMUONRegionalTriggerConfig::CreateCrateIterator() const 
{ 
  /// Return trigger crates iterator
  return fTriggerCrates.CreateIterator(); 
}




