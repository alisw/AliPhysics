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
// Class AliMpRegionalTriggerReader
// --------------------
// The class to read the  regional trigger crate files
// Author: Ch. Finck, Subatech Nantes
//-----------------------------------------------------------------------------

#include "AliMpRegionalTriggerReader.h"
#include "AliMpTriggerCrate.h"
#include "AliMpLocalBoard.h"
#include "AliMUONRegionalTriggerConfig.h"
#include "AliMUONTriggerCrateConfig.h"
#include "AliMpConstants.h"
#include "AliMpFiles.h"
#include "AliMpHelper.h"
#include "AliMpExMap.h"

#include "AliLog.h"

#include <TArrayI.h>
#include <TObjArray.h>
#include <Riostream.h>
#include <TClass.h>
#include <TSystem.h>


/// \cond CLASSIMP
ClassImp(AliMpRegionalTriggerReader)
/// \endcond


//______________________________________________________________________________
AliMpRegionalTriggerReader::AliMpRegionalTriggerReader()
  : TObject()
{
/// constructor
}

//______________________________________________________________________________
AliMpRegionalTriggerReader::~AliMpRegionalTriggerReader()
{
/// Destructor
}

//
// public methods
//

//______________________________________________________________________________
Int_t AliMpRegionalTriggerReader::ReadData(TList& list, istream& in)
{
/// Load the Regional trigger from ASCII data files
/// and fill objects
    
    if (!in) {
      return -1;
    }
     
    AliMpExMap* triggerCrates = 0x0;
    AliMpExMap* localBoardMap = 0x0;
    TObjArray*  localBoardArray = 0x0;
    
    // common trigger crate map for mapping/config
    triggerCrates = static_cast<AliMpExMap*> (list.At(0));

    Bool_t mapping = false;
    // only for mapping
    if (list.GetSize() == 3) 
    {
      mapping = true;
      localBoardMap   = static_cast<AliMpExMap*> (list.At(1));
      localBoardArray = static_cast<TObjArray*>  (list.At(2));
    }
              
    AliMpLocalBoard* board = 0x0;
    AliMpTriggerCrate* crate = 0x0;
    AliMUONTriggerCrateConfig* crateConfig = 0x0;

    Int_t localBoardId = 0;
    TArrayI listInt;
    UShort_t crateId, mask;
    Int_t mode, coincidence;
    Int_t nofBoards;
    char line[80];
   
    // decode file and store in objects
    while (!in.eof())
    {
      in.getline(line,80);
      if (!strlen(line)) break;
      TString crateName(AliMpHelper::Normalize(line));
      
      in.getline(line,80);    
      sscanf(line,"%hx",&crateId);
  
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

      if (mapping)
      {
        crate = (AliMpTriggerCrate*)(triggerCrates->GetValue(crateName.Data()));
        if (!crate) 
        {
          crate = new AliMpTriggerCrate(crateName.Data(), crateId);
          triggerCrates->Add(crateName.Data(), crate);
        }
      }
      else
      {
        crateConfig = (AliMUONTriggerCrateConfig*)(triggerCrates->GetValue(crateName.Data()));
        if (!crateConfig) 
        {
          crateConfig = new AliMUONTriggerCrateConfig(crateName.Data(), crateId, mask, mode, coincidence);
          triggerCrates->Add(crateName.Data(), crateConfig);
        }
      }
      
      Char_t localBoardName[20];
      Int_t slot;
      UInt_t switches;
      
      for ( Int_t i = 0; i < nofBoards; ++i ) 
      {
          in.getline(line,80);
          sscanf(line,"%02d %s %03d %03x",&slot,localBoardName,&localBoardId,&switches);
          if (mapping)
          {
            board = new AliMpLocalBoard(localBoardId, localBoardName, slot); 
            board->SetSwitch(switches);
            board->SetCrate(crateName);
            
            if (localBoardId > AliMpConstants::NofLocalBoards())
              board->SetNotified(false); // copy cards
            
            crate->AddLocalBoard(localBoardId);
          }
          else
          {
            crateConfig->AddLocalBoard(localBoardId);
          }
          
          // add  list of DEs for local board
          listInt.Reset();
          in.getline(line,80);
          if (mapping)
          {
            TString tmp(AliMpHelper::Normalize(line));
            AliMpHelper::DecodeName(tmp,' ',listInt);
            for (Int_t ii = 0; ii < listInt.GetSize(); ++ii) { 
              if ( listInt[ii] ) board->AddDE(listInt[ii]);
            }  
          }
          
          // set copy number and transverse connector
          in.getline(line,80);
          if (mapping)
          {
            TString tmp1 = AliMpHelper::Normalize(line);
            AliMpHelper::DecodeName(tmp1,' ',listInt);
            
            board->SetInputXfrom(listInt[0]);
            board->SetInputXto(listInt[1]);
            
            board->SetInputYfrom(listInt[2]);
            board->SetInputYto(listInt[3]);
            
            board->SetTC(listInt[4]);
            
            // add local board into array
            localBoardArray->AddAt(board,board->GetId());
            localBoardMap->Add(board->GetId(),board);
          
        }
      }
    }
    
    return triggerCrates->GetSize();
  }
