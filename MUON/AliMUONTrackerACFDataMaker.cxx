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

#include "AliMUONTrackerACFDataMaker.h"

#include "AliMUONTrackerData.h"
#include "AliMUONTrackerIO.h"
#include "AliMUON2DMap.h"
#include "AliMUON1DMap.h"
#include "AliLog.h"
#include <TString.h>
#include <TSystem.h>

///\class AliMUONTrackerACFDataMaker
///
/// Producer of AliMUONVTrackerData from ASCII calibration files
///
///\author Laurent Aphecetche, Subatech

///\cond CLASSIMP
ClassImp(AliMUONTrackerACFDataMaker)
///\endcond

//_____________________________________________________________________________
AliMUONTrackerACFDataMaker::AliMUONTrackerACFDataMaker(const char* acfPath,
                                                        const char* type)
: AliMUONVTrackerDataMaker(),
  fIsValid(kTRUE),
  fData(0x0),
  fSource(Form("%s-%s",acfPath,type))
{
    /// Ctor

    static Int_t number(0);
    
    ++number;
    
    AliMUONVStore* store(0x0);
    
    TString stype(type);
    stype.ToUpper();
    Bool_t isSingleEvent(kTRUE);
    TString filename(gSystem->ExpandPathName(acfPath));
    
    if ( stype == "PEDESTALS" )
    {
      fData = new AliMUONTrackerData(Form("PED",number),"Pedestals",2,isSingleEvent);
      fData->SetDimensionName(0,"Mean");
      fData->SetDimensionName(1,"Sigma");
      store = new AliMUON2DMap(kTRUE);
      AliMUONTrackerIO::ReadPedestals(filename.Data(),*store);
    }
    else if ( stype == "GAINS" ) 
    {
      fData = new AliMUONTrackerData(Form("GAIN%d",number),"Gains",5,isSingleEvent);
      fData->SetDimensionName(0,"a1");
      fData->SetDimensionName(1,"a2");
      fData->SetDimensionName(2,"thres");
      fData->SetDimensionName(3,"qual");
      fData->SetDimensionName(4,"sat");
      store = new AliMUON2DMap(kTRUE);
      TString comment;
      AliMUONTrackerIO::ReadGains(filename.Data(),*store,comment);
    }
    else if ( stype == "CAPACITANCES" )
    {
      fData = new AliMUONTrackerData(Form("CAPA%d",number),"Capacitances",2,isSingleEvent);
      fData->SetDimensionName(0,"Capa");
      fData->SetDimensionName(1,"Injection gain");
      store = new AliMUON1DMap(20000);
      AliMUONTrackerIO::ReadCapacitances(filename.Data(),*store);
    }
    
    if (!store)
    {
      fIsValid = kFALSE;
      delete fData;
      fData = 0x0;
      AliError("Could not create store");
      return;
    }
    
    fData->Add(*store);
}

//_____________________________________________________________________________
AliMUONTrackerACFDataMaker::~AliMUONTrackerACFDataMaker()
{
  /// dtor
  delete fData;
}

//_____________________________________________________________________________
Long64_t 
AliMUONTrackerACFDataMaker::Merge(TCollection*)
{
  /// Merge
  AliError("Not implemented. Does it have sense ?");
  return 0;
}

