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

#include <cstdlib>
#include <iostream>
#include "AliMUONTrackerIO.h"

/// \class AliMUONTrackerIO
///
/// Reader class for ASCII calibration files for MUON tracker : 
/// converts those ASCII files into AliMUONVStore (suitable to e.g. feed
/// the OCDB).
///
/// \author Laurent Aphecetche, Subatech

using std::ostringstream;
using std::istringstream;
using std::cout;
using std::endl;
/// \cond CLASSIMP
ClassImp(AliMUONTrackerIO)
/// \endcond

#include "AliDCSValue.h"
#include "AliLog.h"
#include "AliMUONCalibParamND.h"
#include "AliMUONCalibParamNF.h"
#include "AliMUONVStore.h"
#include "AliMpConstants.h"
#include "AliMpDDLStore.h"
#include "AliMpDEManager.h"
#include "AliMpDetElement.h"
#include <Riostream.h>
#include <TClass.h>
#include <TObjString.h>
#include <TSystem.h>
#include <sstream>

//_____________________________________________________________________________
AliMUONTrackerIO::AliMUONTrackerIO()
{
  /// ctor
}

//_____________________________________________________________________________
AliMUONTrackerIO::~AliMUONTrackerIO()
{
  /// dtor
}

//_____________________________________________________________________________
Int_t 
AliMUONTrackerIO::ReadOccupancy(const char* filename,AliMUONVStore& occupancyMap)
{
  /// Read occupancy file created by online DA
  /// and append values to the occupancyMap store.
  /// Expected format of the file is :
  /// busPatchId manuId sumofn nevt
  
  TString sFilename(gSystem->ExpandPathName(filename));
  
  std::ifstream in(sFilename.Data());
  if (!in.good()) 
  {
    return kCannotOpenFile;
  }
  
  TString datastring;
  ostringstream stream;
  char line[1024];
  while ( in.getline(line,1024) )
  	stream << line << "\n";
  
  in.close();
  
  return DecodeOccupancy(stream.str().c_str(),occupancyMap);
  
}

//_____________________________________________________________________________
Int_t 
AliMUONTrackerIO::DecodeOccupancy(const char* data, AliMUONVStore& occupancyMap)
{
  /// Decode occupancy string created append values to the occupancyMap store.
  /// Expected format of the file is :
  /// busPatchId manuId sumofn nevt
 
  if ( ! AliMpDDLStore::Instance(kFALSE) )
  {
    AliErrorClass("Mapping not loaded. Cannot work");
    return kNoMapping;
  }
  
  char line[1024];
  istringstream in(data);
  
  Int_t n(0);
  
  while ( in.getline(line,1024) )
  {
    AliDebugClass(3,Form("line=%s",line));
    if ( line[0] == '/' && line[1] == '/' ) continue;
    std::istringstream sin(line);
    
    Int_t busPatchId, manuId;
    Int_t numberOfEvents;
    Double_t sumn;

    sin >> busPatchId >> manuId >> sumn >> numberOfEvents;
    
    if ( busPatchId == -1 && manuId == -1 && sumn == 0 && numberOfEvents == 0 )
    {
      /// DA could not produce information (because run failed somehow). 
      /// That's OK, but let the world know about it
      return kNoInfoFile;
    }
    
    Int_t detElemId = AliMpDDLStore::Instance()->GetDEfromBus(busPatchId);

    AliMpDetElement* de = AliMpDDLStore::Instance()->GetDetElement(detElemId);

    Int_t numberOfChannelsInManu = -1;
    
    if (de) numberOfChannelsInManu = de->NofChannelsInManu(manuId);

    if ( numberOfChannelsInManu <= 0 ) 
    {
      AliErrorClass(Form("BP %5d DE %5d MANU %5d nchannels=%d",busPatchId,detElemId,manuId,numberOfChannelsInManu));
      continue;      
    }
    
    AliMUONVCalibParam* occupancy = 
    static_cast<AliMUONVCalibParam*>(occupancyMap.FindObject(detElemId,manuId));
    if (occupancy) 
    {
      AliErrorClass(Form("DE %5d MANU %5d is already there ?!",detElemId,manuId));
      continue;
    }
        
    occupancy = new AliMUONCalibParamND(5,1,detElemId,manuId,0);

    occupancyMap.Add(occupancy);
    
    occupancy->SetValueAsDouble(0,0,sumn);
    occupancy->SetValueAsDouble(0,1,sumn); // with only 0 and 1s, sumw = sumw2 = sumn
    occupancy->SetValueAsDouble(0,2,sumn);
    occupancy->SetValueAsInt(0,3,numberOfChannelsInManu);
    occupancy->SetValueAsInt(0,4,numberOfEvents);
    ++n;
  }
  
  return n;
}

//_____________________________________________________________________________
Int_t 
AliMUONTrackerIO::ReadPedestals(const char* filename, AliMUONVStore& pedStore)
{
  /// Read pedestal file (produced by the MUONTRKda.exe program for instance)
  /// and append the read values into the given VStore
  /// To be used when the input is a file (for instance when reading data 
  /// from the OCDB).
  
  TString sFilename(gSystem->ExpandPathName(filename));
  
  std::ifstream in(sFilename.Data());
  if (!in.good()) 
  {
    return kCannotOpenFile;
  }
  
  ostringstream stream;
  char line[1024];
  while ( in.getline(line,1024) )
  	stream << line << "\n";
  
  in.close();

  return DecodePedestals(stream.str().c_str(),pedStore);
  
}

//_____________________________________________________________________________
Int_t 
AliMUONTrackerIO::DecodePedestals(const char* data, AliMUONVStore& pedStore)
{
  /// Read pedestal Data (produced by the MUONTRKda.exe program for instance)
  /// and append the read values into the given VStore
  /// To be used when the input is a TString (for instance when getting data 
  /// from AMORE DB).
  
  char line[1024];
  Int_t busPatchID, manuID, manuChannel;
  Float_t pedMean, pedSigma;
  Int_t n(0);
  istringstream in(data);
  
  while ( in.getline(line,1024) )
  {
    AliDebugClass(3,Form("line=%s",line));
    if ( line[0] == '/' && line[1] == '/' ) continue;
    std::istringstream sin(line);
    sin >> busPatchID >> manuID >> manuChannel >> pedMean >> pedSigma;
    Int_t detElemID = AliMpDDLStore::Instance()->GetDEfromBus(busPatchID);
    
    if ( !AliMpDEManager::IsValidDetElemId(detElemID) )
    {
      AliErrorClass(Form("Got an invalid DE = %d from busPatchId=%d manuId=%d",
                         detElemID,busPatchID,manuID));
      continue;
    }
    
    AliDebugClass(3,Form("BUSPATCH %3d DETELEMID %4d MANU %3d CH %3d MEAN %7.2f SIGMA %7.2f",
                    busPatchID,detElemID,manuID,manuChannel,pedMean,pedSigma));
		    
    AliMUONVCalibParam* ped = 
      static_cast<AliMUONVCalibParam*>(pedStore.FindObject(detElemID,manuID));
    if (!ped) 
    {
      ped = new AliMUONCalibParamNF(2,AliMpConstants::ManuNofChannels(),
                                    detElemID,manuID,
                                    AliMUONVCalibParam::InvalidFloatValue());  
      pedStore.Add(ped);
    }
    ped->SetValueAsFloat(manuChannel,0,pedMean);
    ped->SetValueAsFloat(manuChannel,1,pedSigma);
    ++n;
  }

  return n;
}

//_____________________________________________________________________________
Int_t 
AliMUONTrackerIO::ReadConfig(const char* filename, AliMUONVStore& confStore)
{
  /// Read config file (produced by the MUONTRKda.exe program for instance)
  /// and append the read values into the given VStore
  /// To be used when the input is a file (for instance when reading data 
  /// from the OCDB).
  
  TString sFilename(gSystem->ExpandPathName(filename));
  
  std::ifstream in(sFilename.Data());
  if (!in.good()) 
  {
    return kCannotOpenFile;
  }
  
  ostringstream stream;
  char line[1024];
  while ( in.getline(line,1024) )
  	stream << line << "\n";
  
  in.close();
  
  return DecodeConfig(stream.str().c_str(),confStore);
}

//_____________________________________________________________________________
Int_t 
AliMUONTrackerIO::DecodeConfig(const char* data, AliMUONVStore& confStore)
{
  /// Read config data (produced by the MUONTRKda.exe program for instance)
  /// and append the read values into the given VStore
  /// To be used when the input is a TString (for instance when getting data 
  /// from AMORE DB).

  char line[1024];
  Int_t busPatchID, manuID;
  Int_t n(0);
  istringstream in(data);
  
  while ( in.getline(line,1024) )
  {
    AliDebugClass(3,Form("line=%s",line));
    if ( line[0] == '#' ) 
    {
      continue;
    }
    std::istringstream sin(line);
    sin >> busPatchID >> manuID;

    Int_t detElemID = AliMpDDLStore::Instance()->GetDEfromBus(busPatchID);

    if ( !AliMpDEManager::IsValidDetElemId(detElemID) )
    {
      AliErrorClass(Form("Got an invalid DE = %d from busPatchId=%d manuId=%d",
                         detElemID,busPatchID,manuID));
      continue;
    }
    
    AliMUONVCalibParam* conf = 
    static_cast<AliMUONVCalibParam*>(confStore.FindObject(detElemID,manuID));
    if (!conf) 
    {
      conf = new AliMUONCalibParamNF(1,1,detElemID,manuID,1);
      confStore.Add(conf);
    }
    ++n;
  }
  
  return n;
}

//_____________________________________________________________________________
Int_t 
AliMUONTrackerIO::WriteConfig(ofstream& out, const AliMUONVStore& confStore)
{
  /// Write the conf store as an ASCII file
  /// Note that we are converting (back) the detElemId into a busPatchId
  /// Return the number of lines written
  
  if ( !AliMpDDLStore::Instance() ) 
  {
    cout << "ERROR: mapping not loaded. Cannot work" << endl;
    return 0;
  }
  
  TIter next(confStore.CreateIterator());
  AliMUONVCalibParam* param;
  Int_t n(0);
  
  while ( (param=static_cast<AliMUONVCalibParam*>(next())) )
  {
    Int_t detElemId = param->ID0();
    Int_t manuId = param->ID1();
    
    Int_t busPatchId = AliMpDDLStore::Instance()->GetBusPatchId(detElemId,manuId);
    ++n;
    
    out << busPatchId << " " << manuId << endl;
  }
  return n;
}

