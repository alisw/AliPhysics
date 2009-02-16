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
#include "AliMUONTrackerIO.h"

/// \class AliMUONTrackerIO
///
/// Reader class for ASCII calibration files for MUON tracker : 
/// converts those ASCII files into AliMUONVStore (suitable to e.g. feed
/// the OCDB).
///
/// \author Laurent Aphecetche, Subatech

/// \cond CLASSIMP
ClassImp(AliMUONTrackerIO)
/// \endcond

#include "AliLog.h"
#include "AliMUONCalibParamNF.h"
#include "AliMUONVStore.h"
#include "AliMpConstants.h"
#include "AliMpDDLStore.h"
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
  
  TString datastring;
  datastring.ReadFile(in);
    
  in.close();

  return DecodePedestals(datastring,pedStore);
  
}

//_____________________________________________________________________________
Int_t 
AliMUONTrackerIO::DecodePedestals(TString data, AliMUONVStore& pedStore)
{
  /// Read pedestal Data (produced by the MUONTRKda.exe program for instance)
  /// and append the read values into the given VStore
  /// To be used when the input is a TString (for instance when getting data 
  /// from AMORE DB).
  
  char line[1024];
  Int_t busPatchID, manuID, manuChannel;
  Float_t pedMean, pedSigma;
  Int_t n(0);
  istringstream in(data.Data());
  
  while ( in.getline(line,1024) )
  {
    AliDebugClass(3,Form("line=%s",line));
    if ( line[0] == '/' && line[1] == '/' ) continue;
    std::istringstream sin(line);
    sin >> busPatchID >> manuID >> manuChannel >> pedMean >> pedSigma;
    Int_t detElemID = AliMpDDLStore::Instance()->GetDEfromBus(busPatchID);
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
AliMUONTrackerIO::ReadGains(const char* filename, AliMUONVStore& gainStore,
                            TString& comment)
{
  /// Read gain file (produced by the MUONTRKda.exe program for instance)
  /// and append the read values into the given VStore
  /// To be used when the input is a file (for instance when reading data 
  /// from the OCDB).
  
  comment = "";
  
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
  datastring = TString(stream.str().c_str());
  
  in.close();
  
  return DecodeGains(datastring,gainStore,comment);

}

//_____________________________________________________________________________
Int_t 
AliMUONTrackerIO::DecodeGains(TString data, AliMUONVStore& gainStore,
                            TString& comment)
{
  /// Read gain file (produced by the MUONTRKda.exe program for instance)
  /// and append the read values into the given VStore
  /// To be used when the input is a TString (for instance when getting data 
  /// from AMORE DB).
  
  char line[1024];
  istringstream in(data.Data());
  Int_t busPatchID, manuID, manuChannel;
  Float_t a0, a1;
  Int_t thres;
  UInt_t qual;
  const Int_t kSaturation(3000); // FIXME: how to get this number ?
  Int_t n(0);
  Int_t runNumber(-1);
  Int_t* runs(0x0);
  Int_t* dac(0x0);
  Int_t nDAC(0);
  Int_t iDAC(0);
  
  while ( in.getline(line,1024) )
  {
    if ( strlen(line) < 10 ) continue;
    if ( line[0] == '/' && line[1] == '/' ) 
    {
      TString sline(line);
      if ( sline.Contains("DUMMY") )
      {
        AliDebugClass(1,"Got a dummy file here");
        return kDummyFile;
      }
      if ( sline.Contains("* Run") )
      {
        TObjArray* a = sline.Tokenize(":");
        if ( a->GetLast() >= 1 ) 
        {
          TString s = static_cast<TObjString*>(a->At(1))->String();
          runNumber = s.Atoi();
          AliDebugClass(1,Form("runNumber is %d",runNumber));
        }            
      }
      if ( sline.Contains("DAC values") )
      {
        nDAC = TString(sline(2,sline.Length()-2)).Atoi();
        AliDebugClass(1,Form("# of DAC values = %d",nDAC));
        if ( nDAC > 0 )
        {
          if ( nDAC < 100 ) 
          {
            runs = new Int_t[nDAC];
            dac = new Int_t[nDAC];
            // skip two lines
            in.getline(line,1024);
            in.getline(line,1024);
            // then get run and dac values
            for ( Int_t i = 0; i < nDAC; ++i ) 
            {
              in.getline(line,1024);
              Int_t a,b;
              sscanf(line,"// %d %d",&a,&b);
              runs[iDAC] = a;
              dac[iDAC] = b;
              AliDebugClass(1,Form("RUN %d is DAC %d",runs[iDAC],dac[iDAC]));
              ++iDAC;
            }
          }
          else
          {
            AliErrorClass(Form("Something went wrong, as I get too big nDAC = %d",nDAC));
            nDAC = 0;
            return kFormatError;
          }
        }
      }
      continue;
    }
    
    sscanf(line,"%d %d %d %f %f %d %x",&busPatchID,&manuID,&manuChannel,
           &a0,&a1,&thres,&qual); 
    AliDebugClass(3,Form("line=%s",line));
    Int_t detElemID = AliMpDDLStore::Instance()->GetDEfromBus(busPatchID);
    AliDebugClass(3,Form("BUSPATCH %3d DETELEMID %4d MANU %3d CH %3d A0 %7.2f "
                    "A1 %e THRES %5d QUAL %x",
                    busPatchID,detElemID,manuID,manuChannel,a0,a1,thres,qual));
    if ( qual == 0 ) continue;
    
    AliMUONVCalibParam* gain = 
      static_cast<AliMUONVCalibParam*>(gainStore.FindObject(detElemID,manuID));
    
   if (!gain) 
    {
      gain = new AliMUONCalibParamNF(5,AliMpConstants::ManuNofChannels(),detElemID,manuID,0);
      gainStore.Add(gain);
    }
    gain->SetValueAsFloat(manuChannel,0,a0);
    gain->SetValueAsFloat(manuChannel,1,a1);
    gain->SetValueAsInt(manuChannel,2,thres);
    gain->SetValueAsInt(manuChannel,3,qual);
    gain->SetValueAsInt(manuChannel,4,kSaturation);
    ++n;
  }

  comment = "";
  
  if ( runNumber > 0 )
  {
    comment = Form("RUN %d",runNumber);
  }
  
  for ( Int_t i = 0; i < nDAC; ++i )
  {
    comment += Form(";(RUN %d = DAC %d)",runs[i],dac[i]);
  }
  
  delete[] runs;
  delete[] dac;
  
  return n;
}

//_____________________________________________________________________________
Int_t
AliMUONTrackerIO::ReadCapacitances(const char* file, AliMUONVStore& capaStore)
{
  /// Read capacitance file
  /// and append the read values into the given VStore
  
  ifstream in(gSystem->ExpandPathName(file));
  if (in.bad()) return kCannotOpenFile;
  
  Int_t ngenerated(0);
  
  char line[1024];
  Int_t serialNumber(-1);
  AliMUONVCalibParam* param(0x0);
  
  while ( in.getline(line,1024,'\n') )
  {
    if ( isdigit(line[0]) ) 
    {
      serialNumber = atoi(line);
      param = static_cast<AliMUONVCalibParam*>(capaStore.FindObject(serialNumber));
      if (param)
      {
        AliErrorClass(Form("serialNumber %d appears several times !",serialNumber));
        continue;
      }
      param = new AliMUONCalibParamNF(2,AliMpConstants::ManuNofChannels(),serialNumber,0,1.0);
      Bool_t ok = capaStore.Add(param);
      if (!ok)
      {
        AliErrorClass(Form("Could not set serialNumber=%d",serialNumber));
        continue;
      }      
      continue;
    }
    Int_t channel;
    Float_t capaValue;
    Float_t injectionGain;
    sscanf(line,"%d %f %f",&channel,&capaValue,&injectionGain);
    AliDebugClass(1,Form("SerialNumber %10d Channel %3d Capa %f injectionGain %f",
                    serialNumber,channel,capaValue,injectionGain));
    param->SetValueAsFloat(channel,0,capaValue);
    param->SetValueAsFloat(channel,1,injectionGain);
    ++ngenerated;
  }
  
  in.close();
  
  return ngenerated;
}
