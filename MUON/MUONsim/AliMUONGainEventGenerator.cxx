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

#include "AliMUONGainEventGenerator.h"

//-----------------------------------------------------------------------------
/// \class AliMUONGainEventGenerator
///
/// Generate raw data files that mimics the one we'll get from the 
/// online gain calibration procedure.
///
/// We start from one set of gain values, and one set of pedestal values,
/// to be found in the OCDB 
/// (generated e.g. with AliMUONCDB class)
///
/// We then use those gains and to generate n sets of pedestal values, 
/// that we store in OCDB in runs = [firstRunNumber+1,...firstRunNumber+n-1]
///
/// Then we loop from 0 to n-1, and for each, we generate a pedestal file,
/// using AliMUONPedestalEventGenerator, where the pedestals used are those
/// stored in the previous step.
///
/// Output files are supposed to be processed by the MUONTRKda (gain part)
/// to produce 4 (one per LDC) ascii files (with computed gain values), that
/// in turn will be processed by the Shuttle to put gain values in the OCDB.
/// Those last gains in OCDB can be then compared with the ones generated in 
/// the first step, to check that everything is OK in the whole procedure.
///
/// \author Laurent Aphecetche
//-----------------------------------------------------------------------------

#include "AliMUONCalibrationData.h"
#include "AliMUONVStore.h"
#include "AliMUONVCalibParam.h"
#include <TRandom.h>
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliMUONPedestalEventGenerator.h"
#include "AliLog.h"
#include <TROOT.h>
#include <TSystem.h>
#include "AliMUONDigitizerV3.h"
#include "AliCodeTimer.h"

/// \cond CLASSIMP
ClassImp(AliMUONGainEventGenerator)
/// \endcond

//_____________________________________________________________________________
AliMUONGainEventGenerator::AliMUONGainEventGenerator( Int_t sourceGainRunNumber,
                                                      Int_t sourcePedRunNumber,                                                      
                          Int_t nEventsPerFile, 
                          const char* dateBaseFileName)
  : TTask("AliMUONGainEventGenerator","Generate gain raw data files"),
  fNofEventsPerFile(nEventsPerFile),
  fSourcePedestalRunNumber(sourcePedRunNumber),
  fDateBaseFileName(dateBaseFileName),
  fSourceGains(AliMUONCalibrationData::CreateGains(sourceGainRunNumber)),
  fSourcePedestals(AliMUONCalibrationData::CreatePedestals(fSourcePedestalRunNumber))
{
    /// Ctor
    
    if (!fSourceGains)
    {
      AliFatal(Form("Cannot get gains for run %d",sourceGainRunNumber));
    }
    if (!fSourcePedestals)
    {
      AliFatal(Form("Cannot get pedestals for run %d",sourcePedRunNumber));
    }
}

//_____________________________________________________________________________
AliMUONGainEventGenerator::~AliMUONGainEventGenerator()
{  
  /// dtor
  delete fSourceGains;
  delete fSourcePedestals;
}

//_____________________________________________________________________________
void
AliMUONGainEventGenerator::Exec(Option_t*)
{
  /// Main method
  
  AliCodeTimer::Instance()->Reset();
  
  const Int_t kNInjections = 9;
  Float_t injections[kNInjections] = { 0,  200 , 400,  800, 1200, 1600, 
                                      2000, 2500, 3000 };
  
  for ( Int_t i = 0; i < kNInjections; ++i ) 
  {
    Int_t runNumber = fSourcePedestalRunNumber + i;
    if (i)
    {
      GeneratePedestals(runNumber,injections[i]);
    }
    TString pwd(gSystem->WorkingDirectory());
    TString dir(Form("%s/RUN%d",pwd.Data(),runNumber));
    AliInfo(Form("Creating directory %s",dir.Data()));
    gSystem->MakeDirectory(dir.Data());
    gSystem->ChangeDirectory(dir.Data());
    AliCodeTimerAuto(Form("generation of pedestal for run %d",runNumber),0);
    TString pedfile;
    if ( fDateBaseFileName.Length() > 0 ) pedfile = Form("%s.%d",fDateBaseFileName.Data(),runNumber);
    AliMUONPedestalEventGenerator pgen(runNumber,fNofEventsPerFile,pedfile.Data());
    AliInfo(Form("Generating pedestal events for injection number %d. Please be patient.",i));
    pgen.Exec("");
    gSystem->ChangeDirectory(pwd.Data());
  }
  
  AliCodeTimer::Instance()->Print();
}

//_____________________________________________________________________________
void
AliMUONGainEventGenerator::GeneratePedestals(Int_t runNumber, Float_t injection)
{
  /// Generate "pedestal" values for a given run, by "decalibrating"
  /// charge injection
  
  AliCodeTimerAuto(Form("Run %d injection %7.2f",runNumber,injection),0);
  TIter next(fSourceGains->CreateIterator());
  
  AliMUONVStore* generatedPedestals = fSourcePedestals->Create();
  
  AliMUONVCalibParam* gain;
  AliMUONVCalibParam* ped;
  
  while ( ( gain = static_cast<AliMUONVCalibParam*>(next()) ) )
  {
    ped = static_cast<AliMUONVCalibParam*>(fSourcePedestals->FindObject(gain->ID0(),gain->ID1()));
    
    AliMUONVCalibParam* genPed = static_cast<AliMUONVCalibParam*>(generatedPedestals->FindObject(gain->ID0(),gain->ID1()));
    if (!genPed)
    {
      genPed = static_cast<AliMUONVCalibParam*>(ped->Clone());
      generatedPedestals->Add(genPed);
    }
    
    for ( Int_t i = 0; i < ped->Size(); ++i ) 
    {
      Float_t mean = ped->ValueAsFloat(i,0);
      if ( mean == AliMUONVCalibParam::InvalidFloatValue() ) 
      {
        // non existing channel
        continue;
      }
      Int_t adc = AliMUONDigitizerV3::DecalibrateTrackerDigit(*ped,gain,i,
                                                              injection,kFALSE);

      Float_t res = (ped->ValueAsFloat(i,1)/mean);
            
      genPed->SetValueAsFloat(i,0,adc);
      genPed->SetValueAsFloat(i,1,adc*res);
    }
  }
  
  WriteToCDB(generatedPedestals,runNumber);
  
  delete generatedPedestals;
}

//_____________________________________________________________________________
void 
AliMUONGainEventGenerator::WriteToCDB(TObject* object, Int_t runNumber)
{
  /// Write a given object to OCDB
  
  AliCDBId id("MUON/Calib/Pedestals",runNumber,runNumber);
  AliCDBMetaData md;
  md.SetAliRootVersion(gROOT->GetVersion());
  md.SetComment("Pedestal values generated from AliMUONGainEventGenerator");
  md.SetResponsible("AliMUONGainEventGenerator");
  
  AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE_ROOT/OCDB/SHUTTLE/TestShuttle/TestCDB");
  man->Put(object,id,&md);
}
