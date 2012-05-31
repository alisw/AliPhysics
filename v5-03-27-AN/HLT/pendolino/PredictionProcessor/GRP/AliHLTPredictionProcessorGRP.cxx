// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no         *
//*                  for The ALICE HLT Project.                            *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

//  @file   AliHLTPredictionProcessorGRP.cxx
//  @author Matthias Richter
//  @date   2010-04-12
//  @brief  Prediction processor for the GRP entry
// 

#include "AliHLTPredictionProcessorGRP.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliGRPObject.h"
#include "AliCDBMetaData.h"
#include "AliCDBId.h"
#include <cassert>

//#include <TObjArray.h>
//#include <AliDCSValue.h>
//#include <TTimeStamp.h>


ClassImp(AliHLTPredictionProcessorGRP)

const float kNormAtmosPressure = 1013.0; // norm pressure in mbar

AliHLTPredictionProcessorGRP::AliHLTPredictionProcessorGRP(AliHLTPendolino* pendolino)
  : AliHLTPredictionProcessorInterface("GRP", pendolino)
{
  // constructor
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}


AliHLTPredictionProcessorGRP::~AliHLTPredictionProcessorGRP()
{
  // destructor
}


UInt_t AliHLTPredictionProcessorGRP::makePrediction(Bool_t /*doPrediction*/)
{
  // switch for prediction making

  // not exactly clear what to do here, have to check the interface
  return 0;
}


void AliHLTPredictionProcessorGRP::Initialize(Int_t run, UInt_t startTime, 
					      UInt_t endTime)
{
  // initializes AliHLTPredictionProcessorGRP
  AliHLTPredictionProcessorInterface::Initialize(run, startTime, endTime);

  TString msg("Initialized GRP PredictProc. Run: ");
  msg += GetRunNumber();
  msg += ", start time: ";
  msg += StartTime();
  msg += ", end time: ";
  msg += EndTime();
  msg += ".";	
  Log(msg.Data());
}


UInt_t AliHLTPredictionProcessorGRP::Process(TMap* dcsAliasMap)
{
  // processes the DCS value map
  
  if (!dcsAliasMap) return 9;
  if (dcsAliasMap->GetEntries() == 0 ) return 9;

  Float_t l3Current=0.0;
  Bool_t l3Polarity=0; // positive field according to LHC convention
  Float_t dipoleCurrent=0.0;
  Bool_t dipolePolarity=0; // positive field according to LHC convention
  Float_t cavernAtmosPressure=kNormAtmosPressure;
  Float_t cavernAtmosPressure2=kNormAtmosPressure;
  Float_t surfaceAtmosPressure=kNormAtmosPressure;
  
  //Bool_t bRet = kTRUE;
  const char* key="";

  key="L3Current";
  if (!GetSensorValue(dcsAliasMap, key, &l3Current)) {
    Log(Form("failed to extract %s from alias map", key));
    //bRet=kFALSE;
  }

  key="L3Polarity";
  if (!GetSensorValue(dcsAliasMap, key, &l3Polarity)) {
    Log(Form("failed to extract %s from alias map", key));
    //bRet=kFALSE;
  }

  key="DipoleCurrent";
  if (!GetSensorValue(dcsAliasMap, key, &dipoleCurrent)) {
    Log(Form("failed to extract %s from alias map", key));
    //bRet=kFALSE;
  }

  key="DipolePolarity";
  if (!GetSensorValue(dcsAliasMap, key, &dipolePolarity)) {
    Log(Form("failed to extract %s from alias map", key));
    //bRet=kFALSE;
  }

  key="CavernAtmosPressure";
  if (!GetSensorValue(dcsAliasMap, key, &cavernAtmosPressure)) {
    Log(Form("failed to extract %s from alias map", key));
    //bRet=kFALSE;
  }

  key="CavernAtmosPressure2";
  if (!GetSensorValue(dcsAliasMap, key, &cavernAtmosPressure2)) {
    Log(Form("failed to extract %s from alias map", key));
    //bRet=kFALSE;
  }

  key="SurfaceAtmosPressure";
  if (!GetSensorValue(dcsAliasMap, key, &surfaceAtmosPressure)) {
    Log(Form("failed to extract %s from alias map", key));
    //bRet=kFALSE;
  }

  AliGRPObject* grpObj=NULL;
  AliGRPObject* grpExist=NULL;
  AliCDBMetaData* cdbMetaData=NULL;
  bool bUpdate=false;

  // read existing GRP object
  AliCDBEntry* pEntry=GetFromOCDB("GRP","Data");
  if (pEntry) {
    grpExist=dynamic_cast<AliGRPObject*>(pEntry->GetObject());
    if (grpExist) {
      // possible values for polarities: 0 and 1
      // nominal value of currents is thousands of ampere 
      bUpdate=bUpdate || (TMath::Abs(grpExist->GetL3Polarity()- l3Polarity)>0.5);
      bUpdate=bUpdate || (TMath::Abs(grpExist->GetL3Current(AliGRPObject::kMean) - l3Current)>10.0);
      bUpdate=bUpdate || (TMath::Abs(grpExist->GetDipolePolarity()- dipolePolarity)>0.5);
      bUpdate=bUpdate || (TMath::Abs(grpExist->GetDipoleCurrent(AliGRPObject::kMean) - dipoleCurrent)>10.0);
      bUpdate=bUpdate 
	|| (grpExist->GetCavernAtmosPressure()==NULL)
	|| (TMath::Abs(grpExist->GetCavernAtmosPressure()->GetValue(0) - cavernAtmosPressure)> 1.0);
      bUpdate=bUpdate 
	|| (grpExist->GetCavernAtmosPressure2()==NULL)
	|| (TMath::Abs(grpExist->GetCavernAtmosPressure2()->GetValue(0) - cavernAtmosPressure2) > 1.0);
      bUpdate=bUpdate 
	|| (grpExist->GetSurfaceAtmosPressure()==NULL)
	|| (TMath::Abs(grpExist->GetSurfaceAtmosPressure()->GetValue(0) - surfaceAtmosPressure) > 1.0);


      cdbMetaData=pEntry->GetMetaData();
    } else {
      bUpdate=true;
    }
  } else {
    bUpdate=true;
  }
  // TString msg="AliHLTPredictionProcessorGRP::Process: update ";
  // msg+=bUpdate;
  // msg+="   \nL3 current = "; msg+=l3Current;
  // msg+="   \nL3 polarity = "; msg+=l3Polarity;
  // msg+="   \nDipole current = "; msg+=dipoleCurrent;
  // msg+="   \nDipole polarity = "; msg+=l3Polarity;
  // msg+="   \nCavern pressure = "; msg+=cavernAtmosPressure;
  // msg+="   \nCavern pressure2 = "; msg+=cavernAtmosPressure2;
  // msg+="   \nSurface pressure = "; msg+=surfaceAtmosPressure;
  // Log(msg.Data());

  // don't write if object is there and up-to-date
  if (grpExist && !bUpdate) return 0;

  if (grpExist) {
    // update existing object
    grpObj=grpExist;
  } else {
    // generate GRP object
    // TODO: correct energy has to be extracted somewhere
    grpObj=new AliGRPObject;
    float cmsEnergy=14000;
    if (grpObj) {
      grpObj->SetBeamEnergy(cmsEnergy/0.120); // LHC convention
      grpObj->SetBeamType("p-p");
    }
  }
  if (!grpObj) {
    Log(" *** Failed to generate GRP object");
    return 7; // the return value needs to be checked
  }

  grpObj->SetL3Current(l3Current,(AliGRPObject::Stats)0);
  grpObj->SetDipoleCurrent(dipoleCurrent,(AliGRPObject::Stats)0);  
  grpObj->SetL3Polarity(l3Polarity);  
  grpObj->SetDipolePolarity(dipolePolarity);
  grpObj->SetPolarityConventionLHC();                    // LHC convention +/+ current -> -/- field main components
  grpObj->SetCavernAtmosPressure(CreateSensor("CavernAtmosPressure",
					      cavernAtmosPressure,
					      GetStartTime(),
					      GetEndTime())
				 );
  grpObj->SetCavernAtmosPressure2(CreateSensor("CavernAtmosPressure2",
					       cavernAtmosPressure2,
					       GetStartTime(),
					       GetEndTime())
				  );
  grpObj->SetSurfaceAtmosPressure(CreateSensor("SurfaceAtmosPressure",
					       surfaceAtmosPressure,
					       GetStartTime(),
					       GetEndTime())
				  );

  if (!cdbMetaData) {
    cdbMetaData=new AliCDBMetaData;
    cdbMetaData->SetResponsible("Matthias.Richter@cern.ch");
    cdbMetaData->SetComment(Form("GRP entry for the magnetic field initialization of HLT components, produced by %s", ClassName()));
  }

  // note: 'validityStart' specifies run no w.r.t. current run no -> thus 0
  // create with run specific validity -> kFALSE 
  if (Store("GRP", "Data", grpObj, cdbMetaData, 0, kFALSE)) {
  } else {
    Log(" *** Failed to store GRP object");
    return 7;
  }
  if (!grpExist) {
    if (grpObj) delete grpObj;
    grpObj=NULL;
    if (cdbMetaData) delete cdbMetaData;
    cdbMetaData=NULL;
  }
  
  return 0;
}

TMap* AliHLTPredictionProcessorGRP::produceTestData(TString /*aliasName*/)
{
  // produces test data for AliHLTPredictionProcessorGRP
  TMap* resultMap = 0;

  return resultMap;
}

bool AliHLTPredictionProcessorGRP::CreateInitialGRPEntry(int runno, const char* beamtype, const char* runtype, const char* /*detectorList*/)
{
  // Create the initial GRP entry.
  // The beam type and run type parameters are propagated form the ECS to
  // the run manager, and provided as arguments to the pendolino.
  // The initial GRP entry is created ignoring the magnetic field, this
  // is updated by the pendolino afterwords.

  AliGRPObject* grpObj=new AliGRPObject;
  float cmsEnergy=14000; // check if this can be obtained from somewhere
  grpObj->SetBeamEnergy(cmsEnergy/0.120); // LHC convention
  grpObj->SetBeamType(beamtype);
  grpObj->SetRunType(runtype);
  grpObj->SetL3Current(0.0,(AliGRPObject::Stats)0);
  grpObj->SetDipoleCurrent(0.0,(AliGRPObject::Stats)0);  
  grpObj->SetL3Polarity(0);  
  grpObj->SetDipolePolarity(0);
  grpObj->SetPolarityConventionLHC();
  grpObj->SetTimeStart(time(NULL));   // check if time can be fetched from ECS

  /* the following needs to be checked
  grpObj->SetLHCPeriod(getenv("LHC_PERIOD"));
  grpObj->SetNumberOfDetectors(numOfDetectors);
  grpObj->SetDetectorMask(detectMask);  
  */
  
  AliCDBManager* man = AliCDBManager::Instance();
  if (!man || !man->IsDefaultStorageSet()) {
    return false;
  }

  AliCDBPath cdbPath("GRP/GRP/Data");
  AliCDBId cdbId(cdbPath, runno, runno);
  AliCDBMetaData cdbMetaData;
  cdbMetaData.SetResponsible("matthias.richter@cern.ch");
  cdbMetaData.SetComment(Form("Online GRP entry, produced by AliHLTPredictionProcessorGRP"));
  return man->Put(grpObj, cdbId, &cdbMetaData);
}
