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

/*

*/   
// T0 preprocessor:
// 2) takes data after  pass0 , 
// processes it, and stores either to OCDB .


#include "AliT0PreprocessorOffline.h"
#include "AliT0RecoParam.h"
#include "AliT0CalibTimeEq.h"
#include "AliCDBStorage.h"
#include "AliCDBMetaData.h"
#include "AliCDBManager.h"
#include "AliCTPTimeParams.h"
#include "AliLHCClockPhase.h"
#include "AliT0CalibSeasonTimeShift.h"
#include "AliT0CalibLatency.h"
#include "AliCDBEntry.h"
#include "AliLog.h"

#include <TTimeStamp.h>
#include <TFile.h>
#include <TObjString.h>
#include <TNamed.h>
#include "TClass.h"


ClassImp(AliT0PreprocessorOffline)

//____________________________________________________
AliT0PreprocessorOffline::AliT0PreprocessorOffline():
TNamed("AliT0PreprocessorOffline","AliT0PreprocessorOffline"),
  startRun(0),      
  endRun(0),  
  startTime(0),    
  endTime(0),     
  ocdbStorage("")  
{
  //constructor
}
//____________________________________________________

AliT0PreprocessorOffline::~AliT0PreprocessorOffline()
{
  //destructor

}
//____________________________________________________
void AliT0PreprocessorOffline::Process(TString filePhysName, Int_t ustartRun, Int_t uendRun, TString pocdbStorage)
{
  CalibOffsetChannels(filePhysName, ustartRun, uendRun, pocdbStorage);
  CalibT0sPosition(filePhysName, ustartRun, uendRun, pocdbStorage);
}
//____________________________________________________

void AliT0PreprocessorOffline::CalibOffsetChannels(TString filePhysName, Int_t ustartRun, Int_t uendRun, TString pocdbStorage)
{

  Float_t zero_timecdb[24]={0};
  Float_t *timecdb = zero_timecdb;
  Float_t *cfdvalue = zero_timecdb;
  Int_t badpmt=-1;
  //Processing data from DAQ Physics run
  AliInfo("Processing Time Offset between channels");
  if (pocdbStorage.Length()>0) ocdbStorage=pocdbStorage;
  else
    ocdbStorage="local://"+gSystem->GetFromPipe("pwd")+"/OCDB";
 
  AliT0CalibTimeEq *offline = new AliT0CalibTimeEq();
  Int_t writeok = offline->ComputeOfflineParams(filePhysName.Data(), timecdb, cfdvalue, badpmt);
  printf(" AliT0PreprocessorOffline::CalibOffsetChannels :: writeok %i \n", writeok);
  AliCDBMetaData metaData;
  metaData.SetBeamPeriod(1);
  metaData.SetResponsible("Alla Maevskaya");
  metaData.SetComment("Time equalizing result with slew");
  if (writeok<=0)  {
    AliCDBId* id1=NULL;
    id1=new AliCDBId("T0/Calib/TimeDelay", ustartRun, uendRun );
    AliCDBStorage* gStorage = AliCDBManager::Instance()->GetStorage(ocdbStorage);
    gStorage->Put(offline, (*id1), &metaData);
  }
  else {
    
    AliWarning(Form("writeok = %d data is not OK to be in OCDB",writeok));
      }		  
  
  delete offline;
 
}
//-------------------------------------------------------------------------------------
void AliT0PreprocessorOffline::CalibT0sPosition(TString filePhysName, Int_t ustartRun, Int_t uendRun, TString pocdbStorage)
{
  Float_t zero_timecdb[4]={0};
  Float_t *timecdb = zero_timecdb;
  if (pocdbStorage.Length()>0) ocdbStorage=pocdbStorage;
  else
    ocdbStorage="local://"+gSystem->GetFromPipe("pwd")+"/OCDB";
  
  AliT0CalibSeasonTimeShift *offline = new AliT0CalibSeasonTimeShift();
  Int_t writeok = offline->SetT0Par(filePhysName.Data(), timecdb);
  printf(" AliT0PreprocessorOffline::CalibT0sPosition :: writeok %i \n", writeok);
  AliCDBMetaData metaData;
  metaData.SetBeamPeriod(1);
  metaData.SetResponsible("Alla Maevskaya");
  metaData.SetComment("Time equalizing result with slew");
  
  if (writeok <= 0)  {
    AliCDBId* id1=NULL;
    id1=new AliCDBId("T0/Calib/TimeAdjust", ustartRun, uendRun);
    AliCDBStorage* gStorage = AliCDBManager::Instance()->GetStorage(ocdbStorage);
    gStorage->Put(offline, (*id1), &metaData);
  }

}
