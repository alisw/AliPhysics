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
#include "AliGRPObject.h"

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
  ocdbStorage(0),
  fNewDArun(9999999),
  fStatusDelay(0),
  fStatusAdjust(0)
  
{
  //constructor
}
//____________________________________________________

AliT0PreprocessorOffline::~AliT0PreprocessorOffline()
{
  //destructor

}
//____________________________________________________
void AliT0PreprocessorOffline::Process(TString filePhysName, Int_t ustartRun, Int_t uendRun, AliCDBStorage* pocdbStorage)
{
  if ( ustartRun < fNewDArun) 
    CalibOffsetChannels(filePhysName, ustartRun, uendRun, pocdbStorage);
  CalibT0sPosition(filePhysName, ustartRun, uendRun, pocdbStorage);
}
//____________________________________________________

void AliT0PreprocessorOffline::CalibOffsetChannels(TString filePhysName, Int_t ustartRun, Int_t uendRun, AliCDBStorage* pocdbStorage)
{

  Float_t zero_timecdb[24]={0};
  Float_t zero_cfdvalue[100]={0};
  Float_t *timecdb = zero_timecdb;
  Float_t *cfdvalue = zero_cfdvalue;
  Float_t cfd[24][5];
  for(Int_t i=0; i<24; i++) 
    for (Int_t i0=0; i0<5; i0++)
      cfd[i][i0] = 0;
  Int_t badpmt=-1;
  Float_t meanvertex=0;
  //Processing data from DAQ Physics run
  AliInfo("Processing Time Offset between channels");
  if (pocdbStorage) ocdbStorage=pocdbStorage;
  else {
   TString localStorage = "local://"+gSystem->GetFromPipe("pwd")+"/OCDB"; 
   ocdbStorage= AliCDBManager::Instance()->GetStorage(localStorage.Data());
  }
 AliCDBEntry* entryGRP = AliCDBManager::Instance()->Get("GRP/GRP/Data");
 AliGRPObject* grpData = dynamic_cast<AliGRPObject*>(entryGRP->GetObject());
 AliCDBEntry *entryCalib = AliCDBManager::Instance()->Get("T0/Calib/TimeDelay");
  if(!entryCalib) {
    AliWarning(Form("Cannot find any AliCDBEntry for [Calib, TimeDelay]!"));
  }
  else
    {
      AliT0CalibTimeEq *clb = (AliT0CalibTimeEq*)entryCalib->GetObject();
      //  clb->Print();
      timecdb = clb->GetTimeEq();
      for(Int_t i=0; i<24; i++) 
	for (Int_t i0=0; i0<5; i0++)  cfd[i][i0] = clb->GetCFDvalue(i, i0);
      meanvertex=clb->GetMeanVertex();
    }
  Float_t pedestalLHC15fghij[24] ={1428., 1421.,  1488.,  1488., 1582., 1405.,
				   1392., 1373.,  1591.,  1333., 1449., 1091.,
				   1607., 1350.,  1445.,  1526., 1379., 1695,
				   1475., 1448.,  1573.,  1466., 1334., 1441.};
  
  for(Int_t i=0; i<24; i++) cfdvalue[i]=cfd[i][0];
  for(Int_t i=24; i<48; i++) cfdvalue[i]=cfd[i-24][1];
  cfdvalue[48] = meanvertex;
  for(Int_t i=49; i<51; i++) cfdvalue[i]=cfd[i-49][2];
  for(Int_t i=52; i<76; i++) {
    cfdvalue[i]=cfd[i-52][3];
    TString LHCperiod = grpData->GetLHCPeriod();
    if(cfdvalue[i]<1 && LHCperiod.Contains("15") )  cfdvalue[i] = pedestalLHC15fghij[i-52];
    printf("@@@ %i pedestal %f\n",i,cfdvalue[i]);
  } 
  AliT0CalibTimeEq *offline = new AliT0CalibTimeEq();
  Int_t writeok = offline->ComputeOfflineParams(filePhysName.Data(), timecdb, cfdvalue, badpmt);
  printf(" AliT0PreprocessorOffline::CalibOffsetChannels :: writeok %i \n", writeok);
  AliCDBMetaData metaData;
  metaData.SetBeamPeriod(1);
  metaData.SetResponsible("Alla Maevskaya");
  metaData.SetComment("Time equalizing result with slew");
  fStatusDelay = writeok;
  if (writeok==0)  {
    AliCDBId* id1=NULL;
    id1=new AliCDBId("T0/Calib/TimeDelay", ustartRun, uendRun );
    ocdbStorage->Put(offline, (*id1), &metaData);
  }
  else {
    
    AliWarning(Form("writeok = %d data is not OK to be in OCDB",writeok));
      }		  
  
  delete offline;
 
}
//-------------------------------------------------------------------------------------
void AliT0PreprocessorOffline::CalibT0sPosition(TString filePhysName, Int_t ustartRun, Int_t uendRun, AliCDBStorage* pocdbStorage)
{
  printf(" AliT0PreprocessorOffline::CalibT0sPosition \n");
  Float_t zero_timecdb[4]={0};
  Float_t *timecdb = zero_timecdb;
  if (pocdbStorage) ocdbStorage=pocdbStorage;
  else {
    TString localStorage = "local://"+gSystem->GetFromPipe("pwd")+"/OCDB";
    ocdbStorage=AliCDBManager::Instance()->GetStorage(localStorage.Data());
  }
  
  AliT0CalibSeasonTimeShift *offline = new AliT0CalibSeasonTimeShift();
  Int_t writeok = offline->SetT0Par(filePhysName.Data(), timecdb);
  printf(" AliT0PreprocessorOffline::CalibT0sPosition :: writeok %i \n", writeok);
  AliCDBMetaData metaData;
  metaData.SetBeamPeriod(1);
  metaData.SetResponsible("Alla Maevskaya");
  metaData.SetComment("Time equalizing result with slew");
  fStatusAdjust = writeok;
  if (writeok == 0)  {
    AliCDBId* id1=NULL;
    id1=new AliCDBId("T0/Calib/TimeAdjust", ustartRun, uendRun);
    ocdbStorage->Put(offline, (*id1), &metaData);
  }

}
//_____________________________________________________________________________
Int_t AliT0PreprocessorOffline::GetStatus() const
{
  //
  // Checks the status
  // fStatusPos: errors
  // fStatusNeg: only info
  //
  printf("!!!! GetStatus %i delay %i adjust %i \n",fStatusAdjust+ fStatusDelay, fStatusDelay, fStatusAdjust);
  return   fStatusAdjust+ fStatusDelay;
}
