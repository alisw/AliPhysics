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
  startRun(0),                         // start Run - used to make fast selection in THnSparse
  endRun(0),                           // end   Run - used to make fast selection in THnSparse
  startTime(0),                        // startTime - used to make fast selection in THnSparse
  endTime(0),                          // endTime   - used to make fast selection in THnSparse
  ocdbStorage("")                  // path to the OCDB storage
  
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

  Float_t *timecdb, *cfdvalue;
  for (Int_t i=0; i<24; i++) {  timecdb[i]=999; cfdvalue=0; }
  Int_t badpmt=0;
  //Processing data from DAQ Physics run
  AliInfo("Processing Time Offset between channels");

  if (filePhysName)
    {
   if (pocdbStorage.Length()>0) ocdbStorage=pocdbStorage;
   else
     ocdbStorage="local://"+gSystem->GetFromPipe("pwd")+"/OCDB";
   //      AliCDBManager* man = AliCDBManager::Instance();
   // man->SetDefaultStorage("raw://");
   //   man->SetDefaultStorage("ocdbStorage");
   //     man->SetRun(ustartRun);
   AliCDBEntry *entry = AliCDBManager::Instance()->Get("GRP/CTP/CTPtiming");
   if (!entry) AliFatal("CTP timing parameters are not found in OCDB !");
   AliCTPTimeParams *ctpParams = (AliCTPTimeParams*)entry->GetObject();
   Float_t l1Delay = (Float_t)ctpParams->GetDelayL1L0()*25.0;
   
   AliCDBEntry *entry1 = AliCDBManager::Instance()->Get("GRP/CTP/TimeAlign");
   if (!entry1) AliFatal("CTP time-alignment is not found in OCDB !");
   AliCTPTimeParams *ctpTimeAlign = (AliCTPTimeParams*)entry1->GetObject();
   l1Delay += ((Float_t)ctpTimeAlign->GetDelayL1L0()*25.0);
   
   AliCDBEntry *entry4 = AliCDBManager::Instance()->Get("GRP/Calib/LHCClockPhase");
   if (!entry4) AliFatal("LHC clock-phase shift is not found in OCDB !");
   AliLHCClockPhase *phase = (AliLHCClockPhase*)entry4->GetObject();
   Float_t fGRPdelays = l1Delay - phase->GetMeanPhase();
   AliCDBEntry *entryL = AliCDBManager::Instance()->Get("T0/Calib/Latency");
   AliT0CalibLatency *lat = (AliT0CalibLatency*)entryL->GetObject();
   Float_t fLatencyHPTDC = lat->GetLatencyHPTDC();
   Float_t fLatencyL1 = lat->GetLatencyL1();
   AliCDBEntry *entryCalib = AliCDBManager::Instance()->Get("T0/Calib/TimeDelay");
   if(!entryCalib) {
     AliWarning(Form("Cannot find any AliCDBEntry for [Calib, TimeDelay]!"));
   }
   else
     {
       AliT0CalibTimeEq *clb = (AliT0CalibTimeEq*)entryCalib->GetObject();
       timecdb = clb->GetTimeEq();
       cfdvalue = clb->GetCFDvalue();
       for (Int_t i=0; i<24; i++) {
	 if( cfdvalue[i] < 500 ) cfdvalue[i] =( 1000.*fLatencyHPTDC - 1000.*fLatencyL1 + 1000.*fGRPdelays)/24.4;
	 printf("Calc  mean CFD time %i %f \n",i,cfdvalue[i]);
       }
     }
      AliCDBEntry *entryCalibreco = AliCDBManager::Instance()->Get("T0/Calib/RecoParam");
      if(entryCalibreco) {
	AliT0RecoParam *rpr = (AliT0RecoParam*) entryCalibreco->GetObject();
  	badpmt = rpr->GetRefPoint();
	printf(" bad PMT %i \n", badpmt);
      }
      AliT0CalibTimeEq *offline = new AliT0CalibTimeEq();
      offline->Reset();
      Bool_t writeok = offline->ComputeOfflineParams(filePhysName.Data(), timecdb, cfdvalue, badpmt);
      AliCDBMetaData metaData;
      metaData.SetBeamPeriod(1);
      metaData.SetResponsible("Alla Maevskaya");
      metaData.SetComment("Time equalizing result with slew");
      
      if (writeok)  {
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
	

}
//-------------------------------------------------------------------------------------
void AliT0PreprocessorOffline::CalibT0sPosition(TString filePhysName, Int_t ustartRun, Int_t uendRun, TString pocdbStorage)
{
  if (filePhysName)
    {
      if (pocdbStorage.Length()>0) ocdbStorage=pocdbStorage;
      else
	ocdbStorage="local://"+gSystem->GetFromPipe("pwd")+"/OCDB";
      
      AliT0CalibSeasonTimeShift *offline = new AliT0CalibSeasonTimeShift();
      Bool_t writeok = offline->SetT0Par(filePhysName.Data());
      AliCDBMetaData metaData;
      metaData.SetBeamPeriod(1);
      metaData.SetResponsible("Alla Maevskaya");
      metaData.SetComment("Time equalizing result with slew");
      
      if (writeok)  {
	AliCDBId* id1=NULL;
	id1=new AliCDBId("T0/Calib/TimeAdjust", ustartRun, uendRun);
	AliCDBStorage* gStorage = AliCDBManager::Instance()->GetStorage(ocdbStorage);
	gStorage->Put(offline, (*id1), &metaData);
      }
    }
}
