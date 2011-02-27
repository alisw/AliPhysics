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
#include "AliT0CalibTimeEq.h"
#include "AliCDBStorage.h"
#include "AliCDBMetaData.h"
#include "AliCDBManager.h"

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
TNamed("AliT0PreprocessorOffline","AliT0PreprocessorOffline")
{
  //constructor
}
//____________________________________________________

AliT0PreprocessorOffline::~AliT0PreprocessorOffline()
{
  //destructor

}
//____________________________________________________
//____________________________________________________

void AliT0PreprocessorOffline::CalibOffsetChannels(TString filePhysName, Int_t ustartRun, Int_t uendRun, TString ocdbStorage)
{


  //Processing data from DAQ Physics run
  AliInfo("Processing Time Offset between channels");
  if (filePhysName)
    {
      AliT0CalibTimeEq *offline = new AliT0CalibTimeEq();
      offline->Reset();
      Bool_t writeok = offline->ComputeOfflineParams(filePhysName.Data());
      AliCDBMetaData metaData;
      metaData.SetBeamPeriod(1);
      metaData.SetResponsible("Alla Maevskaya");
      metaData.SetComment("Time equalizing result with slew");
      
      if (writeok)  {
	AliCDBId* id1=NULL;
	id1=new AliCDBId("T0/Calib/TimeDelay", ustartRun, uendRun);
	AliCDBStorage* gStorage = AliCDBManager::Instance()->GetStorage(ocdbStorage);
	gStorage->Put(offline, (*id1), &metaData);
      }
      else {
	
	   AliWarning(Form("writeok = %d not enough data for equalizing",writeok));
      }		  
   
      delete offline;
    }
	

}
