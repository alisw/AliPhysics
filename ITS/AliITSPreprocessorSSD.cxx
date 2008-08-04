#include "AliITSPreprocessorSSD.h"
 
#include "AliCDBMetaData.h"
#include "AliLog.h"
#include "TFile.h"

#include <TTimeStamp.h>
#include <TObjString.h>

#include "AliITSRawStreamSSD.h"
#include "AliITSNoiseSSDv2.h"
#include "AliITSPedestalSSDv2.h"
#include "AliITSBadChannelsSSDv2.h"
#include <Riostream.h>


//
// Author: Enrico Fragiacomo
// Date: 13/10/2006
// 
// SHUTTLE preprocessing class for SSD calibration files

/* $Id$ */

const Int_t AliITSPreprocessorSSD::fgkNumberOfSSD = 1698;

ClassImp(AliITSPreprocessorSSD)

//-----------------------------------------------------------------------
AliITSPreprocessorSSD::AliITSPreprocessorSSD(AliShuttleInterface* shuttle) :
  AliPreprocessor("SSD", shuttle)
{
  // constructor

  AddRunType("ELECTRONICS_CALIBRATION_RUN");
  AddRunType("PEDESTAL");
  AddRunType("PHYSICS");

}

///______________________________________________________________________________________________
void AliITSPreprocessorSSD::Initialize(Int_t run, UInt_t startTime,
	UInt_t endTime)
{
 
  AliPreprocessor::Initialize(run, startTime, endTime);
  
  Log(Form("\n\tRun %d \n\tStartTime %s \n\tEndTime %s", run,
	       TTimeStamp(startTime).AsString(),
	       TTimeStamp(endTime).AsString()));
  
}

//______________________________________________________________________________________________
UInt_t AliITSPreprocessorSSD::Process(TMap* /*dcsAliasMap*/)
{

  // Note. To be modified: dcsAliasMap is not needed but I can not get rid
  // of it unless the base class AliPreprocessor is modified accordingly.

  //  TObjArray calib_array(fgkNumberOfSSD); 
  //TObjArray badch_array(fgkNumberOfSSD); 
  //TObjArray ped_array(fgkNumberOfSSD); 
  //Float_t noise=0, gain=0;
  
  //---------------------------------------
  // initialize the calibration objects
  AliITSNoiseSSDv2 *calib = new AliITSNoiseSSDv2();
  AliITSBadChannelsSSDv2 *badch = new AliITSBadChannelsSSDv2();
  AliITSPedestalSSDv2 *pedel = new AliITSPedestalSSDv2();
  
  TString runType = GetRunType();
  if(runType == "ELECTRONICS_CALIBRATION_RUN") {
    
  }
  else if(runType == "PEDESTAL") {

    TList* list = GetFileSources(kDAQ, "CALIBRATION");
    if (list && list->GetEntries() > 0)
      {
	Log("The following sources produced files with the id CALIBRATION");
	list->Print();
	
	// create iterator over list of LDCs (provides list of TObjString)
	TIter next(list);
	TObjString *ok;
	
	// expect to iterate 3 times (LDC0, LDC1, LDC2)
	while ( (ok = (TObjString*) next()) ) {                               
	  
	  TString key = ok->String();
	  
	  TString fileName = GetFile(kDAQ, "CALIBRATION", key.Data());
	  if (fileName.Length() > 0) {
	    
	    Log(Form("Got the file %s, now we can extract some values.", fileName.Data()));
	    
	    TFile *f = new TFile(fileName.Data());
	    if(!f || !f->IsOpen()){
	    	Log("Error opening file!");
		delete list;
		return 2;
	    }
	    
	    AliITSNoiseSSDv2 *cal; 
	    f->GetObject("AliITSNoiseSSDv2;1", cal); 
	    if(!cal) {
	    	Log("File does not contain expected data for the noise!");
		delete list;
		return 3;
	    }	    
	    AliITSPedestalSSDv2 *ped;
	    f->GetObject("AliITSPedestalSSDv2;1", ped); 
	    if(!ped) {
	    	Log("File does not contain expected data for the pedestals!");
		delete list;
		return 5;
	    }	    
	    AliITSBadChannelsSSDv2 *bad;
	    f->GetObject("AliITSBadChannelsSSDv2;1", bad); 
	    if(!bad) {
	    	Log("File does not contain expected data for bad channels  !");
		delete list;
		return 4;
	    }	    

	    for(Int_t module=0; module<fgkNumberOfSSD; module++) {
	      for(Int_t strip=0; strip<768; strip++) {
		if(cal->GetNoiseP(module,strip)) 
		  calib->AddNoiseP(module,strip,cal->GetNoiseP(module,strip));
		if(cal->GetNoiseN(module,strip)) 
		  calib->AddNoiseN(module,strip,cal->GetNoiseN(module,strip));
		if(ped->GetPedestalP(module,strip)) 
		  pedel->AddPedestalP(module,strip,
				      ped->GetPedestalP(module,strip));
		if(ped->GetPedestalN(module,strip)) 
		  pedel->AddPedestalN(module,strip,
				      ped->GetPedestalN(module,strip));
		if(bad->GetBadChannelP(module,strip)) 
		  badch->AddBadChannelP(module,strip,
					 bad->GetBadChannelP(module,strip));
		if(bad->GetBadChannelN(module,strip)) 
		  badch->AddBadChannelN(module,strip,
					 bad->GetBadChannelN(module,strip));
	      }
	    }

	    f->Close(); delete f;	    
		
	  } else {
	  	Log("GetFile error!");
		delete list;
		return 6;
	  } // if filename
	} // end iteration over LDCs
	
	delete list;
      } else {
      	  Log("GetFileSources error!");
	  if(list) delete list;
	  return 7;
      } // if list
    
    //Now we have to store the final CDB file
    AliCDBMetaData metaData;
    metaData.SetBeamPeriod(0);
    metaData.SetResponsible("Enrico Fragiacomo");
    metaData.SetComment("Fills noise, pedestal and bad channels TObjArray");
    
      if(!Store("Calib", "NoiseSSD", (TObject *)calib, &metaData, 0, 1)) {
	Log("no store");
        return 1;
      }  
      
      if(!Store("Calib", "BadChannelsSSD", (TObject*)badch, &metaData, 0, 1)) {
	Log("no store");
        return 1;
      }  
      
      if(!StoreReferenceData("Calib","PedestalSSD",  (TObject*)pedel, &metaData)) {
	Log("no store");
	return 1;
      }
	 
  } // end if pedestal run
  else {
    Log("Nothing to do");
    return 0;
  }
  
  Log("Database updated");
  return 0; // 0 means success

}

