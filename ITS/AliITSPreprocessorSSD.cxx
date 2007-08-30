#include "AliITSPreprocessorSSD.h"

#include "AliCDBMetaData.h"
#include "AliLog.h"
#include "TFile.h"

#include <TTimeStamp.h>
#include <TObjString.h>

#include "AliITSRawStreamSSD.h"
#include "AliITSNoiseSSD.h"
#include <Riostream.h>


//
// Author: Enrico Fragiacomo
// Date: 13/10/2006
// 
// SHUTTLE preprocessing class for SSD calibration files

/* $Id$ */

ClassImp(AliITSPreprocessorSSD)

//______________________________________________________________________________________________
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

  TObjArray calib_array(1698); 
  //Float_t noise=0, gain=0;

  TString runType = GetRunType();
  if(runType == "ELECTRONICS_CALIBRATION_RUN") {

  }
  else if(runType == "PEDESTAL_RUN") {

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
	    
	    TObjArray *cal; 
	    f->GetObject("TObjArray;1", cal); 
	    f->Close(); delete f;
	    
	    if(!cal) {
	    	Log("File does not contain expected data!");
		delete list;
	    }
	    
	    Int_t nmod = cal->GetEntries();

	    for(Int_t mod=0; mod<nmod; mod++) {

	      AliITSNoiseSSD *calib = (AliITSNoiseSSD*) cal->At(mod);
	      calib_array.AddAt(calib,calib->GetMod());
	    }
		
	  } else {
	  	Log("GetFile error!");
		delete list;
		return 3;
	  } // if filename
	} // end iteration over LDCs
	
	delete list;
      } else {
      	  Log("GetFileSources error!");
	  if(list) delete list;
	  return 4;
      } // if list
    
      //Now we have to store the final CDB file
      AliCDBMetaData metaData;
      metaData.SetBeamPeriod(0);
      metaData.SetResponsible("Enrico Fragiacomo");
      metaData.SetComment("This preprocessor fills the TObjArray of AliITSNoiseSSD");
  
     if(!Store("Calib", "NoiseSSD", &calib_array, &metaData, 0, 1)) {
          Log("no store");
        return 1;
     }  

  } // end if noise_run
  else {
    Log("Nothing to do");
    return 0;
  }
  
  Log("Database updated");
  return 0; // 0 means success

}

