/*
Contact: annalisa.mastroserio@cern.ch
Link: tydes.home.cern.ch/tydes/doc/CalibrationOverview/CalibrationAlgorithms/
Run Type: DAQ_FO_UNIF_SCAN
DA Type: LDC
Number of events needed: Depending on scan type
Input Files: spd_focalib_params, raw data
Output Files: ./calibResults/ScanDCSconfigToFXS/* 
Trigger types used: PHYSICS
*/

////////////////////////////////////////////////////////////////////////////////
// This program can be compiled in two modes.                                 //
//                                                                            //
// 1. With the DAQ DA framework on. This is the default operating mode.       //
// Call this program with the name of the executable followed by the          //
// data files to process.                                                     //
//                                                                            //
// 2. Without the DAQ DA framework on. Define the SPD_DA_OFF environment var. //
// Call this program with the name of the executable followed by the          //
// runNr and the data files to process.                                       //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef SPD_DA_OFF
extern "C" {
#include <daqDA.h>
}
#endif
#include "event.h"
#include "monitor.h"
#include "AliRawReaderDate.h"
#include "AliITSRawStreamSPD.h"
#include "AliITSOnlineSPDfoChip.h"
#include "AliITSOnlineSPDfoInfo.h"
#include "AliITSOnlineSPDfo.h"
#include "AliITSOnlineSPDfoAnalyzer.h"
#include "AliLog.h"
#include <Riostream.h>
#include <fstream>
#include <TROOT.h>
#include <TStopwatch.h>
#include <TPluginManager.h>
#include <TObjArray.h>
#include <TClonesArray.h>
#include <TObjString.h>
#include <TString.h>
#include <TFitter.h>
#include <TFile.h>

int main(int argc, char **argv) {
  if (argc<2) {
    printf("Wrong number of arguments\n");
    return -1;
  }
  
  
  TStopwatch timer;
  timer.Start();
  
  // directory structure, hard coded
  char *saveDirDCSconfigToFXS= "./calibResults/ScanDCSconfigToFXS"; //     may delete content
  char *configFilesDir       = "./configFiles";                     //     may delete content
  char *saveDirIdsToFXS      = "./calibResults/IdsToFXS"; 

  // make sure the directory structure is correct:
  system("mkdir ./calibResults >& /dev/null");
  system("mkdir ./calibResults/ScanDCSconfigToFXS >& /dev/null");
  system("mkdir ./calibResults/IdsToFXS >& /dev/null");
  system("mkdir ./configFiles >& /dev/null");


  // parameters config files
  TString thresholdsFileName = Form("%s/focalib_params.txt",configFilesDir); 
  
  TFitter *fitter = new TFitter(3);
  TVirtualFitter::SetFitter(fitter);
  
  // This line is needed in case of a stand-alone application w/o
  // $ROOTSYS/etc/system.rootrc file
  gROOT->GetPluginManager()->AddHandler("TVirtualStreamerInfo",
					"*",
					"TStreamerInfo",
					"RIO",
					"TStreamerInfo()");
  
  // turn off annoying warning messages
  new AliLog;
  AliLog::Instance()->SetGlobalDebugLevel(-20);
  
  
  
  
// ********* STEP 0: Get configuration files from db (if there are any) , then read parameters*********  
  
  //chip efficiency selection parameters (needed afterwards for the data analysis)
  Int_t status = 0;
#ifndef SPD_DA_OFF
  TString idp = "spd_focalib_params";
  status=daqDA_DB_getFile(idp.Data(),thresholdsFileName.Data());
  if (status) {
    printf("Failed to get config file %s: status=%d. Using default tuning parameters.\n",idp.Data(),status);
    TString rmCmd = Form("rm -f %s",thresholdsFileName.Data());
    system(rmCmd.Data());
  }
#endif
   
  
  
  // ********* STEP 1: Produce FO scan container files (Reference Data). ***********************************
  int startSeg = 1;
  
#ifndef SPD_DA_OFF
  if (getenv("DATE_RUN_NUMBER")==0) {
    printf("DATE_RUN_NUMBER not properly set.\n");
    return -1;
  }
  int runNr = atoi(getenv("DATE_RUN_NUMBER"));
#else
  int runNr = atoi(argv[1]);
  startSeg = 2;
#endif
  
  Int_t evType =-1;
  AliITSOnlineSPDfoInfo *info[20]; Int_t ntriggers[20]; Int_t vDB[20]; Bool_t iseq[20];
  AliITSOnlineSPDfo *fomanager[20];
  TString s = "focalib";
  
  for(Int_t equip =0; equip < 20; equip++ ) {        
    info[equip] = new AliITSOnlineSPDfoInfo();
    info[equip]->SetRunNumber(runNr);
    info[equip]->SetRouter(equip);
    ntriggers[equip] = 0;
    vDB[equip] =0;  
    iseq[equip]=kFALSE;
    fomanager[equip] = new AliITSOnlineSPDfo(s,runNr,equip);
  }
  
  
  // loop over run segments
  for (int segNr=startSeg; segNr<argc; segNr++) {
    
    int status;
    
    // define data source : this is argument 1   
    status=monitorSetDataSource( argv[segNr] );
    if (status!=0) {
      printf("monitorSetDataSource() failed : %s\n",monitorDecodeError(status));
      return -1;
    }
    // declare monitoring program 
    status=monitorDeclareMp("ITS_SPD_CAL");
    if (status!=0) {
      printf("monitorDeclareMp() failed : %s\n",monitorDecodeError(status));
      return -1;
    }
    // define wait event timeout - 1s max 
    monitorSetNowait();
    monitorSetNoWaitNetworkTimeout(1000);
    
    Int_t eventType;
    UInt_t eventNr=0;
      
    // main loop (infinite) 
    for(;;) {
      
      struct eventHeaderStruct *event;
      eventTypeType eventT;
      
      /* check shutdown condition */
#ifndef SPD_DA_OFF
      if (daqDA_checkShutdown()) {break;}
#endif
      /* get next event (blocking call until timeout) */
      status=monitorGetEventDynamic((void **)&event);
      if (status==MON_ERR_EOF) {
	printf ("End of File detected\n");
	break; /* end of monitoring file has been reached */
      }
      if (status!=0) {
	printf("monitorGetEventDynamic() failed : %s\n",monitorDecodeError(status));
	break;
      }
      /* retry if got no event */
      if (event==NULL) {
	continue;
      }
      
      eventT=event->eventType;
      eventType = (Int_t) eventT;
      
      if (eventT == PHYSICS_EVENT) {
	
	eventNr++;
	//if(eventNr%5000 == 0 )printf(" eventNr %d\n",eventNr);
	
	AliRawReader *reader = new AliRawReaderDate((void*)event);
	AliITSRawStreamSPD *str = new AliITSRawStreamSPD(reader);
        
	for (UInt_t eqId=0; eqId<20; eqId++) {
	  
	  reader->Reset();
	  reader->Select("ITSSPD",eqId,eqId);
          
	  if (str->ReadCalibHeader()>0) {
        
            if(evType<0) evType = str->GetFOHtype();
               
            if(!iseq[eqId]){ // create output files           
             fomanager[eqId]->CreateOutputFile();
             fomanager[eqId]->SetFOscanParams(info[eqId]);
             iseq[eqId]=kTRUE;
            }
                     
            if(info[eqId]->GetNumDACindex()<1) {
              Int_t ind =0;
              while(str->GetFOHdacIndex(ind)>0) {
                info[eqId]->AddDACindex(str->GetFOHdacIndex(ind));
                ind++;
              }
            }
            	    
            if(!ntriggers[eqId]) {
	      ntriggers[eqId] = str->GetFOHtriggers();
	      info[eqId]->SetNumTriggers(str->GetFOHtriggers()); 
            }
	    if(!vDB[eqId])      {
	      vDB[eqId]       = str->GetFOHglobalDBversion();
	      info[eqId]->SetDBversion(str->GetFOHglobalDBversion());
	    }	    
	   
            if(!fomanager[eqId]->GetNdacs()) fomanager[eqId]->SetNdacs(str->GetFOHnumDacs());
           
            TArrayS dacvalues(str->GetFOHnumDacs());
            for(Int_t n = 0; n<(Int_t)str->GetFOHnumDacs(); n++) dacvalues.AddAt(str->GetFOHdacValue(n),n);
            
            TArrayS dacs = fomanager[eqId]->CreateDACArray(dacvalues, info[eqId]->GetDACIndexArray());
	    
            for(Int_t ihs =0; ihs < 6; ihs++) { // needed in the header to access the HS and ChipId info (in data it is different)
	      for(Int_t ich =0; ich < 10; ich++){
		if(!str->GetFOHchipPresent(ihs, ich)) continue;
                  info[eqId]->SetActiveChipsAndHS(ihs,ich);
		  Short_t measure[4] = {str->GetFOHMatrixID(),str->GetFOHpixelRow(), str->GetFOHpixelCol(), str->GetFOHchipCount(ihs,ich)}; 
                  fomanager[eqId]->AddMeasurement(dacs,measure,ihs,ich);                                            
	      } // chip loop      
            }// HS loop
	  }// if str->ReadHeader()>0;          
	}//if eqId      
	delete str;
	delete reader;
      }
      
      // free resources 
      free(event);   
    }// infinite loop
    
    
 
    
#ifndef SPD_DA_OFF
    daqDA_progressReport((unsigned int)( ((Float_t)(segNr-startSeg+1))/(argc-startSeg)*50 ));
#else
    printf("progress: %d\n",(unsigned int)( ((Float_t)(segNr-startSeg+1))/(argc-startSeg)*50 ));
#endif
    
  }// loop over run segments
  
   
     TString id[20], files[20];
    for(Int_t ifile =0; ifile < 20; ifile++) {
      if(iseq[ifile]){
      id[ifile] = Form("SPD_ref_fo%02i",ifile);
      files[ifile] = fomanager[ifile]->GetFile()->GetName();
      fomanager[ifile]->WriteToFile();
    }
      delete fomanager[ifile];
    }
    
      
    // ANALYSIS part
    
    for(Int_t iff =0; iff<20 ; iff++){
     if(!iseq[iff]) continue;
          
     AliITSOnlineSPDfoAnalyzer * analyzer = new AliITSOnlineSPDfoAnalyzer(Form("%i_%s%02i.root",runNr,s.Data(),iff));      
     analyzer->ReadParamsFromLocation(configFilesDir);
     analyzer->Process();
     
     TString dcsConfigFileName = Form("%s/dcsConfig_run_%d_eq_%d.txt",saveDirDCSconfigToFXS,runNr,iff);
     ofstream dcsfile;
     dcsfile.open(dcsConfigFileName.Data());
     dcsfile << "[SPD SCAN]\n"; 
     dcsfile << "RunNumber=" << runNr << "\n";
     dcsfile << "Type="<< evType <<"\n";
     dcsfile << "Router=" << iff << "\n";
     dcsfile << "ActualDetConfiguration=" << vDB[iff]<<"\n\n";
     dcsfile << "[DACvalues]\n"; 
     
       for(Int_t hs =0; hs<6; hs++){
	 for(Int_t ichip =0; ichip < 10; ichip++){
	   TArrayI dacs = analyzer->ChooseDACValues(hs,ichip);
	   
	   if(dacs.GetSize() == 0) continue;
	   for(Int_t idac =0; idac < dacs.GetSize() - 1; idac++) { // -1 (the last one is the quality flag)
	     if(dacs.At(idac) >=0 ) {
	     
	       dcsfile << ((analyzer->GetFOHandler())->GetFOscanInfo())->GetDACindex(idac) << ",";
	       dcsfile << iff << ",";
	       dcsfile << hs << ",";
	       dcsfile << ichip << "=" ;
	       dcsfile << dacs.At(idac) << ",";
	       dcsfile << dacs.At(dacs.GetSize() - 1) << "\n";
	     }
	   }
	 }
       }
       dcsfile.close();
    }
    
    
    printf("Preparing DCS config files\n");
    // send a tared file of all the dcsConfig text files
    TString command = Form("cd %s; tar -cf dcsConfig.tar *",saveDirDCSconfigToFXS);
    //printf("\n\n%s\n\n",command.Data());
    system(command.Data());
    TString fileName = Form("%s/dcsConfig.tar",saveDirDCSconfigToFXS);
    TString iddcs = "SPD_dcsConfig";
    
#ifndef SPD_DA_OFF
    status = daqDA_FES_storeFile(fileName.Data(),iddcs.Data());
    if (status!=0) {
      printf("Failed to export file %s , status %d\n",fileName.Data(),status);
      return -1;
    }
#endif
    
    
    printf("Opening id list file\n");
    TString idsFXSFileName = Form("%s/FXSids_run_%d.txt",saveDirIdsToFXS,runNr);
    ofstream idsFXSfile;
    idsFXSfile.open(idsFXSFileName.Data());
    
    // send reference data to FXS
    for (UInt_t eqId=0; eqId<20; eqId++) {
      if(!iseq[eqId]) continue;   
      //printf("Preparing reference data for eq %d\n",eqId);
      
      TString idf = Form("SPD_fo_scan_%d",eqId);
#ifndef SPD_DA_OFF
      status = daqDA_FES_storeFile(files[eqId].Data(),idf.Data());
      if (status!=0) {
	printf("Failed to export file %s , status %d\n",files[eqId].Data(),status);
	return -1;
      }
#endif
      idsFXSfile << Form("%s\n",idf.Data());   
    }
    
    timer.Stop();
    timer.Print();   
    printf("DA finished.\n");
    return 0;  
}
