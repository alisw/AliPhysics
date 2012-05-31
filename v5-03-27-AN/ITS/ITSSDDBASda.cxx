/*
- Contact: - prino@to.infn.it
- Link: - alien:///alice/data/2009/LHC09c_SDD/000079094/raw/09000079094024.10.root
- Run Type: - PEDESTAL_RUN
- DA Type: - LDC
- Number of events needed: > 10
- Input Files: - 
- Output Files: - SDDbase_step1_ddl*c*_sid*.data SDDbase_step2_ddl*c*_sid*.data
- Trigger types used: 
*/


//////////////////////////////////////////////////////////////////////////////
// Detector Algorithm for analysis of SDD baseline runs.                    //
//                                                                          //
// Produces ASCII and ROOT output files with:                               //
// 1. anode quality bit                                                     //
// 1. Baseline values                                                       //
// 2. Raw noise                                                             //
// 3. Common mode coefficients                                              //
// 4. Noise corrected for common mode                                       //
// Files are stored locally.                                                //
// The next DAQ-DA step on Test Pulse run writes files to FXS               //
//                                                                          //
// Author: F. Prino (prino@to.infn.it)                                      //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////



extern "C" {
#include "daqDA.h"
}

#include "event.h"
#include "monitor.h"

#include <stdio.h>
#include <stdlib.h>

// ROOT includes
#include <TFile.h>
#include <TSystem.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TROOT.h>
#include <TPluginManager.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TDatime.h>

// AliRoot includes
#include "AliRawReaderDate.h"
#include "AliITSOnlineSDDBase.h"
#include "AliITSOnlineSDDCMN.h"
#include "AliITSRawStreamSDD.h"
#include "AliITSRawStreamSDDCompressed.h"

#ifdef ALI_AMORE
#include <AmoreDA.h>
#endif

/* Main routine
      Arguments: list of DATE raw data files
*/
int main(int argc, char **argv) {
  // main - Arguments: list of DATE raw data files
  int status = 0;

  // line added to solve IO problems
  gROOT->GetPluginManager()->AddHandler("TVirtualStreamerInfo",
                                        "*",
                                        "TStreamerInfo",
                                        "RIO",
                                        "TStreamerInfo()");

  /* log start of process */
  printf("ITS SDD BASELINE+NOISE algorithm program started\n");  


  /* check that we got some arguments = list of files */
  if (argc<2) {
    printf("Wrong number of arguments\n");
    return -1;
  }


  Int_t maxNEvents=10; // maximum number of events to be analyzed
  const Int_t kTotDDL=24;
  const Int_t kModPerDDL=12;
  const Int_t kSides=2;
  UInt_t amSamplFreq=40;

  gSystem->Exec("rm -f SDDbase_*.data");
  gSystem->Exec("rm -f  SDDbase_step2_LDC.tar");
  
  AliITSOnlineSDDBase **base=new AliITSOnlineSDDBase*[kTotDDL*kModPerDDL*kSides];
  AliITSOnlineSDDCMN **corr=new AliITSOnlineSDDCMN*[kTotDDL*kModPerDDL*kSides];
  TH2F **histo=new TH2F*[kTotDDL*kModPerDDL*kSides];
  Bool_t isFilled[kTotDDL*kModPerDDL*kSides];
  Bool_t writtenoutput=kFALSE;

  Char_t hisnam[20];
  for(Int_t iddl=0; iddl<kTotDDL;iddl++){
    for(Int_t imod=0; imod<kModPerDDL;imod++){
      for(Int_t isid=0;isid<kSides;isid++){
	Int_t index=kSides*(kModPerDDL*iddl+imod)+isid;
	base[index]=new AliITSOnlineSDDBase(iddl,imod,isid);
	sprintf(hisnam,"h%02dc%02ds%d",iddl,imod,isid);
	histo[index]=new TH2F(hisnam,"",256,-0.5,255.5,256,-0.5,255.5);
	isFilled[index]=0;
      }
    }
  }
  
  /* report progress */
  daqDA_progressReport(8);

  Int_t iev;
  Int_t ievPed;
  Int_t ievUsed;
  Int_t nEvToBeSkipped=5;

  for(Int_t iStep=0;iStep<2;iStep++){
    printf("Start Analysis Step %d\n",iStep);
    iev=0;
    ievPed=0;
    ievUsed=0;
    if(iStep==1){
      for(Int_t iddl=0; iddl<kTotDDL;iddl++){
	for(Int_t imod=0; imod<kModPerDDL;imod++){
	  for(Int_t isid=0;isid<kSides;isid++){
	    Int_t index=kSides*(kModPerDDL*iddl+imod)+isid;
	    corr[index]=new AliITSOnlineSDDCMN(iddl,imod,isid);
	    isFilled[index]=0;
	  }
	}
      }
    }

    /* read the data files */
    int n;
    for (n=1;n<argc;n++) {
   
      status=monitorSetDataSource( argv[n] );
      if (status!=0) {
	printf("monitorSetDataSource() failed : %s\n",monitorDecodeError(status));
	return -1;
      }

      /* report progress */
      /* in this example, indexed on the number of files */
      // Progress report inside the event loop as well?
      daqDA_progressReport(10+40*iStep*n/argc);

      /* read the file */
      for(;;) {
	struct eventHeaderStruct *event;
	eventTypeType eventT;

	/* get next event */
	status=monitorGetEventDynamic((void **)&event);
	if (status==MON_ERR_EOF) break; /* end of monitoring file has been reached */
	if (status!=0) {
	  printf("monitorGetEventDynamic() failed : %s\n",monitorDecodeError(status));
	  return -1;
	}

	/* retry if got no event */
	if (event==NULL) {
	  break;
	}

	iev++; 

	/* use event - here, just write event id to result file */
	eventT=event->eventType;
	switch (event->eventType){
      
	  /* START OF RUN */
	case START_OF_RUN:
	  break;
      
	  /* END OF RUN */
	case END_OF_RUN:
	  break;

 	  //      case PHYSICS_EVENT:  // comment this line for test raw data
	  //	break;               // comment this line for test raw data


	case CALIBRATION_EVENT:
	  break;  // uncomment this line for test raw data
	case PHYSICS_EVENT: // uncomment this line for test raw data
	  printf(" Event number = %i ",iev);
	  ievPed++;
	  if(ievPed<=nEvToBeSkipped){
	    printf(" -> SKIP\n");
	    break;
	  }
	  printf("  -> Analyze\n");
	  ievUsed++;
	  AliRawReader *rawReader = new AliRawReaderDate((void*)event);
	  rawReader->Reset();
	  UChar_t cdhAttr=AliITSRawStreamSDD::ReadBlockAttributes(rawReader);
	  amSamplFreq=AliITSRawStreamSDD::ReadAMSamplFreqFromCDH(cdhAttr);
	  AliITSRawStream* s=AliITSRawStreamSDD::CreateRawStreamSDD(rawReader,cdhAttr);
	  if(!writtenoutput){
	    printf("Use %s raw stream, sampling frequency %d MHz\n",s->ClassName(),amSamplFreq);
	    writtenoutput=kTRUE;
	  }


	  for(Int_t iddl=0; iddl<kTotDDL;iddl++){
	    for(Int_t imod=0; imod<kModPerDDL;imod++){
	      for(Int_t isid=0;isid<kSides;isid++){
		Int_t index=kSides*(kModPerDDL*iddl+imod)+isid;
		histo[index]->Reset();
	      }
	    }
	  }

	  while(s->Next()){
	    Int_t iDDL=rawReader->GetDDLID();
	    Int_t iCarlos=s->GetCarlosId();
	    if(s->IsCompletedModule()) continue;
	    if(s->IsCompletedDDL()) continue;
	    if(iDDL>=0 && iDDL<kTotDDL){ 
	      Int_t index=kSides*(kModPerDDL*iDDL+iCarlos)+s->GetChannel(); 
	      histo[index]->Fill(s->GetCoord2(),s->GetCoord1(),s->GetSignal());
	      isFilled[index]=1;
	    }
	  }
	  delete s;
	  delete rawReader;
	  for(Int_t iddl=0; iddl<kTotDDL;iddl++){
	    for(Int_t imod=0; imod<kModPerDDL;imod++){
	      for(Int_t isid=0;isid<kSides;isid++){
		Int_t index=kSides*(kModPerDDL*iddl+imod)+isid;
		if(iStep==0){
		  if(amSamplFreq==20) base[index]->SetLastGoodTB(126);
		  else base[index]->SetLastGoodTB(254);
		  base[index]->AddEvent(histo[index]);
		}
		if(iStep==1){
		  if(amSamplFreq==20) corr[index]->SetLastGoodTB(126);
		  else corr[index]->SetLastGoodTB(254);
		  corr[index]->AddEvent(histo[index]);
		}
	      }
	    }
	  }

	  /* free resources */
	  free(event);
	}
	if(ievUsed>=maxNEvents) break;
      }
    }

    if(iStep==0){
      for(Int_t iddl=0; iddl<kTotDDL;iddl++){
	for(Int_t imod=0; imod<kModPerDDL;imod++){
	  for(Int_t isid=0;isid<kSides;isid++){
	    Int_t index=kSides*(kModPerDDL*iddl+imod)+isid;
	    base[index]->ValidateAnodes();
	    base[index]->WriteToASCII();
	  }
	}
      }
    }
  }

  /* write report */
  TDatime time;
  TObjString timeinfo(Form("%02d%02d%02d%02d%02d%02d",time.GetYear()-2000,time.GetMonth(),time.GetDay(),time.GetHour(),time.GetMinute(),time.GetSecond()));
  printf("Run #%s, received %d calibration events, time %s\n",getenv("DATE_RUN_NUMBER"),ievUsed,timeinfo.GetString().Data());

  /* report progress */
  daqDA_progressReport(90);

  TObjArray* basHistos=new TObjArray();
  TObjArray* noiseHistos=new TObjArray();
  TObjArray* cmnHistos=new TObjArray();
  TObjArray* corrnHistos=new TObjArray();
  TObjArray* statusHistos=new TObjArray();


  Char_t filnam[100],command[150];
  for(Int_t iddl=0; iddl<kTotDDL;iddl++){
    for(Int_t imod=0; imod<kModPerDDL;imod++){
      for(Int_t isid=0;isid<kSides;isid++){
	Int_t index=kSides*(kModPerDDL*iddl+imod)+isid;
	corr[index]->ValidateAnodes();
	corr[index]->WriteToASCII();
	if(isFilled[index]){
	  basHistos->AddLast(corr[index]->GetBaselineAnodeHisto());
	  noiseHistos->AddLast(corr[index]->GetRawNoiseAnodeHisto());
	  cmnHistos->AddLast(corr[index]->GetCMNCoefAnodeHisto());
	  corrnHistos->AddLast(corr[index]->GetCorrNoiseAnodeHisto());
	  statusHistos->AddLast(corr[index]->GetStatusAnodeHisto());
	  sprintf(filnam,"SDDbase_step2_ddl%02dc%02d_sid%d.data",iddl,imod,isid);
	  sprintf(command,"tar -rf SDDbase_step2_LDC.tar %s",filnam);
	  gSystem->Exec(command);
	}
      }
    }
  }

#ifdef ALI_AMORE
  amore::da::AmoreDA amoreDA(amore::da::AmoreDA::kSender);
  Int_t statusamore =0;
  statusamore += amoreDA.Send("TimeInfoPedestal",&timeinfo);
  statusamore += amoreDA.Send("Baselines",basHistos);
  statusamore += amoreDA.Send("RawNoise",noiseHistos);
  statusamore += amoreDA.Send("CommonMode",cmnHistos);
  statusamore += amoreDA.Send("CorrectedNoise",corrnHistos);
  statusamore += amoreDA.Send("NoisyChannels",statusHistos);
  if ( statusamore )
    printf("Warning: Failed to write Arrays in the AMORE database\n");
  else 
    printf("amoreDA.Send() OK\n");
#else
  printf("Warning: SDDBAS DA not compiled with AMORE support\n");
#endif
    
    
  TFile *fh=new TFile("SDDbaseHistos.root","RECREATE");
  basHistos->Write();
  noiseHistos->Write();
  cmnHistos->Write();
  corrnHistos->Write();
  statusHistos->Write();
  fh->Close();

  /* report progress */
  daqDA_progressReport(100);



  return status;
}
