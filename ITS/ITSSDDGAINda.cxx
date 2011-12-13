/*
- Contact: - prino@to.infn.it
- Link: - alien:///alice/data/2009/LHC09c_SDD/000079095/raw/09000079095024.10.root
- Run Type: - PULSER_RUN
- DA Type: - LDC
- Number of events needed: >15
- Input Files: - SDDbase_step1_ddl*c*_sid*.data (output of previous PEDESTAL run)
- Output Files: - SDDbase_ddl*c*_sid*.data
- Trigger types used: 
*/

//////////////////////////////////////////////////////////////////////////////
// Detector Algorithm for analysis of SDD test pulse runs.                  //
//                                                                          //
// Produces ASCII and ROOT output files with:                               //
// 1. anode quality bit                                                     //
// 1. Baseline values                                                       //
// 2. Raw noise                                                             //
// 3. Common mode coefficients                                              //
// 4. Noise corrected for common mode                                       //
// 5. Gain                                                                  //
// Files are written to FXS                                                 //
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
#include <TH1F.h>
#include <TH2F.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TPluginManager.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TDatime.h>

// AliRoot includes
#include "AliRawReaderDate.h"
#include "AliITSOnlineSDDTP.h"
#include "AliITSRawStreamSDD.h"
#include "AliITSRawStreamSDDCompressed.h"

#ifdef ALI_AMORE
#include <AmoreDA.h>
#endif
/* Main routine
      Arguments: list of DATE raw data files
*/
int main(int argc, char **argv) {
  // main -  Arguments: list of DATE raw data files
  int status = 0;

  // line added to solve IO problems
  gROOT->GetPluginManager()->AddHandler("TVirtualStreamerInfo",
                                        "*",
                                        "TStreamerInfo",
                                        "RIO",
                                        "TStreamerInfo()");

  /* log start of process */
  printf("ITS SDD TEST-PULSE algorithm program started\n");  


  /* check that we got some arguments = list of files */
  if (argc<2) {
    printf("Wrong number of arguments\n");
    return -1;
  }


  Int_t maxNEvents=15; // maximum number of events to be analyzed
  const Int_t kTotDDL=24;
  const Int_t kModPerDDL=12;
  const Int_t kSides=2;
  UInt_t amSamplFreq=40;
  UChar_t cdhAttr=0;

  gSystem->Exec("rm -f  SDDbase_LDC.tar");

  AliITSOnlineSDDTP **tpan=new AliITSOnlineSDDTP*[kTotDDL*kModPerDDL*kSides];
  TH2F **histo=new TH2F*[kTotDDL*kModPerDDL*kSides];
  Bool_t isFilled[kTotDDL*kModPerDDL*kSides];
  Bool_t writtenoutput=kFALSE;
  Char_t hisnam[20];
  for(Int_t iddl=0; iddl<kTotDDL;iddl++){
    for(Int_t imod=0; imod<kModPerDDL;imod++){
      for(Int_t isid=0;isid<kSides;isid++){
	Int_t index=kSides*(kModPerDDL*iddl+imod)+isid;
	tpan[index]=new AliITSOnlineSDDTP(iddl,imod,isid,100.);
	sprintf(hisnam,"h%02dc%02ds%d",iddl,imod,isid);
	histo[index]=new TH2F(hisnam,"",256,-0.5,255.5,256,-0.5,255.5);
	isFilled[index]=0;
      }
    }
  }
  
  /* report progress */
  daqDA_progressReport(10);

  Int_t iev=0;
  Int_t ievPul=0;
  Int_t ievUsed=0;
  Int_t nEvToBeSkipped=5;

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
    daqDA_progressReport(10+80*n/argc);

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
	ievPul++;
	if(ievPul<=nEvToBeSkipped){
	  printf(" -> SKIP\n");
	  break;
	}
	printf("  -> Analyze\n");
	ievUsed++;

	AliRawReader *rawReader = new AliRawReaderDate((void*)event);
	rawReader->Reset();
	cdhAttr=AliITSRawStreamSDD::ReadBlockAttributes(rawReader);
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
	      if(amSamplFreq==20) tpan[index]->SetLastGoodTB(126);
	      else tpan[index]->SetLastGoodTB(254);
	      if(isFilled[index]) tpan[index]->AddEvent(histo[index]);    
	    }
	  }
	}
	
	/* free resources */
	free(event);
      }
      if(ievUsed>=maxNEvents) break;
    }    
  }
    
  /* write report */
  TDatime time;
  TObjString timeinfo(Form("%02d%02d%02d%02d%02d%02d",time.GetYear()-2000,time.GetMonth(),time.GetDay(),time.GetHour(),time.GetMinute(),time.GetSecond()));
  printf("Run #%s, received %d calibration events, time %s\n",getenv("DATE_RUN_NUMBER"),ievUsed,timeinfo.GetString().Data());

  /* report progress */
  daqDA_progressReport(90);

  TObjArray* gainHistos=new TObjArray();
  TObjArray* statusHistos=new TObjArray();
  
  Char_t filnam[100],command[120];
  for(Int_t iddl=0; iddl<kTotDDL;iddl++){
    for(Int_t imod=0; imod<kModPerDDL;imod++){
      for(Int_t isid=0;isid<kSides;isid++){
	Int_t index=kSides*(kModPerDDL*iddl+imod)+isid;
	if(isFilled[index]){
	  tpan[index]->ValidateAnodes();
	  tpan[index]->WriteToASCII();
	  gainHistos->AddLast(tpan[index]->GetGainAnodeHisto());
	  statusHistos->AddLast(tpan[index]->GetStatusAnodeHisto());
	  sprintf(filnam,"SDDbase_ddl%02dc%02d_sid%d.data",iddl,imod,isid);
	  sprintf(command,"tar -rf SDDbase_LDC.tar %s",filnam);
	  gSystem->Exec(command);
	}
      }
    }  
  }


  FILE *conffil=fopen("fee.conf","w");
  fprintf(conffil,"%d\n",amSamplFreq);
  fprintf(conffil,"%02X\n",cdhAttr);
  fclose(conffil);
  gSystem->Exec("tar -rf SDDbase_LDC.tar fee.conf");
  status=daqDA_FES_storeFile("./SDDbase_LDC.tar","SDD_Calib");


#ifdef ALI_AMORE
  amore::da::AmoreDA amoreDA(amore::da::AmoreDA::kSender);
  Int_t statusamore =0;
  statusamore += amoreDA.Send("TimeInfoPulser",&timeinfo);
  statusamore += amoreDA.Send("Gain",gainHistos);
  if ( statusamore )
    printf("Warning: Failed to write Arrays in the AMORE database\n");
  else 
    printf("amoreDA.Send() OK\n");
#else
  printf("Warning: SDDGAIN DA not compiled with AMORE support\n");
#endif
    
 
  TFile *fh=new TFile("SDDgainHistos.root","RECREATE");
  gainHistos->Write();
  statusHistos->Write();
  fh->Close();


  /* report progress */
  daqDA_progressReport(100);



  return status;
}
