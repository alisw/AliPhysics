/*

This program reads the DAQ data files passed as argument using the monitoring library.

It computes the average event size and populates local "./result.txt" file with the 
result.

The program reports about its processing progress.

Messages on stdout are exported to DAQ log system.

DA for ZDC standalone pedestal runs

Contact: Chiara.Oppedisano@to.infn.it
Link: 
Run Type: STANDALONE_LASER_RUN
DA Type: LDC
Number of events needed: no constraint (tipically ~10^3)
Input Files: 
Output Files: ZDCLaser.dat
Trigger Types Used: Standalone Trigger

*/
#define PEDDATA_FILE  "ZDCPedestal.dat"
#define MAPDATA_FILE  "ZDCChMapping.dat"
#define LASDATA_FILE  "ZDCLaserCalib.dat"

#include <stdio.h>
#include <stdlib.h>
#include <Riostream.h>

// DATE
#include <event.h>
#include <monitor.h>
#include <daqDA.h>

//ROOT
#include <TRandom.h>
#include <TH1F.h>
#include <TF1.h>
#include <TFile.h>
#include <TFitter.h>

//AliRoot
#include <AliRawReaderDate.h>
#include <AliRawEventHeaderBase.h>
#include <AliZDCRawStream.h>


/* Main routine
      Arguments: list of DATE raw data files
*/
int main(int argc, char **argv) {
  
  TFitter *minuitFit = new TFitter(4);
  TVirtualFitter::SetFitter(minuitFit);

  int status = 0;

  /* log start of process */
  printf("\nZDC LASER program started\n");  

  /* check that we got some arguments = list of files */
  if (argc<2) {
    printf("Wrong number of arguments\n");
    return -1;
  }

  // --- Histograms for LASER runs 
  //     20 signal channels + 2 reference PTMs
  //
  TH1F::AddDirectory(0);
  //
  TH1F *hPMRefChg = new TH1F("hPMRefChg","hPMRefChg", 100,0.,1000.);
  TH1F *hPMRefAhg = new TH1F("hPMRefAhg","hPMRefAhg", 100,0.,1000.);
  //
  TH1F *hPMRefClg = new TH1F("hPMRefClg","hPMRefClg", 100,0.,4000.);
  TH1F *hPMRefAlg = new TH1F("hPMRefAlg","hPMRefAlg", 100,0.,4000.);


  /* open result file */
  FILE *fp=NULL;
  fp=fopen("./result.txt","a");
  if (fp==NULL) {
    printf("Failed to open file\n");
    return -1;
  }
  
  FILE *mapFile4Shuttle;
        
  // *** To analyze LASER events you MUST have a pedestal data file!!!
  // *** -> check if a pedestal run has been analyzied
  int read = 0;
  read = daqDA_DB_getFile(PEDDATA_FILE,PEDDATA_FILE);
  if(read){
    printf("\t ERROR!!! ZDCPedestal.dat file NOT FOUND in DAQ db!!!\n");
    return -1;
  }
  else printf("\t ZDCPedestal.dat file retrieved from DAQ db\n");
  
  FILE *filePed = fopen(PEDDATA_FILE,"r");
  if (filePed==NULL) {
    printf("\t ERROR!!! Can't open ZDCPedestal.dat file!!!\n");
    return -1;
  }

  // 144 = 48 in-time + 48 out-of-time + 48 correlations
  Float_t readValues[2][144], MeanPed[44], MeanPedWidth[44], 
  	MeanPedOOT[44], MeanPedWidthOOT[44];
  // ***************************************************
  //   Unless we have a narrow correlation to fit we
  //	don't fit and store in-time vs. out-of-time
  //	histograms -> mean pedstal subtracted!!!!!!
  // ***************************************************
  //Float_t CorrCoeff0[44], CorrCoeff1[44];
  //
  for(int jj=0; jj<144; jj++){
    for(int ii=0; ii<2; ii++){
       fscanf(filePed,"%f",&readValues[ii][jj]);
    }
    if(jj<48){
      MeanPed[jj] = readValues[0][jj];
      MeanPedWidth[jj] = readValues[1][jj];
      //printf("\t MeanPed[%d] = %1.1f\n",jj, MeanPed[jj]);
    }
    else if(jj>48 && jj<96){
      MeanPedOOT[jj-48] = readValues[0][jj];
      MeanPedWidthOOT[jj-48] = readValues[1][jj];
    }
    /*else if(jj>144){
      CorrCoeff0[jj-96] = readValues[0][jj]; 
      CorrCoeff1[jj-96] = readValues[1][jj];;
    }
    */
  }

  /* report progress */
  daqDA_progressReport(10);


  /* init some counters */
  int nevents_physics=0;
  int nevents_total=0;

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

      // Initalize raw-data reading and decoding
      AliRawReader *reader = new AliRawReaderDate((void*)event);
      reader->Select("ZDC");
      // --- Reading event header
      //UInt_t evtype = reader->GetType();
      //printf("\n\t ZDCLASERda -> ev. type %d\n",evtype);
      //printf("\t ZDCLASERda -> run # %d\n",reader->GetRunNumber());
      //
      AliZDCRawStream *rawStreamZDC = new AliZDCRawStream(reader);
        

      /* use event - here, just write event id to result file */
      eventT=event->eventType;
      
      Int_t ich=0, adcMod[48], adcCh[48], sigCode[48], det[48], sec[48];
      if(eventT==START_OF_DATA){
	  	
	if(!rawStreamZDC->Next()) printf(" \t No raw data found!! \n");
        else{
	  while(rawStreamZDC->Next()){
            if(rawStreamZDC->IsChMapping()){
	      adcMod[ich] = rawStreamZDC->GetADCModFromMap(ich);
	      adcCh[ich] = rawStreamZDC->GetADCChFromMap(ich);
	      sigCode[ich] = rawStreamZDC->GetADCSignFromMap(ich);
	      det[ich] = rawStreamZDC->GetDetectorFromMap(ich);
	      sec[ich] = rawStreamZDC->GetTowerFromMap(ich);
	      ich++;
	    }
	  }
	}
	// --------------------------------------------------------
	// --- Writing ascii data file for the Shuttle preprocessor
        mapFile4Shuttle = fopen(MAPDATA_FILE,"w");
        for(Int_t i=0; i<ich; i++){
	   fprintf(mapFile4Shuttle,"\t%d\t%d\t%d\t%d\t%d\t%d\n",i,
	     adcMod[i],adcCh[i],sigCode[i],det[i],sec[i]);
	   //
	   //printf("ZDCLASERDA.cxx ->  ch.%d mod %d, ch %d, code %d det %d, sec %d\n",
	   //	   i,adcMod[i],adcCh[i],sigCode[i],det[i],sec[i]);
        }
        fclose(mapFile4Shuttle);
      }

      /* use event - here, just write event id to result file */
      eventT=event->eventType;
    
      if(eventT==PHYSICS_EVENT){
        //
	// --- Reading data header
        reader->ReadHeader();
        const AliRawDataHeader* header = reader->GetDataHeader();
        if(header) {
         UChar_t message = header->GetAttributes();
	 if(message & 0x20){ // DEDICATED LASER RUN
	    //printf("\t STANDALONE_LASER_RUN raw data found\n");
	    continue;
	 }
	 else{
	    printf("\t NO STANDALONE_LASER_RUN raw data found\n");
	    return -1;
	 }
  	}
  	else{
  	   printf("\t ATTENTION! No Raw Data Header found!!!\n");
  	   return -1;
  	}

  	if (!rawStreamZDC->Next()) printf(" \t No raw data found!! \n");
	//
	// ----- Setting ch. mapping -----
	for(Int_t jk=0; jk<48; jk++){
	  rawStreamZDC->SetMapADCMod(jk, adcMod[jk]);
	  rawStreamZDC->SetMapADCCh(jk, adcCh[jk]);
	  rawStreamZDC->SetMapADCSig(jk, sigCode[jk]);
	  rawStreamZDC->SetMapDet(jk, det[jk]);
	  rawStreamZDC->SetMapTow(jk, sec[jk]);
	}
  	//
	while(rawStreamZDC->Next()){
  	  Int_t index=-1;
	  // Getting data only for reference PMTs (sector[1]=5)
  	  if((rawStreamZDC->IsADCDataWord()) && (rawStreamZDC->GetSector(1)==5)){
	    index = rawStreamZDC->GetADCChannel();
	    Float_t Pedestal = MeanPed[index];
	    Float_t CorrADC = rawStreamZDC->GetADCValue() - Pedestal;
	    
	    // ==== HIGH GAIN CHAIN
	    if(rawStreamZDC->GetADCGain() == 0){
	      // %%%%% PMRef chain side C
	      if(rawStreamZDC->GetSector(0)==1) hPMRefChg->Fill(CorrADC);
	      // %%%%% PMRef side A
	      else if(rawStreamZDC->GetSector(0)==4) hPMRefAhg->Fill(CorrADC);
	    }
	    // ==== LOW GAIN CHAIN
	    else{
	      // %%%%% PMRef chain side C
	      if(rawStreamZDC->GetSector(0)==1) hPMRefClg->Fill(CorrADC);
	      // %%%%% PMRef side A
	      else if(rawStreamZDC->GetSector(0)==4) hPMRefAlg->Fill(CorrADC);
	    }
  	  }//IsADCDataWord()
  	  //
         }
         //
         nevents_physics++;
         //
	 delete reader;
         delete rawStreamZDC;

      }//(if PHYSICS_EVENT) 
      nevents_total++;

      /* free resources */
      free(event);
    
    }
  }  
  
  /* Analysis of the histograms */
  //
  Int_t maxBinRef[4], nBinRef[4];
  Float_t xMaxRef[4], maxXvalRef[4], xlowRef[4]; 
  Float_t meanRef[2], sigmaRef[2];
  TF1 *funRef[4];
  
  // ~~~~~~~~ PM Ref side C high gain chain ~~~~~~~~
  maxBinRef[0] = hPMRefChg->GetMaximumBin();
  nBinRef[0] = (hPMRefChg->GetXaxis())->GetNbins();
  xMaxRef[0] = (hPMRefChg->GetXaxis())->GetXmax();
  maxXvalRef[0] = maxBinRef[0]*xMaxRef[0]/nBinRef[0];
  // 
  if(maxXvalRef[0]-100.<0.) {xlowRef[0]=0.;}
  else xlowRef[0] = maxXvalRef[0];
  hPMRefChg->Fit("gaus","Q","",xlowRef[0],maxXvalRef[0]+100.);
  funRef[0] = hPMRefChg->GetFunction("gaus");
  meanRef[0] = (Float_t) (funRef[0]->GetParameter(1));
  sigmaRef[0] = (Float_t) (funRef[0]->GetParameter(2));
  
  // ~~~~~~~~ PM Ref side A high gain chain ~~~~~~~~
  maxBinRef[1] = hPMRefAhg->GetMaximumBin();
  nBinRef[1] = (hPMRefAhg->GetXaxis())->GetNbins();
  xMaxRef[1] = (hPMRefAhg->GetXaxis())->GetXmax();
  maxXvalRef[1] = maxBinRef[1]*xMaxRef[1]/nBinRef[1];
  //
  if(maxXvalRef[1]-100.<0.) {xlowRef[1]=0.;}
  else xlowRef[1] = maxXvalRef[1];
  hPMRefAhg->Fit("gaus","Q","",xlowRef[1],maxXvalRef[1]+100.);
  funRef[1] = hPMRefAhg->GetFunction("gaus");
  meanRef[1] = (Float_t) (funRef[1]->GetParameter(1));
  sigmaRef[1] = (Float_t) (funRef[1]->GetParameter(2));
  
  // ~~~~~~~~ PM Ref side C low gain chain ~~~~~~~~
  maxBinRef[2] = hPMRefClg->GetMaximumBin();
  nBinRef[2] = (hPMRefClg->GetXaxis())->GetNbins();
  xMaxRef[2] = (hPMRefClg->GetXaxis())->GetXmax();
  maxXvalRef[2] = maxBinRef[2]*xMaxRef[2]/nBinRef[2];
  //
  if(maxXvalRef[2]-100.<0.) {xlowRef[2]=0.;}
  else xlowRef[2] = maxXvalRef[2];
  hPMRefClg->Fit("gaus","Q","",xlowRef[2],maxXvalRef[2]+100.);
  funRef[2] = hPMRefClg->GetFunction("gaus");
  meanRef[2] = (Float_t) (funRef[2]->GetParameter(1));
  sigmaRef[2] = (Float_t) (funRef[2]->GetParameter(2));
  
  // ~~~~~~~~ PM Ref side A low gain chain ~~~~~~~~
  maxBinRef[3] = hPMRefAlg->GetMaximumBin();
  nBinRef[3] = (hPMRefAlg->GetXaxis())->GetNbins();
  xMaxRef[3] = (hPMRefAlg->GetXaxis())->GetXmax();
  maxXvalRef[3] = maxBinRef[3]*xMaxRef[3]/nBinRef[3];
  //
  if(maxXvalRef[3]-100.<0.) {xlowRef[3]=0.;}
  else xlowRef[3] = maxXvalRef[3];
  hPMRefAlg->Fit("gaus","Q","",xlowRef[3],maxXvalRef[3]+100.);
  funRef[3] = hPMRefAlg->GetFunction("gaus");
  meanRef[3] = (Float_t) (funRef[3]->GetParameter(1));
  sigmaRef[3] = (Float_t) (funRef[3]->GetParameter(2));
  //
  FILE *fileShuttle;
  fileShuttle = fopen(LASDATA_FILE,"w");
  for(Int_t i=0; i<4; i++)  fprintf(fileShuttle,"\t%f\t%f\n",meanRef[i], sigmaRef[i]); 
  //						       
  fclose(fileShuttle);
  //
  delete hPMRefChg;
  delete hPMRefAhg;
  delete hPMRefClg;
  delete hPMRefAlg;

  //delete minuitFit;
  TVirtualFitter::SetFitter(0);

  /* write report */
  fprintf(fp,"Run #%s, received %d physics events out of %d\n",getenv("DATE_RUN_NUMBER"),nevents_physics,nevents_total);

  /* close result file */
  fclose(fp);
  
  /* report progress */
  daqDA_progressReport(90);

  /* store the result file on FES */
  status = daqDA_FES_storeFile(MAPDATA_FILE,MAPDATA_FILE);
  if(status){
    printf("Failed to export file : %d\n",status);
    return -1;
  }
  //
  status = daqDA_FES_storeFile(LASDATA_FILE,LASDATA_FILE);
  if(status){
    printf("Failed to export file : %d\n",status);
    return -1;
  }

  /* report progress */
  daqDA_progressReport(100);


  return status;
}
