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
  // --- Histos for reference PMTs (high gain chains)
  TH1F *hPMRefC = new TH1F("hPMRefC","hPMRefC", 100,0.,1400.);
  TH1F *hPMRefA = new TH1F("hPMRefA","hPMRefA", 100,0.,1400.);
  //
  // --- Histos for detector PMTs (just high gain chain)
  TH1F *hZNC[5], *hZPC[5], *hZNA[5], *hZPA[5];
  char hnamZNC[20], hnamZPC[20], hnamZNA[20], hnamZPA[20];
  for(Int_t j=0; j<5; j++){
    sprintf(hnamZNC,"ZNC-tow%d",j);
    sprintf(hnamZPC,"ZPC-tow%d",j);
    sprintf(hnamZNA,"ZNA-tow%d",j);
    sprintf(hnamZPA,"ZPA-tow%d",j);
    //
    hZNC[j] = new TH1F(hnamZNC, hnamZNC, 100, 0., 1400.);
    hZPC[j] = new TH1F(hnamZPC, hnamZPC, 100, 0., 1400.);
    hZNA[j] = new TH1F(hnamZNA, hnamZNA, 100, 0., 1400.);
    hZPA[j] = new TH1F(hnamZPA, hnamZPA, 100, 0., 1400.);
  }

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

	rawStreamZDC->SetSODReading(kTRUE);
	  	
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

	rawStreamZDC->SetSODReading(kTRUE);

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
	  Int_t detector = rawStreamZDC->GetSector(0);
	  
  	  if(rawStreamZDC->IsADCDataWord() && !(rawStreamZDC->IsUnderflow())
	     && !(rawStreamZDC->IsOverflow()) && detector!=-1){
	    
	    printf("  IsADCWord %d, IsUnderflow %d, IsOverflow %d\n",
	      rawStreamZDC->IsADCDataWord(),rawStreamZDC->IsUnderflow(),rawStreamZDC->IsOverflow());
 
 	    if(rawStreamZDC->GetSector(1)!=5){ // Physics signals
    	      if(detector==1) index = rawStreamZDC->GetSector(1);        // *** ZNC
	      else if(detector==2) index = rawStreamZDC->GetSector(1)+5; // *** ZPC
	      else if(detector==4) index = rawStreamZDC->GetSector(1)+12;// *** ZNA
	      else if(detector==5) index = rawStreamZDC->GetSector(1)+17;// *** ZPA
	    }
	    else{ // Reference PMs
	      index = (detector-1)/3+22;
	    }
	    
	    Float_t Pedestal = MeanPed[index];
	    Float_t CorrADC = rawStreamZDC->GetADCValue() - Pedestal;
	    
	    // **** Detector PMs
	    if(rawStreamZDC->GetSector(1)!=5 && rawStreamZDC->GetADCGain()==0){
	      // ---- side C
	      hZNC[rawStreamZDC->GetSector(1)]->Fill(CorrADC);
	      hZPC[rawStreamZDC->GetSector(1)]->Fill(CorrADC);
	      // ---- side A
	      hZNA[rawStreamZDC->GetSector(1)]->Fill(CorrADC);
	      hZPA[rawStreamZDC->GetSector(1)]->Fill(CorrADC);
	    }
	    // **** Reference PMs
	    if(rawStreamZDC->GetSector(1)==5 && rawStreamZDC->GetADCGain()==0){
	      // ---- PMRef chain side C
	      if(detector==1) hPMRefC->Fill(CorrADC);
	      // ---- PMRef side A
	      else if(detector==4) hPMRefA->Fill(CorrADC);
	    }	  
  	  }//IsADCDataWord()+NOunderflow+NOoverflow
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
  Int_t maxBin[22], nBin[22];
  Float_t xMax[22], maxXval[22], xlow[22]; 
  Float_t mean[22], sigma[22];
  TF1 *fun[4];
  
  for(Int_t k=0; k<5; k++){
    // --- ZNC
    maxBin[k] = hZNC[k]->GetMaximumBin();
    nBin[k] = (hZNC[k]->GetXaxis())->GetNbins();
    xMax[k] = (hZNC[k]->GetXaxis())->GetXmax();
    if(nBin[k]!=0) maxXval[k] = maxBin[k]*xMax[k]/nBin[k];
    //
    if(maxXval[k]-150.<0.) xlow[k]=0.;
    else xlow[k] = maxXval[k]-150.;
    hZNC[k]->Fit("gaus","Q","",xlow[k],maxXval[k]+150.);
    fun[k] = hZNC[k]->GetFunction("gaus");
    mean[k]  = (Float_t) (fun[k]->GetParameter(1));
    sigma[k] = (Float_t) (fun[k]->GetParameter(2));
    // --- ZPC
    maxBin[k+5] = hZPC[k]->GetMaximumBin();
    nBin[k+5] = (hZPC[k]->GetXaxis())->GetNbins();
    xMax[k+5] = (hZPC[k]->GetXaxis())->GetXmax();
    if(nBin[k+5]!=0) maxXval[k+5] = maxBin[k+5]*xMax[k+5]/nBin[k+5];
    //
    if(maxXval[k+5]-150.<0.) xlow[k+5]=0.;
    else xlow[k+5] = maxXval[k+5]-150.;
    hZPC[k]->Fit("gaus","Q","",xlow[k+5],maxXval[k+5]+150.);
    fun[k+5] = hZPC[k]->GetFunction("gaus");
    mean[k+5]  = (Float_t) (fun[k+5]->GetParameter(1));
    sigma[k+5] = (Float_t) (fun[k+5]->GetParameter(2));
    // --- ZNA
    maxBin[k+10] = hZNA[k]->GetMaximumBin();
    nBin[k+10] = (hZNA[k]->GetXaxis())->GetNbins();
    xMax[k+10] = (hZNA[k]->GetXaxis())->GetXmax();
    if(nBin[k+10]!=0) maxXval[k+10] = maxBin[k+10]*xMax[k+10]/nBin[k+10];
    //
    if(maxXval[k+10]-150.<0.) xlow[k+10]=0.;
    else xlow[k+10] = maxXval[k+10]-150.;
    hZNA[k]->Fit("gaus","Q","",xlow[k+10],maxXval[k+10]+150.);
    fun[k+10] = hZNA[k]->GetFunction("gaus");
    mean[k+10]  = (Float_t) (fun[k+10]->GetParameter(1));
    sigma[k+10] = (Float_t) (fun[k+10]->GetParameter(2));
    // --- ZPA
    maxBin[k+15] = hZPA[k]->GetMaximumBin();
    nBin[k+15] = (hZPA[k]->GetXaxis())->GetNbins();
    xMax[k+15] = (hZPA[k]->GetXaxis())->GetXmax();
    if(nBin[k+15]!=0) maxXval[k+15] = maxBin[k+15]*xMax[k+15]/nBin[k+15];
    //
    if(maxXval[k+15]-150.<0.) xlow[k+15]=0.;
    else xlow[k+15] = maxXval[k+15]-150.;
    hZPA[k]->Fit("gaus","Q","",xlow[k+15],maxXval[k+15]+150.);
    fun[k+15] = hZPA[k]->GetFunction("gaus");
    mean[k+15]  = (Float_t) (fun[k+15]->GetParameter(1));
    sigma[k+15] = (Float_t) (fun[k+15]->GetParameter(2));
    
  }
  
  // ~~~~~~~~ PM Ref side C ~~~~~~~~
  maxBin[20] = hPMRefC->GetMaximumBin();
  nBin[20] = (hPMRefC->GetXaxis())->GetNbins();
  xMax[20] = (hPMRefC->GetXaxis())->GetXmax();
  if(nBin[20]!=0) maxXval[20] = maxBin[20]*xMax[20]/nBin[20];
  // 
  if(maxXval[20]-150.<0.) xlow[20]=0.;
  else xlow[20] = maxXval[20];
  hPMRefC->Fit("gaus","Q","",xlow[20],maxXval[20]+150.);
  fun[20] = hPMRefC->GetFunction("gaus");
  mean[20]  = (Float_t) (fun[20]->GetParameter(1));
  sigma[20] = (Float_t) (fun[20]->GetParameter(2));
  
  // ~~~~~~~~ PM Ref side A ~~~~~~~~
  maxBin[21] = hPMRefA->GetMaximumBin();
  nBin[21] = (hPMRefA->GetXaxis())->GetNbins();
  xMax[21] = (hPMRefA->GetXaxis())->GetXmax();
  if(nBin[21]!=0) maxXval[21] = maxBin[21]*xMax[21]/nBin[21];
  //
  if(maxXval[21]-100.<0.) xlow[21]=0.;
  else xlow[21] = maxXval[21];
  hPMRefA->Fit("gaus","Q","",xlow[21],maxXval[21]+100.);
  fun[21] = hPMRefA->GetFunction("gaus");
  mean[21]  = (Float_t) (fun[21]->GetParameter(1));
  sigma[21] = (Float_t) (fun[21]->GetParameter(2));
    
  FILE *fileShuttle;
  fileShuttle = fopen(LASDATA_FILE,"w");
  Int_t det[22]  = {1,1,1,1,1,2,2,2,2,2,4,4,4,4,4,5,5,5,5,5,1,4};
  Int_t quad[22] = {0,1,2,3,4,0,1,2,3,4,0,1,2,3,4,0,1,2,3,4,5,5};
  for(Int_t i=0; i<22; i++){
    fprintf(fileShuttle,"\t%d\t%d\t%f\t%f\n",det[i],quad[i],mean[i], sigma[i]); 
  }
  //						       
  fclose(fileShuttle);
  //
  for(Int_t j=0; j<5; j++){
    delete hZNC[j];
    delete hZPC[j];
    delete hZNA[j];
    delete hZPA[j];
  }
  delete hPMRefC;
  delete hPMRefA;

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
