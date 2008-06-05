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
DA Type: 
Number of events needed: no constraint (tipically ~10^3)
Input Files: 
Output Files: ZDCLaser.dat
Trigger Types Used: Standalone Trigger

*/

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

//AliRoot
#include <AliRawReaderDate.h>
#include <AliZDCRawStream.h>


/* Main routine
      Arguments: list of DATE raw data files
*/
int main(int argc, char **argv) {

  int status = 0;

  /* log start of process */
  printf("ZDC LASER program started\n");  

  /* check that we got some arguments = list of files */
  if (argc<2) {
    printf("Wrong number of arguments\n");
    return -1;
  }

  // --- Histograms for LASER runs 
  //     20 signal channels + 2 reference PTMs
  //
  TH1F::AddDirectory(0);
  TH1F *hZDCsideC[10], *hZDCsideA[10];
  char nhistZDCC[50], nhistZDCA[50];
  for(Int_t j=0; j<10; j++){
     if(j<5){ // ZNs
       sprintf(nhistZDCC,"ZNCtow%d",j);
       sprintf(nhistZDCA,"ZNAtow%d",j);
     }
     else if(j>=5 && j<10){ // ZPs
       sprintf(nhistZDCC,"ZPCtow%d",j);
       sprintf(nhistZDCA,"ZPAtow%d",j);
     }
     hZDCsideC[j] = new TH1F(nhistZDCC, nhistZDCC, 100, 0., 1000.);
     hZDCsideA[j] = new TH1F(nhistZDCA, nhistZDCA, 100, 0., 1000.);
  }
  TH1F *hPMRefsideC = new TH1F("hPMRefsideC","hPMRefsideC", 100,0.,1000.);
  TH1F *hPMRefsideA = new TH1F("hPMRefsideA","hPMRefsideA", 100,0.,1000.);


  /* open result file */
  FILE *fp=NULL;
  fp=fopen("./result.txt","a");
  if (fp==NULL) {
    printf("Failed to open file\n");
    return -1;
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


      /* use event - here, just write event id to result file */
      eventT=event->eventType;
    
      if(eventT==PHYSICS_EVENT){
        //
        // *** To analyze EMD events you MUST have a pedestal data file!!!
        // *** -> check if a pedestal run has been analyzied
        FILE *filePed=NULL;
        filePed=fopen("./ZDCPedestal.dat","r");
        if (filePed==NULL) {
          printf("\t ERROR!!! You MUST have a ZDCPedestal.dat file!!!\n");
          return -1;
        }
        // 144 = 48 in-time + 48 out-of-time + 48 correlations
        Float_t readValues[2][144], MeanPed[44], MeanPedWidth[44], 
   	      MeanPedOOT[44], MeanPedWidthOOT[44],
   	      CorrCoeff0[44], CorrCoeff1[44];
        //
        for(int jj=0; jj<144; jj++){
          for(int ii=0; ii<2; ii++){
             fscanf(filePed,"%f",&readValues[ii][jj]);
          }
	  if(jj<48){
	    MeanPed[jj] = readValues[0][jj];
	    MeanPedWidth[jj] = readValues[1][jj];
	  }
	  else if(jj>48 && jj<96){
	    MeanPedOOT[jj-48] = readValues[0][jj];
	    MeanPedWidthOOT[jj-48] = readValues[1][jj];
	  }
	  else if(jj>144){
	    CorrCoeff0[jj-96] = readValues[0][jj]; 
	    CorrCoeff1[jj-96] = readValues[1][jj];;
	  }
        }
        //
        // Initalize raw-data reading and decoding
        AliRawReader *reader = new AliRawReaderDate((void*)event);
        const AliRawDataHeader* header = reader->GetDataHeader();
        if(header) {
         UChar_t message = header->GetAttributes();
	 if(message & 0x20){ // DEDICATED LASER RUN
	    printf("\t STANDALONE_LASER_RUN raw data found\n");
	    continue;
	 }
	 else{
	    printf("\t NO STANDALONE_LASER_RUN raw data found\n");
	    return -1;
	 }
  	}
  	//Commented until we won't have true Raw Data Header...
  	//else{
  	//   printf("\t ATTENTION! No Raw Data Header found!!!\n");
  	  // return -1;
  	//}
  	//
  	AliZDCRawStream *rawStreamZDC = new AliZDCRawStream(reader);
  	//
  	if (!rawStreamZDC->Next()) printf(" \t No raw data found!! \n");
  	//
	while(rawStreamZDC->Next()){
  	  Int_t index=-1;
	  // Implemented only for HIGH gain chain
  	  if((rawStreamZDC->IsADCDataWord()) && (rawStreamZDC->GetADCGain()==0)){
	    index = rawStreamZDC->GetADCChannel();
	    Float_t Pedestal = MeanPed[index];
	    Float_t CorrADC = rawStreamZDC->GetADCValue() - Pedestal;
	    if(rawStreamZDC->GetSector(0)==1){
	      if(rawStreamZDC->GetSector(1)==5) hPMRefsideC->Fill(CorrADC);
	      else hZDCsideC[rawStreamZDC->GetSector(1)]->Fill(CorrADC);
	    }
	    else if(rawStreamZDC->GetSector(0)==2){
	      hZDCsideC[rawStreamZDC->GetSector(1)+5]->Fill(CorrADC);
	    }
	    else if(rawStreamZDC->GetSector(0)==4){
	      if(rawStreamZDC->GetSector(1)==5) hPMRefsideA->Fill(CorrADC);
	      else hZDCsideA[rawStreamZDC->GetSector(1)]->Fill(CorrADC);
	    }
	    else if(rawStreamZDC->GetSector(0)==5){
	      hZDCsideA[rawStreamZDC->GetSector(1)+5]->Fill(CorrADC);
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
  Int_t maxBinC[10], maxBinA[10], maxBinRef[2];
  Int_t nBinC[10], nBinA[10], nBinRef[2];
  Float_t xMaxC[10], xMaxA[10], xMaxRef[2];
  Float_t maxXvalC[10], maxXvalA[10], maxXvalRef[2];
  Float_t xlowC[10], xlowA[10], xlowRef[10]; 
  TF1 *funA[10], *funC[10], *funRef[2];
  //
  Float_t meanC[10], meanA[10], meanRef[2];
  Float_t sigmaA[10], sigmaC[10], sigmaRef[10];
  //
  for(Int_t k=0; k<10; k++){
     maxBinC[k] = hZDCsideC[k]->GetMaximumBin();
     nBinC[k] = (hZDCsideC[k]->GetXaxis())->GetNbins();
     xMaxC[k] = (hZDCsideC[k]->GetXaxis())->GetXmax();
     maxXvalC[k] = maxBinC[k]*xMaxC[k]/nBinC[k];
     //
     if(maxXvalC[k]-100.<0.) {xlowC[k]=0.;}
     else xlowC[k] = maxXvalC[k];
     hZDCsideC[k]->Fit("gaus","Q","",xlowC[k],maxXvalC[k]+100.);
     funC[k] = hZDCsideC[k]->GetFunction("gaus");
     meanC[k] = (Float_t) (funC[k]->GetParameter(1));
     sigmaC[k] = (Float_t) (funC[k]->GetParameter(2));
     //
     maxBinA[k] = hZDCsideA[k]->GetMaximumBin();
     nBinA[k] = (hZDCsideA[k]->GetXaxis())->GetNbins();
     xMaxA[k] = (hZDCsideA[k]->GetXaxis())->GetXmax();
     maxXvalA[k] = maxBinA[k]*xMaxA[k]/nBinA[k];
     //
     if(maxXvalA[k]-100.<0.) {xlowA[k]=0.;}
     else xlowA[k] = maxXvalA[k];
     hZDCsideA[k]->Fit("gaus","Q","",xlowA[k],maxXvalA[k]+100.);
     funA[k] = hZDCsideC[k]->GetFunction("gaus");
     meanA[k] = (Float_t) (funA[k]->GetParameter(1));
     sigmaA[k] = (Float_t) (funA[k]->GetParameter(2));
  }
  //
  maxBinRef[0] = hPMRefsideC->GetMaximumBin();
  nBinRef[0] = (hPMRefsideC->GetXaxis())->GetNbins();
  xMaxRef[0] = (hPMRefsideC->GetXaxis())->GetXmax();
  maxXvalRef[0] = maxBinRef[0]*xMaxRef[0]/nBinRef[0];
  //
  if(maxXvalRef[0]-100.<0.) {xlowRef[0]=0.;}
  else xlowRef[0] = maxXvalRef[0];
  hPMRefsideC->Fit("gaus","Q","",xlowRef[0],maxXvalRef[0]+100.);
  funRef[0] = hPMRefsideC->GetFunction("gaus");
  meanRef[0] = (Float_t) (funRef[0]->GetParameter(1));
  sigmaRef[0] = (Float_t) (funRef[0]->GetParameter(2));
  //
  maxBinRef[1] = hPMRefsideA->GetMaximumBin();
  nBinRef[1] = (hPMRefsideA->GetXaxis())->GetNbins();
  xMaxRef[1] = (hPMRefsideA->GetXaxis())->GetXmax();
  maxXvalRef[1] = maxBinRef[1]*xMaxRef[1]/nBinRef[1];
  //
  if(maxXvalRef[1]-100.<0.) {xlowRef[1]=0.;}
  else xlowRef[1] = maxXvalRef[1];
  hPMRefsideA->Fit("gaus","Q","",xlowRef[1],maxXvalRef[1]+100.);
  funRef[1] = hPMRefsideA->GetFunction("gaus");
  meanRef[1] = (Float_t) (funRef[1]->GetParameter(1));
  sigmaRef[1] = (Float_t) (funRef[1]->GetParameter(2));
  //
  FILE *fileShuttle;
  const char *fName = "ZDCLaser.dat";
  fileShuttle = fopen(fName,"w");
  for(Int_t i=0; i<10; i++) fprintf(fileShuttle,"\t%f\t%f\n",meanC[i], sigmaC[i]);
  for(Int_t i=0; i<10; i++) fprintf(fileShuttle,"\t%f\t%f\n",meanA[i], sigmaA[i]);
  for(Int_t i=0; i<2; i++)  fprintf(fileShuttle,"\t%f\t%f\n",meanRef[i], sigmaRef[i]); 
  //						       
  fclose(fileShuttle);
  //
  for(Int_t j=0; j<10; j++){
     delete hZDCsideC[j];
     delete hZDCsideA[j];
     delete hPMRefsideC;
     delete hPMRefsideA;
  }

  /* write report */
  fprintf(fp,"Run #%s, received %d physics events out of %d\n",getenv("DATE_RUN_NUMBER"),nevents_physics,nevents_total);

  /* close result file */
  fclose(fp);
  
  /* report progress */
  daqDA_progressReport(90);

  /* store the result file on FES */
  status = daqDA_FES_storeFile(fName,"ZDCLASER_data");
  if(status){
    printf("Failed to export file : %d\n",status);
    return -1;
  }

  /* report progress */
  daqDA_progressReport(100);


  return status;
}
