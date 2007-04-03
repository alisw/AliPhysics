extern "C" {
#include <daqDA.h>
}

#include "event.h"
#include "monitor.h"

#include <Riostream.h>
#include <stdio.h>
#include <stdlib.h>

//AliRoot
#include "AliRawReaderDate.h"
#include "AliRawReader.h"
#include "AliT0digit.h"
#include "AliT0RawReader.h"
#include "AliT0Dqclass.h"

//ROOT
#include "TFile.h"
#include "TKey.h"
#include "TH2S.h"
#include "TObject.h"
#include "TBenchmark.h"
#include "TRandom.h"

#include "TCanvas.h"
#include "TString.h"
#include "TH1.h"
#include "TF1.h"
#include "TSpectrum.h"
#include "TVirtualFitter.h"

void fitv2DA(char *filename = "t0histdate.root");

/* Main routine
      Arguments: 
      1- monitoring data source
*/
int main(int argc, char **argv) {
  int status;
  
  if (argc!=2) {
    printf("Wrong number of arguments\n");
    return -1;
  }


  /* define data source : this is argument 1 */  
  status=monitorSetDataSource( argv[1] );
  if (status!=0) {
    printf("monitorSetDataSource() failed : %s\n",monitorDecodeError(status));
    return -1;
  }


  /* declare monitoring program */
  status=monitorDeclareMp( __FILE__ );
  if (status!=0) {
    printf("monitorDeclareMp() failed : %s\n",monitorDecodeError(status));
    return -1;
  }


  /* define wait event timeout - 1s max */
  monitorSetNowait();
  monitorSetNoWaitNetworkTimeout(1000);
  

  /* log start of process */
  printf("T0 monitoring program started\n");  


  // Allocation of histograms - start
  TH1F * hCFD[24]; TH1F *hLED[24]; TH1F*hQT01[4];TH1F*hQT02[4];
  TH1F*hQTD[4]; TH1F*hADC[24]; TH2*hQTCFD[4];
  Char_t buf1[10], buf2[10], buf3[10], buf4[10], buf5[10], buf6[10],buf7[10];
  for (Int_t ic=0; ic<24; ic++)
    {
      sprintf(buf1,"CFD%i",ic+1);
      hCFD[ic]= new TH1F(buf1,"CFD",6000,0,6000);
      sprintf(buf2,"LED%i",ic+1);
      hLED[ic]= new TH1F(buf2,"LED",2000,0,2000);
      //LED-CFD
      sprintf(buf6,"ADC%i",ic+1);
      hADC[ic]= new TH1F(buf6,"ADC",6000,0,6000);
    }
  
  for (Int_t iq=0; iq<4; iq++)
    {
      //QT01 - QT04
      sprintf(buf3,"QT0%i",iq+1);
      hQT01[iq]= new TH1F(buf3,"QT01",6000,0,6000);
      sprintf(buf4,"QT1_%i",iq+1);
      //QT11 - QT14 
      hQT02[iq]= new TH1F(buf4,"QT02",6000,0,6000);
      sprintf(buf5,"QTD_%i",iq+1);
      //QT11-QT01 ....
      hQTD[iq]= new TH1F(buf5,"QTdiff",4500,1500,6000);

      sprintf(buf7,"QTCFD_%i",iq+1);
      hQTCFD[iq]= new TH2F(buf7,"QT vs CFD",500,0,6000,500,0,5000);
  }

  TH1F *hORA= new TH1F("hORA"," T0 A ",1000, 1000,2000);
  TH1F *hORC= new TH1F("hORC"," T0 C ",1000, 1000,2000);
  TH1F*hEffCFD= new TH1F("hEffCFD","Effeciency",8,0.25,4.25);

  //  TH2F*hQTCFD= new TH2F("hQTCFD","QT vs CFD",500,0.5,6000.5,500,0.5,5000.5);
  TH2F*hQTLED= new TH2F("hQTLED","QT vs LED",500,0.5,6000.5,500,0.5,5000.5);
  TH2F*hLEDCFD= new TH2F("hLEDCFD","LEd vs CFD",500,-0.5,10000.5,500,-0.5,10000.5);
  // Allocation of histograms - end

  Int_t iev=0;
  /* main loop (infinite) */
  for(;;) {
    struct eventHeaderStruct *event;
    eventTypeType eventT;
  
    /* check shutdown condition */
    if (daqDA_checkShutdown()) {break;}
    
    /* get next event (blocking call until timeout) */
    status=monitorGetEventDynamic((void **)&event);
    if (status==(int)MON_ERR_EOF) {
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

    iev++; 

    /* use event - here, just write event id to result file */
    eventT=event->eventType;
    
    if (eventT==PHYSICS_EVENT) {
      printf(" event number = %i \n",iev);

      // Initalize raw-data reading and decoding
      AliRawReader *reader = new AliRawReaderDate((void*)event);
      // Enable the following two lines in case of real-data
      //    reader->LoadEquipmentIdsMap("T0map.txt");
      //    reader->RequireHeader(kFALSE);
      AliT0RawReader *start = new AliT0RawReader(reader);

      // Read raw data
      Int_t allData[110][5];
      start->Next();
      for (Int_t i=0; i<110; i++) {
	allData[i][0]= start->GetData(i,0);
	if (allData[i][0] != 0) cout<<"event "<<event<<" i "<<iev<<" "<<allData[i][0] - allData[0][0]<<endl;
      }

      // Fill the histograms
      for(Int_t ik=0; ik<24; ik++) { 
	for (Int_t iHit=0; iHit<5; iHit++) {
	  if(allData[ik+1][iHit] != 0 ) {
	    //	   cout<<"event "<<event<<" iHit "<< iHit<< allData[ik+1][0]<<" "<< allData[0][0]<< endl;
	    hCFD[ik]->Fill((allData[ik+1][iHit] - allData[0][0])/1000.); 
	  }
	  if(allData[ik+25][iHit] != 0 )    hLED[ik]->Fill(allData[ik+25][iHit] - allData[0][0]);
	  if(allData[ik+25][iHit] != 0 || allData[ik+1][iHit] !=0)   
	    hADC[ik]->Fill(allData[ik+25][iHit] - allData[ik+1][iHit]);
	}
      }
      hLEDCFD->Fill(allData[9][0],allData[1][0]);
      
      hQT01[1]->Fill(allData[18][0] - allData[0][0]);
      hQT01[2]->Fill(allData[19][0] - allData[0][0]);
      hQT01[3]->Fill(allData[20][0] - allData[0][0]);
      hQT02[0]->Fill(allData[21][0] - allData[0][0]);
      hQT02[0]->Fill(allData[22][0] - allData[0][0]);
      hQT02[0]->Fill(allData[23][0] - allData[0][0]);
      hQT02[0]->Fill(allData[24][0] - allData[0][0]);

      
      hQTD[0]->Fill(allData[21][0] - allData[17][0]);
      hQTLED->Fill(allData[9][0] - allData[0][0], allData[21][0]-allData[17][0]);
      hQTCFD[0]->Fill(allData[21][0]-allData[17][0], allData[1][0] - allData[0][0]);
		      
      hQTD[1]->Fill(allData[22][0]-allData[18][0]);
      hQTCFD[1]->Fill(allData[22][0]-allData[18][0], allData[2][0] - allData[0][0]);
      
      hQTD[2]->Fill(allData[23][0]-allData[19][0]);
      hQTCFD[2]->Fill(allData[23][0]-allData[19][0], allData[3][0] - allData[0][0]);
      
      hQTD[3]->Fill(allData[24][0]-allData[20][0]);
      hQTCFD[3]->Fill(allData[24][0]-allData[20][0], allData[4][0] - allData[0][0]);
      

      hORA->Fill(allData[25][0] - allData[0][0]);
      hORA->Fill(allData[26][0] - allData[0][0]);
      // End of fill histograms

    }



    /* free resources */
    free(event);
    
    /* exit when last event received, no need to wait for TERM signal */
    if (eventT==END_OF_RUN) {
      printf("EOR event detected\n");
      break;
    }
  }

  // write a file with the histograms
  Char_t filehist[20]; 
  sprintf(filehist,"t0histdate.root");
  TFile *hist = new TFile(filehist,"RECREATE");
  
  hLEDCFD->Write();

  for (Int_t i=0; i<24; i++)
    {
      hCFD[i]->Write();
      hLED[i]->Write();
      hADC[i]->Write();
    }
  for (Int_t i=0; i<4; i++)
    {
      hQT01[i]->Write();
      hQT02[i]->Write();
      hQTD[i]->Write();
      hQTCFD[i]->Write();
    }

  hEffCFD->Write();
  hORA->Write();
  hORC->Write();
  hQTLED->Write();
  hist->Close();

  // Fit the histograms and write the output file
  fitv2DA();  

  return status;
}
  
void fitv2DA(char *filename)
{
	TFile *file = TFile::Open(filename);
	FILE *dafile;
	dafile = fopen("da.txt","w");
	fprintf(dafile,"Const\t\tMean\t\tSigma\t\tLBorder\tRBorder\n");

//	TFile *fithist = new TFile("fitedhist.root","RECREATE");

	char histoname[10];
	char canvname[10];
	Float_t p[3] = {0.,0.,0.};
	Float_t LBordX = 0.;
	Float_t RBordX = 0.;
/*
	TIter next(file.GetListOfKeys());
   	TKey *key;
      	while ((key = (TKey*)nextkey())) {
	      if(key->GetName()==)
	}
*/
	Int_t npeaks = 10;
	cout<<"npeaks = "<<npeaks<<endl;
	TH1F *hist[24];
	TCanvas *c[24];
	for(Int_t i=0;i<24;i++){///////////#######Cos tu sie dzieje....
		sprintf(histoname,"CFD%d",i+1);
		sprintf(canvname,"C%d",i+1);
		TH1F *histtemp = (TH1F*)file->Get(histoname);
		hist[i] = histtemp;
	}	
	cout<<"Wczytany plik z wykresami..."<<endl;
//	TCanvas *c[8]; 
	for(Int_t ii=0;ii<24;ii++){
		c[ii] = new TCanvas(canvname,canvname,800,600);
	p[0] = 0.;
	p[1] = 0.;
	p[2] = 0.;
	LBordX = 0.;
	RBordX = 0.;

	Double_t max=0.0; 
//	Double_t peak,peakmax=0.0;
	Double_t tabmax[2] = {0.0, 0.0};

	TSpectrum *s = new TSpectrum(2*npeaks);
	Int_t nfound = s->Search(hist[ii],2,"goff",0.05);
	cout<<"Found "<<nfound<<" peaks\n";
	if(nfound!=0){
	Float_t *xpeak = s->GetPositionX();
	for(Int_t k=0;k<nfound;k++)
	{
		cout<<"Petla po pikach..."<<endl;
		Float_t xp = xpeak[k];
		Int_t xbin = hist[ii]->GetXaxis()->FindBin(xp);
		Float_t yp = hist[ii]->GetBinContent(xbin);
		if(yp>max) {
			max = yp;
			tabmax[0] = xp;
			tabmax[1] = yp;
		}
		cout<<"xbin = "<<xbin<<"\txpeak = "<<xpeak[k]<<"\typeak = "<<yp<<endl;
	}
	
	Double_t MinBordX, MaxBordX, MaxBordY;

	MinBordX = tabmax[0]-tabmax[0]*0.05;
	MaxBordX = tabmax[0]+tabmax[0]*0.05;
	MaxBordY = tabmax[1]+tabmax[1]*0.5;

        Double_t border;
        border = tabmax[1]*0.1;
        
	Double_t LBordY = 99999.;
	Double_t RBordY = 99999.;
        LBordX = tabmax[0];
	RBordX = tabmax[0];
        do{
                LBordY = hist[ii]->GetBinContent(LBordX);
		LBordX -=1;
        }while(LBordY>border);
	do{
		RBordX +=1;
		RBordY = hist[ii]->GetBinContent(RBordX);
	}while(RBordY>border);
	
	
	hist[ii]->GetXaxis()->SetRangeUser(MinBordX,MaxBordX);
	hist[ii]->GetYaxis()->SetRangeUser(0,MaxBordY);
	TF1 *fit = new TF1("fit","gaus",LBordX,RBordX);
	hist[ii]->Fit("fit","IR");
	Double_t maxx=0.;
//	Double_t chi2 = fit->GetChisquare();
//	Float_t p[3];
	
	for(Int_t j=0;j<3;j++){
		p[j] = fit->GetParameter(j);
	}

	max = fit->GetMaximum();
	maxx = fit->GetMaximumX();
//	hist[ii]->Write();
	cout<<"\nFitMaxY = "<<max<<"\tFitMaxX "<<maxx<<"\tPeakMax*10% "<<border<<endl;
	cout<<"LBordY "<<LBordY<<"\tLBordX "<<LBordX+1<<"\tRBordY "<<RBordY<<"\tRBordX "<<RBordX<<endl;
	cout<<"Parameters: Const - MaxY - "<<p[0]<<"\tMean - "<<p[1]<<"\tSigma - "<<p[2]<<endl;
	cout<<"Borders: Min - "<<MinBordX<<"\tMax - "<<MaxBordX<<endl;
	if(p[0]==0.&&p[1]==0.&&p[2]==0.)	RBordX = 0;
	}
//	else{
//	return;
//	}
	AliT0Dqclass *daqpar = new AliT0Dqclass();
	cout<<"Alit0dqclass..."<<endl;
	daqpar->SetTime(0,p[0]);
	cout<<"First parameter..."<<endl;
	daqpar->SetTime(1,p[1]);
	daqpar->SetTime(2,p[2]);
	cout<<"gauss parameters finished..."<<endl;
	daqpar->SetTime(3,LBordX+1);
	cout<<"first border..."<<endl;
	daqpar->SetTime(4,RBordX);
	cout<<"second border..."<<endl;
	TFile *DA = new TFile("DAfile.root","RECREATE");
	DA->cd();
	cout<<"file"<<endl;
	daqpar->Write("Time");
	cout<<"file written"<<endl;
	DA->Close();
	cout<<"file closed"<<endl;
	delete DA;
	cout<<"file deleted"<<endl;

	fprintf(dafile,"%lf\t%lf\t%lf\t%d\t%d\n",p[0],p[1],p[2],LBordX+1,RBordX);
	}	
	

	fclose(dafile);
}
