#include "TH1F.h"
#include "TH2F.h"
#include "TMap.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TGLabel.h"
#include "TGFileDialog.h"
#include "TGrid.h"
#include "TGFileDialog.h"
#include "TControlBar.h"
#include <Riostream.h>

#include "AliCDBManager.h"
#include "AliRawReader.h"
#include "AliRawReaderRoot.h"
#include "AliT0LookUpValue.h"
#include "AliT0LookUpKey.h"
#include "AliT0Parameters.h"
#include "AliT0RawReader.h"
#include "AliT0Calibrator.h"


void calibLaser1st(Int_t run, Int_t nruns)
{

  // reading laser data for set of MIPs 
  // filling histograms

  /* 176424
 Mips RUN 0.5 176418 0.6 176419 0.7 176420 0.8 176421 0.9 176422 1 176423 1.2 176424 1.4 176425 1.6 176426 1.8 176427 2.0 176428  3 176429 4 176430 5 176431 6 176432 7 176433 8 176434 9 176435 10 176436  15 176437 20 176438 25 176439 
  */
   TString histFile;
  
  TH1F *hChannel[107];  TH1I *hQTC[24];  TH1I *h1CFDminLED[24];
  Int_t allData[110][5];
  Int_t time[24], adc[24];
  
  Char_t  buf0[20],buf1[20], buf2[20], buf3[20], buf4[20];
  
  TString names[107], type;
  AliT0LookUpKey* lookkey= new AliT0LookUpKey();
  AliT0LookUpValue*  lookvalue= new AliT0LookUpValue();
  AliCDBManager * man = AliCDBManager::Instance();
  man->SetDefaultStorage("raw://");
  man->SetRun(run);
  AliT0Parameters *fParam = AliT0Parameters::Instance();
  fParam->Init();
  AliT0Calibrator *fCalib = new AliT0Calibrator();
  TMap *lookup = fParam->GetMapLookup();
  TMapIter *iter = new TMapIter(lookup);
  
  for( Int_t iline=0; iline<107; iline++)
    {
      lookvalue = ( AliT0LookUpValue*) iter->Next();
      lookkey = (AliT0LookUpKey*) lookup->GetValue((TObject*)lookvalue);
      if(lookkey){
	Int_t key=lookkey->GetKey();
	names[key]=lookkey->GetChannelName();
	if(!names[key].Contains("CFD"))
	  hChannel[key] = new TH1F(names[key].Data(),names[key].Data(),
				 5000, 0, 20000 );
	else
	  hChannel[key] = new TH1F(names[key].Data(),names[key].Data(),
				   1000, 2000, 3000 );
       }
      else
	{printf(" no such value %i \n", iline);}
      
    } 
  
  
  for(Int_t ic=0; ic<24; ic++) {
    {
      sprintf(buf1,"QTC%i",ic+1);
      sprintf(buf4,"LEDminCFD%i",ic+1);
      hQTC[ic] = new TH1I(buf1,"QTC", 2500, 0, 10000);
      h1CFDminLED[ic] = new TH1I(buf4,"LED - CFD", 500, 0, 1000);
      
      
    }
  }
  
  TH1F* hMPDA= new TH1F("hMPDA","MPD A",1000, 0,10000);
  TH1F* hMPDC= new TH1F("hMPDC","MPD C",1000, 0,10000);
   
  Int_t event=0;
  //FILENAMES TO INCLUDE
  const char *fFileName="";
  
  TGrid::Connect("alien://");
  for ( Int_t indexfile=0; indexfile < nruns; indexfile++ ) 
    {
      fFileName=Form("alien:///alice/data/2012/LHC12a_T0/000%i/raw/12000%i005.10.root", run+indexfile, run+indexfile);
      printf(" File %s   (%i)\n ",fFileName,indexfile);
      
      AliRawReader *reader = new AliRawReaderRoot(fFileName);
      if(!reader) continue;
      reader->LoadEquipmentIdsMap("T0map.txt");
      reader->RequireHeader(kTRUE);
      AliT0RawReader *start = new AliT0RawReader(reader);
      
      for (Int_t i0=0; i0<107; i0++ )
	for (Int_t j0=0; j0<5; j0++)  allData[i0][j0]=0; 	
      
      while (reader->NextEvent()) {
	start->Next();
	for (Int_t i0=0; i0<107; i0++)
	  for (Int_t j0=0; j0<5; j0++) allData[i0][j0]=0; 	
	for (Int_t i=0; i<107; i++) 
	  for (Int_t iHit=0; iHit<5; iHit++) 
	    allData[i][iHit]= start->GetData(i,iHit);
	
	if(event%1000 == 0) 
	  printf("Event:%d\n",event);
	
	//     if(event > 10000) break;
	
	if (allData[53][0]>0 && allData[54][0]) 
	  hMPDA->Fill(allData[53][0] - allData[54][0] );
	
	if (allData[105][0]>0 && allData[106][0]) 
	  hMPDC->Fill(allData[105][0] - allData[106][0] );
	
	for (Int_t it = 0; it<24; it=it+2)
	  {
	    Int_t cc=it/2;
	    time[cc]=0;
	    adc[cc]=0;
	    if(allData[cc+1][0] != 0 && allData[it+25][0] != 0 && allData[it+26][0] !=0)
	      hQTC[cc]->Fill(allData[it+25][0]-allData[it+26][0]);
	      
	    if(allData[cc+1][0] != 0 && allData[cc+13][0]!=0 ) 
	      {
		adc[cc] = allData[cc+13][0]-allData[cc+1][0];
		h1CFDminLED[cc]->Fill(adc[cc]);
	      }
	  }
	for (Int_t it = 24; it<48; it=it+2)
	  {
	    Int_t cc=(Int_t)(it/2);
	    time[cc]=0;
	    adc[cc]=0;
	    if( allData[cc+45][0]<2000 || allData[cc+45][0]>3000) continue;
	    if(  allData[it+57][0] > 0 && allData[it+58][0] >0)
	      {
		hQTC[cc]->Fill(allData[it+57][0]-allData[it+58][0]);
	      }
	    if(allData[cc+57][0] != 0 && allData[cc+45][0]!=0 ) 
	      {
		  adc[cc] = allData[cc+57][0]-allData[cc+45][0];
		  h1CFDminLED[cc]->Fill( adc[cc]);
		  //CFD calib//
	      }
	  }
	
	
	// only 1sh hit collected
	for(Int_t ik=0; ik<107; ik++) {
	  if(  allData[ik][0] >100 ) {
	    hChannel[ik] -> Fill(allData[ik][0]);
	  }
	}
	  event++;
	  
      } //event
      
      histFile =Form("caliblaser/histCalib%i.root",run+indexfile);
      cout<<" hists will be in file "<<histFile<<endl; 
      TFile *hist = new TFile(histFile.Data(),"RECREATE");
      hist->cd();
      for(Int_t ik=0; ik<107; ik++)	hChannel[ik] ->Write();
      
      for (Int_t i=0; i<24; i++)
	  {
	    hQTC[i]->Write();
	    h1CFDminLED[i]->Write();
	  }
      
      for(Int_t ik=0; ik<107; ik++)	hChannel[ik] ->Reset();
      
      for (Int_t i=0; i<24; i++)
	{
	  hQTC[i]->Reset();
	  h1CFDminLED[i]->Reset();
	}
	
      hMPDC->Write();
      hMPDA->Write();
      
    } //rfile
      
}
