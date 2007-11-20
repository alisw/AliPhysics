/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id:  */

//____________________________________________________________________
//                                                                          
// T0 - T0. 
//
// This class privides GIU service for reading RAW data from Laser
// during electronics test 
//                                                       
#include "TH1F.h"
#include "TH2F.h"
#include "TMap.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "AliT0RawReader.h"
#include "TGLabel.h"
#include <iostream.h>

#include "AliT0CalibLaserData.h"

#include "AliCDBManager.h"
#include "AliRawReader.h"
#include "AliRawReaderDate.h"
#include "AliT0LookUpValue.h"
#include "AliT0Parameters.h"
#include "AliT0RawReader.h"

ClassImp(AliT0CalibLaserData)

AliT0CalibLaserData::AliT0CalibLaserData() : TObject(),
	                                 fRunNumber(905)
{
//
}
/*
//________________________________________________________________

AliT0CalibLaserData::AliT0CalibLaserData(const AliT0CalibLaserData& calibda) : TObject(),
	                                 fRunNumber(905)
{
//copy constructor
  
}
//________________________________________________________________

AliT0CalibLaserData &AliT0CalibLaserData::operator =(const AliT0CalibLaserData& calibda)
{
// assignment operator

  return *this;
}
//________________________________________________________________
AliT0CalibLaserData::~AliT0CalibLaserData()
{
  //
}
*/
//________________________________________________________________

void AliT0CalibLaserData::ReadHistSize(Int_t rNumber)
{
  //build GUI frame for reading:
  // - run number
  // - histograms rates

    fRunNumber = rNumber;
    
    TGMainFrame* fMain = new TGMainFrame(0,1500,1500);
    fMain->SetLayoutManager( new TGMatrixLayout(fMain,7,7) );
 
    fMain->AddFrame( new TGLabel(fMain, " Histogram") );
    fMain->AddFrame( new TGLabel(fMain, "X min") );
    fMain->AddFrame( new TGLabel(fMain, "X max") );
    fMain->AddFrame( new TGLabel(fMain, "X N# channels") );
  
    fMain->AddFrame( new TGLabel(fMain, "Y min") );
    fMain->AddFrame( new TGLabel(fMain, "Y max") );
    fMain->AddFrame( new TGLabel(fMain, "Y N# channels") );

    fMain->AddFrame( new TGLabel(fMain, "QTC" ) );
    fEntries[0] = new TGNumberEntry(fMain, 0);
    fMain->AddFrame(fEntries[0]);
    fEntries[1] = new TGNumberEntry(fMain, 10000);
    fMain->AddFrame(fEntries[1]);
    fEntries[2] = new TGNumberEntry(fMain, 10000);
    fMain->AddFrame(fEntries[2]);
    fMain->AddFrame( new TGLabel(fMain, " ") );
    fMain->AddFrame( new TGLabel(fMain, " ") );
    fMain->AddFrame( new TGLabel(fMain, " ") );
  
    fMain->AddFrame( new TGLabel(fMain, "LED - CFD" ) );
    fEntries[3] = new TGNumberEntry(fMain, 0);
    fMain->AddFrame(fEntries[3]);
    fEntries[4] = new TGNumberEntry(fMain, 10000);
    fMain->AddFrame(fEntries[4]);
    fEntries[5] = new TGNumberEntry(fMain, 10000);
    fMain->AddFrame(fEntries[5]);
    fMain->AddFrame( new TGLabel(fMain, " ") );
    fMain->AddFrame( new TGLabel(fMain, " ") );
    fMain->AddFrame( new TGLabel(fMain, " ") );
    fMain->AddFrame( new TGLabel(fMain, "CFD vs QTC " ) );
 // QTC axis X
    fEntries[6] = new TGNumberEntry(fMain, 1000.);
    fMain->AddFrame(fEntries[6]);
    fEntries[7] = new TGNumberEntry(fMain, 8000.5);
    fMain->AddFrame(fEntries[7]);
    fEntries[8] = new TGNumberEntry(fMain, 700);
    fMain->AddFrame(fEntries[8]);
// CFD axis Y 
   fEntries[9] = new TGNumberEntry(fMain, 12000.5);
    fMain->AddFrame(fEntries[9]);
    fEntries[10] = new TGNumberEntry(fMain, 18000.5);
    fMain->AddFrame(fEntries[10]);
    fEntries[11] = new TGNumberEntry(fMain, 2000);
    fMain->AddFrame(fEntries[11]);
//
    fMain->AddFrame( new TGLabel(fMain, "CFD vs LED-CFD " ) );
//LED-CFD axis X
    fEntries[12] = new TGNumberEntry(fMain, 0);
    fMain->AddFrame(fEntries[12]);
    fEntries[13] = new TGNumberEntry(fMain, 500);
    fMain->AddFrame(fEntries[13]);
    fEntries[14] = new TGNumberEntry(fMain, 500);
    fMain->AddFrame(fEntries[14]);
// CFD axis Y
    fEntries[15] = new TGNumberEntry(fMain, 14600);
    fMain->AddFrame(fEntries[15]);
    fEntries[16] = new TGNumberEntry(fMain, 14700);
    fMain->AddFrame(fEntries[16]);
    fEntries[17] = new TGNumberEntry(fMain, 100);
    fMain->AddFrame(fEntries[17]);
 
   fMain->AddFrame( new TGLabel(fMain, " Number of run") );
    fEntries[18] = new TGNumberEntry(fMain, 905);
     fMain->AddFrame(fEntries[18]);
    
    for ( int i=0; i<19; i++ ) fEntries[i]->SetWidth(70);
    //    printf( "Max Length %d\n", fEntries[0]->GetMaxWidth() );

    TGTextButton *fOk = new TGTextButton(fMain, "OK");
    fOk->Connect("Clicked()","AliT0CalibLaserData",this,"DoOk()");
    //    fOk->SetCommand(".q");
    fMain->AddFrame(fOk);
    
    fMain->MapSubwindows();
    fMain->Resize();
    fMain->SetWindowName("Dialog");
    fMain->MapWindow();
}
    
void AliT0CalibLaserData::DoOk()
{
 
    //    delete fMain;
    fRunNumber = Int_t (fEntries[18]->GetNumber());
    printf(" RUN NUMBER %i\n",fRunNumber);
    for( int i=0; i<18; i++ ) 
      fHistLimits[i] = fEntries[i]->GetNumber();
 
    ReadData();
}

void AliT0CalibLaserData::ReadData()
{
  // reading RAW data from test LCS
  // filling histograms
  // fillinf tree


      TH1F *hChannel[105];  TH1F *hQTC[12];  
      TH2F *hCFDvsQTC[12]; TH2F *hCFDvsLED[12]; TH1F *h1CFDminLED[12];
      TH1F *hmpd[12];
      Int_t allData[110][5];
      Int_t numberOfHits[105];

      Char_t  buf1[20], buf2[20], buf3[20], buf4[20], buf7[20];

      TTree* digitsTree = new TTree("testData","Tree of test data Digits");
      TBranch *b[106];

       Int_t channels[106];

      TString names[106], type;
      AliT0LookUpKey* lookkey= new AliT0LookUpKey();
      AliT0LookUpValue*  lookvalue= new AliT0LookUpValue();
      AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT");
      AliCDBManager::Instance()->SetRun(0);
      AliT0Parameters *fParam = AliT0Parameters::Instance();
      fParam->Init();
      TMap *lookup = fParam->GetMapLookup();
      TMapIter *iter = new TMapIter(lookup);

      for( Int_t iline=0; iline<106; iline++)
	{
	  lookvalue = ( AliT0LookUpValue*) iter->Next();
	   lookkey = (AliT0LookUpKey*) lookup->GetValue((TObject*)lookvalue);
	  if(lookkey){
	    Int_t key=lookkey->GetKey();
	    names[key]=lookkey->GetChannelName();
	    //	cout<<lookkey->GetKey()<<" "<<lookkey->GetChannelName()<<" trm "<<lookvalue->GetTRM()<<" tdc "<<lookvalue->GetTDC()<<" chain  "<<lookvalue->GetChain()<<" channel "<<lookvalue->GetChannel()<<endl;
	    hChannel[key] = new TH1F(names[key].Data(),names[key].Data(),30000,0,30000);
	    //	hitsname="xHits" + names[key];
	    //	hNumHits[key] = new TH1F(hitsname.Data(),hitsname.Data(),50,-0.25,24.25);
	    type =names[key] + "/I";
	    b[key]=digitsTree->Branch(names[key].Data(),&channels[key], type);
	   }
	  else
	    {printf(" no such value %i \n", iline);}

	} 
      for(Int_t ic=0; ic<12; ic++) {
	{
	  sprintf(buf1,"QTC%i",ic+1);
	  sprintf(buf2,"CFDvsQTC%i",ic+1);
	  sprintf(buf3,"CFDvsLED%i",ic+1);
	  sprintf(buf4,"LEDminCFD%i",ic+1);
	  sprintf(buf7,"mpd%i",ic+1);

	  hQTC[ic] = new TH1F(buf1,"QTC",(Int_t)fHistLimits[2],fHistLimits[0],fHistLimits[1]);
	  h1CFDminLED[ic] = new TH1F(buf4,"LED - CFD",(Int_t)fHistLimits[5],fHistLimits[3],fHistLimits[4]);

	  hmpd[ic] = new TH1F(buf7,"mpd",20000,-10000.0,10000.0);
	  hCFDvsQTC[ic] = new TH2F(buf2,"CFD vs	QTC",
	  	                 (Int_t)fHistLimits[8],fHistLimits[6],fHistLimits[7],
		                 (Int_t)fHistLimits[11],fHistLimits[9],fHistLimits[10]);
	  hCFDvsLED[ic] = new TH2F(buf3,"CFD vs LED-CFD",
	                         (Int_t)fHistLimits[14],fHistLimits[12],fHistLimits[13],
	                         (Int_t)fHistLimits[17],fHistLimits[15],fHistLimits[16]);


	}


      }

      TH1F*hEffCFD= new TH1F("hEffCFD","Effeciency",50,-0.25,24.25);
      TH1F*hEffLED= new TH1F("hEffLED","Effeciency",50,-0.25,24.25);
      TH1F*hEffQT0= new TH1F("hEffQT0","Effeciency",50,-0.25,24.25);
      TH1F*hEffQT1= new TH1F("hEffQT1","Effeciency",50,-0.25,24.25);

      
      Char_t filename[13];
      sprintf(filename,"t0%i.raw",fRunNumber);
       AliRawReader *reader = new AliRawReaderDate(filename);
     // AliRawReader *reader = new AliRawReaderFile();
       reader->LoadEquipmentIdsMap("T0map.txt");
      //    reader->RequireHeader(kFALSE);
       reader->RequireHeader(kTRUE);
       AliT0RawReader *start = new AliT0RawReader(reader);
       //  start->SetNumberOfTRM(1);
       for (Int_t i0=0; i0<105; i0++)
	 {
	   for (Int_t j0=0; j0<5; j0++) allData[i0][j0]=0; 	
	   numberOfHits[i0]=0;
	 }
       Int_t event=0;

       while (reader->NextEvent()) {
	 start->Next();
	 for (Int_t i=0; i<105; i++) {
	   for (Int_t iHit=0; iHit<5; iHit++) 
	    {
	      allData[i][iHit]= start->GetData(i,iHit);
	    }
	}
      
       if(event%1000 == 0) printf("Event:%d\n",event);

       //          if(event > 200000) break;

	    for (Int_t it = 0; it<24; it=it+2)
	      {
		for (Int_t iHit=0; iHit<5; iHit++)
		  {
		if(allData[it+25][iHit] != 0 && allData[it+26][iHit] !=0)
		  {
		    Int_t cc=it/2;
		    hQTC[cc]->Fill(allData[it+25][iHit]-allData[it+26][iHit]);
		    hmpd[cc]->Fill(allData[it+26][iHit]-allData[it+25][iHit]);
		    if(allData[cc+1][iHit] != 0 ) hCFDvsQTC[cc]->Fill(allData[it+25][iHit]-allData[it+26][iHit],allData[cc+1][iHit]-allData[0][0]);
		    if(allData[cc+1][iHit] != 0 && allData[cc+13][iHit]!=0 ) 
		      {
			hCFDvsLED[cc]->Fill(allData[cc+13][iHit]-allData[cc+1][iHit],allData[cc+1][iHit]-allData[0][0]);  
			h1CFDminLED[cc]->Fill(allData[cc+13][iHit]-allData[cc+1][iHit]);
		      }
		  }
	      }
	  }

	for (Int_t iHit=0; iHit<5; iHit++) 
	  {

	    for(Int_t ik=1; ik<105; ik++)
	      { 
		channels[ik] = -100;
		if((allData[ik][iHit] - allData[0][0]) > 0 ) 
		  {
		    numberOfHits[ik]++;
		    hChannel[ik] -> Fill(allData[ik][iHit] - allData[0][0]);
		    channels[ik] = allData[ik][iHit] - allData[0][0];
		  }
	      }
	      digitsTree->Fill();   
	  }

	event++;

      } //event


     if (event>1)
	{
	  printf("efficiency for %i events \n",event);
	  for (Int_t i0=1; i0<13;  i0++)
	    {
	      printf("%s  %f  %s  %f %s  %f  %s  %f \n",
		     names[i0].Data(), Float_t(numberOfHits[i0])/Float_t(event),
		     names[i0+13].Data(),Float_t(numberOfHits[i0+13])/Float_t(event),
		     names[i0+57].Data(),Float_t(numberOfHits[i0+57])/Float_t(event),
		     names[i0+69].Data(),Float_t(numberOfHits[i0+69])/Float_t(event));

	      hEffCFD->Fill(i0,Float_t(numberOfHits[i0])  / Float_t(event));
	      hEffLED->Fill(i0,Float_t(numberOfHits[i0+13]) / Float_t(event));
	      hEffCFD->Fill(i0+12,Float_t(numberOfHits[i0+57]) /Float_t(event));
	      hEffLED->Fill(i0+12,Float_t(numberOfHits[i0+69]) /Float_t(event));
	    }

	  for (Int_t i0=0; i0<24;  i0=i0+2)
	    {
	      hEffQT1->Fill(i0, Float_t (numberOfHits[i0+25]) / Float_t(event));
	      hEffQT0->Fill(i0, Float_t (numberOfHits[i0]+26) / Float_t(event));
	      hEffQT1->Fill((i0+12), Float_t (numberOfHits[i0]+81) /  Float_t(event));
	      hEffQT0->Fill((i0+12), Float_t (numberOfHits[i0]+82) /  Float_t(event));

	      printf("%s  %f  %s  %f %s  %f  %s  %f \n",
		     names[i0+25].Data(), Float_t(numberOfHits[i0+25])/Float_t(event),
		     names[i0+26].Data(),Float_t(numberOfHits[i0+26])/Float_t(event),
		     names[i0+81].Data(),Float_t(numberOfHits[i0+81])/Float_t(event),
		     names[i0+82].Data(),Float_t(numberOfHits[i0+82])/Float_t(event));
	    }
	}	      


      Char_t filehist[40]; 
     sprintf(filehist,"t0tree%i.root",fRunNumber);
     //  sprintf(filehist,"test.root",runNumber);
     TFile *hist = new TFile(filehist,"RECREATE");
       hist->cd();
      //  digitsTree->Write("",TObject::kOverwrite);

      hEffCFD->Write();
      hEffLED->Write();
      hEffQT0->Write();
      hEffQT1->Write();

      for(Int_t ik=0; ik<105; ik++)	hChannel[ik] ->Write();

      for (Int_t i=0; i<12; i++)
	{
	  hQTC[i]->Write();
	  hmpd[i]->Write();
	  hCFDvsQTC[i]->Write();
	  hCFDvsLED[i]->Write();
	  h1CFDminLED[i]->Write();
	}
      
}
