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
#include "TGFileDialog.h"
#include <iostream>

#include "AliT0CalibLaserData.h"

#include "AliCDBManager.h"
#include "AliRawReader.h"
#include "AliRawReaderRoot.h"
#include "AliT0LookUpValue.h"
#include "AliT0LookUpKey.h"
#include "AliT0Parameters.h"
#include "AliT0RawReader.h"

ClassImp(AliT0CalibLaserData)
  //const char *fFileName;
AliT0CalibLaserData::AliT0CalibLaserData() : TObject(),
					     fTEntry(0),
					     fFileName(" ")
  
{
//
for ( int i=0; i<30; i++ ) { fEntries[i] = NULL; fHistLimits[i] = 0.0;}
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

void AliT0CalibLaserData::ReadHistSize()
{
  //build GUI frame for reading:
  // - run number
  // - histograms rates

    
    TGMainFrame* fMain = new TGMainFrame(0,1500,1500);
    fMain->SetLayoutManager( new TGMatrixLayout(fMain,10,7) );
 
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
    fEntries[2] = new TGNumberEntry(fMain, 2500);
    fMain->AddFrame(fEntries[2]);
    fMain->AddFrame( new TGLabel(fMain, " ") );
    fMain->AddFrame( new TGLabel(fMain, " ") );
    fMain->AddFrame( new TGLabel(fMain, " ") );
  
    fMain->AddFrame( new TGLabel(fMain, "LED - CFD" ) );
    fEntries[3] = new TGNumberEntry(fMain, 0);
    fMain->AddFrame(fEntries[3]);
    fEntries[4] = new TGNumberEntry(fMain, 1000);
    fMain->AddFrame(fEntries[4]);
    fEntries[5] = new TGNumberEntry(fMain, 1000);
    fMain->AddFrame(fEntries[5]);
    fMain->AddFrame( new TGLabel(fMain, " ") );
    fMain->AddFrame( new TGLabel(fMain, " ") );
    fMain->AddFrame( new TGLabel(fMain, " ") );
    fMain->AddFrame( new TGLabel(fMain, "CFD vs QTC " ) );
 // QTC axis X
    fEntries[6] = new TGNumberEntry(fMain, 0.);
    fMain->AddFrame(fEntries[6]);
    fEntries[7] = new TGNumberEntry(fMain, 8000);
    fMain->AddFrame(fEntries[7]);
    fEntries[8] = new TGNumberEntry(fMain, 800);
    fMain->AddFrame(fEntries[8]);
// CFD axis Y 
   fEntries[9] = new TGNumberEntry(fMain, 0);
    fMain->AddFrame(fEntries[9]);
    fEntries[10] = new TGNumberEntry(fMain, 5000);
    fMain->AddFrame(fEntries[10]);
    fEntries[11] = new TGNumberEntry(fMain, 5000);
    fMain->AddFrame(fEntries[11]);
//
    fMain->AddFrame( new TGLabel(fMain, "CFD vs LED-CFD " ) );
//LED-CFD axis X
    fEntries[12] = new TGNumberEntry(fMain, 0);
    fMain->AddFrame(fEntries[12]);
    fEntries[13] = new TGNumberEntry(fMain, 1000);
    fMain->AddFrame(fEntries[13]);
    fEntries[14] = new TGNumberEntry(fMain, 1000);
    fMain->AddFrame(fEntries[14]);
// CFD axis Y
    fEntries[15] = new TGNumberEntry(fMain, 1000);
    fMain->AddFrame(fEntries[15]);
    fEntries[16] = new TGNumberEntry(fMain, 5000);
    fMain->AddFrame(fEntries[16]);
    fEntries[17] = new TGNumberEntry(fMain, 4000);
    fMain->AddFrame(fEntries[17]);
 
   fMain->AddFrame( new TGLabel(fMain, "CFD C " ) );
//CFD side C
    fEntries[18] = new TGNumberEntry(fMain, 0);
    fMain->AddFrame(fEntries[18]);
    fEntries[19] = new TGNumberEntry(fMain, 5000);
    fMain->AddFrame(fEntries[19]);
    fEntries[20] = new TGNumberEntry(fMain, 5000);
    fMain->AddFrame(fEntries[20]);
 
    fMain->AddFrame( new TGLabel(fMain, "") );
    fMain->AddFrame( new TGLabel(fMain, "") );
    fMain->AddFrame( new TGLabel(fMain, "") );

    fMain->AddFrame( new TGLabel(fMain, "CFD A " ) );
//CFD side A
    fEntries[21] = new TGNumberEntry(fMain,0);
    fMain->AddFrame(fEntries[21]);
    fEntries[22] = new TGNumberEntry(fMain, 5000);
    fMain->AddFrame(fEntries[22]);
    fEntries[23] = new TGNumberEntry(fMain, 5000);
    fMain->AddFrame(fEntries[23]);

    fMain->AddFrame( new TGLabel(fMain, "") );
    fMain->AddFrame( new TGLabel(fMain, "") );
    fMain->AddFrame( new TGLabel(fMain, "") );


    fMain->AddFrame( new TGLabel(fMain, "LED C " ) );
//LED axis X
    fEntries[24] = new TGNumberEntry(fMain, 0);
    fMain->AddFrame(fEntries[24]);
    fEntries[25] = new TGNumberEntry(fMain, 5000);
    fMain->AddFrame(fEntries[25]);
    fEntries[26] = new TGNumberEntry(fMain, 5000);
    fMain->AddFrame(fEntries[26]);

   
   for ( int i=0; i<27; i++ ) fEntries[i]->SetWidth(70);
   fMain->AddFrame( new TGLabel(fMain, " File name") );
   fTEntry = new TGTextEntry(fMain,"");
   fMain->AddFrame(fTEntry);
   fTEntry->SetWidth(80);
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
  
  fFileName = (fTEntry->GetText());
  //OpenFile();
  printf(" DoOK >>  File %s\n",fFileName);
  for( int i=0; i<27; i++ ) 
    fHistLimits[i] = fEntries[i]->GetNumber();
  
  ReadData();
}

void AliT0CalibLaserData::ReadData()
{
  // reading RAW data from test LCS
  // filling histograms
  // fillinf tree
  
  
  TH1I *hChannel[105];  TH1I *hQTC[24];  
  TH2F *hCFDvsQTC[24]; TH2F *hCFDvsLED[24]; TH1I *h1CFDminLED[24];

  Int_t allData[110][50];
  Int_t numberOfHits[105];
  
  Char_t  buf1[20], buf2[20], buf3[20], buf4[20], buf7[20];
  
  Int_t channels[106];
  
  TString names[106], type;
  AliT0LookUpKey* lookkey;//= new AliT0LookUpKey();
  AliT0LookUpValue*  lookvalue;//= new AliT0LookUpValue();
  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
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
	//	if(names[key].Contains("QT0"))
	//  hChannel[key] = new TH1F(names[key].Data(),names[key].Data(),1000,3000,5000);
	//	else
	  hChannel[key] = new TH1I(names[key].Data(),names[key].Data(),Int_t(fHistLimits[20]),fHistLimits[18],fHistLimits[19]);
	if(key >0 && key<13)
	  hChannel[key] = new TH1I(names[key].Data(),names[key].Data(),Int_t(fHistLimits[20]),fHistLimits[18],fHistLimits[19]);
	if(key >13 && key<25)
	  hChannel[key] = new TH1I(names[key].Data(),names[key].Data(),Int_t(fHistLimits[26]),fHistLimits[24],fHistLimits[25]);
	
	if(key >57 && key<69)
	  hChannel[key] = new TH1I(names[key].Data(),names[key].Data(),Int_t (fHistLimits[23]),fHistLimits[21],fHistLimits[22]);
	//	hitsname="xHits" + names[key];
	//	hNumHits[key] = new TH1F(hitsname.Data(),hitsname.Data(),50,-0.25,24.25);
      }
      else
	{printf(" no such value %i \n", iline);}
      
    } 
  for(Int_t ic=0; ic<24; ic++) {
    {
      sprintf(buf1,"QTC%i",ic+1);
      sprintf(buf2,"CFDvsQTC%i",ic+1);
      sprintf(buf3,"CFDvsLED%i",ic+1);
      sprintf(buf4,"LEDminCFD%i",ic+1);
      sprintf(buf7,"mpd%i",ic+1);
      
      hQTC[ic] = new TH1I(buf1,"QTC",(Int_t)fHistLimits[2],fHistLimits[0],fHistLimits[1]);
      h1CFDminLED[ic] = new TH1I(buf4,"LED - CFD",(Int_t)fHistLimits[5],fHistLimits[3],fHistLimits[4]);
      
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

      
     // Char_t filename[80];
  //     sprintf(filename,"t0%i.raw",fRunNumber);
  AliRawReader *reader = new AliRawReaderRoot(fFileName);
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
	  //	  if( allData[i][iHit]>0)	  cout<<i<<" "<<iHit<<" "<<allData[i][iHit]<<endl;
	}
    }
    
     if(event%1000 == 0) 
    printf("Event:%d\n",event);
    
     //   if(event > 100000) break;
    
    for (Int_t it = 0; it<24; it=it+2)
      {
	Int_t cc=it/2;
	for (Int_t iHit=0; iHit<5; iHit++)
	  {
	    if(allData[it+25][iHit] != 0 && allData[it+26][iHit] !=0)
	      {
		hQTC[cc]->Fill(allData[it+25][iHit]-allData[it+26][iHit]);
		if(allData[cc+1][iHit] != 0 ) 
		  hCFDvsQTC[cc]->Fill(allData[it+25][iHit]-allData[it+26][iHit],
				      allData[cc+1][iHit]-allData[0][0]+5000.);
	      }
	    if(allData[cc+1][iHit] != 0 && allData[cc+13][iHit]!=0 ) 
	      {
		hCFDvsLED[cc]->Fill(allData[cc+13][iHit]-allData[cc+1][iHit],
				    allData[cc+1][iHit]-allData[0][0]+5000.);  
		h1CFDminLED[cc]->Fill(allData[cc+13][iHit]-allData[cc+1][iHit]);
	      }
	    
	  }
      }
    for (Int_t it = 24; it<48; it=it+2)
      {
	Int_t cc=(Int_t)(it/2);
	for (Int_t iHit=0; iHit<5; iHit++)
	  {
	    if(allData[it+57][iHit] != 0 && allData[it+58][iHit] !=0)
	      {
		hQTC[cc]->Fill(allData[it+57][iHit]-allData[it+58][iHit]);
		//		    hmpd[cc]->Fill(allData[it+26][iHit]-allData[it+25][iHit]);
		if(allData[cc+1][iHit] != 0 ) 
		  hCFDvsQTC[cc]->Fill(allData[it+57][iHit]-allData[it+58][iHit],
				      allData[cc+45][iHit]-allData[0][0]+5000);
	      }
	    if(allData[cc+57][iHit] != 0 && allData[cc+45][iHit]!=0 ) 
	      {
		hCFDvsLED[cc]->Fill(allData[cc+57][iHit]-allData[cc+45][iHit],
				    allData[cc+1][iHit]-allData[0][0]+5000);  
		h1CFDminLED[cc]->Fill(allData[cc+57][iHit]-allData[cc+45][iHit]);
	      }
	      
	}
      }
    
    
    for (Int_t iHit=0; iHit<5; iHit++) 
      {
	
	for(Int_t ik=1; ik<105; ik++)
	  { 
	    channels[ik] = -100;
	    if((allData[ik][iHit] - allData[0][0] +5000) != 0 &&   //!!!!! Uncomment it !!!!! and comment next line
	       //      if((allData[ik][iHit] - allData[1][0] ) != 0 &&     
	         allData[ik][iHit] >0 ) 
	      {  
		numberOfHits[ik]++;
		//	    	hChannel[ik] -> Fill(allData[ik][iHit] - allData[1][0]);            //Comment this line !!!
		hChannel[ik] -> Fill(allData[ik][iHit] - allData[0][0]+5000);
		//	cout<<" zpis'>> "<<iHit<<" "<<ik<<" "<<allData[ik][iHit] - allData[0][0]<<endl;
		//	hChannel[ik] -> Fill(allData[0][0]-allData[ik][iHit] );
		//		channels[ik] = allData[ik][iHit] - allData[0][0];
	      }
	  }
	//	      digitsTree->Fill();   
      }

    event++;
    
  } //event


  if (event>1)
    {
      printf("efficiency for %i events \n",event);
      for (Int_t i0=1; i0<13;  i0++)
	{
	  printf("%s  %f  %s  %f %s  %f  %s  %f \n ",
		 names[i0].Data(), Float_t(numberOfHits[i0])/Float_t(event),
		 names[i0+12].Data(),Float_t(numberOfHits[i0+12])/Float_t(event),
		 names[i0+56].Data(),Float_t(numberOfHits[i0+56])/Float_t(event),
		 names[i0+68].Data(),Float_t(numberOfHits[i0+68])/Float_t(event));
	  
	  hEffCFD->Fill(i0,Float_t(numberOfHits[i0])  / Float_t(event));
	  hEffLED->Fill(i0,Float_t(numberOfHits[i0+12]) / Float_t(event));
	  hEffCFD->Fill(i0+12,Float_t(numberOfHits[i0+56]) /Float_t(event));
	  hEffLED->Fill(i0+12,Float_t(numberOfHits[i0+68]) /Float_t(event));
	}
      printf("\n");      
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
  

  Char_t filehist[50]; 
  // sprintf(filehist,"t0treeDA%s",fFileName);
   sprintf(filehist,"t0tree%s",fFileName);
  printf("\n Wrote data in %s !!\n",filehist);
     TFile *hist = new TFile(filehist,"RECREATE");
     hist->cd();
     //  digitsTree->Write("",TObject::kOverwrite);
     
     hEffCFD->Write();
     hEffLED->Write();
     hEffQT0->Write();
      hEffQT1->Write();
      
      for(Int_t ik=0; ik<105; ik++)	hChannel[ik] ->Write();
      
      for (Int_t i=0; i<24; i++)
	{
	  hQTC[i]->Write();
	  hCFDvsQTC[i]->Write();
	  hCFDvsLED[i]->Write();
	  h1CFDminLED[i]->Write();
	}
      
}

void AliT0CalibLaserData::OpenFile()
{

const char *ft[]={"T0 raw files","*.root","All files","*",0,0};
  TString dir(".");
  TGFileInfo fi; fi.fFileTypes=ft; fi.fIniDir=StrDup(dir);
  new TGFileDialog(gClient->GetRoot(), 0x0, kFDOpen, &fi);
  if(!fi.fFilename) return;
  // fFileName =*( fi.fFilename);
  fFileName = fi.fFilename;
  printf(" AliT0CalibLaserData::OpenFile %s %s\n",fi.fFilename, fFileName );
  //  if(gFile){ gFile->Close(); gFile=0;}

  // gFile=TFile::Open(fi.fFilename);

}

