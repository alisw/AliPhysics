#include "TH1F.h"
#include "TH2F.h"
#include "TMap.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
/*
#include "/scratch/alla/alice/AliRoot/STEER/AliCDBManager.h"
#include "/scratch/alla/alice/AliRoot/RAW/AliRawReader.h"
#include "/scratch/alla/alice/AliRoot/T0/AliT0LookUpValue.h"
#include "/scratch/alla/alice/AliRoot/T0/AliT0Parameters.h"
#include "/scratch/alla/alice/AliRoot/T0/AliT0RawReader.h"
*/
void readLaserData(Int_t runNumber=905)
{

  TH1F *hChannel[105];  TH1F *hNumHits[105];   TH1F *hQTC[12];  
  TH2F *hCFD_QTC[12]; TH2F *hCFD_LED[12]; TH1F *h1CFD_LED[12];
  TH1F *hQTCc[12];TH1F *hmpd[12];
  Int_t allData[110][5];
  Int_t numberOfHits[105];
  
  Char_t  buf1[10], buf2[10], buf3[10], buf4[10], buf5[10], buf6[10], buf7[10];

  TTree* digitsTree = new TTree("testData","Tree of test data Digits");
  TBranch *b[106];

   Int_t channels[106];
 
  TString names[106], type;
  AliT0LookUpKey* lookkey= new AliT0LookUpKey();
  AliT0LookUpValue*  lookvalue= new AliT0LookUpValue();
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
	//	cout<<lookkey->GetKey()<<" "<<lookkey->GetChannelName()<<" trm "<<lookvalue->GetTRM()<<" tdc "<<lookvalue->GetTDC()<<" chain  "<<lookvalue->GetChain()<<" channel "<<lookvalue->GetChannel()<<endl;
	hChannel[key] = new TH1F(names[key].Data(),names[key].Data(),30000,0,30000);
	//	hitsname="xHits" + names[key];
	//	hNumHits[key] = new TH1F(hitsname.Data(),hitsname.Data(),50,-0.25,24.25);
	type =names[key] + "/I";
	b[key]=digitsTree->Branch(names[key].Data(),&channels[key], type);
       }
      else
	{cout<<iline<<" no such value "<<endl;}
	
    } 
  for(Int_t ic=0; ic<12; ic++) {
    {
      sprintf(buf1,"QTC%i",ic+1);
      sprintf(buf2,"CFD_QTC%i",ic+1);
      sprintf(buf3,"CFD_LED%i",ic+1);
      sprintf(buf4,"LED-CFD%i",ic+1);
      sprintf(buf5,"56__%i",ic+1);
      sprintf(buf6,"55__%i",ic+1);
      sprintf(buf7,"mpd%i",ic+1);
      
      
      hQTC[ic] = new TH1F(buf1,"QTC",10000,0,10000.0);
      hQTCc[ic] = new TH1F(buf6,"QTCsmall",10000,0.0,10000.0);
      hmpd[ic] = new TH1F(buf7,"mpd",20000,-10000.0,10000.0);
      //     hCFD_QTC[ic] = new TH2F(buf2,"CFD_QTC",7000,1000.5,8000.5,2000,12000.5,18000.5);
      hCFD_QTC[ic] = new TH2F(buf2,"CFD_QTC",700,1000.5,8000.5,2000,12000.5,18000.5);
      hCFD_LED[ic] = new TH2F(buf3,"CFD_LED",500,0.0,500.0,100,14600.0,14700.0);
      h1CFD_LED[ic] = new TH1F(buf4,"CFD_LED",1000,0.0,1000.0);


    }
    

  }
  //   cout<<" hist created "<<endl;
   
  TH1F*hEffCFD= new TH1F("hEffCFD","Effeciency",50,-0.25,24.25);
  TH1F*hEffLED= new TH1F("hEffLED","Effeciency",50,-0.25,24.25);
  TH1F*hEffQT0= new TH1F("hEffQT0","Effeciency",50,-0.25,24.25);
  TH1F*hEffQT1= new TH1F("hEffQT1","Effeciency",50,-0.25,24.25);
  
  
  Char_t filename[13];
  sprintf(filename,"t0%i.raw",runNumber);
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

   //    if(event > 200000) break;
    //  cout<<"!!!!  Event Number "<< event-2<<endl;
   
  	for (Int_t it = 0; it<24; it=it+2)
	  {
	    for (Int_t iHit=0; iHit<5; iHit++)
	      {
	    if(allData[it+25][iHit] != 0 && allData[it+26][iHit] !=0)
	      {
		Int_t cc=it/2;
		// if( allData[56][0]-allData[0][0] > 0) 
		hQTC[cc]->Fill(allData[it+25][iHit]-allData[it+26][iHit]);
		if( allData[55][0]-allData[0][0] > 0) hQTCc[cc]->Fill(allData[it+26][iHit]-allData[it+25][iHit]);
		hmpd[cc]->Fill(allData[it+26][iHit]-allData[it+25][iHit]);
		if(allData[cc+1][iHit] != 0 ) hCFD_QTC[cc]->Fill(allData[it+25][iHit]-allData[it+26][iHit],allData[cc+1][iHit]-allData[0][0]);
		if(allData[cc+1][iHit] != 0 && allData[cc+13][iHit]!=0 ) 
		  {
		    hCFD_LED[cc]->Fill(allData[cc+13][iHit]-allData[cc+1][iHit],allData[cc+1][iHit]-allData[0][0]);  
		    h1CFD_LED[cc]->Fill(allData[cc+13][iHit]-allData[cc+1][iHit]);
		  }
		//  cout<<allData[cc+1][iHit]-allData[0][0]<<" "<<cc<<endl;
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
      cout<<"efficiency for "<<event<<" events"<<endl;
      for (Int_t i0=1; i0<13;  i0++)
	{

	  cout<<names[i0].Data()<<" "<<Float_t(numberOfHits[i0])/Float_t(event)<<" ";
	  cout<<names[i0+13].Data()<<" "<<Float_t(numberOfHits[i0])/Float_t(event)<<" ";
	  cout<<names[i0+57].Data()<<" "<<Float_t(numberOfHits[i0])/Float_t(event)<<" ";
	  cout<<names[i0+69].Data()<<" "<<Float_t(numberOfHits[i0])/Float_t(event)<<endl;
	  
	  hEffCFD->Fill(i0,Float_t(numberOfHits[i0])  / Float_t(event));
	  hEffLED->Fill(i0,Float_t(numberOfHits[i0+13]) / Float_t(event));
	  hEffCFD->Fill(i0+12,Float_t(numberOfHits[i0+57]) /Float_t(event));
	  hEffLED->Fill(i0+12,Float_t(numberOfHits[i0+69]) /Float_t(event));
	}
      cout<<endl;
  
      for (Int_t i0=0; i0<24;  i0=i0+2)
	{
	  hEffQT1->Fill(i0, Float_t (numberOfHits[i0+25]) / Float_t(event));
	  hEffQT0->Fill(i0, Float_t (numberOfHits[i0]+26) / Float_t(event));
	  hEffQT1->Fill((i0+12), Float_t (numberOfHits[i0]+81) /  Float_t(event));
	  hEffQT0->Fill((i0+12), Float_t (numberOfHits[i0]+82) /  Float_t(event));
	  cout<<names[i0+25].Data()<<" "<<Float_t(numberOfHits[i0+25])/Float_t(event)<<" ";
	  cout<<names[i0+26].Data()<<" "<<Float_t(numberOfHits[i0+26])/Float_t(event)<<" ";
	  cout<<names[i0+81].Data()<<" "<<Float_t(numberOfHits[i0]+81)/Float_t(event)<<" ";
	  cout<<names[i0+82].Data()<<" "<<Float_t(numberOfHits[i0]+82)/Float_t(event)<<endl;


	}
    }	      
      
     
  Char_t filehist[40]; 
 sprintf(filehist,"t0tree%i.root",runNumber);
 //  sprintf(filehist,"test.root",runNumber);
 TFile *hist = new TFile(filehist,"RECREATE");
  cout<<" writing hist in file "<<filehist<<endl;

    digitsTree->Write("",TObject::kOverwrite);


  hEffCFD->Write();
  hEffLED->Write();
  hEffQT0->Write();
  hEffQT1->Write();

  for(Int_t ik=0; ik<105; ik++)	hChannel[ik] ->Write();
   
  for (i=0; i<12; i++)
    {
      hQTC[i]->Write();
      hQTCc[i]->Write();
      hmpd[i]->Write();
      hCFD_QTC[i]->Write();
      hCFD_LED[i]->Write();
      h1CFD_LED[i]->Write();
    }

  cout<<" hist in file"<<endl;
 
}
 
