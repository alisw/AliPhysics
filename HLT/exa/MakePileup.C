//$Id$

/**
  Macro for making pileup events in pp. Remember to have the MC
  optione switched on and that the resulting digit files will 
  not be RLE so use the PileUp option when tracking.

  Compile and run it by:
  root>.L MakePileup.C++
  root>MakePileup(path,triggerevent);
  
  Remember to have a "digitfile" link to the root 
  file containing the pp events in "path". 

  Authors: Roland Bramm & Anders Vestbo & Constantin Loizides
*/

#ifndef __CINT__
#include "AliHLTLogger.h"
#include "AliHLTFileHandler.h"
#include "AliHLTDigitData.h"
#include "AliHLTTransform.h"
#include <TNtuple.h>
#include <TRandom.h>
#include <TSystem.h>
#include <stdio.h>
#include <iostream.h>
#include <time.h>

void QSort(void **dPt,Int_t first,Int_t last);
Int_t CompareDigits(void *dPt1,void *dPt2);
#endif

/* Nice script I am using for the generation of pileup events:

--------------------------------------------------------
#!/bin/bash

from=$1
to=$2

dir=/tmp/pp

mkdir -p $dir
ln -s -f /data1/AliRoot/PythiaPP_1000Events_TPConly_ev0-375_digits.root $dir/digitfile.root
cd ~/l3/level3code/exa

for i in `seq $from $to`; do 
 echo $i;
 pdir="pileup-100-$i";

 aliroot -b runMakePileup.C\(\"$dir\",0,375,\"$pdir\",101\);
done
--------------------------------------------------------

  Path should point to the path in which there is a link to the root digitsfile.
  Startev and Endev mark the beginning and ending of valid event numbers in that file.
  Pileupdir is the directory in which the binary pileuped data will be written
  Npiles specifies the number of pileups wanted.
  TriggerEvent specifies the triggerevent in the linear mode or if -1 is given the random mode.
  Good_tpc_ev_filename is the rootfile name where the offline good-tracks are stored (see fill.C)
*/

void MakePileup(Char_t *path,Int_t startev,Int_t endev,Char_t *pileupdir="pileup",Int_t npiles=25,Int_t triggerevent=-1,Char_t *good_tpc_ev_filename=0)
{
  Int_t ntotalev=endev-startev+1;

  if(ntotalev<npiles){
    cerr << "Total number of events in root file must be at least equal to the number of pp events to be piledup." <<endl;
    return;
  }
  if(npiles<1 || npiles % 2 != 1 ){
    cerr << "npiles must be at least one and a odd number." <<endl;
    return;
  }

  Int_t const minTriggerEvents=5;
  //Int_t RefOffset[25] = {30,60,90,120,150,180,210,240,270,300,330,360,-30,-60,-90,-120,-150,-180,-210,-240,-270,-300,-330,-360,-390};
  Int_t tshift=30*25/npiles;
  Int_t *Offset=new Int_t[npiles];
  Offset[0]=0;
  for(Int_t i=1;i<=(npiles-1)/2;i++)
    Offset[i] = i*tshift;
  for(Int_t i=1;i<=(npiles-1)/2;i++)
    Offset[i+(npiles-1)/2] = -i*tshift;

  Int_t *EventNumbers=new Int_t[npiles];
  if(triggerevent<0){
    Int_t evdiff=endev-startev;
    time_t t;
    Int_t iseed = (abs((12*1000+1)*time(&t)))%900000000;
    TRandom *rand=new TRandom(iseed);

    while(triggerevent<0){
      Double_t dummy=evdiff*rand->Rndm();

      if(good_tpc_ev_filename){
	TFile goodfile=TFile(good_tpc_ev_filename,"READ");
	if(!goodfile.IsOpen()){
	  cerr << "Could not open good tpc tracks file "<<good_tpc_ev_filename<<endl;
	  return;
	}
	Char_t gooddummy[1000];
	sprintf(gooddummy,"good-pptracks-%d",(Int_t)dummy+startev);
	TNtuple *ngtuppel = (TNtuple*)goodfile.Get(gooddummy);
	if(!ngtuppel) continue;
	Int_t ngoodevents=(Int_t)ngtuppel->GetEntries();
	if(ngoodevents>=minTriggerEvents) triggerevent=(Int_t)dummy+startev;
	goodfile.Close();
      } else
	triggerevent=(Int_t)dummy+startev;
    } //found a trigger event
    for(Int_t i=1;i<npiles;i++){
      Int_t ev=triggerevent;
      while(ev==triggerevent){
	Double_t dummy=evdiff*rand->Rndm()+startev;
	ev=(Int_t)dummy;
      }
      EventNumbers[i]=ev;
    }
    delete rand;
  } else {
    Int_t ev=startev;
    for(Int_t i=1;i<npiles;i++){
      if(ev==triggerevent) ev++;
      EventNumbers[i]=ev;
      ev++;
    }
  }
  EventNumbers[0]=triggerevent;
  cout << "Events: ";for(Int_t i=0;i<npiles;i++) cout << EventNumbers[i] << " "; cout << endl;
  cout << "Offsets: ";for(Int_t i=0;i<npiles;i++) cout << Offset[i] << " ";cout << endl;

  /////////////
  //exit(1);
  /////////////
	     
  UInt_t size;
  Int_t patch = -1; //All of slice.

  Int_t NEvents = npiles; //for compatibility with old version
  AliHLTFileHandler **hand=new AliHLTFileHandler*[NEvents];
  AliHLTDigitRowData **data=new AliHLTDigitRowData*[NEvents];
  AliHLTDigitData **rowData=new AliHLTDigitData*[NEvents];
  AliHLTDigitData *test;
  Char_t Carry[1024];
  Int_t DigitsTot = 0;
  Int_t tot_dig=0;
  Int_t MeVDigit;
  Int_t DigitsPerRow = 0;
  Byte_t *NData;

  const Int_t maxdigits=50000;
  AliHLTDigitData *temp[maxdigits];
  AliHLTDigitData *temp2[maxdigits];
  
  //Create 1 filehander per event:
  sprintf(Carry,"%s/digitfile.root",path);
  for(Int_t event=0; event<NEvents; event++)
    {
      hand[event] = new AliHLTFileHandler();
      if(!hand[event]->SetAliInput(Carry))
	cerr<<" Error opening file :"<<Carry<<endl;
    }
  
  //Looping over slices and the over events
  for(Int_t slice=0 ; slice<36 ; slice++){
    cout<<"\nLoading slice "<<slice<<" "<<endl;
    for(Int_t event = 0; event < NEvents ; event++){
      Int_t eventno=EventNumbers[event];
      hand[event]->Init(slice,patch);  //give merge option to FileHandler
      data[event] = hand[event]->AliAltroDigits2Memory(size,eventno,kTRUE);
      cout << "."<<flush;
    }
    cout<<" done"<<endl;

    //Find out how many digits to allocate per patch
    DigitsTot=0;
    Int_t sign = slice < 18 ? -1 : 1;
    for(Int_t row = 0 ; row < AliHLTTransform::GetNRows(patch) ; row++){
      for(Int_t event = 0; event < NEvents ; event++){
	rowData[event] = (AliHLTDigitData *)data[event]->fDigitData;
	for(UInt_t dig=0; dig<data[event]->fNDigit; dig++)
	  {
	    Int_t time = rowData[event][dig].fTime + sign*Offset[event];
	    if(time < AliHLTTransform::GetNTimeBins() && time >= 0)
	      DigitsTot++;
	  }

	hand[event]->UpdateRowPointer(data[event]);
      }
    }

    //cout << "Try to allocate : " << DigitsTot << endl;
    Int_t AllDigitsSize = sizeof(AliHLTDigitData) * DigitsTot + sizeof(AliHLTDigitRowData) * AliHLTTransform::GetNRows(patch);
    NData = new Byte_t[AllDigitsSize];
    memset(NData,0,AllDigitsSize);
    AliHLTDigitRowData *AllRowData = (AliHLTDigitRowData*)NData;
    //cout << "Allocated " << endl;
    
    //Reset the data pointers, because they changed when doing UpdateRowPointer.
    for(Int_t event=0; event<NEvents; event++)
      data[event] = (AliHLTDigitRowData*)hand[event]->GetDataPointer(size);

    //Pileup the event    
    tot_dig=0;
    for(Int_t row = 0 ; row < AliHLTTransform::GetNRows(patch) ; row++){
      DigitsPerRow = 0;
      for(Int_t event = 0; event < NEvents ; event++){
      	rowData[event] = (AliHLTDigitData *)data[event]->fDigitData;
	for(UInt_t dig=0; dig<data[event]->fNDigit; dig++)
	  {
	    Int_t time = rowData[event][dig].fTime + sign*Offset[event];
	    if(time < AliHLTTransform::GetNTimeBins() && time >= 0)
	      DigitsPerRow++;
	  }
      }
      //cout << "Found : " << DigitsPerRow << " in all Events, Row :" << row << endl;

      //Copy the digits to a temporary storage, 
      //because we have to sort them before storing them :
      MeVDigit = 0;
      for(Int_t event = 0; event < NEvents ; event++){
	rowData[event] = (AliHLTDigitData *)data[event]->fDigitData;
	for(UInt_t digit = 0; digit < data[event]->fNDigit ; digit++){

	  Int_t time = rowData[event][digit].fTime + sign*Offset[event];
	  if(time >= AliHLTTransform::GetNTimeBins() || time < 0)
	    continue;
	  
	  temp[MeVDigit] = &rowData[event][digit];
	  temp[MeVDigit]->fTime = temp[MeVDigit]->fTime + sign*Offset[event];
	  //cout << MeVDigit << " : " << digit << endl;
	  MeVDigit++;
	}
	hand[event]->UpdateRowPointer(data[event]);
      }

      //Sort the digits:
      QSort((void**)temp,0,MeVDigit);

      //Now, loop and check for overlaps:
      Int_t final_count=0;
      for(Int_t c=0; c<MeVDigit; c++)
	{

	  if(final_count>0 && temp[c]->fTime == temp[c-1]->fTime && temp[c]->fPad == temp[c-1]->fPad)
	    {
	      //this one is overlapping with the previous one:
	      temp2[final_count-1]->fCharge += temp[c]->fCharge;
	      if(temp2[final_count-1]->fCharge>=1024)
		temp2[final_count-1]->fCharge=1023; //Saturation
	      
	      Int_t testev=( (temp[c]->fTrackID[0])>>22 ) & 0x3ff;
	      if( testev == triggerevent) 
		{
		  //This digit comes from trigger event, so store the mcid
		  temp2[final_count-1]->fTrackID[0] = temp[c]->fTrackID[0];
		  temp2[final_count-1]->fTrackID[1] = temp[c]->fTrackID[1];
		  temp2[final_count-1]->fTrackID[2] = temp[c]->fTrackID[2];
		}
	    }
	  else
	    {
	      temp2[final_count] = temp[c];
	      final_count++;
	    }
	}

      tot_dig+=final_count;
      AllRowData->fRow = row;
      AllRowData->fNDigit = final_count;
      test = (AliHLTDigitData *) AllRowData->fDigitData;

      //and copy them to the new event:
      for(Int_t c=0; c<final_count; c++)
	{
	  test[c].fCharge = temp2[c]->fCharge;
	  test[c].fPad = temp2[c]->fPad;
	  test[c].fTime = temp2[c]->fTime;


	  Int_t testev=( (temp2[c]->fTrackID[0])>>22 ) & 0x3ff;

	  //Save the MCID
	  //Only store the mc corresponding to the trigger event:
	  if ( testev == triggerevent)
	    {
	      test[c].fTrackID[0] = ( (temp2[c]->fTrackID[0]) & 0x3fffff ) - 128;
	      test[c].fTrackID[1] = ( (temp2[c]->fTrackID[1]) & 0x3fffff ) - 128;
	      test[c].fTrackID[2] = ( (temp2[c]->fTrackID[2]) & 0x3fffff ) - 128;
	    }
	  else //Digit not corresponding to the triggerevent, so mark it as noise.
	    {
	      test[c].fTrackID[0]=-1;
	      test[c].fTrackID[1]=-1;
	      test[c].fTrackID[2]=-1;
	    }
	  if(c>0 && test[c].fTime == test[c-1].fTime && test[c].fPad == test[c-1].fPad)
	    cerr<<"Overlap at row "<<row<<" pad "<<(int)test[c].fPad<<" time "<<(int)test[c].fTime<<endl;
	  //cout<<"Copying back, pad "<<(int)test[c].fPad<<" time "<<(int)test[c].fTime<<" charge "<<(int)test[c].fCharge<<endl;
	}
      
      if(MeVDigit!=DigitsPerRow)
	cerr<<endl<<"Error: "<<MeVDigit<<" "<<DigitsPerRow<<endl;
      if(final_count > MeVDigit)
	cerr<<"Error; final_count "<<final_count<<" MeVDigit "<<MeVDigit<<endl;

      Byte_t *tmp = (Byte_t *)AllRowData;
      Int_t UpdateSize = sizeof(AliHLTDigitRowData) + sizeof(AliHLTDigitData)*final_count;
      tmp += UpdateSize;
      AllRowData = (AliHLTDigitRowData *) tmp; 
    }//end looping over row

    if(tot_dig>DigitsTot)
      cerr<<endl<<"Mismatching digitcount "<<tot_dig<<" "<<DigitsTot<<endl;

    //cout << "Clearing data structures. " << endl;
    for(Int_t event = 0; event < NEvents ; event++)
      {
	data[event]=0;
	hand[event]->Free();
	hand[event]->FreeDigitsTree();
      }
    
    AliHLTFileHandler *Out = new AliHLTFileHandler();
    AllRowData = (AliHLTDigitRowData*)NData;
    if(slice==0){
      if(pileupdir[0]=='/')
	sprintf(Carry,"mkdir -p %s/",pileupdir);
      else
	sprintf(Carry,"mkdir -p %s/%s/",path,pileupdir);
      gSystem->Exec(Carry);
    }
    if(pileupdir[0]=='/')
      sprintf(Carry,"%s/digits_%d_%d_%d.raw",pileupdir,triggerevent,slice,patch);
    else 
      sprintf(Carry,"%s/%s/digits_%d_%d_%d.raw",path,pileupdir,triggerevent,slice,patch);


    if(!Out->SetBinaryOutput(Carry))
      {
	cerr<<"Cannot open file: "<<Carry<<endl;
	return;
      }
    cout << "Writing to file: " << Carry <<endl;
    Out->Memory2Binary(AliHLTTransform::GetNRows(patch),AllRowData);
    Out->CloseBinaryOutput();
    delete Out;
    delete [] NData;
  }
  
  for(Int_t event = 0; event < NEvents ; event++)
    delete hand[event];

  delete[] hand;
  delete[] rowData;
  delete[] data;
}

void QSort(void **dPt,Int_t first,Int_t last)
{
  //General sorting routine. only sorting the pointers.
  AliHLTDigitData **a = (AliHLTDigitData**)dPt;
  
  AliHLTDigitData *tmp;
  int i; 
  int j;
  
  while (last - first > 1) {
    i = first;
    j = last;
    for (;;) {
      while (++i < last && CompareDigits((void*)a[i], (void*)a[first]) < 0)
	;
      while (--j > first && CompareDigits((void*)a[j], (void*)a[first]) > 0)
	;
      if (i >= j)
	break;

      tmp  = a[i];
      a[i] = a[j];
      a[j] = tmp;
    }
    if (j == first) {
      ++first;
      continue;
    }
    tmp = a[first];
    a[first] = a[j];
    a[j] = tmp;
    if (j - first < last - (j + 1)) {
      QSort((void**)a, first, j);
      first = j + 1;   // QSort(j + 1, last);
    } else {
      QSort((void**)a, j + 1, last);
      last = j;        // QSort(first, j);
    }
  }
}

Int_t CompareDigits(void *dPt1,void *dPt2)
{
  AliHLTDigitData *a = (AliHLTDigitData*)dPt1;
  AliHLTDigitData *b = (AliHLTDigitData*)dPt2;
  if(a->fPad==b->fPad && a->fTime == b->fTime) return 0;
  
  if(a->fPad<b->fPad) return -1;
  if(a->fPad==b->fPad && a->fTime<b->fTime) return -1;
  
  return 1;
}
