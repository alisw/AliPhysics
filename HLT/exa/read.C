// $Id$

/* Different small routines for reading and testing 
   the data.*/


#ifndef __CINT__
#include "AliL3FileHandler.h"
#include "AliL3DigitData.h"
#include "AliL3Transform.h"
#include "AliL3Logger.h"
#include <stdio.h>
#include <iostream.h>
#endif


void read(Char_t *path="./",Int_t min=0,Int_t max=35)
{
  AliL3Transform::Init(path);

  for(Int_t slice=0; slice<35; slice++)
    {
      Char_t fname[256];
      sprintf(fname,"%s/digits_%d_0.raw",path,slice);
      AliL3FileHandler *file = new AliL3FileHandler();
      if(!file->SetBinaryInput(fname))
	{
	  cerr<<"Error opening file "<<fname<<endl;
	  return;
	}

      Int_t row[2]={0,175};
      file->Init(slice,0,row);

      UInt_t size;
      AliL3DigitRowData *data = file->CompBinary2Memory(size);
      
      for(Int_t r=0; r<175; r++)
	{
	  UInt_t padrow=data->fRow;
	  AliL3DigitData *dPt = (AliL3DigitData*)data->fDigitData;
	  cout<<"padrow "<<padrow<<" ndigits "<<data->fNDigit<<endl;
	  
	  for(Int_t d=0; d<(Int_t)data->fNDigit; d++)
	    {
	      if(d>0 && dPt[d].fPad == dPt[d-1].fPad && dPt[d].fTime == dPt[d-1].fTime)
		cout<<"Slice "<<slice<<" padrow "<<padrow<<" pad "<<(int)dPt[d].fPad<<" time "
		    <<(int)dPt[d].fTime<<" charge "<<(int)dPt[d].fCharge<<endl;
	    }
	  
	  file->UpdateRowPointer(data);
	}
      
      file->CloseBinaryInput();
      delete file;
    }
}

void read_ali(Char_t *fname, Int_t sl=0, Int_t sh=35)
{
  //need galice file or alirunfile.root link
  if(AliL3Transform::Init(fname,kTRUE))
  {
    cout << "created temp init file!" << endl;
  }

  AliL3FileHandler *fileHandler = new AliL3FileHandler();

  if(!fileHandler->SetAliInput(fname))
    {
      cerr<<"Error opening file "<<fname<<endl;
      return;
    }

  Int_t event=0;
  UInt_t nrow=0;

  for(Int_t slice=sl; slice<=sl; slice++){
    for(Int_t patch=0;patch<AliL3Transform::GetNPatches();patch++){

      cerr<<"reading slice: "<<slice<<" patch: "<<patch<<endl;

      fileHandler->Free();
      fileHandler->Init(slice,patch);      
      AliL3DigitRowData *data=fileHandler->AliDigits2Memory(nrow,event);
      if(!data) cerr << "Obscure error while reading data." << endl;
      cerr<<" found "<< nrow << " rows" <<endl;
    }      
  }
  fileHandler->CloseAliInput();
}

void read_pp(Char_t *path="./",Int_t min=0,Int_t max=35,Int_t ev=0)
{
  AliL3Transform::Init(path);

  for(Int_t slice=min; slice<max; slice++)
    {
      Char_t fname[256];
      sprintf(fname,"%s/digits_%d_%d_-1.raw",path,ev,slice);
      AliL3FileHandler *file = new AliL3FileHandler();
      if(!file->SetBinaryInput(fname))
	{
	  cerr<<"Error opening file "<<fname<<endl;
	  return;
	}

      file->Init(slice,-1);

      UInt_t size;
      AliL3DigitRowData *data;
      data=(AliL3DigitRowData*)file->Allocate(); //size from binary input
      file->Binary2Memory(size,data);

      for(Int_t r=AliL3Transform::GetFirstRow(-1); r<AliL3Transform::GetLastRow(-1); r++)
	{

	  UInt_t padrow=data->fRow;
	  AliL3DigitData *dPt = (AliL3DigitData*)data->fDigitData;
	  cout<<r<<" "<<"padrow "<<padrow<<" ndigits "<<data->fNDigit<<endl;
	  
	  for(Int_t d=0; d<(Int_t)data->fNDigit; d++)
	    {
	      //if((int)dPt[d]->fTime>1023)
	      cout<<"Slice "<<slice<<" padrow "<<padrow<<" pad "<<(int)dPt[d].fPad<<" time "
		  <<(int)dPt[d].fTime<<" charge "<<(int)dPt[d].fCharge<<" mc "<<(int)dPt[d].fTrackID[0]<<endl;
	    }
	  
	  file->UpdateRowPointer(data);
	}
      
      file->CloseBinaryInput();
      delete file;
    }
}

void read_event_tree(Char_t *rootfile,Int_t startev=0)
{
  AliL3FileHandler *handler = new AliL3FileHandler();
  if(!handler->SetAliInput(rootfile)){
    cerr<<" Error opening file: "<<rootfile<<endl;
    return;
  }
  
  //Looping over slices and the over events
  Int_t ev=startev;
  UInt_t size=0;
  Bool_t cont=kTRUE;
  while(cont){
    cout<<"\nEvent " << ev << " Loading slices "<<flush;
    for(Int_t slice=0 ; slice<1 ; slice++){
      handler->Init(slice,-1);
      if(!handler->AliAltroDigits2Memory(size,ev)){
	cout << " not a valid event, breaking!!!" <<endl;
	cont=kFALSE;
	break;
      }
      cout << "."<<flush;
    }
    if(cont){
      handler->Free();
      handler->FreeDigitsTree();
      cout<<" done"<<endl;
      ev++;
    }
  }

  delete handler;
}
