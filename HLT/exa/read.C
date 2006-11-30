// $Id$

/* Different small routines for reading and testing 
   the data.*/


#ifndef __CINT__
#include "AliHLTFileHandler.h"
#include "AliHLTDigitData.h"
#include "AliHLTTransform.h"
#include "AliHLTLogger.h"
#include <stdio.h>
#include <iostream.h>
#endif


void read(Char_t *path="./",Int_t min=0,Int_t max=35)
{
  AliHLTTransform::Init(path);

  for(Int_t slice=0; slice<35; slice++)
    {
      Char_t fname[256];
      sprintf(fname,"%s/digits_%d_0.raw",path,slice);
      AliHLTFileHandler *file = new AliHLTFileHandler();
      if(!file->SetBinaryInput(fname))
	{
	  cerr<<"Error opening file "<<fname<<endl;
	  return;
	}

      Int_t row[2]={0,175};
      file->Init(slice,0,row);

      UInt_t size;
      AliHLTDigitRowData *data = file->CompBinary2Memory(size);
      
      for(Int_t r=0; r<175; r++)
	{
	  UInt_t padrow=data->fRow;
	  AliHLTDigitData *dPt = (AliHLTDigitData*)data->fDigitData;
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
  if(AliHLTTransform::Init(fname,kTRUE))
  {
    cout << "created temp init file!" << endl;
  }

  AliHLTFileHandler *fileHandler = new AliHLTFileHandler();

  if(!fileHandler->SetAliInput(fname))
    {
      cerr<<"Error opening file "<<fname<<endl;
      return;
    }

  Int_t event=0;
  UInt_t nrow=0;

  for(Int_t slice=sl; slice<=sl; slice++){
    for(Int_t patch=0;patch<AliHLTTransform::GetNPatches();patch++){

      cerr<<"reading slice: "<<slice<<" patch: "<<patch<<endl;

      fileHandler->Free();
      fileHandler->Init(slice,patch);      
      AliHLTDigitRowData *data=fileHandler->AliDigits2Memory(nrow,event);
      if(!data) cerr << "Obscure error while reading data." << endl;
      cerr<<" found "<< nrow << " rows" <<endl;
    }      
  }
  fileHandler->CloseAliInput();
}

void read_pp(Char_t *path="./",Int_t min=0,Int_t max=35,Int_t ev=0)
{
  AliHLTTransform::Init(path);

  for(Int_t slice=min; slice<max; slice++)
    {
      Char_t fname[256];
      sprintf(fname,"%s/digits_%d_%d_-1.raw",path,ev,slice);
      AliHLTFileHandler *file = new AliHLTFileHandler();
      if(!file->SetBinaryInput(fname))
	{
	  cerr<<"Error opening file "<<fname<<endl;
	  return;
	}

      file->Init(slice,-1);

      UInt_t size;
      AliHLTDigitRowData *data;
      data=(AliHLTDigitRowData*)file->Allocate(); //size from binary input
      file->Binary2Memory(size,data);

      for(Int_t r=AliHLTTransform::GetFirstRow(-1); r<AliHLTTransform::GetLastRow(-1); r++)
	{

	  UInt_t padrow=data->fRow;
	  AliHLTDigitData *dPt = (AliHLTDigitData*)data->fDigitData;
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
  AliHLTFileHandler *handler = new AliHLTFileHandler();
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
