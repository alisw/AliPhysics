/**************************************************************************
 * Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
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
/* $Id$ */

// Storing digits in a binary file
// according to the DDL mapping
// To be used in Alice Data Challenges
// This class is used by AliTPCDDL.C macro
// Author: D.Favretto

#include <Riostream.h>
#include <TObjArray.h>
#include "AliTPCBuffer.h"
#include "AliSimDigits.h"

//#include "TFile.h"
//#include "TTree.h"

ClassImp(AliTPCBuffer)
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//___________________________________________________________
  AliTPCBuffer::AliTPCBuffer():TObject(),
    fVerbose(0),
    fNumberOfDigits(0),
    f()
{
  //
  // default
  //
}
//____________________________________________________________
  AliTPCBuffer::AliTPCBuffer(const char* fileName):TObject(),
    fVerbose(0),
    fNumberOfDigits(0),
    f()
{
  // Constructor
#ifndef __DECCXX
  f.open(fileName,ios::binary|ios::out);
#else
  f.open(fileName,ios::out);
#endif
  // fout=new TFile(fileName,"recreate");
  // tree=new TTree("tree","Values");

  remove("TPCdigits.txt");
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
AliTPCBuffer::~AliTPCBuffer(){
  // The destructor closes the IO stream
  f.close();
  //delete tree;
  //delete fout;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
AliTPCBuffer::AliTPCBuffer(const AliTPCBuffer &source):TObject(source),
    fVerbose(0),
    fNumberOfDigits(0),
    f()
{
  // Copy Constructor
  this->fVerbose=source.fVerbose;
  return;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
AliTPCBuffer& AliTPCBuffer::operator=(const AliTPCBuffer &source){
  //Assigment operator
  this->fVerbose=source.fVerbose;
  return *this;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
void AliTPCBuffer::WriteRow(Int_t eth,AliSimDigits *digrow,Int_t minPad,Int_t maxPad,Int_t flag,Int_t sec,Int_t SubSec,Int_t row){
  //flag=0 the whole row is written to the root file
  //flag=1 only value in the range [minPad,MaxPasd] are written to the root file
  //flag=2 complementary case of 1
  Int_t Pad;
  Int_t Dig;
  Int_t Time;
  tree->Branch("sec",&sec,"sec/I");
  tree->Branch("SubSec",&SubSec,"SubSec/I");
  tree->Branch("row",&row,"row/I");
  tree->Branch("Pad",&Pad,"Pad/I");
  tree->Branch("Dig",&Dig,"Dig/I");
  tree->Branch("Time",&Time,"Time/I");
  digrow->First();
  do{
    Dig=digrow->CurrentDigit(); //adc
    Time=digrow->CurrentRow(); //time
    Pad =digrow->CurrentColumn(); // pad 
    //    cout<<"Sec "<<sec<<" row "<<row<<" Pad "<<Pad<<" Dig "<<Dig<<" Time "<<Time<<endl; 
    if(Dig>eth){
      switch (flag){
      case 0:{
	tree->Fill();
	fNumberOfDigits++;
      	break;
      }//end case 0
      case 1:{
	  if((Pad>=minPad)&&(Pad<=maxPad)){
	    tree->Fill();
	    fNumberOfDigits++;
	  }
	break;
      }//end case 1
      case 2:{
	if((Pad<minPad)||(Pad>maxPad)){
	  tree->Fill();
	  fNumberOfDigits++;
	}
	break;
      }//end case 2
      };//end switch
    }//end if
  }while (digrow->Next());
  tree->Write();
  return;
}
*/
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void AliTPCBuffer::WriteRowBinary(Int_t eth,AliSimDigits *digrow,Int_t minPad,Int_t maxPad,Int_t flag,Int_t sec,Int_t SubSec,Int_t row){
  //It writes TPC digits as par the flag specifications. Being called by AliTPCDDL.C
  //flag=0 the whole row is written into the file
  //flag=1 only value in the range [minPad,MaxPasd] are written into the file
  //flag=2 complementary case of 1

  struct DataPad{
    Int_t Sec;
    Int_t SubSec;
    Int_t Row;
    Int_t Pad;
    Int_t Dig;
    Int_t Time;
  };
  DataPad data;
  data.Sec=sec;
  data.SubSec=SubSec;
  data.Row=row;
  digrow->First();
  Int_t padID=-1;
  Int_t ddlNumber=0;
  ofstream ftxt;
  if (fVerbose==2){
    ftxt.open("TPCdigits.txt",ios::app);
    if(sec<36)
      ddlNumber=sec*2+SubSec;
    else
      ddlNumber=72+(sec-36)*4+SubSec;
  }//end if
  do{
    data.Dig=digrow->CurrentDigit();    //adc
    data.Time=digrow->CurrentRow();     //time
    data.Pad =digrow->CurrentColumn();  // pad 
    if(fVerbose==2)
      if (padID!=data.Pad){
	ftxt<<"S:"<<data.Sec<<" DDL:"<<ddlNumber<<" R:"<<data.Row<<" P:"<<data.Pad<<endl;
	padID=data.Pad;
      }//end if
    if(data.Dig>eth){
      switch (flag){
      case 0:{
	fNumberOfDigits++;
	f.write((char*)(&data),sizeof(data));
	if(fVerbose==2)
	  ftxt<<"A:"<<data.Dig<<" T:"<<data.Time<<endl; 
	break;
      }//end case 0
      case 1:{
	if((data.Pad>=minPad)&&(data.Pad<=maxPad)){
	  f.write((char*)(&data),sizeof(data));
	  if(fVerbose==2)
	    ftxt<<"A:"<<data.Dig<<" T:"<<data.Time<<endl; 
	  fNumberOfDigits++;
	}
	break;
      }//end case 1
      case 2:{
	if((data.Pad<minPad)||(data.Pad>maxPad)){
	  f.write((char*)(&data),sizeof(data));
	  if(fVerbose==2)
	    ftxt<<"A:"<<data.Dig<<" T:"<<data.Time<<endl; 
	  fNumberOfDigits++;
	}
	break;
      }//end case 2
      };//end switch
    }//end if
  }while (digrow->Next());
  if (fVerbose==2)
    ftxt.close();
  return;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
