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

#include "Riostream.h"
#include "TObjArray.h"
#include "AliTPCBuffer.h"
#include "AliSimDigits.h"

//#include "TFile.h"
//#include "TTree.h"

ClassImp(AliTPCBuffer)
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
AliTPCBuffer::AliTPCBuffer(const char* fileName){
  f.open("AliTPCDDL.dat",ios::binary|ios::out);
  // fout=new TFile(fileName,"recreate");
  // tree=new TTree("tree","Values");
  NumberOfDigits=0;
  fVerbose=0;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
AliTPCBuffer::~AliTPCBuffer(){
  f.close();
  //delete tree;
  //delete fout;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
AliTPCBuffer::AliTPCBuffer(const AliTPCBuffer &source){
  // Copy Constructor
  return;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
AliTPCBuffer& AliTPCBuffer::operator=(const AliTPCBuffer &source){
  //Assigment operator
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
	NumberOfDigits++;
      	break;
      }//end case 0
      case 1:{
	  if((Pad>=minPad)&&(Pad<=maxPad)){
	    tree->Fill();
	    NumberOfDigits++;
	  }
	break;
      }//end case 1
      case 2:{
	if((Pad<minPad)||(Pad>maxPad)){
	  tree->Fill();
	  NumberOfDigits++;
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
  //flag=0 the whole row is written intto the file
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
  do{
    data.Dig=digrow->CurrentDigit(); //adc
    data.Time=digrow->CurrentRow(); //time
    data.Pad =digrow->CurrentColumn(); // pad 
    if(data.Dig>eth){
      switch (flag){
      case 0:{
	NumberOfDigits++;
	f.write((char*)(&data),sizeof(data));
	break;
      }//end case 0
      case 1:{
	if((data.Pad>=minPad)&&(data.Pad<=maxPad)){
	  f.write((char*)(&data),sizeof(data));
	  NumberOfDigits++;
	}
	break;
      }//end case 1
      case 2:{
	if((data.Pad<minPad)||(data.Pad>maxPad)){
	  f.write((char*)(&data),sizeof(data));
	  NumberOfDigits++;
	}
	break;
      }//end case 2
      };//end switch
    }//end if
  }while (digrow->Next());
  return;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
