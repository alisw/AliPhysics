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
/* $Id:*/
////////////////////////////////////////////////
//  Huffman classes for set:TPC               //
////////////////////////////////////////////////
//This file contains two classes and it implements 
//the Huffman algorithm for creating tables
//used in the compression phase.
//The class AliTPCHNode represents a node of the Huffman tree, while
//the class AliTPCHTable represents a compression table

#include <TObjArray.h>
#include <TMath.h>
#include "AliAltroBuffer.h"
#include "AliTPCHNode.h"
#include "AliTPCHTable.h"

ClassImp(AliTPCHTable)
  
AliTPCHTable::AliTPCHTable(){
  //Constructor
  fCodeLen=0;
  fCode=0;
  fHNodes=0;
  fNnodes=0;
  fNum=0;
  fVerbose=0;
}

//////////////////////////////////////////////////////////////////////////////

AliTPCHTable::AliTPCHTable(Int_t size){
  //Initialization
  fSize=size;
  fCodeLen = new UChar_t[fSize]; 
  fCode = new Double_t[fSize]; 
  fHNodes = new TObjArray;
  fNnodes=0;
  fNum=0;
  fVerbose=0;
  fSym= new Short_t[fSize];
  for (Short_t i=0;i<fSize;i++) {
    fSym[i]=i;
  }//end for
  ClearTable(); 
}

//////////////////////////////////////////////////////////////////////////////

AliTPCHTable::AliTPCHTable(const AliTPCHTable &source)
  :TObject(source){
  //Copy Constructor 
  if(&source == this) return;
  this->fSize = source.fSize;
  this->fCodeLen = source.fCodeLen;
  this->fCode = source.fCode;
  this->fSym = source.fSym;
  this->fHNodes = source.fHNodes;
  this->fNnodes = source.fNnodes;
  this->fNum=source.fNum;
  this->fVerbose=source.fVerbose;
  return;
}

//////////////////////////////////////////////////////////////////////////////

AliTPCHTable& AliTPCHTable::operator=(const AliTPCHTable &source) {
  //Assignment operator
  if(&source == this) return *this;
  this->fSize = source.fSize;
  this->fCodeLen = source.fCodeLen;
  this->fCode = source.fCode;
  this->fSym = source.fSym;
  this->fHNodes = source.fHNodes;
  this->fNnodes = source.fNnodes;
  this->fNum=source.fNum;
  this->fVerbose=source.fVerbose;
  return *this;
}

//////////////////////////////////////////////////////////////////////////////

AliTPCHTable::~AliTPCHTable(){
  //HTable destructor
  if(fVerbose)
    cout<<"HTable destructor !\n";
  if (fCodeLen) delete[] fCodeLen;
  if (fCode) delete [] fCode;
  if (fHNodes) {
    fHNodes->Delete(); //Clear out the collection and deletes the removed objects
    delete fHNodes;
  }//end if
}

//////////////////////////////////////////////////////////////////////////////

void AliTPCHTable::SetCodeLen(UChar_t len,Int_t val){
  //Sets codelength of "val" to the variable "len"
  fCodeLen[val]=len;
  return;
}

//////////////////////////////////////////////////////////////////////////////

void AliTPCHTable::SetCode(Double_t code,Int_t val){
  //Sets the binary code of the variable "val"
  fCode[val]=code;
  return;
}

//////////////////////////////////////////////////////////////////////////////

void AliTPCHTable::PrintTable()const{
  //This method prints a table
  cout<<"Table for Huffman coding\n";
  cout<<"  Symbol|   Code   |   Length \n";
  for (Int_t i=0;i<fSize;i++){
    if (fCodeLen[i]){
      cout.width(6);cout<<fSym[i];
      cout.width(3);cout<<"|";
      cout.width(6);cout<<hex<<(UInt_t)fCode[i]<<dec;
      cout.width(5);cout<<"|";
      cout.width(6);cout<<(UInt_t)fCodeLen[i]<<endl;  
    }//end if
  }//end for
}

//////////////////////////////////////////////////////////////////////////////

Bool_t AliTPCHTable::SpanTree(AliTPCHNode *start, UInt_t code, UChar_t len){
  //Hoffman codes are generated spanning the Huffman tree
  //In an Huffman tree any internal node has always two children
  AliTPCHNode * visited;
  visited = start;
  Int_t idx=visited->GetSymbol();
  if (!visited->GetLeft()) {
    fCode[idx] = code; 
    fCodeLen[idx] = len;
    return kTRUE;
  }//end if
  SpanTree(visited->GetLeft(), code << 1, len + 1);
  SpanTree(visited->GetRight(), code << 1 | 0x0001, len + 1);
  return kTRUE;
}

//////////////////////////////////////////////////////////////////////////////

void AliTPCHTable::ResetHNodes(){
  //Reset number of HNodes and the HNodes array 
  if (fHNodes)  fHNodes->Clear();
  if (fNnodes)  fNnodes=0;
}

//////////////////////////////////////////////////////////////////////////////

void AliTPCHTable::ClearTable(){
  //Clear the table
  memset(fCodeLen,0,sizeof(UChar_t)*fSize);
  memset(fCode,0,sizeof(Double_t)*fSize);
}

//////////////////////////////////////////////////////////////////////////////

Int_t  AliTPCHTable::GetFrequencies(const char *fname){
  //It fills the "fCode" array with the frequencies of the symbols read from the file
  AliAltroBuffer buff(fname,0);
  UInt_t numberOfWords=0;
  Int_t val;
  while((val=buff.GetNext())!=-1){
    fCode[val]++;
    fNum++;
    numberOfWords++;
  }
  cout<<"Total number of words: "<<numberOfWords<<endl;
  //Print out the frequencies 
  /*
    for (Int_t i=0;i<fSize;i++){
    if (fCode[i])cout<<"Symbol: "<<i<<" Freq: "<<fCode[i]<<endl;
    }
    cout<<endl;
  */
  return 0;
}

Int_t AliTPCHTable::SetValFrequency(Int_t Val,Double_t Value){
  //This method sets to "Value" the frequency of the symbol "Val"
  fCode[Val]=Value;
  fNum=1;
  return 0;
}

//////////////////////////////////////////////////////////////////////////////

Int_t AliTPCHTable::SetFrequency(Int_t Val){
  //It increments by one the frequency of the symbol "Val" whose frequency is 
  //stored in the fCode array
  fCode[Val]++;
  fNum++;
  return 0;
}

//////////////////////////////////////////////////////////////////////////////

Int_t AliTPCHTable::NormalizeFrequencies(){
  //This method normalized the frequencies
  //Frequencies normalization
  Double_t sum=0.;
  for (Int_t i=0; i< fSize; i++) {
    sum+=fCode[i];
  }//end for 
  if (fVerbose){
    cout<<"Frequency sum: "<<sum<<endl;
  }//end if
  if(sum!=0.){
    for (Int_t i=0; i< fSize; i++) {
      fCode[i]/=sum;
      if ((fCode[i]!=0.) && (fCode[i]<10e-20))cout<<"Frequency value very small !!! "<<fCode[i]<<endl;
    }//end for 
  }
  return 0;
}
//////////////////////////////////////////////////////////////////////////////
Int_t AliTPCHTable::StoreFrequencies(const char *fname)const{
  //It stores the frequencies in a text file
  ofstream ftxt(fname);  
  for (Int_t i=0;i<fSize;i++){
    ftxt<<fCode[i]<<endl;
  }
  ftxt.close();
  return 0;
}

//////////////////////////////////////////////////////////////////////////////

void AliTPCHTable::CompleteTable(Int_t k){
  //According to the kind of table (0..4) it associates a dummy frequency (1) to 
  //every symbols whose real frequency is zero, in a given range 0..max
  Int_t max;
  UInt_t val;
  switch(k){
  case 0:
    max=fSize;
    val=1;
    break;
  case 1:
    max=445;
    val=1;
    break;
  default:
    max=fSize;
    val=1;
    break;
  }//end switch

  for(Int_t i=0;i<max;i++){
    if(fCode[i]==0.0)fCode[i]=val;
  }
  return;
}

//////////////////////////////////////////////////////////////////////////////

Double_t  AliTPCHTable::GetEntropy()const{
  //This method calculates the value of the entropy 
  Double_t entropy=0;
  Double_t prob=0;
  for (Int_t i=0;i<fSize;i++){
    if (fCode[i]){
      prob=fCode[i]/(Double_t)fNum;
      entropy+=prob*(TMath::Log2(prob));
    }
  }
  return -entropy;
}

//////////////////////////////////////////////////////////////////////////////

Int_t AliTPCHTable::BuildHTable(){
  //It builds a Huffman tree 
  if(GetWordsNumber()){
    for (Int_t i=0; i< fSize; i++) {
      if (fCode[i] > 0.){
	fNnodes++;
	//cout<< "Symbol:"<<i<<" Freq:"<<fCode[i]<<endl;
	fHNodes->Add(new AliTPCHNode(i,fCode[i]));
      }//end if
    }//end for
    Int_t nentries=fHNodes->GetEntriesFast();  
    //cout<<"Number of symbols: "<<nentries<<endl;
    Int_t nindex=nentries-1;  
    while (nindex > 0){
      fHNodes->Sort(nindex+1);
      AliTPCHNode *aux = new AliTPCHNode(0,0);
      AliTPCHNode *node= (AliTPCHNode*)fHNodes->UncheckedAt(nindex-1);
      AliTPCHNode *node1= (AliTPCHNode*)fHNodes->UncheckedAt(nindex);
      
      aux->SetLeft(node);
      aux->SetRight(node1);
      aux->SetFrequency(node->GetFrequency() + node1->GetFrequency());
      fHNodes->RemoveAt(nindex-1);
      fHNodes->AddAt(aux,nindex-1);
      nindex--;
    }//end while
    ClearTable();  
    AliTPCHNode *start= (AliTPCHNode*)fHNodes->UncheckedAt(0);
    SpanTree(start,0,0);
  }//end if
  else{
    cout<<"Table contains 0 elements\n";
  }//end else
  return 0;
}
//////////////////////////////////////////////////////////////////////////////
