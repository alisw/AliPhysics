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

////////////////////////////////////////////////
//  Huffman classes for set:TPC               //
////////////////////////////////////////////////

#include <TObjArray.h>
#include "Riostream.h"
#include "TMath.h"
#include "AliTPCHuffman.h"
#include "AliTPCBuffer160.h"

ClassImp(AliTPCHNode)

AliTPCHNode::AliTPCHNode(){
  // constructor
  fLeft=0;
  fRight=0;
}

//////////////////////////////////////////////////////////////////////////////

AliTPCHNode::AliTPCHNode(Int_t sym, Double_t freq){
  // standard constructor
  fSymbol=sym;
  fFrequency=freq;
  fLeft=0;
  fRight=0;
}

//////////////////////////////////////////////////////////////////////////////

AliTPCHNode::AliTPCHNode(const AliTPCHNode &source){
  //     Copy Constructor 
  if(&source == this) return;
  this->fSymbol = source.fSymbol;
  this->fFrequency = source.fFrequency;
  this->fLeft = source.fLeft;
  this->fRight = source.fRight;
  return;
}

//////////////////////////////////////////////////////////////////////////////

AliTPCHNode& AliTPCHNode::operator=(const AliTPCHNode &source){
  //    Assignment operator
  if(&source == this) return *this;
  this->fSymbol = source.fSymbol;
  this->fFrequency = source.fFrequency;
  this->fLeft = source.fLeft;
  this->fRight = source.fRight;
  return *this;
}

//////////////////////////////////////////////////////////////////////////////

Int_t AliTPCHNode::Compare(const TObject *obj)const{
  // function called by Sort method of TObjArray
  AliTPCHNode *node=(AliTPCHNode *)obj;
  Double_t f=fFrequency;
  Double_t fo=node->fFrequency;
  if (f<fo) return 1;
  else if (f>fo) return -1;
  else return 0;
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

ClassImp(AliTPCHTable)
  
AliTPCHTable::AliTPCHTable(){
  // constructor
  fCodeLen=0;
  fCode=0;
  fHNodes=0;
  fNnodes=0;
  fNum=0;
  fVerbose=0;
}

//////////////////////////////////////////////////////////////////////////////

AliTPCHTable::AliTPCHTable(Int_t size){
  //initialise
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

AliTPCHTable::AliTPCHTable(const AliTPCHTable &source){
  //     Copy Constructor 
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
  //    Assignment operator
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
  // HTable
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
  fCodeLen[val]=len;
  return;
}

//////////////////////////////////////////////////////////////////////////////

void AliTPCHTable::SetCode(Double_t code,Int_t val){
  fCode[val]=code;
  return;
}

//////////////////////////////////////////////////////////////////////////////

void AliTPCHTable::PrintTable()const{
  cout<<"Table for Huffman coding\n";
  cout<<"  Symbol|   Code   |   Length \n";
  for (Int_t i=0;i<fSize;i++){
    if (fCodeLen[i]){
      cout.width(6);cout<<fSym[i];
      cout.width(3);cout<<"|";
      cout.width(6);cout<<hex<<(ULong_t)fCode[i]<<dec;
      cout.width(5);cout<<"|";
      cout.width(6);cout<<(ULong_t)fCodeLen[i]<<endl;  
    }//end if
  }//end for
}

//////////////////////////////////////////////////////////////////////////////

Bool_t AliTPCHTable::SpanTree(AliTPCHNode *start, ULong_t code, UChar_t len){
  // span tree
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
  // Reset number of HNodes and the HNodes array 
  if (fHNodes)  fHNodes->Clear();
  if (fNnodes)  fNnodes=0;
}

//////////////////////////////////////////////////////////////////////////////

void AliTPCHTable::ClearTable(){
  // Clear the table
  memset(fCodeLen,0,sizeof(UChar_t)*fSize);
  memset(fCode,0,sizeof(Double_t)*fSize);
}

//////////////////////////////////////////////////////////////////////////////

Int_t  AliTPCHTable::GetFrequencies(const char *fname){
  AliTPCBuffer160 buff(fname,0);
  ULong_t NumberOfWords=0;
  Int_t Val;
  while((Val=buff.GetNext())!=-1){
    fCode[Val]++;
    fNum++;
    NumberOfWords++;
  }
  cout<<"Total number of words: "<<NumberOfWords<<endl;
  //Print out the frequencies 
  /*
    for (Int_t i=0;i<fSize;i++){
    if (fCode[i])cout<<"Symbol: "<<i<<" Freq: "<<fCode[i]<<endl;
    }
    cout<<endl;
  */
  return 0;
}

Int_t AliTPCHTable::SetValFrequency(const Int_t Val,Double_t Value){
  fCode[Val]=Value;
  fNum=1;
  return 0;
}

//////////////////////////////////////////////////////////////////////////////

Int_t AliTPCHTable::SetFrequency(const Int_t Val){
  fCode[Val]++;
  fNum++;
  return 0;
}

//////////////////////////////////////////////////////////////////////////////

Int_t AliTPCHTable::StoreFrequencies(const char *fname){
  ofstream ftxt(fname);  
  for (Int_t i=0;i<fSize;i++){
    ftxt<<(ULong_t)fCode[i]<<endl;
  }
  ftxt.close();
  return 0;
}

//////////////////////////////////////////////////////////////////////////////

void AliTPCHTable::CompleteTable(Int_t k){
  Int_t max;
  ULong_t val;
  switch(k){
  case 0:
    max=700;
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
  // build Htable
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
