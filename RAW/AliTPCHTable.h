/* Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////
//  Huffman Table associated classes for set:TPC //
///////////////////////////////////////////////////

 
#ifndef AliTPCHTABLE_H
#define AliTPCHTABLE_H

#include <TObject.h>

class TObjArray;
class AliTPCHNode;

class AliTPCHTable: public TObject{ 
 public:
  AliTPCHTable(); 
  AliTPCHTable(Int_t size);
  virtual   ~AliTPCHTable();
  AliTPCHTable(const AliTPCHTable &source); // copy constructor
  AliTPCHTable& operator=(const AliTPCHTable &source); // ass. op.
  
  Int_t      Size()const {return fSize;}
  UChar_t*   CodeLen()const {return fCodeLen;}
  Double_t*  Code()const {return fCode;}
  Short_t*   Sym()const {return fSym;}
  void       SetCodeLen(UChar_t len,Int_t val);
  void       SetCode(Double_t code,Int_t val);
  TObjArray* HNodes()const {return fHNodes;}
  void       PrintTable()const;
  //This method builds the Huffman tree starting from the frequencies that are 
  //strored temporary in fCode array
  Int_t      BuildHTable();
  //This method returns the number of words stored in the fSym array
  UInt_t    GetWordsNumber()const{return fNum;}
  //This method increase by one the frequency of each value that is present
  //in the specified file
  Int_t      GetFrequencies(const char* fname);
  //This method increase by one the frequency of a given value
  Int_t      SetFrequency(Int_t Val);
  //This method stores the frequency of the symbol in a text file
  Int_t      StoreFrequencies(const char *fname)const;
  void       CompleteTable(Int_t k);
  Double_t   GetEntropy()const;
  void       SetVerbose(Int_t val){fVerbose=val;}
  //Method to set directly a frequency 
  Int_t      SetValFrequency(Int_t Val,Double_t Value);
  Int_t      NormalizeFrequencies();
 private:
  //This method executes the pre-order visit of an Huffman tree and calculates the 
  //codeword for each leaf
  Bool_t     SpanTree(AliTPCHNode*start, UInt_t code, UChar_t len);
  void       ResetHNodes();  //Reset the array fHNodes but not delete the removed objects
  void       ClearTable();   //Reset the table
  Int_t       fSize;         //size of the arrays fCodelen and fCode
  UChar_t     *fCodeLen;     //![fSize] number of bits array
  Double_t    *fCode;        //![fSize] coded symbols array
  
  Short_t    *fSym;          //![fSize] array of input symbols
  TObjArray  *fHNodes;       // array of nodes
  Int_t       fNnodes;       // number of nodes
  UInt_t      fNum;          // number of words
  Int_t       fVerbose;      // if fVerbose== 0 no output messages; fVerbose!=0 output messages are printed out 
  ClassDef(AliTPCHTable,1)   //Huffman Table object for set:TPC
};

#endif
