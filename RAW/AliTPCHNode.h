/* Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////
//  Huffman Table associated classes for set:TPC //
///////////////////////////////////////////////////

 
#ifndef AliTPCHNODE_H
#define AliTPCHNODE_H

#include <TObject.h>

class AliTPCHNode: public TObject  {
 public:
  AliTPCHNode(); //default constructor
  AliTPCHNode(Int_t symbol, Double_t freq);
  virtual ~AliTPCHNode() {}
  AliTPCHNode(const AliTPCHNode &source); // copy constructor
  AliTPCHNode& operator=(const AliTPCHNode &source); // ass. op.

  Bool_t  IsSortable() const{return kTRUE;}
  Int_t   Compare(const TObject *obj) const;
  void    SetLeft(AliTPCHNode* point){fLeft=point;}
  void    SetRight(AliTPCHNode* point){fRight=point;}
  AliTPCHNode* GetRight()const{return fRight;}
  AliTPCHNode* GetLeft()const{return fLeft;}
  void     SetSymbol(Int_t sym){fSymbol=sym;}
  void     SetFrequency(Double_t freq){fFrequency=freq;}
  Double_t GetFrequency()const{return fFrequency;}
  Int_t    GetSymbol()const{return fSymbol;}

 private:
  Int_t         fSymbol;       // Symbols
  Double_t      fFrequency;    // Frequency of the Symbol
  AliTPCHNode   *fLeft;        // Pointer to the left son
  AliTPCHNode   *fRight;       // Pointer to the right son
  ClassDef(AliTPCHNode,1)     
};

#endif
