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

#include "AliTPCHNode.h"

ClassImp(AliTPCHNode)

AliTPCHNode::AliTPCHNode(){
  //Constructor
  fLeft=0;
  fRight=0;
}

//////////////////////////////////////////////////////////////////////////////

AliTPCHNode::AliTPCHNode(Int_t sym, Double_t freq){
  //Standard constructor
  fSymbol=sym;
  fFrequency=freq;
  fLeft=0;
  fRight=0;
}

//////////////////////////////////////////////////////////////////////////////

AliTPCHNode::AliTPCHNode(const AliTPCHNode &source)
  :TObject(source){
  //Copy Constructor 
  if(&source == this) return;
  this->fSymbol = source.fSymbol;
  this->fFrequency = source.fFrequency;
  this->fLeft = source.fLeft;
  this->fRight = source.fRight;
  return;
}

//////////////////////////////////////////////////////////////////////////////

AliTPCHNode& AliTPCHNode::operator=(const AliTPCHNode &source){
  //Assignment operator
  if(&source == this) return *this;
  this->fSymbol = source.fSymbol;
  this->fFrequency = source.fFrequency;
  this->fLeft = source.fLeft;
  this->fRight = source.fRight;
  return *this;
}

//////////////////////////////////////////////////////////////////////////////

Int_t AliTPCHNode::Compare(const TObject *obj)const{
  //Function called by Sort method of TObjArray
  AliTPCHNode *node=(AliTPCHNode *)obj;
  Double_t f=fFrequency;
  Double_t fo=node->fFrequency;
  if (f<fo) return 1;
  else if (f>fo) return -1;
  else return 0;
}

//////////////////////////////////////////////////////////////////////////////
