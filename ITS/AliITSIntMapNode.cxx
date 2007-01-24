//////////////////////////////////////////////////////////////////////
// Author: Henrik Tydesjo                                           //
// Implementation of the nodes to put in the integer map            //
// (AliITSIntMap).                                                  //
//////////////////////////////////////////////////////////////////////   

#include "AliITSIntMapNode.h"

AliITSIntMapNode::AliITSIntMapNode():
  fKey(0),
  fVal(0),
  fNext(0)
{}

AliITSIntMapNode::AliITSIntMapNode(Int_t key, Int_t val, AliITSIntMapNode* next):
  fKey(key),
  fVal(val),
  fNext(next)
{}

AliITSIntMapNode::AliITSIntMapNode(const AliITSIntMapNode& obj):
  fKey(obj.fKey),
  fVal(obj.fVal),
  fNext(obj.fNext)
{
  // copy constructor
}

AliITSIntMapNode::~AliITSIntMapNode() 
{}

AliITSIntMapNode& AliITSIntMapNode::operator=(const AliITSIntMapNode& obj) 
{
  // assignment operator
  if (this!=&obj) {
    fKey = obj.fKey;
    fVal = obj.fVal;
    fNext = obj.fNext;
  }
  return *this;
}

