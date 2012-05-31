//////////////////////////////////////////////////////////////////////
// Author: Henrik Tydesjo                                           //
// Implementation of the nodes to put in the integer map            //
// (AliITSIntMap).                                                  //
//////////////////////////////////////////////////////////////////////   

#include "AliITSIntMapNode.h"

AliITSIntMapNode::AliITSIntMapNode():
  fKey(0),
  fVal(0),
  fLeft(NULL),
  fRight(NULL)
{}

AliITSIntMapNode::AliITSIntMapNode(Int_t key, Int_t val, AliITSIntMapNode* left, AliITSIntMapNode* right):
  fKey(key),
  fVal(val),
  fLeft(left),
  fRight(right)
{}

AliITSIntMapNode::AliITSIntMapNode(const AliITSIntMapNode& obj):
  fKey(obj.fKey),
  fVal(obj.fVal),
  fLeft(obj.fLeft),
  fRight(obj.fRight)
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
    fLeft = obj.fLeft;
    fRight = obj.fRight;
  }
  return *this;
}

