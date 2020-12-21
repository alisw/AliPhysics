/* $Id$ */
///////////////////////////////////////////////////////////////////////////
// Class AliRDHFNode
//
// Base node class for binary decision tree. One topological cut is applied
// for each node.
//
// Author: M. Cai, cai.mengke@cern.ch
///////////////////////////////////////////////////////////////////////////

#include <Riostream.h>

#include "AliRDHFDTNode.h"

using std::cout;
using std::endl;

/// \cond CLASSIMP
ClassImp(AliRDHFDTNode);
/// \endcond

// Empty constructor
AliRDHFDTNode::AliRDHFDTNode():
  AliVCuts("RDHFDTNode","Node"),
  fLNodeInd(-1),
  fRNodeInd(-1),
  fMNodeInd(-1),
  fSelector(-1),
  fCutValue(0.),
  fCutType(kNo),
  fNodeType(kNull),
  fPurity(0.),
  fResponse(0.)
{
	// ~
}

// Standard constructor
AliRDHFDTNode::AliRDHFDTNode( Int_t LNodeInd, Int_t RNodeInd, Int_t MNodeInd, Int_t selector, 
							  Double_t cutValue, ECutType cutType, ENodeType nodeType, 
							  Double_t purity, Double_t Response):
  AliVCuts("RDHFDTNode","Node"),
  fLNodeInd(LNodeInd),
  fRNodeInd(RNodeInd),
  fMNodeInd(MNodeInd),
  fSelector(selector),
  fCutValue(cutValue),
  fCutType(cutType),
  fNodeType(nodeType),
  fPurity(purity),
  fResponse(Response)
{
	// ~
}

// Copy constructor
AliRDHFDTNode::AliRDHFDTNode(const AliRDHFDTNode& source):
  AliVCuts(source),
  fLNodeInd(source.fLNodeInd),
  fRNodeInd(source.fRNodeInd),
  fMNodeInd(source.fMNodeInd),
  fSelector(source.fSelector),
  fCutValue(source.fCutValue),
  fCutType(source.fCutType),
  fNodeType(source.fNodeType),
  fPurity(source.fPurity),
  fResponse(source.fResponse)
{
	// ~
}

// Assignment constructor
AliRDHFDTNode &AliRDHFDTNode::operator=(const AliRDHFDTNode& source)
{
  if(&source == this) return *this;
  AliVCuts::operator=(source);
  fLNodeInd=source.fLNodeInd;
  fRNodeInd=source.fRNodeInd;
  fMNodeInd=source.fMNodeInd;
  fSelector=source.fSelector;
  fCutValue=source.fCutValue;
  fCutType=source.fCutType;
  fNodeType=source.fNodeType;
  fPurity=source.fPurity;
  fResponse=source.fResponse;
  //~ PrintAll();
  return *this;
}

AliRDHFDTNode::~AliRDHFDTNode()
{
   //~ if (fLeft  != NULL) delete fLeft;
   //~ if (fRight != NULL) delete fRight;
}
   
//_______________________________________________________________________

Bool_t AliRDHFDTNode::CheckNodeIntegrity()
{
   if(GetCutType()==kNo||GetNodeType()==kNull){
	   std::cout << "ERROR: Invalid node/cut type!" <<endl;
	   return kFALSE;
   }
   else if(GetNodeType()==kRoot&&(GetMNodeInd()!=-1)){
	   std::cout << "ERROR: Invalid root node!" <<endl;
	   return kFALSE;
   }
   else if(GetNodeType()==kMed&&!(GetMNodeInd()>=0&&GetLNodeInd()>0&&GetRNodeInd()>0)){
	   std::cout << "ERROR: Invalid internal node!" <<endl;
	   return kFALSE;
   }
   else if((GetNodeType()==kBkg||GetNodeType()==kSignal)&&!(GetMNodeInd()>=0&&GetLNodeInd()<0&&GetRNodeInd()<0)){
	   std::cout << "ERROR: Invalid top node!" <<endl;
	   return kFALSE;
   }
   return kTRUE; // So far no problem (not guaranteed)
}

//_______________________________________________________________________

Int_t AliRDHFDTNode::Decision( const std::vector<Double_t>& inputValues )
{
   Int_t result(0);
   if(!CheckNodeIntegrity()) return result;
   if(inputValues[GetSelector()] > GetCutValue()) result = 1;
   else result = -1;
   if (GetCutType() == kGT) return result;
   else return -result;
}

void AliRDHFDTNode::PrintAll()
{
	std::cout<<"Node: L("<<GetLNodeInd()<<"), R("<<GetRNodeInd()<<"), M("<<GetMNodeInd()<<"), Cut(v_"<<GetSelector()<<"_"<<GetCutType()<<"_"<< GetCutValue()<<"), Type("<<GetNodeType()<<")"<<endl;
	return;
}

