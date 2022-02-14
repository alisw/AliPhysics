/* $Id$ */
////////////////////////////////////////////////////////////////////////////
// Class AliRDHFDecisionTree
//
// Base decision tree class. Support multiple of variables for classification
// applied in reconstructed HF vertices. 
//
// Author: M. Cai, cai.mengke@cern.ch
////////////////////////////////////////////////////////////////////////////

#include <Riostream.h>

#include "AliRDHFDecisionTree.h"

using std::cout;
using std::endl;
//~ using namespace AliRDHFDecisionTree;

/// \cond CLASSIMP
ClassImp(AliRDHFDecisionTree);
/// \endcond

//--------------------------------------------------------------------------
AliRDHFDecisionTree::AliRDHFDecisionTree() :
 AliVCuts("DecisionTree","DT"),
 fNodes("AliRDHFDTNode",20),
 fRootNode(0x0),
 fmaxDepth(0)
{
  //
  // Empty Constructor
  //
}
//--------------------------------------------------------------------------
AliRDHFDecisionTree::AliRDHFDecisionTree(AliRDHFDTNode *node) :
 AliVCuts("DecisionTree","DT"),
 fNodes("AliRDHFDTNode",20),
 fRootNode(0x0),
 fmaxDepth(0)
{
  //
  // Default Constructor
  //
  fRootNode = AddRootNode(node);
}
//--------------------------------------------------------------------------
AliRDHFDecisionTree::AliRDHFDecisionTree(const AliRDHFDecisionTree &source) :
  AliVCuts(source),
  fNodes(*((TClonesArray*)source.fNodes.Clone())),
  fRootNode(source.fRootNode),
  fmaxDepth(source.fmaxDepth)
{
  //
  // Copy constructor
  //

}
//--------------------------------------------------------------------------
AliRDHFDecisionTree::~AliRDHFDecisionTree() {
  //
  // Default Destructor
  //
  fNodes.Delete();
}

//--------------------------------------------------------------------------
AliRDHFDTNode *AliRDHFDecisionTree::AddNode(AliRDHFDTNode *node)
{
	if(node->GetNodeType()==AliRDHFDTNode::kRoot) return AddRootNode(node); 
	else if(!GetRootNode()){
		std::cout << "WARNING: Root node not set yet! Please add the root node first." << endl;
		return 0x0;
	}
	else return new(fNodes[fNodes.GetEntriesFast()]) AliRDHFDTNode(*node);
}
//--------------------------------------------------------------------------
AliRDHFDTNode *AliRDHFDecisionTree::GetRootNode()
{
	if(fRootNode) return fRootNode;
	else{
		std::cout << "ERROR: Root node was not set yet!" <<endl;
		return 0;
	}
}
//--------------------------------------------------------------------------
// External call for the decision
Int_t AliRDHFDecisionTree::Decision(const std::vector<Double_t>& inputValues )
{
	//~ if(GetRootNode()!=GetNode(0)){
		//~ std::cout<<"ERROR: Root node was not at the first element in Array!"<<endl;
		//~ return 0;
	//~ }
	// Start from the root node
	//~ else
	return Decision(inputValues, 0);
}
//--------------------------------------------------------------------------
// Internal recursive all for the decision
Int_t AliRDHFDecisionTree::Decision( const std::vector<Double_t>& inputValues, Int_t currentNodeInd)
{
	AliRDHFDTNode *currentnode = GetNode(currentNodeInd);
	if(currentnode->GetNodeType()==AliRDHFDTNode::kNull) return 0;
	else if(currentnode->GetNodeType()==AliRDHFDTNode::kSignal) return 1;
	else if(currentnode->GetNodeType()==AliRDHFDTNode::kBkg) return -1;
	// Internal node
	Int_t result = currentnode->Decision(inputValues);
	if(result==1) return Decision(inputValues,currentnode->GetRNodeInd()); // Go right
	else if(result==-1) return Decision(inputValues,currentnode->GetLNodeInd()); // Go left
	else{
		std::cout<<"ERROR: Some problem with the node links..."<<endl;
		return 0;
	}
}
//--------------------------------------------------------------------------
Bool_t AliRDHFDecisionTree::VerifyNodes()
{
	Int_t nBadIntNodes(0), nBadLinkNodes(0);
	if(GetNNodes()<3){
		std::cout<<"WARNING: No. of nodes less than 3..."<<endl;
		return kFALSE;
	}
    else for(Int_t i=0;i<GetNNodes();i++){
		AliRDHFDTNode *currentnode = GetNode(i);
		// Integrity check
		if(!currentnode->CheckNodeIntegrity()) nBadIntNodes++;
		// Like check
		if(currentnode->GetNodeType()==AliRDHFDTNode::kRoot){
			AliRDHFDTNode *lnode = GetNode(currentnode->GetLNodeInd());
			AliRDHFDTNode *rnode = GetNode(currentnode->GetRNodeInd());
			if(!(lnode->GetMNodeInd()==i&&rnode->GetMNodeInd()==i)) nBadLinkNodes++;
		}
		if(currentnode->GetNodeType()==AliRDHFDTNode::kMed){
			AliRDHFDTNode *mnode = GetNode(currentnode->GetMNodeInd());
			AliRDHFDTNode *lnode = GetNode(currentnode->GetLNodeInd());
			AliRDHFDTNode *rnode = GetNode(currentnode->GetRNodeInd());
			if(!(lnode->GetMNodeInd()==i&&rnode->GetMNodeInd()==i&&(mnode->GetLNodeInd()==i||mnode->GetRNodeInd()==i))) nBadLinkNodes++;
		}
		if(currentnode->GetNodeType()==AliRDHFDTNode::kSignal){
			AliRDHFDTNode *mnode = GetNode(currentnode->GetMNodeInd());
			if(mnode->GetLNodeInd()!=i) nBadLinkNodes++;
		}
		if(currentnode->GetNodeType()==AliRDHFDTNode::kBkg){
			AliRDHFDTNode *mnode = GetNode(currentnode->GetMNodeInd());
			if(mnode->GetRNodeInd()!=i) nBadLinkNodes++;
		}
	}
	std::cout<<"INFO: Node Check - " <<nBadIntNodes<<" nodes have bad integrity, "<<nBadLinkNodes<<" nodes have bad link."<<endl;
	return nBadIntNodes+nBadLinkNodes>0?kFALSE:kTRUE;
}
//--------------------------------------------------------------------------
Int_t AliRDHFDecisionTree::TraversalNodes(Int_t startNodeInd)
{
	if(!GetNode(startNodeInd)->CheckNodeIntegrity()) return startNodeInd;
	else if(GetNode(startNodeInd)->GetNodeType()==AliRDHFDTNode::kSignal||GetNode(startNodeInd)->GetNodeType()==AliRDHFDTNode::kBkg) return -1;
	else{
		Int_t result = TraversalNodes(GetNode(startNodeInd)->GetLNodeInd());
		if(result>=0) return result;
		else return TraversalNodes(GetNode(startNodeInd)->GetRNodeInd());
	}
	//return -1 if all nodes are fine
}
// DFT, not implemented yet...
//--------------------------------------------------------------------------
void AliRDHFDecisionTree::PrintAll()
{
	std::cout<<"	Decision Tree:"<<endl;
	for(Int_t i=0;i<GetNNodes();i++){
		std::cout<<i<<": "; GetNode(i)->PrintAll();
	}
	return;
}
