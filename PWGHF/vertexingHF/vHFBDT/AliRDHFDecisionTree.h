#ifndef ALIRDHFDECISIONTREE_H
#define ALIRDHFDECISIONTREE_H

///*************************************************************************
/// Class AliRDHFDecisionTree
///
/// Base decision tree class. Support multiple of variables for classification
/// applied in reconstructed HF vertices. 
///
/// Author: M. Cai, cai.mengke@cern.ch
///*************************************************************************

#include <TString.h>
#include <TClonesArray.h>

#include "AliVCuts.h"
#include "AliRDHFDTNode.h"

class TClonesArray;
class AliRDHFDTNode;

class AliRDHFDecisionTree : public AliVCuts
{
 public:
 
    AliRDHFDecisionTree();
    AliRDHFDecisionTree(AliRDHFDTNode *node);
    AliRDHFDecisionTree(const AliRDHFDecisionTree& source);  
    
    virtual ~AliRDHFDecisionTree();
       
    AliRDHFDTNode *GetNode(Int_t i) { return (AliRDHFDTNode *)fNodes.UncheckedAt(i); }
    AliRDHFDTNode *GetRootNode();
    
    Int_t GetNNodes() {return fNodes.GetEntriesFast();}
    void PrintAll();
    
    AliRDHFDTNode *AddRootNode(AliRDHFDTNode *node) {fRootNode = new(fNodes[0]) AliRDHFDTNode(*node); return fRootNode; }
    AliRDHFDTNode *AddNode(AliRDHFDTNode *node);
    Bool_t VerifyNodes();
    Int_t TraversalNodes(Int_t startNodeInd);
    
    Int_t Decision( const std::vector<Double_t>& inputValues);
    
    // Bypass the abstract classes in VCuts
   Bool_t IsSelected(TObject* obj) { return kFALSE; }
   UInt_t GetSelectionMask(const TObject* obj) { return 0; }
   TObject *GetStatistics(Option_t *) const { return 0; }
    
 protected:
    AliRDHFDecisionTree& operator=(const AliRDHFDecisionTree& source);
    Int_t Decision( const std::vector<Double_t>& inputValues, Int_t currentNodeInd);
 
    TClonesArray fNodes;
    AliRDHFDTNode *fRootNode; // We should keep the root node as the first element of fNodes
    UInt_t fmaxDepth;
    
    /// \cond CLASSIMP
    ClassDef(AliRDHFDecisionTree,1);
    /// \endcond
};

#endif
