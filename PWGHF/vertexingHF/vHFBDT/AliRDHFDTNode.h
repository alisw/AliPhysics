#ifndef ALIRDHFDTNODE_H
#define ALIRDHFDTNODE_H

///*************************************************************************
/// Class AliRDHFNode
///
/// Base node class for binary decision tree. One topological cut is applied
/// for each node.
///
/// Author: M. Cai, cai.mengke@cern.ch
///*************************************************************************

#include "AliVCuts.h"
#include <TObject.h>

class TObject;

class AliRDHFDTNode : public AliVCuts {
	
 public:
   enum ENodeType {kNull, kRoot, kMed, kSignal, kBkg};
   enum ECutType {kNo, kGT, kLT};

   AliRDHFDTNode();
   // More reasonable BDT node which can be accessed via pointer
   AliRDHFDTNode( Int_t LNodeInd, Int_t RNodeInd, Int_t MNodeInd, Int_t selector, 
				  Double_t cutValue, ECutType cutType, ENodeType nodeType, Double_t purity, Double_t Response);
   
   // Copy			  
   AliRDHFDTNode(const AliRDHFDTNode& source);
   AliRDHFDTNode& operator=(const AliRDHFDTNode& source);	  

   virtual ~AliRDHFDTNode();
   
   virtual Int_t Decision( const std::vector<Double_t>& inputValues );
   Int_t GetRNodeInd() {return fRNodeInd; };
   Int_t GetLNodeInd() {return fLNodeInd; };
   Int_t GetMNodeInd() {return fMNodeInd; }; 

   Double_t GetPurity() 	{ return fPurity; 	} 
   ENodeType    GetNodeType() 	{ return fNodeType; }
   ECutType    GetCutType() 	{ return fCutType; }
   Double_t GetCutValue() 	{ return fCutValue;	}
   Double_t GetResponse() 	{ return fResponse;	}
   Int_t 	GetSelector() 	{ return fSelector;	}
   Bool_t CheckNodeIntegrity();
   void PrintAll();
   
   void SetLNodeInd(Int_t ind) {fLNodeInd = ind;}
   void SetRNodeInd(Int_t ind) {fRNodeInd = ind;}
   void SetMNodeInd(Int_t ind) {fMNodeInd = ind;}
   void SetSelector(Int_t sel) 	{fSelector = sel;}
   void SetCutValue(Double_t val) {fCutValue = val;}
   void SetCutType(ECutType type) {fCutType = type;}
   void SetNodeType(ENodeType type) {fNodeType = type;}
   void SetPurity(Double_t pur) {fPurity = pur;}
   void SetResponse(Double_t res) {fResponse = res;}
   
   // Bypass the abstract classes in VCuts
   Bool_t IsSelected(TObject* obj) { return kFALSE; }
   UInt_t GetSelectionMask(const TObject* obj) { return 0; }
   TObject *GetStatistics(Option_t *) const { return 0; }
   
 protected:
 
   Int_t 	fLNodeInd;
   Int_t	fRNodeInd;
   Int_t 	fMNodeInd;
   Int_t 	fSelector;
   Double_t fCutValue;
   ECutType 	fCutType;
   ENodeType  	fNodeType;
   Double_t fPurity;
   Double_t fResponse; // Temporarily kept for the regression tree
   
   /// \cond CLASSIMP
   ClassDef(AliRDHFDTNode,1);
   /// \endcond
}; 
#endif
