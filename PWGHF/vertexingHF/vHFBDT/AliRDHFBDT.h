#ifndef ALIRDHFBDT_H
#define ALIRDHFBDT_H
/* $Id$ */

///*************************************************************************
/// Class AliRDHFBDT
///
/// Boosted decision tree class. Provide simple BDT I/O and functionality
/// as classifier of S/B. Applied in reconstructed HF vertices. 
///
/// Author: M. Cai, cai.mengke@cern.ch
///*************************************************************************

#include <Riostream.h>
#include <TString.h>
#include <TClonesArray.h>
#include <TArrayD.h>

#include "AliVCuts.h"
#include "AliRDHFDTNode.h"
#include "AliRDHFDecisionTree.h"

using std::cout;
using std::endl;

class AliRDHFDTNode;
class AliRDHFDecisionTree;

class AliRDHFBDT : public AliVCuts
{
 public:
   
   AliRDHFBDT();
   AliRDHFBDT(const TString& strVarName);
   AliRDHFBDT(const AliRDHFBDT& source);  
   
   ~AliRDHFBDT();
   Int_t GetNVars() {return fNVars;}
   Int_t GetNTrees() {return fNTrees;}
   TString GetVarName(Int_t i) {return fVarName[i];}
   TString GetDesc() const {return fDescStr;}
      
   void SetNVars(Int_t i) {fNVars = i; fVarName.resize(i);}
   void SetVarName(Int_t i, TString str){
	   if(i>=fNVars||i<0) {std::cout<<"ERROR: Input index out of bound."<<endl;return;}
	   else fVarName[i]=str;
   }
   void RemoveDesc() {fDescStr="|||BoostedDecisionTree|||";}
   void AddDescLine(const TString str) {fDescStr.Append("\n");fDescStr.Append(str);}
   
   void PrintAll();
   void Description() {std::cout<<fDescStr<<endl;}
   AliRDHFDecisionTree *GetDecisionTree(Int_t i);
   Double_t GetBoostWeight(Int_t i);
   Double_t GetResponse( const std::vector<Double_t>& inputValues );
   Bool_t CompareVarName( const TString& inputVarName );
   
   AliRDHFDecisionTree *AddDecisionTree(AliRDHFDecisionTree *tree, Double_t weight);
   //~ Int_t RemoveDecisionTree();
   
   
	// Bypass the abstract classes in VCuts
   Bool_t IsSelected(TObject* obj) { return kFALSE; }
   UInt_t GetSelectionMask(const TObject* obj) { return 0; }
   TObject *GetStatistics(Option_t *) const { return 0; }

 protected:
 
   TClonesArray fDecisionTrees;
   TArrayD fDTWeights;
   
   Int_t fNTrees;
   TString fDescStr;
   Int_t fNVars;
   std::vector<TString> fVarName;
   
   /// \cond CLASSIMP
   ClassDef(AliRDHFBDT,1);  /// 
   /// \endcond
};

#endif
