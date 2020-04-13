/* $Id$ */

///////////////////////////////////////////////////////////////////////////
// Class AliRDHFBDT
//
// Boosted decision tree class. Provide simple BDT I/O and functionality
// as classifier of S/B. Applied in reconstructed HF vertices. 
//
// Author: M. Cai, cai.mengke@cern.ch
///////////////////////////////////////////////////////////////////////////

#include <Riostream.h>
#include <vector>

#include "AliRDHFBDT.h"

using std::cout;
using std::endl;

/// \cond CLASSIMP
ClassImp(AliRDHFBDT);
/// \endcond
//--------------------------------------------------------------------------
AliRDHFBDT::AliRDHFBDT() :
 AliVCuts("RDHFBDT","BDT"),
 fDecisionTrees("AliRDHFDecisionTree",20),
 fDTWeights(0),
 fNTrees(0),
 fDescStr("BDT"),
 fNVars(0),
 fVarName(0)
{
  //
  // Default Constructor
  //
  RemoveDesc();
}
//--------------------------------------------------------------------------
AliRDHFBDT::AliRDHFBDT(const TString& strVarName) :
 AliVCuts("RDHFBDT","BDT"),
 fDecisionTrees("BoostedDecisionTree",20),
 fDTWeights(0),
 fNTrees(0),
 fDescStr("BDT"),
 fNVars(0),
 fVarName(0)
{
  //
  // Copy constructor
  //
  RemoveDesc();
  TString tmpstr(strVarName);
  TObjArray *strarr= (TObjArray*)tmpstr.Tokenize(":");
  SetNVars(strarr->GetEntriesFast());
  for(Int_t i = 0; i < fNVars; i++) SetVarName(i,strarr->At(i)->GetName());
  //~ fVarName.push_back(strarr->At(i)->GetName());
}
//--------------------------------------------------------------------------
AliRDHFBDT::AliRDHFBDT(const AliRDHFBDT &source):
  AliVCuts(source),
  fDecisionTrees(*((TClonesArray*)source.fDecisionTrees.Clone())),
  fDTWeights(0),
  fNTrees(source.fNTrees),
  fDescStr(source.fDescStr),
  fNVars(source.fNVars),
  fVarName(0)
{
  //
  // assignment operator
  //
  for(Int_t i = 0; i < fNVars; i++) fVarName[i] = source.fVarName[i];
  fDTWeights=source.fDTWeights;
}
//--------------------------------------------------------------------------
AliRDHFBDT::~AliRDHFBDT() {
  //
  // Default Destructor
  //
  fDecisionTrees.Delete();
}
//---------------------------------------------------------------------------
Double_t AliRDHFBDT::GetResponse( const std::vector<Double_t>& inputValues )
{
	Double_t res = 0;
    Double_t norm = 0;
    for (Int_t i=0; i<fDecisionTrees.GetEntriesFast(); i++){
       AliRDHFDecisionTree *current = GetDecisionTree(i);
       Int_t decision = current->Decision(inputValues);
       if(!(decision==1||decision==-1)){
		   std::cout << "ERROR: Problem with the "<<i<<"-th tree's output." <<endl;
		   return -99;
	   }
       else{
           res += GetBoostWeight(i) * decision;
           norm += GetBoostWeight(i);
       }
    }
	return res /= norm;
}
//---------------------------------------------------------------------------
Bool_t AliRDHFBDT::CompareVarName( const TString& inputVarName )
{
	Bool_t result(kTRUE);
	TObjArray *inputVarNameArr = inputVarName.Tokenize(":");
	if(inputVarNameArr->GetEntriesFast()!=GetNVars()){
		std::cout<<"ERROR: Mismatch no. of input variables."<<endl;
		return result;
	}
	for(Int_t i=0; i<GetNVars(); i++)
		if(GetVarName(i)!=inputVarNameArr->At(i)->GetName()){
			std::cout<<"ERROR: Variable name mismatch in the "<<i+1<<"-th input variable (should be "<<GetVarName(i) <<" instead of "<<inputVarNameArr->At(i)->GetName()<<endl;
			result=kFALSE;
	}
	return result;
}
//---------------------------------------------------------------------------
AliRDHFDecisionTree *AliRDHFBDT::AddDecisionTree(AliRDHFDecisionTree *tree, Double_t weight)
{
	if(!tree->VerifyNodes()){
		std::cout << "ERROR: Added tree has problem..." << endl;
		return 0x0;
	}
	else{
		fNTrees++;
		AliRDHFDecisionTree *newtree = new(fDecisionTrees[fDecisionTrees.GetEntriesFast()]) AliRDHFDecisionTree(*tree);
		fDTWeights.Set(GetNTrees());
		fDTWeights.AddAt(weight,GetNTrees()-1);
		return newtree;
	}
}
//---------------------------------------------------------------------------
AliRDHFDecisionTree *AliRDHFBDT::GetDecisionTree(Int_t i)
{
	if(i>=GetNTrees()||i<0){
		std::cout<<"ERROR: Input index is out of bound, we have " << GetNTrees() << " trees." <<endl;
		return 0x0;
	}
	else return (AliRDHFDecisionTree *)fDecisionTrees.UncheckedAt(i);
}
//---------------------------------------------------------------------------
Double_t AliRDHFBDT::GetBoostWeight(Int_t i)
{
	if(i>=GetNTrees()||i<0){
		std::cout<<"ERROR: Input index is out of bound, we have " << GetNTrees() << " trees." <<endl;
		return 0;
	}
	else return fDTWeights.At(i);
}
//---------------------------------------------------------------------------
void AliRDHFBDT::PrintAll()
{
	std::cout<<"		BDT:"<<endl;
	for(Int_t i=0;i<GetNTrees();i++){
		std::cout<<"	"<<i<<"-Tree with weight:  " << GetBoostWeight(i) <<endl;
		GetDecisionTree(i)->PrintAll();
	}
	return;
}
