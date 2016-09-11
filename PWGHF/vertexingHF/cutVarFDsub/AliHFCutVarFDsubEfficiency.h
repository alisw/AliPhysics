#ifndef ALIHFCUTVARFDSUBEFFICIENCY_H
#define ALIHFCUTVARFDSUBEFFICIENCY_H
/// \class AliHFCutVarFDsubEfficiency
/// \brief Eficiency determination for the cut variation feed down method analysis
///
///
///
///
/// \author Felix Reidt <felix.reidt@cern.ch>, CERN
/// \author Fabrizio Grosa <grosa@to.infn.it>, INFN Torino
/// \date Aug 17, 2015

#include "TObject.h"
#include "THnSparse.h"
#include "TF1.h"
#include "TH1F.h"

class TList;
class AliHFCutVarFDsubCutSet;

class AliHFCutVarFDsubEfficiency : public TObject {
protected:
  THnSparseF* fGenLevel;           //!>! MC generator level information
  THnSparseF* fAfterCuts;          //!>! MC information after cuts
  AliHFCutVarFDsubCutSet* fCutSet; //!>! Cut set for which the efficiency will be calculated
  TList* fAxes;                    //!>! List of axes present in the THnSparses
  UInt_t fTHnType;                 //
  Bool_t fPtWeight;                // Flag to activate the pT reweights
  TF1* fFuncWeights;               //!>! Pt reweights function
  Double_t fEfficiency;            /// Efficiency value  
  Double_t fEfficiencyError;       /// Efficiency error

  TH1F* ProjectMCSparse(THnSparseF* thns, UInt_t thnType, TF1* funcWeights, Bool_t reflectedaxes);
  
  AliHFCutVarFDsubEfficiency(const AliHFCutVarFDsubEfficiency& eff); /// Copy constructor
  AliHFCutVarFDsubEfficiency operator=(const AliHFCutVarFDsubEfficiency& eff); /// Assignment operator
  
public:
  AliHFCutVarFDsubEfficiency(); /// Default constructor
  AliHFCutVarFDsubEfficiency(THnSparseF* genLevel, THnSparseF* afterCuts, AliHFCutVarFDsubCutSet* cutSet, TList* axes, Bool_t ptWeight, TF1* funcWeights); /// Constructor

  Double_t GetEfficiency() { return fEfficiency; }
  Double_t GetEfficiencyError() { return fEfficiencyError; }

  /// \cond CLASSDEF
  ClassDef(AliHFCutVarFDsubEfficiency, 2);
  /// \endcond
};
#endif //ALIHFCUTVARFDSUBEFFICIENCY_H
