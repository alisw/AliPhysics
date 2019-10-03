#include "AliHFCutVarFDsubEfficiency.h"

#include <iostream>

#include "TList.h"
#include "THnSparse.h"
#include "TH1F.h"

#include "AliHFCutVarFDsubAxis.h"
#include "AliHFCutVarFDsubCut.h"
#include "AliHFCutVarFDsubCutSet.h"

using std::cerr;
using std::endl;

/// \cond CLASSIMP
ClassImp(AliHFCutVarFDsubEfficiency);
/// \endcond


AliHFCutVarFDsubEfficiency::AliHFCutVarFDsubEfficiency()
  : TObject()
  , fGenLevel(0x0)
  , fAfterCuts(0x0)
  , fCutSet(0x0)
  , fAxes(0x0)
  , fPtWeight(kFALSE)
  , fFuncWeights(0x0)
  , fEfficiency(-1.)
  , fEfficiencyError(-1.)
{
  /// Default constructor
}


AliHFCutVarFDsubEfficiency::AliHFCutVarFDsubEfficiency(THnSparseF* genLevel, THnSparseF* afterCuts, AliHFCutVarFDsubCutSet* cutSet, TList* axes, Bool_t ptWeight, TF1* funcWeights)
  : TObject()
  , fGenLevel(genLevel)
  , fAfterCuts(afterCuts)
  , fCutSet(cutSet)
  , fAxes(axes)
  , fPtWeight(ptWeight)
  , fFuncWeights(funcWeights)
  , fEfficiency(-1.)
  , fEfficiencyError(-1.)
{
  /// Constructor

  // combine the two THnSparses to an array, as every below will be applied to both the same way
  THnSparseF* thns[] = { fGenLevel, fAfterCuts };
  UInt_t thnTypes[]  = { AliHFCutVarFDsubAxis::kMCgenLevel, AliHFCutVarFDsubAxis::kMCafterCuts };
  UInt_t nTHns = sizeof(thnTypes)/sizeof(UInt_t);

  Double_t counts[] = { 0., 0. }; // denominator, nominator of the efficiency equation

  // release all axes
  for (UInt_t iTHn=0; iTHn<nTHns; ++iTHn) {
    for (Int_t iAxis=0; iAxis<thns[iTHn]->GetNdimensions(); ++iAxis) {
      thns[iTHn]->GetAxis(iAxis)->SetRange(-1, -1);
    }
  }

  // apply the cuts
  Int_t SymmCutAxisCounter=0;
  for (Int_t iCut=0; iCut<fCutSet->GetEntries(); ++iCut) {
    AliHFCutVarFDsubCut* cut = fCutSet->GetCut(iCut);
    AliHFCutVarFDsubAxis* axis = (AliHFCutVarFDsubAxis*)fAxes->At(cut->fAxisId);
    if(axis->IsCutSymmetric() && cut->fHigh != -cut->fLow)
      SymmCutAxisCounter++;
  }
  for (UInt_t iTHn=0; iTHn<nTHns; ++iTHn) {
    TH1F* t = 0x0;
    if(iTHn==0 || SymmCutAxisCounter==0) {
      // if none of the axis has a two-region cut symmetric with respect to zero just project
      // or if it is the MC gen level THnSparse (no cuts applied in this case)
      t = (TH1F*)ProjectMCSparse(thns[iTHn],thnTypes[iTHn],funcWeights,kFALSE);
    }
    else { // else get the histograms corresponding to the two symmetric regions and sum them
      TH1F* t1 = (TH1F*)ProjectMCSparse(thns[iTHn],thnTypes[iTHn],funcWeights,kFALSE);
      TH1F* t2 = (TH1F*)ProjectMCSparse(thns[iTHn],thnTypes[iTHn],funcWeights,kTRUE);
      t = (TH1F*)t1->Clone();
      t->Add(t1,t2,1.,1.);
      delete t1;
      t1=0x0;
      delete t2;
      t2=0x0;
    }
    AliHFCutVarFDsubCut* cut = fCutSet->GetCut(0);
    Int_t binMin = t->FindBin(cut->fLow*1.0001);
    Int_t binMax = t->FindBin(cut->fHigh*0.9999);
    counts[iTHn] = t->Integral(binMin, binMax);
    delete t;
    t = 0x0;
  }
  
  TH1F *hTempNum=new TH1F("hTempNum","hTempNum",1,0,1);
  hTempNum->SetBinContent(1,counts[1]);
  hTempNum->SetBinError(1,TMath::Sqrt(counts[1]));

  TH1F *hTempDeNum=new TH1F("hTempDeNum","hTempDeNum",1,0,1);
  hTempDeNum->SetBinContent(1,counts[0]);
  hTempDeNum->SetBinError(1,TMath::Sqrt(counts[0]));
  
  TH1F *hTempEff=(TH1F*)hTempNum->Clone("hTempEff");
  hTempEff->Divide(hTempNum,hTempDeNum,1.,1.,"B");

  fEfficiency      = hTempEff->GetBinContent(1);
  fEfficiencyError = hTempEff->GetBinError(1);

  delete hTempNum;
  hTempNum = 0x0;
  delete hTempDeNum;
  hTempDeNum = 0x0;
  delete hTempEff;
  hTempEff = 0x0;
}

TH1F* AliHFCutVarFDsubEfficiency::ProjectMCSparse(THnSparseF* thns, UInt_t thnType, TF1* funcWeights, Bool_t reflectedaxes)
{
   // apply the cuts
  for (Int_t iCut=1; iCut<fCutSet->GetEntries(); ++iCut) {
    AliHFCutVarFDsubCut* cut = fCutSet->GetCut(iCut);
    AliHFCutVarFDsubAxis* axis = (AliHFCutVarFDsubAxis*)fAxes->At(cut->fAxisId);
    UInt_t axisNo = axis->GetAxisNo(thnType);
    Bool_t isCutSymm = axis->IsCutSymmetric();
    if (axisNo<(UInt_t)-1) {
      TAxis* ax = thns->GetAxis(axisNo);
      Int_t binMin;
      Int_t binMax;
      if(reflectedaxes && isCutSymm) {
        binMin = ax->FindBin(-cut->fHigh*1.0001);
        binMax = ax->FindBin(-cut->fLow*0.9999);
      }
      else {
        binMin = ax->FindBin(cut->fLow*1.0001);
        binMax = ax->FindBin(cut->fHigh*0.9999);
      }
      ax->SetRange(binMin, binMax);
    }
  }

  AliHFCutVarFDsubCut* cut = fCutSet->GetCut(0);
  AliHFCutVarFDsubAxis* axis = (AliHFCutVarFDsubAxis*)fAxes->At(cut->fAxisId);
  UInt_t axisNo = axis->GetAxisNo(thnType);
  TH1F* t=0x0;
  if (axisNo<(UInt_t)-1) {
    t = (TH1F*)thns->Projection(axisNo);
    if(fPtWeight) {
      for(Int_t iBin=0; iBin<t->GetNbinsX(); iBin++) {
        t->SetBinContent(iBin+1,t->GetBinContent(iBin+1)*fFuncWeights->Eval(t->GetBinCenter(iBin+1)));
      }
    }
  }
  else {
    cerr << "Couldn't apply first cut to the THn!" << endl;
    return 0x0;
  }
  
  return t;
}

AliHFCutVarFDsubEfficiency::AliHFCutVarFDsubEfficiency(const AliHFCutVarFDsubEfficiency& eff)
  : TObject()
  , fGenLevel(eff.fGenLevel)
  , fAfterCuts(eff.fAfterCuts)
  , fCutSet(eff.fCutSet)
  , fAxes(eff.fAxes)
  , fEfficiency(eff.fEfficiency)
  , fEfficiencyError(eff.fEfficiencyError)
{
  /// Copy constructor
}


AliHFCutVarFDsubEfficiency AliHFCutVarFDsubEfficiency::operator=(const AliHFCutVarFDsubEfficiency& eff)
{
  /// Assignment operator
  if (this != &eff) {
    fGenLevel = eff.fGenLevel;
    fAfterCuts = eff.fAfterCuts;
    fCutSet = eff.fCutSet;
    fAxes = eff.fAxes;
    fEfficiency = eff.fEfficiency;
    fEfficiencyError = eff.fEfficiencyError;
  }
  return *this;
}
