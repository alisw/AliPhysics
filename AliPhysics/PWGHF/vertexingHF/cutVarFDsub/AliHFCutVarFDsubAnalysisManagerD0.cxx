#include "AliHFCutVarFDsubAnalysisManagerD0.h"

#include <iostream>
#include <vector>

#include "TH1F.h"
#include "TList.h"
#include "TFile.h"

#include "AliNormalizationCounter.h"

#include "AliHFCutVarFDsubAxis.h"
#include "AliHFCutVarFDsubCut.h"
#include "AliHFCutVarFDsubCutSet.h"


/// \cond CLASSIMP
ClassImp(AliHFCutVarFDsubAnalysisManagerD0);
/// \endcond

using std::cerr;
using std::endl;
using std::vector;

AliHFCutVarFDsubAnalysisManagerD0::AliHFCutVarFDsubAnalysisManagerD0()
  : AliHFCutVarFDsubAnalysisManager()
  , fNevents(0.)
{
  /// Default constructor
}


//_________________________________________________________________________________________________
AliHFCutVarFDsubAnalysisManagerD0::~AliHFCutVarFDsubAnalysisManagerD0() {
  /// Destructor
}


//_________________________________________________________________________________________________
Int_t AliHFCutVarFDsubAnalysisManagerD0::GetTHnSparses(const TString strFileData,
                                                       const TString strFilePrompt,
                                                       const TString strContNamePrompt,
                                                       const TString strFileFD,
                                                       const TString strContNameFD,
                                                       Bool_t MConly) {
  /// Retrieve the input data
  //------------------------------------------------------------------------------------------------
  // Real data
  TFile* fDat=0x0;
  if (!MConly) fDat = TFile::Open(strFileData.Data()  ,"READ");
  else         fDat = TFile::Open(strFilePrompt.Data(),"READ");
  if (!fDat) { cerr << "File " << strFileData << " not found!" << endl; return 1; }
  TDirectory* ddirData=(TDirectory*)fDat->Get("PWG3_D2H_D0InvMass");
  if (!ddirData) { cerr << "Directory PWG3_D2H_D0InvMass not found!" << endl; return 2; }
  fData=(THnSparseF*)ddirData->Get("fCutsDataFD");
  if (!fData) { cerr << "THnSparseF fCutsDataFD not found!" << endl; return 3; }
  AliNormalizationCounter* countev = 0x0;
  countev=(AliNormalizationCounter*)ddirData->Get("normalizationCounter0100");
  if (!countev) countev=(AliNormalizationCounter*)ddirData->Get("normalizationCounter010");

  if (!countev) { cerr << "NormalizationCounter* not found!" << endl; return 4; }
  fNevents=countev->GetNEventsForNorm();

  //------------------------------------------------------------------------------------------------
  // MC prompt D
  TFile* fMC = TFile::Open(strFilePrompt.Data(),"READ"); //to open the the root file
  if (!fMC) { cerr << "File " << strFilePrompt << " not found!" << endl; return 5; }
  TString strDir="PWG3_D2H_CFtaskD0toKpi";
  strDir+=strContNamePrompt;
  TDirectory* dirInput=(TDirectory*)fMC->Get(strDir.Data()); // inside the directory
  if (!dirInput) { cerr << "Directory " << strDir << " not found!" << endl; return 6; }
  fMCafterCuts[kPrompt] = (THnSparseF*)dirInput->Get("fCutsMCFD");
  if (!fMCafterCuts[kPrompt]) { cerr << "THnSparseF fCutsMCFD not found!" << endl; return 7; }
  fMCafterCuts[kPrompt]->SetName(Form("%sPrompt",fMCafterCuts[kPrompt]->GetName()));
  fMCgenLevel[kPrompt]=(THnSparseF*)dirInput->Get("fPtMCGenStep");
  if (!fMCgenLevel[kPrompt]) { cerr << "THnSparseF fPtMCGenStep not found!" << endl; return 8; }
  fMCgenLevel[kPrompt]->SetName("fPtMCGenStepPrompt");

  //------------------------------------------------------------------------------------------------
  // MC D from feed-down
  fMC = TFile::Open(strFileFD.Data(),"READ"); //to open the the root file
  if (!fMC) { cerr << "File " << strFileFD << " not found!" << endl; return 9; }
  strDir="PWG3_D2H_CFtaskD0toKpiKeepDfromBOnly";
  strDir+=strContNameFD;
  dirInput=(TDirectory*)fMC->Get(strDir.Data()); // inside the directory
  if (!dirInput) { cerr << "Directory " << strDir << " not found!" << endl; return 10; }
  fMCafterCuts[kFD] = (THnSparseF*)dirInput->Get("fCutsMCFD");
  if (!fMCafterCuts[kFD]) { cerr << "THnSparseF fCutsMCFD not found!" << endl; return 11; }
  fMCafterCuts[kFD]->SetName(Form("%sFD",fMCafterCuts[kFD]->GetName()));
  fMCgenLevel[kFD] = (THnSparseF*)dirInput->Get("fPtMCGenStep");
  if (!fMCgenLevel[kFD]) { cerr << "THnSparseF fPtMCGenStep not found!" << endl; return 12; }
  fMCgenLevel[kFD]->SetName("fPtMCGenStepFD");

  return 0;
}


//_________________________________________________________________________________________________
void AliHFCutVarFDsubAnalysisManagerD0::GetCuts() {
  /// Preparte the cuts

  if (fCuts) {
    delete fCuts;
    fCuts = 0x0;
  }
  fCuts = new TList();
  fCuts->SetOwner();
  TList* ptBin = new TList();
  ptBin->SetOwner();
  // - 1-2 GeV/c in pt
  // Prompt enhanced
  AliHFCutVarFDsubCutSet* cutSet = new AliHFCutVarFDsubCutSet(1.72, 2.05, 1.86484, 0.008);
  cutSet->AddCut(new AliHFCutVarFDsubCut(0, 1    , 2  )); // pt-bin
  cutSet->AddCut(new AliHFCutVarFDsubCut(1, 0    , 7  )); // NormDecLengthXY
  cutSet->AddCut(new AliHFCutVarFDsubCut(2, 0.998, 1.1)); // CosPointXY
  ptBin->Add((TObject*)cutSet);
  // Mixed
  cutSet = new AliHFCutVarFDsubCutSet(1.72, 2.05, 1.86484, 0.008);
  cutSet->AddCut(new AliHFCutVarFDsubCut(0, 1   ,  2  )); // pt-bin
  cutSet->AddCut(new AliHFCutVarFDsubCut(1, 7   , 12  )); // NormDecLengthXY
  cutSet->AddCut(new AliHFCutVarFDsubCut(2, 0.99,  1.1)); // CosPointXY
  ptBin->Add((TObject*)cutSet);
  // FD enhanced
  cutSet = new AliHFCutVarFDsubCutSet(1.72, 2.05, 1.86484, 0.008);
  cutSet->AddCut(new AliHFCutVarFDsubCut(0,  1   ,  2    )); // pt-bin
  cutSet->AddCut(new AliHFCutVarFDsubCut(1, 12   , 99    )); // NormDecLengthXY
  cutSet->AddCut(new AliHFCutVarFDsubCut(2,  0.99,  0.996)); // CosPointXY
  ptBin->Add((TObject*)cutSet);
  fCuts->Add((TObject*)ptBin);
  // - 2-3 GeV/c in pt
  // Prompt enhanced
  ptBin = new TList();
  ptBin->SetOwner();
  cutSet = new AliHFCutVarFDsubCutSet(1.72, 2.05, 1.86484, 0.008);
  cutSet->AddCut(new AliHFCutVarFDsubCut(0, 2    , 5  )); // pt-bin
  cutSet->AddCut(new AliHFCutVarFDsubCut(1, 0    , 7  )); // NormDecLengthXY
  cutSet->AddCut(new AliHFCutVarFDsubCut(2, 0.998, 1.1)); // CosPointXY
  ptBin->Add((TObject*)cutSet);
  // Mixed
  cutSet = new AliHFCutVarFDsubCutSet(1.72, 2.05, 1.86484, 0.008);
  cutSet->AddCut(new AliHFCutVarFDsubCut(0, 2   ,  5  )); // pt-bin
  cutSet->AddCut(new AliHFCutVarFDsubCut(1, 7   , 12  )); // NormDecLengthXY
  cutSet->AddCut(new AliHFCutVarFDsubCut(2, 0.99,  1.1)); // CosPointXY
  ptBin->Add((TObject*)cutSet);
  // FD enhanced
  cutSet = new AliHFCutVarFDsubCutSet(1.72, 2.05, 1.86484, 0.008);
  cutSet->AddCut(new AliHFCutVarFDsubCut(0,  2   ,  5    )); // pt-bin
  cutSet->AddCut(new AliHFCutVarFDsubCut(1, 12   , 99    )); // NormDecLengthXY
  cutSet->AddCut(new AliHFCutVarFDsubCut(2,  0.99,  0.996)); // CosPointXY
  ptBin->Add((TObject*)cutSet);
  fCuts->Add((TObject*)ptBin);
}


//_________________________________________________________________________________________________
void AliHFCutVarFDsubAnalysisManagerD0::GetAxes(Int_t version/*=0*/) {
  if (fAxes) {
    delete fAxes;
    fAxes = 0x0;
  }
  fAxes = new TList();
  fAxes->SetOwner();
  switch (version) {
  case 0:
  default:
    // standard version
    // -1 means axis is not existent (or not used)
    //
    // On the generator level, one doesn't want to cut on the variation variables, but one would
    // like to select e.g. pt bins. In order to be to use the same cut for this, variables on
    // which one doesn't what to cut at generator level have the corresponding axes set to -1
    fAxes->Add((TObject*)new AliHFCutVarFDsubAxis(         1,          0, 0, "pt"));
    fAxes->Add((TObject*)new AliHFCutVarFDsubAxis(         2, (UInt_t)-1, 1, "NormDecLengthXY"));
    fAxes->Add((TObject*)new AliHFCutVarFDsubAxis(         3, (UInt_t)-1, 2, "CosPointXY"));
    fAxes->Add((TObject*)new AliHFCutVarFDsubAxis((UInt_t)-1,          1, 3, "nProngs"));
    fAxes->Add((TObject*)new AliHFCutVarFDsubAxis((UInt_t)-1,          2, 4, "Mother_pt"));
    break;
  }
}


//_________________________________________________________________________________________________
AliHFCutVarFDsubAnalysisManagerD0::AliHFCutVarFDsubAnalysisManagerD0(const AliHFCutVarFDsubAnalysisManagerD0& am)
  : AliHFCutVarFDsubAnalysisManager(am)
  , fNevents(am.fNevents)
{
  /// Copy constructor

  // TODO: do proper deep copies
  Printf("Do not use this copy constructor!");
}
