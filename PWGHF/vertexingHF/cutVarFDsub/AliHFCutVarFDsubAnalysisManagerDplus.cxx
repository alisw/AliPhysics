#include "AliHFCutVarFDsubAnalysisManagerDplus.h"

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
ClassImp(AliHFCutVarFDsubAnalysisManagerDplus);
/// \endcond

using std::cerr;
using std::cout;
using std::endl;
using std::vector;

AliHFCutVarFDsubAnalysisManagerDplus::AliHFCutVarFDsubAnalysisManagerDplus()
  : AliHFCutVarFDsubAnalysisManager()
  , fNevents(0.)
  , fPID(kTRUE)
  , fPIDAxis(3)
{
  /// Default constructor
}


//_________________________________________________________________________________________________
AliHFCutVarFDsubAnalysisManagerDplus::~AliHFCutVarFDsubAnalysisManagerDplus() {
  /// Destructor
}


//_________________________________________________________________________________________________
Int_t AliHFCutVarFDsubAnalysisManagerDplus::GetTHnSparses(const TString strFileMC,
                                                          const TString strFileData,
                                                          const TString strDirMC,
                                                          const TString strDirData,
                                                          const TString strListMC,
                                                          const TString strListData,
                                                          Bool_t MConly,
                                                          const TString sparseData,
                                                          const TString sparseMCpromptreco,
                                                          const TString sparseMCfeeddownreco,
                                                          const TString sparseMCpromptgen,
                                                          const TString sparseMCFDgen) {
  /// Retrieve the input data
  //------------------------------------------------------------------------------------------------
  //MC THnSparses
  TFile *infileMC = TFile::Open(strFileMC.Data(),"READ");//open MC file
  if (!infileMC) {cerr << "File " << strFileMC << " not found!" << endl; return 1;}
  TDirectoryFile *DplusDirMC = (TDirectoryFile*)infileMC->Get(strDirMC.Data());//open MC directory
  if (!DplusDirMC) {cerr << "Directory " << strDirMC << " not found!" << endl; return 2;}
  TList *ListMC = (TList*)DplusDirMC->Get(strListMC.Data());//get MC list
  if (!ListMC) {cerr << "List " << strListMC << " not found!" << endl; return 3;}

  //get THnSparses
  fMCafterCuts[kPrompt] = (THnSparseF*)ListMC->FindObject(sparseMCpromptreco);
  if (!fMCafterCuts[kPrompt]) { cerr << "THnSparseF " <<sparseMCpromptreco << " not found!" << endl; return 4; }
  fMCafterCuts[kFD] = (THnSparseF*)ListMC->FindObject(sparseMCpromptreco);
  if (!fMCafterCuts[kFD]) { cerr << "THnSparseF "<<sparseMCpromptreco<<" not found!" << endl; return 5; }
  fMCgenLevel[kPrompt] = (THnSparseF*)ListMC->FindObject(sparseMCpromptgen);
  if (!fMCgenLevel[kPrompt]) { cerr << "THnSparseF "<<sparseMCpromptgen<<" not found!" << endl; return 6; }
  fMCgenLevel[kFD] = (THnSparseF*)ListMC->FindObject(sparseMCFDgen);
  if (!fMCgenLevel[kFD]) { cerr << "THnSparseF " <<sparseMCFDgen<<" not found!" << endl; return 7; }
  //-------------------------------------------------------------------------------------------------
  //Real Data
  TString normobj=strListData;
  normobj.ReplaceAll("coutputDplus","coutputDplusNorm");

  if(MConly) {
    //get THnSparse from MC file
    fData = (THnSparseF*)ListMC->FindObject(sparseData);
    if (!fData) { cerr << "THnSparseF "<< sparseData <<" not found!" << endl; return 8; }

    //Dummy number of events in case of MC
    fNevents = 1;
  }
  else {
    TFile *infileData = TFile::Open(strFileData.Data(),"READ");//open Data file
    if (!infileData) {cerr << "File " << strFileData << " not found!" << endl; return 10;}
    TDirectoryFile *DplusDirData = (TDirectoryFile*)infileData->Get(strDirData.Data());//open Data directory
    if (!DplusDirData) {cerr << "Directory " << strDirData << " not found!" << endl; return 11;}
    TList *ListData = (TList*)DplusDirData->Get(strListData.Data());//get Data list
    if (!ListData) {cerr << "List " << strListData << " not found!" << endl; return 12;}

    //get THnSparse from Data file
    fData = (THnSparseF*)ListData->FindObject(sparseData);
    if (!fData) { cerr << "THnSparseF "<<sparseData <<" not found!" << endl; return 8; }

    //Number of events from normalization counter
    AliNormalizationCounter* nc=(AliNormalizationCounter*)DplusDirData->Get(normobj.Data());
    fNevents = nc->GetNEventsForNorm();
  }

  return 0;
}


//_________________________________________________________________________________________________
void AliHFCutVarFDsubAnalysisManagerDplus::GetCuts(Double_t*** cutslowset,
                                                   Double_t*** cutshighset,
                                                   Double_t** means,
                                                   Double_t** sigmas,
                                                   Int_t Rebin,
                                                   Int_t fsig,
                                                   Int_t fbkg,
                                                   Double_t min,
                                                   Double_t max,
                                                   Int_t nSets,
                                                   Int_t nPtBins,
                                                   Int_t nCutVariables) {
  /// Prepare the cuts
  if (fCuts) {
    delete fCuts;
    fCuts = 0x0;
  }
  fCuts= new TList();
  fCuts->SetOwner(0);

  TList *ptBin = 0x0;
  AliHFCutVarFDsubCutSet **cutset = new AliHFCutVarFDsubCutSet*[nSets];

  for(Int_t iPt=0; iPt<nPtBins; iPt++) {
    ptBin = new TList();
    ptBin->SetOwner(0);
    for(Int_t iSet=0; iSet<nSets; iSet++) {
      cutset[iSet] = new AliHFCutVarFDsubCutSet(1.72,2.05,means[iSet][iPt],sigmas[iSet][iPt],Rebin,fbkg,fsig,0);
      for(Int_t iCutVar=0; iCutVar<nCutVariables; iCutVar++) {
        cutset[iSet]->AddCut(new AliHFCutVarFDsubCut(iCutVar,cutslowset[iSet][iPt][iCutVar],cutshighset[iSet][iPt][iCutVar]));
      }
      if(fPID) {
        cutset[iSet]->AddCut(new AliHFCutVarFDsubCut(nCutVariables,1.51,2.49));
      }
      ptBin->Add((TObject*)cutset[iSet]);
    }
    fCuts->Add((TObject*)ptBin);
  }
}


//_________________________________________________________________________________________________
void AliHFCutVarFDsubAnalysisManagerDplus::GetAxes(UInt_t* dataAxesNo, UInt_t* MCGenAxesNo, UInt_t* MCCutAxesNo, TString* axesName, Int_t nAxes, Bool_t *isCutSymm=0x0) {
  if (fAxes) {
    delete fAxes;
    fAxes = 0x0;
  }
  fAxes = new TList();
  fAxes->SetOwner();

  // -1 means axis is not existent (or not used)
  //
  // On the generator level, one doesn't want to cut on the variation variables, but one would
  // like to select e.g. pt bins. In order to be to use the same cut for this, variables on
  // which one doesn't what to cut at generator level have the corresponding axes set to -1
  for(Int_t iAxis=0; iAxis<nAxes; iAxis++) {
    if(dataAxesNo[iAxis]==fPIDAxis)
      continue; /// PID axis set with SetPID();
    if(isCutSymm)
      fAxes->Add((TObject*)new AliHFCutVarFDsubAxis(dataAxesNo[iAxis],MCGenAxesNo[iAxis],MCCutAxesNo[iAxis],axesName[iAxis],isCutSymm[iAxis]));
    else
      fAxes->Add((TObject*)new AliHFCutVarFDsubAxis(dataAxesNo[iAxis],MCGenAxesNo[iAxis],MCCutAxesNo[iAxis],axesName[iAxis],kFALSE));
  }
  if(fPID) {
    fAxes->Add((TObject*)new AliHFCutVarFDsubAxis(fPIDAxis,(UInt_t)-1,fPIDAxis,"PID",kFALSE));
  }
}

//_________________________________________________________________________________________________
TH1F* AliHFCutVarFDsubAnalysisManagerDplus::CalculateCrossSection(TString AccFilePath,
                                                                  TString GenLimAccHistoName,
                                                                  TString GenAccHistoName,
                                                                  Int_t PromptOrFD,
                                                                  Double_t BR,
                                                                  Double_t sigma) {

  if(!fCorrYieldPrompt || !fCorrYieldFD) {
    cerr << "Minimisation needed before calculate cross section" << endl;
    return 0x0;
  }
  else if(AccFilePath=="") {
    cerr << "The Acceptance Factor needed to calculate the cross section" << endl;
    return 0x0;
  }
  else {
    TFile *AccFile = TFile::Open(AccFilePath.Data(),"READ");
    if (!AccFile) {
      cerr << "File " << AccFilePath << " not found!" << endl;
      return 0x0;
    }
    else {
      TH1F* hPtGenLimAcc = (TH1F*)AccFile->Get(GenLimAccHistoName.Data());
      TH1F* hPtGenAcc = (TH1F*)AccFile->Get(GenAccHistoName.Data());
      hPtGenAcc->SetDirectory(0);
      hPtGenLimAcc->SetDirectory(0);
      AccFile->Close();

      TAxis* PtAxis = (TAxis*)fCorrYieldPrompt->GetXaxis();
      TArrayD* PtBinsArray = (TArrayD*)PtAxis->GetXbins();
      Double_t* PtLims = (Double_t*)PtBinsArray->GetArray();
      Int_t nPtBins = PtBinsArray->GetSize()-1;

      TH1F* hPtGenLimAccReb = (TH1F*)hPtGenLimAcc->Rebin(nPtBins, "hPtGenLimAccReb", PtLims);
      TH1F* hPtGenAccReb = (TH1F*)hPtGenAcc->Rebin(nPtBins,"hPtGenAccReb", PtLims);
      hPtGenLimAccReb->Sumw2();
      hPtGenAccReb->Sumw2();

      TH1F* hPtAcc = new TH1F("hPtAcc","",nPtBins,PtLims);
      hPtAcc->Divide(hPtGenAccReb,hPtGenLimAccReb,1.,1.,"B");

      TString crosssectype="Prompt";
      if(PromptOrFD==kFD) crosssectype="FD";
      TH1F* hCrossSec = new TH1F(Form("hCrossSec%s",crosssectype.Data()),"",nPtBins,PtLims);

      for(Int_t iBin=0; iBin<nPtBins; iBin++) {
        Double_t CrossSec;
        Double_t CrossSecError;

        if(PromptOrFD==kPrompt) {
          //cross section in mubarn
          CrossSec = 1./hPtAcc->GetBinContent(iBin+1)*fCorrYieldPrompt->GetBinContent(iBin+1)*sigma*1000000./(2*(PtLims[iBin+1]-PtLims[iBin])*fNevents*BR);
          CrossSecError = TMath::Sqrt((fCorrYieldPrompt->GetBinError(iBin+1)/fCorrYieldPrompt->GetBinContent(iBin+1))*(fCorrYieldPrompt->GetBinError(iBin+1)/fCorrYieldPrompt->GetBinContent(iBin+1))+(hPtAcc->GetBinError(iBin+1)/hPtAcc->GetBinContent(iBin+1))*(hPtAcc->GetBinError(iBin+1)/hPtAcc->GetBinContent(iBin+1)))*CrossSec;
        }
        else {
          CrossSec = 1./hPtAcc->GetBinContent(iBin+1)*fCorrYieldFD->GetBinContent(iBin+1)*sigma*1000000./(2*(PtLims[iBin+1]-PtLims[iBin])*fNevents*BR);
          CrossSecError = TMath::Sqrt((fCorrYieldFD->GetBinError(iBin+1)/fCorrYieldFD->GetBinContent(iBin+1))*(fCorrYieldFD->GetBinError(iBin+1)/fCorrYieldFD->GetBinContent(iBin+1))+(hPtAcc->GetBinError(iBin+1)/hPtAcc->GetBinContent(iBin+1))*(hPtAcc->GetBinError(iBin+1)/hPtAcc->GetBinContent(iBin+1)))*CrossSec;
        }
        hCrossSec->SetBinContent(iBin+1,CrossSec);
        hCrossSec->SetBinError(iBin+1,CrossSecError);
      }

      hCrossSec->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
      hCrossSec->GetYaxis()->SetTitle("#frac{d#sigma}{d#it{p}_{T}} (#mub c/GeV)");
      hCrossSec->GetYaxis()->SetTitleOffset(1.1);

      delete hPtAcc;
      hPtAcc = 0x0;

      return hCrossSec;
    }
  }
}


//_________________________________________________________________________________________________
AliHFCutVarFDsubAnalysisManagerDplus::AliHFCutVarFDsubAnalysisManagerDplus(const AliHFCutVarFDsubAnalysisManagerDplus& am)
  : AliHFCutVarFDsubAnalysisManager(am)
  , fNevents(am.fNevents)
  , fPID(am.fPID)
  , fPIDAxis(am.fPIDAxis)
{
  /// Copy constructor

  // TODO: do proper deep copies
  Printf("Do not use this copy constructor!");
}


AliHFCutVarFDsubAnalysisManagerDplus& AliHFCutVarFDsubAnalysisManagerDplus::operator=(const AliHFCutVarFDsubAnalysisManagerDplus& am)
{
  /// Assignment operator

  // TODO: do proper deep copies
  Printf("Do not use this assignment operator!");

  if (this != &am) {
    AliHFCutVarFDsubAnalysisManager::operator=(am);
    fNevents = am.fNevents;
    fPID = am.fPID;
    fPIDAxis = am.fPIDAxis;
  }
  return *this;
}
