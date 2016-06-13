/*
 ***********************************************************
 Manager for event plane corrections framework
Contact: Jaap Onderwaater, j.onderwaater@gsi.de, jacobus.onderwaater@cern.ch
Instructions in AddTask_EPcorrectionsExample.C
2014/12/10
 *********************************************************
 */

//#include "AliSysInfo.h"
#include <iostream>

#include <TROOT.h>
#include <TTimeStamp.h>
#include <TStopwatch.h>
#include <TChain.h>
#include <AliInputEventHandler.h>
#include <AliESDInputHandler.h>
#include <AliAODInputHandler.h>
#include <AliAnalysisManager.h>
#include <AliCentrality.h>
#include <AliESDEvent.h>
#include <TList.h>
#include <TGraphErrors.h>
#include <AliLog.h>
#include "AliQnCorrectionsCutsSet.h"
#include "AliQnCorrectionsManager.h"
#include "AliQnCorrectionsHistos.h"
#include "AliAnalysisTaskQnVectorAnalysis.h"
#include "AliQnCorrectionsQnVector.h"

// make a change

ClassImp(AliAnalysisTaskQnVectorAnalysis)

/* names for the different TProfile */
TString namesQnTrackDetectors[AliAnalysisTaskQnVectorAnalysis::nTrackDetectors] = {"TPC","SPD"};
TString namesQnEPDetectors[AliAnalysisTaskQnVectorAnalysis::nEPDetectors] = {"V0A","V0C","T0A","T0C","FMDA","FMDC"/*,"rawFMDA","rawFMDC"*/};
TString namesQnComponents[AliAnalysisTaskQnVectorAnalysis::kNcorrelationComponents] = {"XX","XY","YX","YY"};

/* centrality binning */
const Int_t nQnCentBins = 11;
Double_t centQnBinning[nQnCentBins+1] = {0.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0};
Double_t centQnBinningmid[nQnCentBins] = {2.5, 7.5, 15.0, 25.0, 35.0, 45.0, 55.0, 65.0, 75.0, 85.0, 95.0};

//_________________________________________________________________________________
AliAnalysisTaskQnVectorAnalysis::AliAnalysisTaskQnVectorAnalysis() :
        AliQnCorrectionsFillEventTask(),
  fEventQAList(0x0),
  fEventCuts(NULL),
  fEventPlaneHistos(0x0),
  fVn(),
  fNDetectorResolutions(0),
  fDetectorResolution(),
  fDetectorResolutionContributors(),
  fDetectorResolutionCorrelations(),
  fTrackDetectorNameInFile(),
  fEPDetectorNameInFile(),
  fCentralityVariable(-1),
  fExpectedCorrectionPass("rec"),
  fAlternativeCorrectionPass("rec")
{
  //
  // Default constructor
  //
}

//_________________________________________________________________________________
AliAnalysisTaskQnVectorAnalysis::AliAnalysisTaskQnVectorAnalysis(const char* name) :
        AliQnCorrectionsFillEventTask(name),
  fEventQAList(0x0),
  fEventCuts(0x0),
  fEventPlaneHistos(0x0),
  fVn(),
  fNDetectorResolutions(0),
  fDetectorResolution(),
  fDetectorResolutionContributors(),
  fDetectorResolutionCorrelations(),
  fTrackDetectorNameInFile(),
  fEPDetectorNameInFile(),
  fCentralityVariable(-1),
  fExpectedCorrectionPass("rec"),
  fAlternativeCorrectionPass("rec")
{
  //
  // Constructor
  //

  fEventQAList = new TList();
  fEventQAList->SetName("EventQA");
  fEventQAList->SetOwner(kTRUE);

  fEventPlaneHistos = new AliQnCorrectionsHistos();

  /* initialize TProfile storage */
  for(Int_t i=0; i < nTrackDetectors*nEPDetectors; i++) {
    for(Int_t j=0; j < kNharmonics; j++) {
      for(Int_t k=0; k < kNcorrelationComponents; k++) {
        fVn[i][j][k]=NULL;
      }
    }
  }

  for (Int_t i=0; i < nTrackDetectors*nEPDetectors; i++) {
    for (Int_t j = 0; j < kNresolutionComponents; j++) {
      fDetectorResolutionContributors[i][j] = -1;
    }
    for (Int_t iCorr = 0; iCorr < nCorrelationPerDetector; iCorr++) {
      for(Int_t h=0; h < kNharmonics; h++) {
        fDetectorResolutionCorrelations[i][iCorr][h] = NULL;
      }
    }
    for(Int_t h=0; h < kNharmonics; h++) {
      fDetectorResolution[i][h] = NULL;
    }
  }

  /* the index in the input TList structure for the data of the different detectors */
  fTrackDetectorNameInFile[kTPC] = "TPC";
  fTrackDetectorNameInFile[kSPD] = "SPD";
  fEPDetectorNameInFile[kVZEROA] = "VZEROA";
  fEPDetectorNameInFile[kVZEROC] = "VZEROC";
  fEPDetectorNameInFile[kTZEROA] = "TZEROA";
  fEPDetectorNameInFile[kTZEROC] = "TZEROC";
  fEPDetectorNameInFile[kFMDA] = "FMDA";
  fEPDetectorNameInFile[kFMDC] = "FMDC";
//  fEPDetectorNameInFile[kRawFMDA] = "FMDAraw";
//  fEPDetectorNameInFile[kRawFMDC] = "FMDCraw";

  /* the detector resolution configurations */
  /* this is only valid in C++11
  fDetectorResolutionContributors = {
      {kTPC,kVZEROC,kVZEROA},
      {kTPC,kFMDA,kFMDC},
      {kSPD,kVZEROC,kVZEROA},
      {kSPD,kFMDA,kFMDC}
  }; */
  /* for the time being we do it in this awful way */
  const Int_t nDetectorResolutions = 10;
  fNDetectorResolutions = nDetectorResolutions;
  Int_t config[nDetectorResolutions][kNresolutionComponents] = {
      {kTPC,kVZEROC,kVZEROA},
      {kTPC,kTZEROC,kTZEROA},
      {kTPC,kFMDA,kFMDC},
//      {kTPC,kRawFMDA,kRawFMDC},
      {kTPC,kVZEROC,kTZEROA},
      {kTPC,kTZEROC,kVZEROA},
      {kSPD,kVZEROC,kVZEROA},
      {kSPD,kTZEROC,kTZEROA},
      {kSPD,kFMDA,kFMDC},
//      {kSPD,kRawFMDA,kRawFMDC},
      {kSPD,kVZEROC,kTZEROA},
      {kSPD,kTZEROC,kVZEROA},
  };
  for (Int_t ixConfig = 0; ixConfig < fNDetectorResolutions; ixConfig++) {
    for (Int_t i = 0; i < kNresolutionComponents; i++) {
      fDetectorResolutionContributors[ixConfig][i] = config[ixConfig][i];
    }
  }

  /* create the needed TProfile for each of the track-EP detector combination
   * and for each of the v_n and for the different correlation namesQnComponents */
  for(Int_t iTrkDetector=0; iTrkDetector < nTrackDetectors; iTrkDetector++) {
    for (Int_t iEPDetector=0; iEPDetector < nEPDetectors; iEPDetector++) {
      for(Int_t h=0; h<kNharmonics; h++) {
        for(Int_t corrComp =0; corrComp< kNcorrelationComponents; corrComp++){
          fVn[iTrkDetector*nEPDetectors+iEPDetector][h][corrComp] =
              new TProfile(Form("vn_%sx%s_%s_h%d", namesQnTrackDetectors[iTrkDetector].Data(), namesQnEPDetectors[iEPDetector].Data(), namesQnComponents[corrComp].Data(),h+1),
                  Form("vn_%sx%s_%s_h%d", namesQnTrackDetectors[iTrkDetector].Data(), namesQnEPDetectors[iEPDetector].Data(), namesQnComponents[corrComp].Data(),h+1),
                  nQnCentBins,
                  centQnBinning);
        }
      }
    }
  }

  /* create the needed correlation TProfile for each the detector resolution and desired additional detector configuration */
  for(Int_t ixDetectorConfig = 0; ixDetectorConfig < nTrackDetectors*nEPDetectors; ixDetectorConfig++) {
    /* filter out not initialized detector configurations */
    if(fDetectorResolutionContributors[ixDetectorConfig][0]==-1) continue;
    for (Int_t iCorr = 0; iCorr < nCorrelationPerDetector; iCorr++) {
      TString detectorOneName;
      TString detectorTwoName;
      TString correlationComponent;
      switch (iCorr) {
      case kABXX: {
        detectorOneName = namesQnTrackDetectors[fDetectorResolutionContributors[ixDetectorConfig][0]];
        detectorTwoName = namesQnEPDetectors[fDetectorResolutionContributors[ixDetectorConfig][1]];
        correlationComponent = namesQnComponents[kXX];
        break;
      }
      case kABYY: {
        detectorOneName = namesQnTrackDetectors[fDetectorResolutionContributors[ixDetectorConfig][0]];
        detectorTwoName = namesQnEPDetectors[fDetectorResolutionContributors[ixDetectorConfig][1]];
        correlationComponent = namesQnComponents[kYY];
        break;
      }
      case kACXX: {
        detectorOneName = namesQnTrackDetectors[fDetectorResolutionContributors[ixDetectorConfig][0]];
        detectorTwoName = namesQnEPDetectors[fDetectorResolutionContributors[ixDetectorConfig][2]];
        correlationComponent = namesQnComponents[kXX];
        break;
      }
      case kACYY: {
        detectorOneName = namesQnTrackDetectors[fDetectorResolutionContributors[ixDetectorConfig][0]];
        detectorTwoName = namesQnEPDetectors[fDetectorResolutionContributors[ixDetectorConfig][2]];
        correlationComponent = namesQnComponents[kYY];
        break;
      }
      case kBCXX: {
        detectorOneName = namesQnEPDetectors[fDetectorResolutionContributors[ixDetectorConfig][1]];
        detectorTwoName = namesQnEPDetectors[fDetectorResolutionContributors[ixDetectorConfig][2]];
        correlationComponent = namesQnComponents[kXX];
        break;
      }
      case kBCYY: {
        detectorOneName = namesQnEPDetectors[fDetectorResolutionContributors[ixDetectorConfig][1]];
        detectorTwoName = namesQnEPDetectors[fDetectorResolutionContributors[ixDetectorConfig][2]];
        correlationComponent = namesQnComponents[kYY];
        break;
      }
      }
      for(Int_t h=0; h < kNharmonics; h++) {
        TString profileName = Form("corr_%s%d_%sx%s_%s_h%d",
            namesQnTrackDetectors[fDetectorResolutionContributors[ixDetectorConfig][0]].Data(),
            ixDetectorConfig,
            detectorOneName.Data(),
            detectorTwoName.Data(),
            correlationComponent.Data(),
            h+1);
        fDetectorResolutionCorrelations[ixDetectorConfig][iCorr][h] = new TProfile(profileName.Data(),profileName.Data(), nQnCentBins, centQnBinning);
      }
    }
  }

  /* create the sub-event structures to store the cross correlations */

  DefineInput(0,TChain::Class());
  DefineInput(1,TList::Class());
  DefineOutput(1, TList::Class());// Event QA histograms
}

AliAnalysisTaskQnVectorAnalysis::~AliAnalysisTaskQnVectorAnalysis() {
  /* clean up everything before leaving */
  for(Int_t iTrkDetector=0; iTrkDetector < nTrackDetectors; iTrkDetector++) {
    for (Int_t iEPDetector=0; iEPDetector < nEPDetectors; iEPDetector++) {
      for(Int_t h=0; h<kNharmonics; h++) {
        for(Int_t corrComp =0; corrComp< kNcorrelationComponents; corrComp++){
          delete fVn[iTrkDetector*nEPDetectors+iEPDetector][h][corrComp];
        }
      }
    }
  }
  for(Int_t ixDetectorConfig = 0; ixDetectorConfig < nTrackDetectors*nEPDetectors; ixDetectorConfig++) {
    /* filter out not initialized detector configurations */
    if(fDetectorResolutionContributors[ixDetectorConfig][0]==-1) continue;
    for (Int_t iCorr = 0; iCorr < nCorrelationPerDetector; iCorr++) {
      for(Int_t h=0; h < kNharmonics; h++) {
        delete fDetectorResolutionCorrelations[ixDetectorConfig][iCorr][h];
      }
    }
    for(Int_t h=0; h < kNharmonics; h++) {
      // delete fDetectorResolution[ixDetectorConfig][h];
      // delete fNewDetectorResolution[ixDetectorConfig][h];
    }
  }
}
//_________________________________________________________________________________
void AliAnalysisTaskQnVectorAnalysis::UserCreateOutputObjects()
{
  //
  // Add all histogram manager histogram lists to the output TList
  //

  PostData(1, fEventQAList);

}



//________________________________________________________________________________________________________
void AliAnalysisTaskQnVectorAnalysis::UserExec(Option_t *){
  //
  // Main loop. Called for every event
  //

#define MAXNOOFDATAVARIABLES 2048

  fEvent = InputEvent();

  /* Get the Qn vectors list */
  TList* qnlist = dynamic_cast<TList*>(GetInputData(1));
  if(!qnlist) return;

  Float_t values[MAXNOOFDATAVARIABLES] = {-999.};
  fDataBank = values;

  FillEventData();

  fEventPlaneHistos->FillHistClass("Event_NoCuts", values);
  if(!IsEventSelected(values)) return;
  fEventPlaneHistos->FillHistClass("Event_Analysis", values);

  TList* newTrkQvecList[nTrackDetectors] = {NULL};
  TList* newEPQvecList[nEPDetectors] = {NULL};
  AliQnCorrectionsQnVector* newTrk_qvec[nTrackDetectors] = {NULL};
  AliQnCorrectionsQnVector* newEP_qvec[nEPDetectors] = {NULL};


  /* get data structures for the different track detectors */
  for (Int_t iTrk = 0; iTrk < nTrackDetectors; iTrk++) {
    newTrkQvecList[iTrk] = dynamic_cast<TList*> (qnlist->FindObject(Form("%s", fTrackDetectorNameInFile[iTrk].Data())));
    if (newTrkQvecList[iTrk] == NULL) continue; /* the detector was not there */
    newTrk_qvec[iTrk] = (AliQnCorrectionsQnVector*) newTrkQvecList[iTrk]->FindObject(fExpectedCorrectionPass.Data());
    if (newTrk_qvec[iTrk] == NULL) {
      newTrk_qvec[iTrk] = (AliQnCorrectionsQnVector*) newTrkQvecList[iTrk]->FindObject(fAlternativeCorrectionPass.Data());
      if (newTrk_qvec[iTrk] == NULL) continue; /* neither the expected nor the alternative were there */
    }
  }

  /* and now for the EP detectors */
  for (Int_t iEP = 0; iEP < nEPDetectors; iEP++) {
    newEPQvecList[iEP] = dynamic_cast<TList*> (qnlist->FindObject(Form("%s",fEPDetectorNameInFile[iEP].Data())));
    if (newEPQvecList[iEP] == NULL) continue; /* the detector was not there */
    newEP_qvec[iEP] = (AliQnCorrectionsQnVector*) newEPQvecList[iEP]->FindObject(fExpectedCorrectionPass.Data());
    if (newEP_qvec[iEP] == NULL) {
      newEP_qvec[iEP] = (AliQnCorrectionsQnVector*) newEPQvecList[iEP]->FindObject(fAlternativeCorrectionPass.Data());
      if (newEP_qvec[iEP] == NULL) continue; /* neither the expected nor the alternative were there */
    }
  }

  /* now fill the Vn profiles with the proper data */
  /* TODO: correct for the normalization */
  for(Int_t iTrkDetector=0; iTrkDetector < nTrackDetectors; iTrkDetector++) {
    /* sanity check */
    if ((newTrkQvecList[iTrkDetector] != NULL) &&
        (newTrk_qvec[iTrkDetector] != NULL) &&
        (newTrk_qvec[iTrkDetector]->IsGoodQuality())) {
      for (Int_t iEPDetector=0; iEPDetector < nEPDetectors; iEPDetector++) {
        /*sanity check */
        if ((newEPQvecList[iEPDetector] != NULL) &&
            (newEP_qvec[iEPDetector] != NULL) &&
            (newEP_qvec[iEPDetector]->IsGoodQuality())) {
          for(Int_t h=0; h < kNharmonics; h++) {
            fVn[iTrkDetector*nEPDetectors+iEPDetector][h][0]->Fill(values[fCentralityVariable],
                    newTrk_qvec[iTrkDetector]->Qx(h+1) * newEP_qvec[iEPDetector]->QxNorm(h+1));
            fVn[iTrkDetector*nEPDetectors+iEPDetector][h][1]->Fill(values[fCentralityVariable],
                    newTrk_qvec[iTrkDetector]->Qx(h+1) * newEP_qvec[iEPDetector]->QyNorm(h+1));
            fVn[iTrkDetector*nEPDetectors+iEPDetector][h][2]->Fill(values[fCentralityVariable],
                    newTrk_qvec[iTrkDetector]->Qy(h+1) * newEP_qvec[iEPDetector]->QxNorm(h+1));
            fVn[iTrkDetector*nEPDetectors+iEPDetector][h][3]->Fill(values[fCentralityVariable],
                    newTrk_qvec[iTrkDetector]->Qy(h+1) * newEP_qvec[iEPDetector]->QyNorm(h+1));
          }
        }
      }
    }
  }

  /* now fill the correlation profiles needed for detector resolution */
  for (Int_t ix = 0; ix < fNDetectorResolutions; ix++) {
    /* sanity checks */
    if(fDetectorResolutionContributors[ix][0] == -1) continue;

    /* sanity checks */
    if((newTrkQvecList[fDetectorResolutionContributors[ix][0]] != NULL) &&
        (newTrk_qvec[fDetectorResolutionContributors[ix][0]] != NULL) &&
        (newEPQvecList[fDetectorResolutionContributors[ix][1]] != NULL) &&
        (newEP_qvec[fDetectorResolutionContributors[ix][1]] != NULL) &&
        (newEPQvecList[fDetectorResolutionContributors[ix][2]] != NULL) &&
        (newEP_qvec[fDetectorResolutionContributors[ix][2]] != NULL) &&
        (newTrk_qvec[fDetectorResolutionContributors[ix][0]]->GetN() != 0) &&
        (newEP_qvec[fDetectorResolutionContributors[ix][1]]->GetN() != 0) &&
        (newEP_qvec[fDetectorResolutionContributors[ix][2]]->GetN() != 0)) {

      for(Int_t h=0; h < kNharmonics; h++) {
        for (Int_t iCorr = 0; iCorr < nCorrelationPerDetector; iCorr++) {
          Float_t detectorOneValue;
          Float_t detectorTwoValue;
          switch (iCorr) {
          case kABXX: {
            detectorOneValue = newTrk_qvec[fDetectorResolutionContributors[ix][0]]->QxNorm(h+1);
            detectorTwoValue = newEP_qvec[fDetectorResolutionContributors[ix][1]]->QxNorm(h+1);
            break;
          }
          case kABYY: {
            detectorOneValue = newTrk_qvec[fDetectorResolutionContributors[ix][0]]->QyNorm(h+1);
            detectorTwoValue = newEP_qvec[fDetectorResolutionContributors[ix][1]]->QyNorm(h+1);
            break;
          }
          case kACXX: {
            detectorOneValue = newTrk_qvec[fDetectorResolutionContributors[ix][0]]->QxNorm(h+1);
            detectorTwoValue = newEP_qvec[fDetectorResolutionContributors[ix][2]]->QxNorm(h+1);
            break;
          }
          case kACYY: {
            detectorOneValue = newTrk_qvec[fDetectorResolutionContributors[ix][0]]->QyNorm(h+1);
            detectorTwoValue = newEP_qvec[fDetectorResolutionContributors[ix][2]]->QyNorm(h+1);
            break;
          }
          case kBCXX: {
            detectorOneValue = newEP_qvec[fDetectorResolutionContributors[ix][1]]->QxNorm(h+1);
            detectorTwoValue = newEP_qvec[fDetectorResolutionContributors[ix][2]]->QxNorm(h+1);
            break;
          }
          case kBCYY: {
            detectorOneValue = newEP_qvec[fDetectorResolutionContributors[ix][1]]->QyNorm(h+1);
            detectorTwoValue = newEP_qvec[fDetectorResolutionContributors[ix][2]]->QyNorm(h+1);
            break;
          }
          }
          fDetectorResolutionCorrelations[ix][iCorr][h]->Fill(values[fCentralityVariable],detectorOneValue*detectorTwoValue);
        }
      }
    }
  }
}  // end loop over events


//__________________________________________________________________
void AliAnalysisTaskQnVectorAnalysis::FinishTaskOutput()
{
  //
  // Finish Task
  //

  THashList* hList = (THashList*) fEventPlaneHistos->HistList();
  for(Int_t i=0; i<hList->GetEntries(); ++i) {
    THashList* list = (THashList*)hList->At(i);
    fEventQAList->Add(list);
  }

  /* TODO: correct vn with detector resolution */
  for(Int_t iTrkDetector=0; iTrkDetector < nTrackDetectors; iTrkDetector++) {
    for (Int_t iEPDetector=0; iEPDetector < nEPDetectors; iEPDetector++) {
      for(Int_t h=0; h < kNharmonics; h++) {
        for(Int_t corrComp =0; corrComp< kNcorrelationComponents; corrComp++){
          if (fVn[iTrkDetector*nEPDetectors+iEPDetector][h][corrComp]) {
            if (fVn[iTrkDetector*nEPDetectors+iEPDetector][h][corrComp]->GetEntries() > 0)
              fEventQAList->Add(fVn[iTrkDetector*nEPDetectors+iEPDetector][h][corrComp]);
          }
        }
      }
    }
  }

  for (Int_t ix = 0; ix < fNDetectorResolutions; ix++) {
    for (Int_t iCorr = 0; iCorr < nCorrelationPerDetector; iCorr++) {
      for(Int_t h=0; h < kNharmonics; h++) {
        fEventQAList->Add(fDetectorResolutionCorrelations[ix][iCorr][h]);
      }
    }
  }

  PostData(1, fEventQAList);
}



//__________________________________________________________________
Bool_t AliAnalysisTaskQnVectorAnalysis::IsEventSelected(Float_t* values) {
  if(!fEventCuts) return kTRUE;
  return fEventCuts->IsSelected(values);
}


