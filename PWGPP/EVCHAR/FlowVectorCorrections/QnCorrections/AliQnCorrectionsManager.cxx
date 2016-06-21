/**************************************************************************************************
 *                                                                                                *
 * Package:       FlowVectorCorrections                                                           *
 * Authors:       Jaap Onderwaater, GSI, jacobus.onderwaater@cern.ch                              *
 *                Ilya Selyuzhenkov, GSI, ilya.selyuzhenkov@gmail.com                             *
 *                Víctor González, UCM, victor.gonzalez@cern.ch                                   *
 *                Contributors are mentioned in the code where appropriate.                       *
 * Development:   2012-2016                                                                       *
 *                                                                                                *
 * This file is part of FlowVectorCorrections, a software package that corrects Q-vector          *
 * measurements for effects of nonuniform detector acceptance. The corrections in this package    *
 * are based on publication:                                                                      *
 *                                                                                                *
 *  [1] "Effects of non-uniform acceptance in anisotropic flow measurements"                      *
 *  Ilya Selyuzhenkov and Sergei Voloshin                                                         *
 *  Phys. Rev. C 77, 034904 (2008)                                                                *
 *                                                                                                *
 * The procedure proposed in [1] is extended with the following steps:                            *
 * (*) alignment correction between subevents                                                     *
 * (*) possibility to extract the twist and rescaling corrections                                 *
 *      for the case of three detector subevents                                                  *
 *      (currently limited to the case of two “hit-only” and one “tracking” detectors)            *
 * (*) (optional) channel equalization                                                            *
 * (*) flow vector width equalization                                                             *
 *                                                                                                *
 * FlowVectorCorrections is distributed under the terms of the GNU General Public License (GPL)   *
 * (https://en.wikipedia.org/wiki/GNU_General_Public_License)                                     *
 * either version 3 of the License, or (at your option) any later version.                        *
 *                                                                                                *
 **************************************************************************************************/
/// \file AliQnCorrectionsManager.cxx
/// \brief Implementation of the class AliQnCorrectionsManager

#include <TFile.h>
#include <TList.h>
#include <TKey.h>
#include "AliQnCorrectionsManager.h"
#include "AliLog.h"

#include <iostream>
#include <iomanip>
#include <fstream>

using std::cout;
using std::endl;
using std::setw;

/// \cond CLASSIMP
ClassImp(AliQnCorrectionsManager);
/// \endcond

const Int_t AliQnCorrectionsManager::nMaxNoOfDetectors = 32;
const Int_t AliQnCorrectionsManager::nMaxNoOfDataVariables = 2048;
const char *AliQnCorrectionsManager::szCalibrationHistogramsKeyName = "CalibrationHistograms";
const char *AliQnCorrectionsManager::szCalibrationQAHistogramsKeyName = "CalibrationQAHistograms";
const char *AliQnCorrectionsManager::szCalibrationNveQAHistogramsKeyName = "CalibrationQANveHistograms";
const char *AliQnCorrectionsManager::szDummyProcessListName = "dummyprocess";
const char *AliQnCorrectionsManager::szAllProcessesListName = "all data";

/// Default constructor.
/// The class owns the detectors and will be destroyed with it
AliQnCorrectionsManager::AliQnCorrectionsManager() :
    TObject(), fDetectorsSet(), fProcessListName(szDummyProcessListName) {

  fDetectorsSet.SetOwner(kTRUE);
  fDetectorsIdMap = NULL;
  fDataContainer = NULL;
  fCalibrationHistogramsList = NULL;
  fSupportHistogramsList = NULL;
  fQAHistogramsList = NULL;
  fNveQAHistogramsList = NULL;
  fQnVectorTree = NULL;
  fQnVectorList = NULL;
  fFillOutputHistograms = kFALSE;
  fFillQAHistograms = kFALSE;
  fFillNveQAHistograms = kFALSE;
  fFillQnVectorTree = kFALSE;
  fProcessesNames = NULL;
}

/// Default destructor
/// Deletes the memory taken
AliQnCorrectionsManager::~AliQnCorrectionsManager() {

  if (fDetectorsIdMap != NULL) delete [] fDetectorsIdMap;
  if (fDataContainer != NULL) delete [] fDataContainer;
  if (fCalibrationHistogramsList != NULL) delete fCalibrationHistogramsList;
  if (fProcessesNames != NULL) delete fProcessesNames;
}

/// Sets the base list that will own the input calibration histograms
/// \param calibrationFile the file
void AliQnCorrectionsManager::SetCalibrationHistogramsList(TFile *calibrationFile) {
  if (calibrationFile) {
    if (calibrationFile->GetListOfKeys()->GetEntries() > 0) {
      fCalibrationHistogramsList = (TList*)((TKey*)calibrationFile->GetListOfKeys()->FindObject(szCalibrationHistogramsKeyName))->ReadObj()->Clone();
      if (fCalibrationHistogramsList != NULL) {
        AliInfo(Form("Stored calibration list %s from file %s",
            fCalibrationHistogramsList->GetName(),
            calibrationFile->GetName()));
        /* we need the histograms ownership once we go to the GRID */
        fCalibrationHistogramsList->SetOwner(kTRUE);
      }
    }
  }
}



/// Adds a new detector
/// Checks for an already added detector and for a detector id
/// out of range. If so, gives a runtime error to inform of misuse.
/// \param detector the new detector to incorporate to the framework
void AliQnCorrectionsManager::AddDetector(AliQnCorrectionsDetector *detector) {
  if (detector->GetId() < nMaxNoOfDetectors) {
    if (fDetectorsSet.FindObject(detector->GetName())) {
      AliFatal(Form("You are trying to add twice %s detector with detector Id %d. FIX IT, PLEASE.",
          detector->GetName(),
          detector->GetId()));
      return;
    }
    fDetectorsSet.Add(detector);
    detector->AttachCorrectionsManager(this);
  }
  else {
    AliFatal(Form("You are trying to add %s detector with detector Id %d " \
        "while the highest Id supported is %d. FIX IT, PLEASE.",
        detector->GetName(),
        detector->GetId(),
        nMaxNoOfDetectors - 1));
    return;
  }
}

/// Searches for a concrete detector by name
/// \param name the name of the detector to find
/// \return pointer to the found detector (NULL if not found)
AliQnCorrectionsDetector *AliQnCorrectionsManager::FindDetector(const char *name) const {
  return (AliQnCorrectionsDetector *) fDetectorsSet.FindObject(name);
}

/// Searches for a concrete detector by detector id
/// \param id the id of the detector to find
/// \return pointer to the found detector (NULL if not found)
AliQnCorrectionsDetector *AliQnCorrectionsManager::FindDetector(Int_t id) const {
  AliQnCorrectionsDetector *detector = NULL;
  for (Int_t ixDetector = 0; ixDetector < fDetectorsSet.GetEntries(); ixDetector++) {
    detector = (AliQnCorrectionsDetector *) fDetectorsSet.At(ixDetector);
    if (detector->GetId() != id) {
      detector = NULL;
      continue;
    }
    else {
      break;
    }
  }
  return detector;
}


/// Searches for a concrete detector configuration by name
/// \param name the name of the detector configuration to find
/// \return pointer to the found detector configuration (NULL if not found)
AliQnCorrectionsDetectorConfigurationBase *AliQnCorrectionsManager::FindDetectorConfiguration(const char *name) const {
  AliQnCorrectionsDetectorConfigurationBase *detectorConfiguration;
  for (Int_t ixDetector = 0; ixDetector < fDetectorsSet.GetEntries(); ixDetector++) {
    detectorConfiguration = ((AliQnCorrectionsDetector *) fDetectorsSet.At(ixDetector))->FindDetectorConfiguration(name);
    if (detectorConfiguration != NULL) {
      return detectorConfiguration;
    }
  }
  return NULL;
}

/// Initializes the correction framework
void AliQnCorrectionsManager::InitializeQnCorrectionsFramework() {

  /* the data bank */
  fDataContainer = new Float_t[nMaxNoOfDataVariables];

  /* let's build the detectors map */
  fDetectorsIdMap = new AliQnCorrectionsDetector *[nMaxNoOfDetectors];
  AliQnCorrectionsDetector *detector = NULL;
  for (Int_t ixDetector = 0; ixDetector < fDetectorsSet.GetEntries(); ixDetector++) {
    detector = (AliQnCorrectionsDetector *) fDetectorsSet.At(ixDetector);
    fDetectorsIdMap[detector->GetId()] = detector;
  }


  /* create the support data structures */
  for (Int_t ixDetector = 0; ixDetector < fDetectorsSet.GetEntries(); ixDetector++) {
    ((AliQnCorrectionsDetector *) fDetectorsSet.At(ixDetector))->CreateSupportDataStructures();
  }

  /* build the support histograms list */
  fSupportHistogramsList = new TList();
  fSupportHistogramsList->SetName(szCalibrationHistogramsKeyName);
  fSupportHistogramsList->SetOwner(kTRUE);

  /* and the QA histograms list if needed */
  if (GetShouldFillQAHistograms()) {
    fQAHistogramsList = new TList();
    fQAHistogramsList->SetName(szCalibrationQAHistogramsKeyName);
    fQAHistogramsList->SetOwner(kTRUE);
    fNveQAHistogramsList = new TList();
    fNveQAHistogramsList->SetName(szCalibrationNveQAHistogramsKeyName);
    fNveQAHistogramsList->SetOwner(kTRUE);
  }


  /* build the support and QA histograms lists for the list of concurrent processes */
  if (fProcessesNames != NULL && fProcessesNames->GetEntries() != 0) {
    for (Int_t i = 0; i < fProcessesNames->GetEntries(); i++) {
      /* the support histgrams list */
      TList *newList = new TList();
      newList->SetName(((TObjString *) fProcessesNames->At(i))->GetName());
      newList->SetOwner(kTRUE);
      fSupportHistogramsList->Add(newList);

      /* the QA histograms list if needed */
      TList *newQAList = NULL;
      if (GetShouldFillQAHistograms()) {
        newQAList = new TList();
        newQAList->SetName(((TObjString *) fProcessesNames->At(i))->GetName());
        newQAList->SetOwner(kTRUE);
        fQAHistogramsList->Add(newQAList);
      }

      /* the non validated entries QA histograms list if needed */
      TList *newNveQAList = NULL;
      if (GetShouldFillNveQAHistograms()) {
        newNveQAList = new TList();
        newNveQAList->SetName(((TObjString *) fProcessesNames->At(i))->GetName());
        newNveQAList->SetOwner(kTRUE);
        fNveQAHistogramsList->Add(newNveQAList);
      }

      /* leave the selected process list name for a latter time */
      if (!fProcessListName.EqualTo(fProcessesNames->At(i)->GetName())) {
        /* build the support histograms list associated to the process */
        for (Int_t ixDetector = 0; ixDetector < fDetectorsSet.GetEntries(); ixDetector++) {
          ((AliQnCorrectionsDetector *) fDetectorsSet.At(ixDetector))->CreateSupportHistograms(newList);
        }
        /* the QA histograms list if needed */
        if (GetShouldFillQAHistograms()) {
          /* and pass it to the detectors for QA histograms creation */
          for (Int_t ixDetector = 0; ixDetector < fDetectorsSet.GetEntries(); ixDetector++) {
            ((AliQnCorrectionsDetector *) fDetectorsSet.At(ixDetector))->CreateQAHistograms(newQAList);
          }
        }
        /* the non validated entries QA histograms list if needed */
        if (GetShouldFillNveQAHistograms()) {
          /* and pass it to the detectors for non validated entries QA histograms creation */
          for (Int_t ixDetector = 0; ixDetector < fDetectorsSet.GetEntries(); ixDetector++) {
            ((AliQnCorrectionsDetector *) fDetectorsSet.At(ixDetector))->CreateNveQAHistograms(newNveQAList);
          }
        }
      }
    }
  }

  /* build the support and QA histograms lists associated to this process */
  /* and pass it to the detectors for support histograms creation */
  if (fProcessListName.Length() != 0) {
    /* let's see first whether we have the current process name within the processes names list */
    TList *processList;
    TList *processQAList;
    TList *processNveQAList;
    if (fProcessesNames != NULL && fProcessesNames->GetEntries() != 0 && fSupportHistogramsList->FindObject(fProcessListName) != NULL) {
      processList = (TList *) fSupportHistogramsList->FindObject(fProcessListName);
      /* the QA histograms list if needed */
      if (GetShouldFillQAHistograms()) {
        processQAList = (TList *) fQAHistogramsList->FindObject(fProcessListName);
      }
      /* the non validated entries QA histograms list if needed */
      if (GetShouldFillNveQAHistograms()) {
        processNveQAList = (TList *) fNveQAHistogramsList->FindObject(fProcessListName);
      }
    }
    else {
      processList = new TList();
      processList->SetName((const char *) fProcessListName);
      processList->SetOwner(kTRUE);
      /* we add it but probably temporarily */
      fSupportHistogramsList->Add(processList);
      /* the QA histograms list if needed */
      if (GetShouldFillQAHistograms()) {
        processQAList = new TList();
        processQAList->SetName((const char *) fProcessListName);
        processQAList->SetOwner(kTRUE);
        /* we add it but probably temporarily */
        fQAHistogramsList->Add(processQAList);
      }
      /* the non validated QA histograms list if needed */
      if (GetShouldFillNveQAHistograms()) {
        processNveQAList = new TList();
        processNveQAList->SetName((const char *) fProcessListName);
        processNveQAList->SetOwner(kTRUE);
        /* we add it but probably temporarily */
        fNveQAHistogramsList->Add(processNveQAList);
      }
    }
    /* now transfer the order to the defined detectors */
    /* so, we always create the histograms to use the latest ones */
    Bool_t retvalue = kTRUE;
    for (Int_t ixDetector = 0; ixDetector < fDetectorsSet.GetEntries(); ixDetector++) {
      retvalue = retvalue && ((AliQnCorrectionsDetector *) fDetectorsSet.At(ixDetector))->CreateSupportHistograms(processList);
      if (!retvalue)
        break;
    }
    if (!retvalue) {
      AliFatal("Failed to build the necessary support histograms.");
    }
    /* the QA histograms list if needed */
    if (GetShouldFillQAHistograms()) {
      /* pass it to the detectors for QA histograms creation */
      for (Int_t ixDetector = 0; ixDetector < fDetectorsSet.GetEntries(); ixDetector++) {
        ((AliQnCorrectionsDetector *) fDetectorsSet.At(ixDetector))->CreateQAHistograms(processQAList);
      }
    }
    /* the non validated QA histograms list if needed */
    if (GetShouldFillNveQAHistograms()) {
      /* pass it to the detectors for non validated entries QA histograms creation */
      for (Int_t ixDetector = 0; ixDetector < fDetectorsSet.GetEntries(); ixDetector++) {
        ((AliQnCorrectionsDetector *) fDetectorsSet.At(ixDetector))->CreateNveQAHistograms(processNveQAList);
      }
    }
  }
  else {
    AliFatal("The process label is missing.");
  }

  /* now get the process list on the calibration histograms list if any */
  /* and pass it to the detectors for input calibration histograms attachment, */
  /* we accept no calibration histograms list and no process list in case we */
  /* are in calibrating phase. We should TODO this. */
  if (fCalibrationHistogramsList != NULL) {
    TList *processList = (TList *)fCalibrationHistogramsList->FindObject((const char *)fProcessListName);
    if (processList != NULL) {
      AliInfo(Form("Assigned process list %s as the calibration histograms list",
          processList->GetName()));
      /* now transfer the order to the defined detectors */
      for (Int_t ixDetector = 0; ixDetector < fDetectorsSet.GetEntries(); ixDetector++) {
        ((AliQnCorrectionsDetector *) fDetectorsSet.At(ixDetector))->AttachCorrectionInputs(processList);
      }
    }
  }

  /* build the Qn vectors list */
  fQnVectorList = new TList();
  /* the list does not own the Qn vectors */
  fQnVectorList->SetOwner(kFALSE);
  /* pass it to the detectors for Qn vector creation and attachment */
  for (Int_t ixDetector = 0; ixDetector < fDetectorsSet.GetEntries(); ixDetector++) {
    ((AliQnCorrectionsDetector *) fDetectorsSet.At(ixDetector))->IncludeQnVectors(fQnVectorList);
  }
}

/// Set the name of the list that should be considered as assigned to the current process
/// If the stored process list name is the default one and the support histograms are
/// already created, change the list name and store the new name and get the new process
/// list on the calibration histograms list if any and pass it to the detectors for input
/// calibration histograms attachment.
/// Changing process list name on the fly during a  running process is not supported.
/// If the list of concurrent processes names is not empty, the new process name should be
/// in the list. If not a run time error is raised.
/// \param name the name of the list
void AliQnCorrectionsManager::SetCurrentProcessListName(const char *name) {
  if (fProcessListName.EqualTo(szDummyProcessListName)) {
    if (fSupportHistogramsList != NULL) {
      /* check the list of concurrent processes */
      if (fProcessesNames != NULL && fProcessesNames->GetEntries() != 0) {
        /* the new process name should be in the list of processes names */
        if (fSupportHistogramsList->FindObject(name) != NULL) {
          /* now we have to substitute the provisional process name list with the temporal one but renamed */
          TList *previousempty = (TList*) fSupportHistogramsList->FindObject(name);
          Int_t finalindex = fSupportHistogramsList->IndexOf(previousempty);
          fSupportHistogramsList->RemoveAt(finalindex);
          delete previousempty;
          TList *previoustemp = (TList *) fSupportHistogramsList->FindObject((const char *)fProcessListName);
          fSupportHistogramsList->Remove(previoustemp);
          previoustemp->SetName(name);
          fSupportHistogramsList->AddAt(previoustemp, finalindex);
        }
        else {
          /* nop! we raise an execution error */
          AliFatal(Form("The name of the process you want to run: %s, is not in the list of concurrent processes", name));
        }
        /* the QA histograms list if needed */
        if (GetShouldFillQAHistograms()) {
          /* the new process name should be in the list of processes names */
          if (fQAHistogramsList->FindObject(name) != NULL) {
            /* now we have to substitute the provisional process name list with the temporal one but renamed */
            TList *previousempty = (TList*) fQAHistogramsList->FindObject(name);
            Int_t finalindex = fQAHistogramsList->IndexOf(previousempty);
            fQAHistogramsList->RemoveAt(finalindex);
            delete previousempty;
            TList *previoustemp = (TList *) fQAHistogramsList->FindObject((const char *)fProcessListName);
            fQAHistogramsList->Remove(previoustemp);
            previoustemp->SetName(name);
            fQAHistogramsList->AddAt(previoustemp, finalindex);
          }
          else {
            /* nop! we raise an execution error */
            AliFatal(Form("The name of the process you want to run: %s, is not in the QA list of concurrent processes", name));
          }
        }
        /* the non validated entries QA histograms list if needed */
        if (GetShouldFillNveQAHistograms()) {
          /* the new process name should be in the list of processes names */
          if (fNveQAHistogramsList->FindObject(name) != NULL) {
            /* now we have to substitute the provisional process name list with the temporal one but renamed */
            TList *previousempty = (TList*) fNveQAHistogramsList->FindObject(name);
            Int_t finalindex = fNveQAHistogramsList->IndexOf(previousempty);
            fNveQAHistogramsList->RemoveAt(finalindex);
            delete previousempty;
            TList *previoustemp = (TList *) fNveQAHistogramsList->FindObject((const char *)fProcessListName);
            fNveQAHistogramsList->Remove(previoustemp);
            previoustemp->SetName(name);
            fNveQAHistogramsList->AddAt(previoustemp, finalindex);
          }
          else {
            /* nop! we raise an execution error */
            AliFatal(Form("The name of the process you want to run: %s, is not in the non validated entries QA list of concurrent processes", name));
          }
        }
      }
      else {
        TList *processList = (TList *) fSupportHistogramsList->FindObject((const char *)fProcessListName);
        processList->SetName(name);
        /* the QA histograms list if needed */
        if (GetShouldFillQAHistograms()) {
          TList *processQAList = (TList *) fQAHistogramsList->FindObject((const char *)fProcessListName);
          processQAList->SetName(name);
        }
        /* the non validated entries QA histograms list if needed */
        if (GetShouldFillNveQAHistograms()) {
          TList *processNveQAList = (TList *) fNveQAHistogramsList->FindObject((const char *)fProcessListName);
          processNveQAList->SetName(name);
        }
      }

      /* now get the process list on the calibration histograms list if any */
      /* and pass it to the detectors for input calibration histograms attachment, */
      /* we accept no calibration histograms list and no process list in case we */
      /* are in calibrating phase. We should TODO this. */
      fProcessListName = name;
      if (fCalibrationHistogramsList != NULL) {
        TList *processList = (TList *)fCalibrationHistogramsList->FindObject((const char *)fProcessListName);
        if (processList != NULL) {
          AliInfo(Form("Assigned process list %s as the calibration histograms list",
              processList->GetName()));
          /* now transfer the order to the defined detectors */
          for (Int_t ixDetector = 0; ixDetector < fDetectorsSet.GetEntries(); ixDetector++) {
            ((AliQnCorrectionsDetector *) fDetectorsSet.At(ixDetector))->AttachCorrectionInputs(processList);
          }
        }
      }
      /* build the Qn vectors list  now that all histograms are loaded */
      if (fQnVectorList == NULL) {
        /* first we build it if it isn't already there */
        fQnVectorList = new TList();
        /* the list does not own the Qn vectors */
        fQnVectorList->SetOwner(kFALSE);
      }
      /* pass it to the detectors for Qn vector creation and attachment */
      for (Int_t ixDetector = 0; ixDetector < fDetectorsSet.GetEntries(); ixDetector++) {
        ((AliQnCorrectionsDetector *) fDetectorsSet.At(ixDetector))->IncludeQnVectors(fQnVectorList);
      }
    }
    else {
      /* histograms list not yet created so, we just change the name */
      fProcessListName = name;
    }
  }
  else {
    AliFatal(Form("Changing process list name on the fly is not supported." \
        " Current name: %s, new name: %s", (const char *)fProcessListName, name));
  }

  /* now that we have everything let's print the configuration before we start */
  PrintFrameworkConfiguration();
}

/// Produce an understandable picture of current correction configuration
void AliQnCorrectionsManager::PrintFrameworkConfiguration() const {
  /* first get the list of detector configurations */
  TList *detectorList = new TList();
  detectorList->SetOwner(kTRUE);
  detectorList->SetName("Detector configurations list");
  /* pass it to the detectors for detector configurations name inclusion */
  for (Int_t ixDetector = 0; ixDetector < fDetectorsSet.GetEntries(); ixDetector++) {
    ((AliQnCorrectionsDetector *) fDetectorsSet.At(ixDetector))->FillDetectorConfigurationNameList(detectorList);
  }

  /* now the list of input correction steps */
  /* First we get an overall list of correction instances that we don't own*/
  TList *inputCorrections = new TList();
  inputCorrections->SetOwner(kFALSE);
  for (Int_t ixDetector = 0; ixDetector < fDetectorsSet.GetEntries(); ixDetector++) {
    ((AliQnCorrectionsDetector *) fDetectorsSet.At(ixDetector))->FillOverallInputCorrectionStepList(inputCorrections);
  }
  /* and now we build the correction step names list */
  TList *inputStepList = new TList();
  /* this one we own its items */
  inputStepList->SetOwner(kTRUE);
  inputStepList->SetName("Input data correction steps");
  for (Int_t i = 0; i < inputCorrections->GetEntries(); i++) {
    /* we got them in execution order which we keep */
    inputStepList->Add(new TObjString(inputCorrections->At(i)->GetName()));
  }
  /* enough for now */
  delete inputCorrections;

  /* now the list of Qn vector correction steps */
  /* First we get an overall list of correction instances that we don't own*/
  TList *vectorCorrections = new TList();
  vectorCorrections->SetOwner(kFALSE);
  for (Int_t ixDetector = 0; ixDetector < fDetectorsSet.GetEntries(); ixDetector++) {
    ((AliQnCorrectionsDetector *) fDetectorsSet.At(ixDetector))->FillOverallQnVectorCorrectionStepList(vectorCorrections);
  }
  /* and now we build the correction step names list */
  TList *vectorStepList = new TList();
  /* this one we own its items */
  vectorStepList->SetOwner(kTRUE);
  vectorStepList->SetName("Qn vector correction steps");
  for (Int_t i = 0; i < vectorCorrections->GetEntries(); i++) {
    /* we got them in execution order which we keep */
    vectorStepList->Add(new TObjString(vectorCorrections->At(i)->GetName()));
  }
  /* enough for now */
  delete vectorCorrections;

  /* and finally the list of correction steps applied to each detector configuration */
  TList *detectorCorrectionsList = new TList(); detectorCorrectionsList->SetOwner(kTRUE); detectorCorrectionsList->SetName("Assigned corrections");
  TList *detectorCalibratingCorrectionsList = new TList(); detectorCalibratingCorrectionsList->SetOwner(kTRUE); detectorCalibratingCorrectionsList->SetName("Calibrating corrections");
  TList *detectorApplyingCorrectionsList = new TList(); detectorApplyingCorrectionsList->SetOwner(kTRUE); detectorApplyingCorrectionsList->SetName("Applying corrections");
  /* pass it to the detectors for detector configuration correction steps inclusion */
  for (Int_t ixDetector = 0; ixDetector < fDetectorsSet.GetEntries(); ixDetector++) {
    ((AliQnCorrectionsDetector *) fDetectorsSet.At(ixDetector))->
        ReportOnCorrections(detectorCorrectionsList, detectorCalibratingCorrectionsList, detectorApplyingCorrectionsList);
  }

  /* let in a first try print everything in a quick way */
  if (kFALSE) {
    detectorList->Print("",-1);
    inputStepList->Print("",-1);
    vectorStepList->Print("",-1);
    detectorCorrectionsList->Print("",-1);
    detectorCalibratingCorrectionsList->Print("",-1);
    detectorApplyingCorrectionsList->Print("",-1);
  }

  /* get the steps involved and the current one */
  Int_t nNoOfSteps = inputStepList->GetEntries() + vectorStepList->GetEntries();
  Int_t nCurrentStep = 0;
  for (Int_t i = 0; i < detectorApplyingCorrectionsList->GetEntries(); i++) {
    if (nCurrentStep < ((TList*) detectorApplyingCorrectionsList->At(i))->GetEntries()) {
      nCurrentStep = ((TList*) detectorApplyingCorrectionsList->At(i))->GetEntries();
    }
  }

  /* let's take some parameters */
  size_t correctionFieldSize = 0;
  size_t detectorFieldSize = 0;
  size_t margin = 3;

  /* the input data correction steps */
  for (Int_t i = 0; i < inputStepList->GetEntries(); i++) {
    if (strlen(inputStepList->At(i)->GetName()) > correctionFieldSize) {
      correctionFieldSize = strlen(inputStepList->At(i)->GetName());
    }
  }
  /* the Qn vector correction steps */
  for (Int_t i = 0; i < vectorStepList->GetEntries(); i++) {
    if (strlen(vectorStepList->At(i)->GetName()) > correctionFieldSize) {
      correctionFieldSize = strlen(vectorStepList->At(i)->GetName());
    }
  }
  /* add pre-margin */
  correctionFieldSize += margin;

  /* the detector configurations */
  for (Int_t i = 0; i < detectorList->GetEntries(); i++) {
    if (strlen(detectorList->At(i)->GetName()) > detectorFieldSize) {
      detectorFieldSize = strlen(detectorList->At(i)->GetName());
    }
  }
  /* add pre-margin */
  detectorFieldSize += margin;

  TString correctionFieldLine('-',correctionFieldSize + margin);
  TString detectorsFieldLine('-', detectorList->GetEntries() * (detectorFieldSize + margin));
  TString correctionFieldSpace(' ',correctionFieldSize + margin);
  TString detectorsFieldSpace(' ', detectorList->GetEntries() * (detectorFieldSize + margin));

  TString line;
  Int_t textAnchor;
  /* the header */
  cout << correctionFieldLine << detectorsFieldLine << endl;
  line = Form("FLOW VECTOR FRAMEWORK (v2.0) - PASS %d/%d", nCurrentStep, nNoOfSteps);
  textAnchor = detectorsFieldLine.Length() / 2 + line.Length() / 2;
  cout << setw(correctionFieldSize+margin) << "|" << setw(textAnchor) << line << setw(detectorsFieldLine.Length() - textAnchor) << "|" << endl;
  cout << setw(correctionFieldSize+margin) << "-" << detectorsFieldLine << endl;
  /* detectors header */
  line = "--FLOW VECTORS--";
  textAnchor = detectorsFieldLine.Length() / 2 + line.Length() / 2;
  cout << setw(correctionFieldSize+margin) << "|" << setw(textAnchor) << line << setw(detectorsFieldLine.Length() - textAnchor) << "|" << endl;
  cout << setw(correctionFieldSize+margin) << "|";
  for (Int_t i = 0; i < detectorList->GetEntries(); i++) {
    cout << setw(detectorFieldSize) << detectorList->At(i)->GetName() << setw(margin) << "|";
  }
  cout << endl;

  /* now the correction steps */
  cout << setw(correctionFieldSize) << "CORRECTIONS" << setw(margin) << "-" << detectorsFieldLine << endl;
  /* first the input correction steps */
  for (Int_t ixInCorr = 0; ixInCorr < inputStepList->GetEntries(); ixInCorr++) {
    cout << setw(correctionFieldSize) << inputStepList->At(ixInCorr)->GetName() << setw(margin) << "|";
    for (Int_t ixDet = 0; ixDet < detectorList->GetEntries(); ixDet++) {
      TString detOut = "-";
      if (detectorCorrectionsList->FindObject(detectorList->At(ixDet)->GetName()) != NULL) {
        if (((TList*) detectorCorrectionsList->FindObject(detectorList->At(ixDet)->GetName()))->FindObject(inputStepList->At(ixInCorr)->GetName())) {
          detOut = "0";
          if (detectorApplyingCorrectionsList->FindObject(detectorList->At(ixDet)->GetName()) != NULL) {
            if (((TList*) detectorApplyingCorrectionsList->FindObject(detectorList->At(ixDet)->GetName()))->FindObject(inputStepList->At(ixInCorr)->GetName())) {
              detOut = "x";
            }
          }
        }
      }
      cout << setw(detectorFieldSize) << detOut << setw(margin) << "|";
    }
    cout << endl;
  }
  /* now the Qn vector correction steps */
  for (Int_t ixInCorr = 0; ixInCorr < vectorStepList->GetEntries(); ixInCorr++) {
    cout << setw(correctionFieldSize) << vectorStepList->At(ixInCorr)->GetName() << setw(margin) << "|";
    for (Int_t ixDet = 0; ixDet < detectorList->GetEntries(); ixDet++) {
      TString detOut = "-";
      if (detectorCorrectionsList->FindObject(detectorList->At(ixDet)->GetName()) != NULL) {
        if (((TList*) detectorCorrectionsList->FindObject(detectorList->At(ixDet)->GetName()))->FindObject(vectorStepList->At(ixInCorr)->GetName())) {
          detOut = "0";
          if (detectorApplyingCorrectionsList->FindObject(detectorList->At(ixDet)->GetName()) != NULL) {
            if (((TList*) detectorApplyingCorrectionsList->FindObject(detectorList->At(ixDet)->GetName()))->FindObject(vectorStepList->At(ixInCorr)->GetName())) {
              detOut = "x";
            }
          }
        }
      }
      cout << setw(detectorFieldSize) << detOut << setw(margin) << "|";
    }
    cout << endl;
  }

  /* now the calibration histograms */
  cout << setw(correctionFieldSize) << "FILL HISTS" << setw(margin) << "-" << detectorsFieldLine << endl;
  /* first the input correction steps */
  for (Int_t ixInCorr = 0; ixInCorr < inputStepList->GetEntries(); ixInCorr++) {
    cout << setw(correctionFieldSize) << inputStepList->At(ixInCorr)->GetName() << setw(margin) << "|";
    for (Int_t ixDet = 0; ixDet < detectorList->GetEntries(); ixDet++) {
      TString detOut = "-";
      if (detectorCorrectionsList->FindObject(detectorList->At(ixDet)->GetName()) != NULL) {
        if (((TList*) detectorCorrectionsList->FindObject(detectorList->At(ixDet)->GetName()))->FindObject(inputStepList->At(ixInCorr)->GetName())) {
          detOut = "0";
          if (detectorCalibratingCorrectionsList->FindObject(detectorList->At(ixDet)->GetName()) != NULL) {
            if (((TList*) detectorCalibratingCorrectionsList->FindObject(detectorList->At(ixDet)->GetName()))->FindObject(inputStepList->At(ixInCorr)->GetName())) {
              detOut = "x";
            }
          }
        }
      }
      cout << setw(detectorFieldSize) << detOut << setw(margin) << "|";
    }
    cout << endl;
  }
  /* now the Qn vector correction steps */
  for (Int_t ixInCorr = 0; ixInCorr < vectorStepList->GetEntries(); ixInCorr++) {
    cout << setw(correctionFieldSize) << vectorStepList->At(ixInCorr)->GetName() << setw(margin) << "|";
    for (Int_t ixDet = 0; ixDet < detectorList->GetEntries(); ixDet++) {
      TString detOut = "-";
      if (detectorCorrectionsList->FindObject(detectorList->At(ixDet)->GetName()) != NULL) {
        if (((TList*) detectorCorrectionsList->FindObject(detectorList->At(ixDet)->GetName()))->FindObject(vectorStepList->At(ixInCorr)->GetName())) {
          detOut = "0";
          if (detectorCalibratingCorrectionsList->FindObject(detectorList->At(ixDet)->GetName()) != NULL) {
            if (((TList*) detectorCalibratingCorrectionsList->FindObject(detectorList->At(ixDet)->GetName()))->FindObject(vectorStepList->At(ixInCorr)->GetName())) {
              detOut = "x";
            }
          }
        }
      }
      cout << setw(detectorFieldSize) << detOut << setw(margin) << "|";
    }
    cout << endl;
  }

  /* finally the legend */
  cout << correctionFieldLine << detectorsFieldLine << endl;
  line = "x: this pass      0: future pass        -: N/A";
  textAnchor = detectorsFieldLine.Length() / 2 + line.Length() / 2;
  cout << setw(correctionFieldSize) << "Legend" << setw(margin) << "|" << setw(textAnchor) << line << setw(detectorsFieldLine.Length() - textAnchor) << "|" << endl;
  cout << correctionFieldLine << detectorsFieldLine << endl;
  cout << endl;

  /* back to clean */
  delete detectorList;
  delete inputStepList;
  delete vectorStepList;
  delete detectorCorrectionsList;
  delete detectorCalibratingCorrectionsList;
  delete detectorApplyingCorrectionsList;
}


/// Produce the final output and release the framework.
/// Produce the all data lists that collect data from all concurrent processes.
void AliQnCorrectionsManager::FinalizeQnCorrectionsFramework() {

  TList *processList = (TList *) fSupportHistogramsList->FindObject((const char *)fProcessListName);
  fSupportHistogramsList->Add(processList->Clone(szAllProcessesListName));
  /* the QA histograms list if needed */
  if (GetShouldFillQAHistograms()) {
    TList *processQAList = (TList *) fQAHistogramsList->FindObject((const char *)fProcessListName);
    fQAHistogramsList->Add(processQAList->Clone(szAllProcessesListName));
  }
  /* the non validated entries QA histograms list if needed */
  if (GetShouldFillNveQAHistograms()) {
    TList *processNveQAList = (TList *) fQAHistogramsList->FindObject((const char *)fProcessListName);
    fNveQAHistogramsList->Add(processNveQAList->Clone(szAllProcessesListName));
  }
}




