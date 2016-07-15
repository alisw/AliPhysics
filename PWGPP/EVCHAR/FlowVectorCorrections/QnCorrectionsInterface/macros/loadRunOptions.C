/**************************************************************************
 * Copyright(c) 2013-2016, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include <stdio.h>
#include <iostream>
#include <fstream>
#include "TSystem.h"
#include "TROOT.h"

using std::cout;
using std::endl;

#include "runAnalysis.H"


Bool_t bOptionsLoaded = kFALSE;
void CleanOptions();

/* tasks to execute: keep the same order than in teh runoptions file */
const Int_t nNoOfTasks = 7;
const char *taskNameInFile[nNoOfTasks] = {
    "Tender", "CDB", "PhysicsSelection", "PIDResponse", "PIDCombinedTask", "CentralityTask", "MultiplicityTask"
};
Bool_t *taskVar[nNoOfTasks] = {
    &bUseTender, &bUseCDB, &bUsePhysicsSelection, &bUsePIDResponse, &bRunPIDCombinedTask, &bUseCentralityTask, &bUseMultiplicityTask
};

const Int_t nNoOfDetectors = 7;
const char *detectorNameInFile[nNoOfDetectors] = {
    "TPC", "SPD", "VZERO", "TZERO", "FMD", "rawFMD", "ZDC"
};
Bool_t *detectorVar[nNoOfDetectors] = {
    &bUseTPC, &bUseSPD, &bUseVZERO, &bUseTZERO, &bUseFMD, &bUseRawFMD, &bUseZDC
};

/// \brief get the condition of "Use name: yes/no"
/// The search for "Use name:" is started on the first read line and on.
/// The file name is left after reading the line that contains "Use name:"
/// or the last read line if not found.
/// \param optionsfile the options file
/// \param name the name to insert to search for its use
/// \return 1 if use = yes, 0 if use = no, -1 if not found or error
Int_t getUse(const ifstream &optionsfile, const char *name) {
  TString currline;
  currline.ReadLine(optionsfile);
  while (currline.BeginsWith("#") || currline.IsWhitespace()) currline.ReadLine(optionsfile);
  if (currline.BeginsWith(Form("Use %s: ",name))) {
    currline.Remove(0,strlen(Form("Use %s: ",name)));
    if (currline.Contains("yes"))
      return 1;
    else if (currline.Contains("no"))
      return 0;
    else
      { printf("ERROR: wrong Use %s option in options file\n", name); return -1; }
  }
  else
    { printf("ERROR: wrong Use %s option in options file\n", name); return -1; }
}

/// \brief Load the run options for the current task
/// The run options as present in the input file are taken over the global variables
/// \param verb flag for verbose output
/// \param filename the options filename
Bool_t loadRunOptions(Bool_t verb,const char *filename) {
  ifstream optionsfile;
  TString currline;

  if (bOptionsLoaded) CleanOptions();

  listOfRuns.SetOwner(kTRUE);
  optionsfile.open(filename);
  currline.ReadLine(optionsfile);
  while (currline.BeginsWith("#") || currline.IsWhitespace()) currline.ReadLine(optionsfile);

  /* Load run options */
  if (!currline.EqualTo("Run options:")) { printf("ERROR: wrong run options file %s\n", filename); return kFALSE; }
  /* run option */
  currline.ReadLine(optionsfile);
  while (currline.BeginsWith("#") || currline.IsWhitespace()) currline.ReadLine(optionsfile);
  if (currline.EqualTo("grid")) {
    bGRIDPlugin = kTRUE;
    bTrainScope = kFALSE;
  }
  else if (currline.EqualTo("local")) {
    bGRIDPlugin = kFALSE;
    bTrainScope = kFALSE;
  }
  else if (currline.EqualTo("train")) {
    bGRIDPlugin = kFALSE;
    bTrainScope = kTRUE;
  }
  else
    { printf("ERROR: wrong run option in options file %s\n", filename); return kFALSE; }
  printf("  Running in %s\n", bTrainScope ? "train" : (bGRIDPlugin ? "grid" : "local"));

  /* MC option */
  currline.ReadLine(optionsfile);
  while (currline.BeginsWith("#") || currline.IsWhitespace()) currline.ReadLine(optionsfile);
  if (currline.EqualTo("real"))
    bMC = kFALSE;
  else if (currline.EqualTo("MC"))
    bMC = kTRUE;
  else
    { printf("ERROR: wrong MC option in options file %s\n", filename); return kFALSE; }
  printf("  with %s data\n", bMC ? "simulated MC" : "real");

  /* SW versions */
  currline.ReadLine(optionsfile);
  while (currline.BeginsWith("#") || currline.IsWhitespace()) currline.ReadLine(optionsfile);
  if (currline.BeginsWith("AliPhysicsVersion: ")) {
    szAliPhysicsVersion = currline.Remove(0,strlen("AliPhysicsVersion: "));
  }
  else
    { printf("ERROR: wrong Aliphysics SW version in options file %s\n", filename); return kFALSE; }
  if (bGRIDPlugin)
    cout << "    Grid AliPhysics version: " << szAliPhysicsVersion << endl;

  /* tasks to execute */
  for (Int_t itask = 0; itask < nNoOfTasks; itask++) {
    Int_t use = getUse(optionsfile,taskNameInFile[itask]);
    switch (use) {
    case 0: *taskVar[itask] = kFALSE; break;
    case 1: *taskVar[itask] = kTRUE; break;
    default: return kFALSE;
    }
    printf("  Use %s: %s\n", taskNameInFile[itask], *taskVar[itask] ? "yes" : "no");
  }


  /* Execution conditions */
  currline.ReadLine(optionsfile);
  while (currline.BeginsWith("#") || currline.IsWhitespace()) currline.ReadLine(optionsfile);
  if (currline.BeginsWith("Use multiplicity: ")) {
    currline.Remove(0,strlen("Use multiplicity: "));
    if (currline.Contains("yes"))
      bUseMultiplicity = kTRUE;
    else if (currline.Contains("no"))
      bUseMultiplicity = kFALSE;
    else
      { printf("ERROR: wrong Use multiplicity option in options file %s\n", filename); return kFALSE; }
  }
  else
    { printf("ERROR: wrong Use multiplicity option in options file %s\n", filename); return kFALSE; }
  printf("  Use multiplicity: %s\n", bUseMultiplicity ? "yes" : "no");

  currline.ReadLine(optionsfile);
  while (currline.BeginsWith("#") || currline.IsWhitespace()) currline.ReadLine(optionsfile);
  if (currline.BeginsWith("Is 2015 dataset")) {
    currline.Remove(0,strlen("Is 2015 dataset"));
    if (currline.Contains("yes"))
      b2015DataSet = kTRUE;
    else if (currline.Contains("no"))
      b2015DataSet = kFALSE;
    else
      { printf("ERROR: wrong Is 2015 dataset option in options file %s\n", filename); return kFALSE; }
  }
  else
    { printf("ERROR: wrong Is 2015 dataset option in options file %s\n", filename); return kFALSE; }
  printf("  Is 2015 dataset: %s\n", b2015DataSet ? "yes" : "no");

  /* ESD input format */
  Int_t useESDinput = getUse(optionsfile, "ESD");
  switch (useESDinput) {
  case 0: bUseESD = kFALSE; break;
  case 1: bUseESD = kTRUE; break;
  default: return kFALSE;
  }
  printf(" Use ESD: %s\n", bUseESD ? "yes" : "no");

  /* AOD input format */
  Int_t useAODinput = getUse(optionsfile, "AOD");
  switch (useAODinput) {
  case 0: bUseAOD = kFALSE; break;
  case 1: bUseAOD = kTRUE; break;
  default: return kFALSE;
  }
  printf(" Use AOD: %s\n", bUseAOD ? "yes" : "no");

  /* task level cuts */
  currline.ReadLine(optionsfile);
  while (currline.BeginsWith("#") || currline.IsWhitespace()) currline.ReadLine(optionsfile);
  if (currline.EqualTo("Task level:")) {
    printf(" Task cuts: \n");
    currline.ReadLine(optionsfile);
    while (currline.BeginsWith("#") || currline.IsWhitespace()) currline.ReadLine(optionsfile);
    while(!currline.EqualTo("end")) {
      /* centrality cuts */
      if (currline.BeginsWith("Centrality: ")) {
        currline.Remove(0, strlen("Centrality: "));
        if (verb) cout << "    Centrality cut string: " << currline << endl;
        Double_t min;
        Double_t max;
        sscanf(currline.Data(), "%lf-%lf", &min, &max);
        centralityMin = min;
        centralityMax = max;
        printf ("      centrality cut: %.1f-%.1f\n",centralityMin,centralityMax);
        currline.ReadLine(optionsfile);
        while (currline.BeginsWith("#") || currline.IsWhitespace()) currline.ReadLine(optionsfile);
      } /* end centrality cuts */

      /* z Vertex cuts */
      if (currline.BeginsWith("zvertex: ")) {
        currline.Remove(0, strlen("zvertex: "));
        if (verb) cout << "    z vertex cut string: " << currline << endl;
        Double_t min;
        Double_t max;
        sscanf(currline.Data(), "%lf:%lf", &min, &max);
        zvertexMin = min;
        zvertexMax = max;
        printf ("      z vertex cut: %.1f:%.1f\n",zvertexMin,zvertexMax);
        currline.ReadLine(optionsfile);
        while (currline.BeginsWith("#") || currline.IsWhitespace()) currline.ReadLine(optionsfile);
      } /* end centrality cuts */
    }
  }
  else
    { printf("ERROR: wrong Task cuts section in options file %s\n", filename); return kFALSE; }

  /* now the detector selected */
  currline.ReadLine(optionsfile);
  while (currline.BeginsWith("#") || currline.IsWhitespace()) currline.ReadLine(optionsfile);
  if (currline.EqualTo("Detectors:")) {
    printf(" Detectors: \n");
    for (Int_t idetector = 0; idetector < nNoOfDetectors; idetector++) {
      Int_t use = getUse(optionsfile,detectorNameInFile[idetector]);
      switch (use) {
      case 0: *detectorVar[idetector] = kFALSE; break;
      case 1: *detectorVar[idetector] = kTRUE; break;
      default: return kFALSE;
      }
      printf("  Use %s: %s\n", detectorNameInFile[idetector], *detectorVar[idetector] ? "yes" : "no");
    }
  }
  else
    { printf("ERROR: wrong Detectors section in options file %s\n", filename); return kFALSE; }

  /* now the corrections file */
  currline.ReadLine(optionsfile);
  while (currline.BeginsWith("#") || currline.IsWhitespace()) currline.ReadLine(optionsfile);
  if (!currline.EqualTo("Corrections file:")) { printf("ERROR: wrong corrections file location in options file %s\n", filename); return kFALSE; }
  printf(" Corrections file:\n");
  /* source option */
  currline.ReadLine(optionsfile);
  while (currline.BeginsWith("#") || currline.IsWhitespace()) currline.ReadLine(optionsfile);
  if (currline.BeginsWith("source: ")) {
    currline.Remove(0,strlen("source: "));
    if (currline.Contains("local"))
      szCorrectionsSource = "local";
    else if (currline.Contains("aliensingle"))
      szCorrectionsSource = "aliensingle";
    else if (currline.Contains("alienmultiple"))
      szCorrectionsSource = "alienmultiple";
    else if (currline.Contains("OADBsingle"))
      szCorrectionsSource = "OADBsingle";
    else if (currline.Contains("OADBmultiple"))
      szCorrectionsSource = "OADBmultiple";
    else
      { printf("ERROR: wrong corrections file source in options file %s\n", filename); return kFALSE; }
  }
  else
    { printf("ERROR: wrong corrections file source in options file %s\n", filename); return kFALSE; }
  printf("  source: %s\n", (const char *) szCorrectionsSource);
  /* path option */
  currline.ReadLine(optionsfile);
  while (currline.BeginsWith("#") || currline.IsWhitespace()) currline.ReadLine(optionsfile);
  if (currline.BeginsWith("path: ")) {
    currline.Remove(0,strlen("path: "));
    szCorrectionsFilePath = currline;
  }
  else
    { printf("ERROR: wrong corrections file path in options file %s\n", filename); return kFALSE; }
  printf("  path: %s\n", (const char *) szCorrectionsFilePath);

  /* file name option */
  currline.ReadLine(optionsfile);
  while (currline.BeginsWith("#") || currline.IsWhitespace()) currline.ReadLine(optionsfile);
  if (currline.BeginsWith("filename: ")) {
    currline.Remove(0,strlen("filename: "));
    szCorrectionsFileName = currline;
  }
  else
    { printf("ERROR: wrong corrections file name in options file %s\n", filename); return kFALSE; }
  printf("  filename: %s\n", (const char *) szCorrectionsFileName);

  /* the Qn vector analysis task */
  Int_t useQnTask = getUse(optionsfile, "QnVectorAnalysisTask");
  switch (useQnTask) {
  case 0: bRunQnVectorAnalysisTask = kFALSE; break;
  case 1: bRunQnVectorAnalysisTask = kTRUE; break;
  default: return kFALSE;
  }
  printf(" Use QnVectorAnalysisTask: %s\n", bRunQnVectorAnalysisTask ? "yes" : "no");

  /* the expected correction step on the Qn vector */
  currline.ReadLine(optionsfile);
  while (currline.BeginsWith("#") || currline.IsWhitespace()) currline.ReadLine(optionsfile);
  if (currline.BeginsWith("Expected correction step: ")) {
    currline.Remove(0,strlen("Expected correction step: "));
    if (currline.EqualTo("raw")) szCorrectionPass = "raw";
    else if (currline.EqualTo("plain")) szCorrectionPass = "plain";
    else if (currline.EqualTo("rec")) szCorrectionPass = "rec";
    else if (currline.EqualTo("align")) szCorrectionPass = "align";
    else if (currline.EqualTo("twist")) szCorrectionPass = "twist";
    else if (currline.EqualTo("scale")) szCorrectionPass = "scale";
    else if (currline.EqualTo("latest")) szCorrectionPass = "latest";
    else
      { printf("ERROR: wrong Expected correction step %s in options file %s\n", currline.Data(), filename); return kFALSE; }
  }
  else
    { printf("ERROR: wrong Expected correction step option in options file %s\n", filename); return kFALSE; }
  printf("  Expected correction step: %s\n", szCorrectionPass.Data());

  /* the alternative expected correction step on the Qn vector */
  currline.ReadLine(optionsfile);
  while (currline.BeginsWith("#") || currline.IsWhitespace()) currline.ReadLine(optionsfile);
  if (currline.BeginsWith("Alternative correction step: ")) {
    currline.Remove(0,strlen("Alternative correction step: "));
    if (currline.EqualTo("raw")) szAltCorrectionPass = "raw";
    else if (currline.EqualTo("plain")) szAltCorrectionPass = "plain";
    else if (currline.EqualTo("rec")) szAltCorrectionPass = "rec";
    else if (currline.EqualTo("align")) szAltCorrectionPass = "align";
    else if (currline.EqualTo("twist")) szAltCorrectionPass = "twist";
    else if (currline.EqualTo("scale")) szAltCorrectionPass = "scale";
    else if (currline.EqualTo("latest")) szAltCorrectionPass = "latest";
    else
      { printf("ERROR: wrong Alternative correction step %s in options file %s\n", currline.Data(), filename); return kFALSE; }
  }
  else
    { printf("ERROR: wrong Alternative correction step option in options file %s\n", filename); return kFALSE; }
  printf("  Alternative correction step: %s\n", szAltCorrectionPass.Data());

  /* closing the options file */
  if (verb) printf(" Closing the options file\n");
  optionsfile.close();
  if (verb) printf(" Options file closed\n");

  /* now get the data files location */
  if (bGRIDPlugin || bTrainScope) {
    TString szDataLocFile;
    ifstream datalocfile;

    if (bMC) { szDataLocFile = "GRIDMCdata.txt"; }
    else { szDataLocFile = "GRIDrealdata.txt"; }

    if (verb) printf(" Opening the data location file: %s\n", szDataLocFile.Data());
    datalocfile.open(szDataLocFile.Data());
    if (datalocfile.fail()) { cout << "ERROR: cannot open location file " << szDataLocFile << endl; return kFALSE; }

    /* the data pattern */
    currline.ReadLine(datalocfile);
    while (currline.BeginsWith("#") || currline.IsWhitespace()) currline.ReadLine(datalocfile);
    if (currline.BeginsWith("Data pattern: ")) {
      szDataPattern = currline.Remove(0,strlen("Data pattern: "));
    }
    else
      { cout << "ERROR: wrong Data pattern in data location file " << szDataLocFile << endl; return kFALSE; }
    cout << " GRID data pattern: " << szDataPattern << endl;
    /* the grid data dir */
    currline.ReadLine(datalocfile);
    while (currline.BeginsWith("#") || currline.IsWhitespace()) currline.ReadLine(datalocfile);
    if (currline.BeginsWith("Grid data dir: ")) {
      szDataDir = currline.Remove(0,strlen("Grid data dir: "));
    }
    else
      { cout << "ERROR: wrong GRID Data dir in data location file " << szDataLocFile << endl; return kFALSE; }
    cout << " GRID data dir: " << szDataDir << endl;
    /* the number of input files */
    currline.ReadLine(datalocfile);
    while (currline.BeginsWith("#") || currline.IsWhitespace()) currline.ReadLine(datalocfile);
    if (currline.BeginsWith("Number of input files: ")) {
      currline.Remove(0,strlen("Number of input files: "));
      nNoOfInputFiles = currline.Atoi();
    }
    else
      { cout << "ERROR: wrong number of input files in data location file " << szDataLocFile << endl; return kFALSE; }
    cout << "  GRID number of input files: " << nNoOfInputFiles << endl;
    /* the run prefix */
    currline.ReadLine(datalocfile);
    while (currline.BeginsWith("#") || currline.IsWhitespace()) currline.ReadLine(datalocfile);
    if (currline.BeginsWith("Run prefix: ")) {
      szRunPrefix = currline.Remove(0,strlen("Run prefix: "));
    }
    else
      { cout << "ERROR: wrong Run prefix in data location file " << szDataLocFile << endl; return kFALSE; }
    if (szRunPrefix.IsWhitespace())
      printf("  NO run prefix\n");
    else
      cout << "  Run prefix: " << szRunPrefix << endl;

    /* now the list of runs to process */
    TObjString *pszRunNo;
    TObjString *pszActiveRunNo;
    currline.ReadLine(datalocfile);
    while (!currline.IsWhitespace()) {
      if (!currline.BeginsWith("#")) {
        pszRunNo = new TObjString(currline);
        if (pszRunNo->GetString().BeginsWith(">")) {
          /* the run is included in the current analysis */
          pszActiveRunNo = new TObjString(pszRunNo->GetString().Remove(0,1));
          listOfActiveRuns.Add(pszActiveRunNo);
          listOfRuns.Add(pszActiveRunNo->Clone());
          delete pszRunNo;
        }
        else {
          /* the run is included as a potential concurrent run */
          listOfRuns.Add(pszRunNo);
        }
      }
      currline.ReadLine(datalocfile);
    }
    printf("  List of concurrent runs: ");
    for (Int_t i=0;i<listOfRuns.GetEntriesFast();i++) {
      if (i%8 == 0) {
        printf("\n    ");
      }
      printf("%s, ", ((TObjString*) listOfRuns.At(i))->GetString().Data());
    }
    printf("\n");
    printf("  List of runs to analyze: ");
    for (Int_t i=0;i<listOfActiveRuns.GetEntriesFast();i++) {
      if (i%8 == 0) {
        printf("\n    ");
      }
      printf("%s, ", ((TObjString*) listOfActiveRuns.At(i))->GetString().Data());
    }
    printf("\n");
    if (verb) printf(" Closing the data location file: %s\n", szDataLocFile.Data());
    datalocfile.close();
  }
  else {
    if (bMC) {
      szLocalFileList = "filelist_mc.txt";
    }
    else
      szLocalFileList = "filelist.txt";
    printf("  Local data list file: %s\n", (const char*) szLocalFileList);
  }
  bOptionsLoaded = kTRUE;
  return kTRUE;
}

void CleanOptions() {
  listOfRuns.Clear();
  listOfActiveRuns.Clear();
}
