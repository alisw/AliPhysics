/* Copyright(c) 2013-2016, ALICE Experiment at CERN, All rights reserved. *
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

#ifndef RUNANALYSIS_H_
#define RUNANALYSIS_H_

Bool_t bGRIDPlugin;
Bool_t bTrainScope;
Bool_t bMC;

/* tasks to use */
/* only valid if not under train scope */
Bool_t bUseTender;
Bool_t bUseCDB;
Bool_t bUsePhysicsSelection;
Bool_t bUsePIDResponse;
Bool_t bRunPIDCombinedTask;
Bool_t bUseCentralityTask;
Bool_t bUseMultiplicityTask;

/* temporal flag to use multiplicity instead of centrality and to inhibit detectors for 2015 dataset */
Bool_t bUseMultiplicity;
Bool_t b2015DataSet;

/* data input */
Bool_t bUseESD;
Bool_t bUseAOD;


/* Task level cuts */
Double_t centralityMin;
Double_t centralityMax;
Double_t zvertexMin;
Double_t zvertexMax;
Double_t bUseOnlyCentCalibEvents;


/* Running conditions */
TString szAliPhysicsVersion;
TString szDataPattern;
TString szDataDir;
TString szRunPrefix;
Int_t nNoOfInputFiles;
TObjArray listOfRuns;
TObjArray listOfActiveRuns;
TString szCorrectionsSource;
TString szCorrectionsFilePath;
TString szCorrectionsFileName;

/* detectors */
Bool_t bUseTPC;
Bool_t bUseSPD;
Bool_t bUseVZERO;
Bool_t bUseTZERO;
Bool_t bUseFMD;
Bool_t bUseRawFMD;
Bool_t bUseZDC;

/* run the Qn vector analysis task */
Bool_t bRunQnVectorAnalysisTask;
TString szCorrectionPass;
TString szAltCorrectionPass;

TString szLocalFileList;

Bool_t loadRunOptions(Bool_t verb = kFALSE, const char *path = ".");

#endif /* RUNANALYSIS_H_ */
