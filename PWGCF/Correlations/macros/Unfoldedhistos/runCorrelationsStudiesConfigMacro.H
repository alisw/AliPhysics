#ifndef RUNCORRELATIONSTUDIESCONFIG_H_
#define RUNCORRELATIONSTUDIESCONFIG_H_

#include <TString.h>
#include <TObjArray.h>

Bool_t bASCIIoutput;

Bool_t bGRIDPlugin;
Bool_t bTrainScope;
Bool_t bMC;
Bool_t bMConlyTruth;

Bool_t isHeavyIon;

/* tasks to use */
/* only valid if not under train scope */
Bool_t bUsePhysicsSelection;
Bool_t bUsePIDResponse;
Bool_t bUseCentralityTask;
Bool_t bUseMultiplicityTask;

/* data input */
Bool_t bUseESD;
Bool_t bUseAOD;
TString szpass;

/* Task level cuts */
Double_t centralityMin;
Double_t centralityMax;
Double_t zvertexMin;
Double_t zvertexMax;
Double_t bUseOnlyCentCalibEvents;

/* Running conditions */
TString szAliPhysicsVersion;
TString szGridWorkingDir;
TString szDataPattern;
TString szDataDir;
TString szRunPrefix;
Int_t nNoOfTestFiles;
Int_t nNoOfInputFiles;
TObjArray listOfActiveRuns;

/* local file list */
TString szLocalFileList;

void runCorrelationsStudiesConfigMacro();

#endif /* RUNCORRELATIONSTUDIESCONFIG_H_ */
