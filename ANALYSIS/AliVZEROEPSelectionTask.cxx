/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
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

//*****************************************************
//   Class AliVZEROEPSelectionTask
//   author: Cvetan Cheshkov
//   30/01/2012
//   This analysis task reads the OADB and
//   provides the parameters needed to flatten
//   the VZERO event plane in AliEventplane
//*****************************************************

#include "AliVZEROEPSelectionTask.h"

#include <TList.h>
#include <TProfile.h>
#include <TFile.h>
#include <TString.h>
#include <TDirectory.h>

#include "AliLog.h"
#include "AliVEvent.h"
#include "AliAnalysisManager.h"
#include "AliOADBContainer.h"
#include "AliEventplane.h"
#include "AliCentrality.h"

ClassImp(AliVZEROEPSelectionTask)

//________________________________________________________________________
AliVZEROEPSelectionTask::AliVZEROEPSelectionTask():
AliAnalysisTaskSE(),
  fRunNumber(-1),
  fUserParams(kFALSE),
  fUseVZEROCentrality(kFALSE),
  fVZEROEPContainer(0)
{   
  // Default constructor
  // Initialize pointers
  AliInfo("VZERO Event Plane Selection enabled.");
  for(Int_t i = 0; i < 11; ++i) fX2In[i] = fY2In[i] = fX2Y2In[i] = fCos8PsiIn[i] = NULL;
}   

//________________________________________________________________________
AliVZEROEPSelectionTask::AliVZEROEPSelectionTask(const char *name):
  AliAnalysisTaskSE(name),
  fRunNumber(-1),
  fUserParams(kFALSE),
  fUseVZEROCentrality(kFALSE),
  fVZEROEPContainer(0)
{
  // Default constructor
  // Initialize pointers
  AliInfo("Event Plane Selection enabled.");
  for(Int_t i = 0; i < 11; ++i) fX2In[i] = fY2In[i] = fX2Y2In[i] = fCos8PsiIn[i] = NULL;
}
 
//________________________________________________________________________
AliVZEROEPSelectionTask::~AliVZEROEPSelectionTask()
{
  // Destructor
  // ...
  if (fUserParams) {
    for(Int_t i = 0; i < 11; ++i) {
      delete fX2In[i];
      fX2In[i] = NULL;
      delete fY2In[i];
      fY2In[i] = NULL;
      delete fX2Y2In[i];
      fX2Y2In[i] = NULL;
      delete fCos8PsiIn[i];
      fCos8PsiIn[i] = NULL;
    }
  }
  if (fVZEROEPContainer){
    delete fVZEROEPContainer;
    fVZEROEPContainer = NULL;
  }
}  

//________________________________________________________________________
void AliVZEROEPSelectionTask::UserCreateOutputObjects()
{  
  // Create the output containers (none in this case)
  // Open the OADB file
  
  if(!fUserParams) {
    TString oadbFileName = Form("%s/COMMON/EVENTPLANE/data/vzero.root", AliAnalysisManager::GetOADBPath());
    TFile *fOADB = TFile::Open(oadbFileName); 
    if(!fOADB->IsOpen()) AliFatal(Form("Cannot open OADB file %s", oadbFileName.Data()));

    AliInfo("Using Standard OADB");
    AliOADBContainer *cont = (AliOADBContainer*)fOADB->Get("vzeroEP");
    if (!cont) AliFatal("Cannot fetch OADB container for VZERO EP selection");
    fVZEROEPContainer = new AliOADBContainer(*cont);
    fOADB->Close();
    delete fOADB;
  }
}

//________________________________________________________________________
void AliVZEROEPSelectionTask::UserExec(Option_t */*option*/)
{ 
  // Execute analysis for current event:
  // Fill the flatenning parameters in
  // AliEventplane object

  AliVEvent* event = InputEvent();
  if (!(fRunNumber == event->GetRunNumber())) {
    fRunNumber = event->GetRunNumber();
    SetParamsFromOADB();
  }

  AliCentrality *centrality = event->GetCentrality();
  Float_t percentile = (fUseVZEROCentrality) ? centrality->GetCentralityPercentile("V0M") : centrality->GetCentralityPercentile("CL1");
  AliEventplane *esdEP = event->GetEventplane();
  if(esdEP) SetEventplaneParams(esdEP,percentile);
}

//________________________________________________________________________
void AliVZEROEPSelectionTask::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  // Nothing here
}

//________________________________________________________________________
void AliVZEROEPSelectionTask::SetEventplaneParams(AliEventplane *esdEP,Float_t percentile)
{
  // Read the OADB histograms and
  // prepare parameters used in order to
  // flatten the event-plane
  if(!esdEP)
    AliFatal("No event plane received");

  if (percentile < 0 || percentile > 100) {
    for(Int_t ring = 0; ring < 11; ++ring) esdEP->SetVZEROEPParams(ring,0.,0.,1.,1.,0.,0.,0.);
    return;
  }

  for(Int_t ring = 0; ring < 11; ++ring) {
    Int_t ibin = fX2In[ring]->FindBin(percentile);
    if (fX2In[ring]->GetBinEntries(ibin) == 0) {
      esdEP->SetVZEROEPParams(ring,0.,0.,1.,1.,0.,0.,0.);
      continue;
    }
    Double_t meanX2 = fX2In[ring]->GetBinContent(ibin);
    Double_t meanY2 = fY2In[ring]->GetBinContent(ibin);
    Double_t sigmaX2 = fX2In[ring]->GetBinError(ibin);
    Double_t sigmaY2 = fY2In[ring]->GetBinError(ibin);
    Double_t rho = (fX2Y2In[ring]->GetBinContent(ibin)-meanX2*meanY2)/sigmaX2/sigmaY2;
  
    Double_t b = rho*sigmaX2*sigmaY2*
      TMath::Sqrt(2.*(sigmaX2*sigmaX2+sigmaY2*sigmaY2-2.*sigmaX2*sigmaY2*TMath::Sqrt(1.-rho*rho))/
		  ((sigmaX2*sigmaX2-sigmaY2*sigmaY2)*(sigmaX2*sigmaX2-sigmaY2*sigmaY2)+
		   4.*sigmaX2*sigmaX2*sigmaY2*sigmaY2*rho*rho));
    Double_t aPlus = TMath::Sqrt(2.*sigmaX2*sigmaX2-b*b);
    Double_t aMinus= TMath::Sqrt(2.*sigmaY2*sigmaY2-b*b);

    Double_t lambdaPlus = b/aPlus;
    Double_t lambdaMinus = b/aMinus;

    Double_t cos8Psi = fCos8PsiIn[ring]->GetBinContent(ibin);
    esdEP->SetVZEROEPParams(ring,meanX2,meanY2,aPlus,aMinus,lambdaPlus,lambdaMinus,cos8Psi);
  }
}

//__________________________________________________________________________
void AliVZEROEPSelectionTask::SetParamsFromOADB() 
{
  if(!fUserParams) {
    TList *list = (TList*)fVZEROEPContainer->GetObject(fRunNumber, "Default");
    if (!list) AliFatal(Form("Cannot find VZERO OADB list for run %d", fRunNumber));
    SetHistograms(list);
  }
  else
    AliInfo("Using custom VZERO event-plane params");
}

//__________________________________________________________________________
void AliVZEROEPSelectionTask::SetUserParams(const char* inFileName, const char* listName)
{
  
  fUserParams = kTRUE;
  
  TFile f(inFileName);
  TList* list = (TList*)f.Get(listName);
  if (!list) AliFatal(Form("Cannot find list %s in file %s", listName, inFileName));
  SetHistograms(list);
  f.Close();
} 

//__________________________________________________________________________
void AliVZEROEPSelectionTask::SetHistograms(TList *list)
{
  // Set the flatenning parameters
  // histograms from a given list

  for(Int_t i = 0; i < 11; ++i) {
    if (fX2In[i]) delete fX2In[i];
    fX2In[i] = (TProfile*)list->FindObject(Form("fX2_%d",i))->Clone(Form("fX2In_%d",i));
    fX2In[i]->SetDirectory(0);
    if (fY2In[i]) delete fY2In[i];
    fY2In[i] = (TProfile*)list->FindObject(Form("fY2_%d",i))->Clone(Form("fY2In_%d",i));
    fY2In[i]->SetDirectory(0);
    if (fX2Y2In[i]) delete fX2Y2In[i];
    fX2Y2In[i] = (TProfile*)list->FindObject(Form("fX2Y2_%d",i))->Clone(Form("fX2Y2In_%d",i));
    fX2Y2In[i]->SetDirectory(0);
    if (fCos8PsiIn[i]) delete fCos8PsiIn[i];
    fCos8PsiIn[i] = (TProfile*)list->FindObject(Form("fCos8Psi_%d",i))->Clone(Form("fCos8PsiIn_%d",i));
    fCos8PsiIn[i]->SetDirectory(0);
  }
}
