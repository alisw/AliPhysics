#if (!defined(__CINT__) && !defined(__CLING__)) || defined __MAKECINT__
#include <TMacro.h>
#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>

#include "AliAnalysisManager.h"
#include "AliESDtrackCuts.h"
#include "AliPWG4HighPtTrackQA.h"
#include "AliVEvent.h"
#endif

void AddTaskPWG4HighPtTrackQApPb(const char *prodType = "LHC13b");
void AddTaskPWG4HighPtTrackQAAll(const char *prodType = "LHC10h", Bool_t isPbPb = kTRUE, Int_t iAODanalysis = 0);
void AddTaskPWG4HighPtTrackQAAll2011(const char *prodType = "LHC11h", Bool_t isPbPb = kTRUE, Int_t iAODanalysis = 0);
void AddTaskPWG4HighPtTrackQAAllReduced(const char *prodType = "LHC11h", Bool_t isPbPb = kTRUE, Int_t iAODanalysis = 0);
void AddTaskPWG4HighPtTrackQALHC11hLTS(const char *prodType = "LHC10h", Bool_t isPbPb = kTRUE, Int_t iAODanalysis = 0);
void AddTaskPWG4HighPtTrackQAAllReduced2011(const char *prodType = "LHC10h", Bool_t isPbPb = kTRUE, Int_t iAODanalysis = 0);
void AddTaskPWG4HighPtTrackQAppRun2(const char *prodtype = "LHC16h", Bool_t usetriggers = kFALSE);
void AddTaskPWG4HighPtTrackQAAOD(const char *prodType = "LHC11h", Bool_t isPbPb = kTRUE, Int_t iAODanalysis = 1, Int_t filterBit = 768, Bool_t useSPDTrackletVsClusterBG = kFALSE, Bool_t usetriggers = kFALSE);
AliPWG4HighPtTrackQA *ConfigureTaskPWG4HighPtTrackQA(const char *prodType = "LHC10e14", Bool_t isPbPb = kTRUE, Int_t iAODanalysis = 0, Int_t centClass = 0, Int_t trackType = 0, Int_t cuts = 0, UInt_t iPhysicsSelectionFlag = AliVEvent::kMB, TString emcaltrigger= "");

void AddTaskPWG4HighPtTrackQA(TString prodType = "LHC10h", Int_t iAODanalysis = 0, Bool_t isPbPb = kTRUE, Bool_t bReduced = kTRUE, Int_t filterBit = 272, Bool_t doEfficiency = kFALSE, Bool_t useSPDTrackletVsClusterBG = kFALSE, Bool_t usetriggers = kFALSE)
{
  if (iAODanalysis == 0)
  { //run on ESDs
    if (prodType.EqualTo("LHC10h") || prodType.EqualTo("LHC11a"))
    {
      if (bReduced)
        AddTaskPWG4HighPtTrackQAAllReduced(prodType.Data(), isPbPb, iAODanalysis);
      else
        AddTaskPWG4HighPtTrackQAAll(prodType.Data(), isPbPb, iAODanalysis);
    }
    else if (prodType.Contains("LHC12") || prodType.Contains("LHC13"))
    {
      AddTaskPWG4HighPtTrackQApPb();
    }
    else if (prodType.Contains("LHC16") || prodType.Contains("LHC17"))
    {
      AddTaskPWG4HighPtTrackQAppRun2(prodType, usetriggers); 
    }
    else
    {
      if (bReduced)
        AddTaskPWG4HighPtTrackQAAllReduced2011(prodType.Data(), isPbPb, iAODanalysis);
      else
        AddTaskPWG4HighPtTrackQAAll2011(prodType.Data(), isPbPb, iAODanalysis);
    }
  }
  else if (iAODanalysis == 1)
  { //run on AODs
    if (doEfficiency)
    {
      TMacro addmacro(gSystem->ExpandPathName("$ALICE_PHYSICS/PWGJE/macros/AddTaskHybridTrackEfficiency.C"));
      addmacro.Exec(Form("\"%s\", %s, AliVEvent::kMB, kTRUE, kFALSE", prodType.Data(), isPbPb ? "kTRUE" : "kFALSE"));
    }
    AddTaskPWG4HighPtTrackQAAOD(prodType.Data(), isPbPb, iAODanalysis, filterBit, useSPDTrackletVsClusterBG, usetriggers);
  }
}

void AddTaskPWG4HighPtTrackQApPb(const char *prodType)
{
  AliPWG4HighPtTrackQA *taskTrackQA05cent10 = ConfigureTaskPWG4HighPtTrackQA(prodType, kFALSE, 0, 10, 0, 5, AliVEvent::kINT7);
  AliPWG4HighPtTrackQA *taskTrackQA75cent10 = ConfigureTaskPWG4HighPtTrackQA(prodType, kFALSE, 0, 10, 7, 5, AliVEvent::kINT7);

  if (!strcmp(prodType, "LHC13d") || !strcmp(prodType, "LHC13e") || !strcmp(prodType, "LHC13f"))
  {
    AliPWG4HighPtTrackQA *taskTrackQA05cent10 = ConfigureTaskPWG4HighPtTrackQA(prodType, kFALSE, 0, 10, 0, 5, AliVEvent::kEMCEJE);
    AliPWG4HighPtTrackQA *taskTrackQA75cent10 = ConfigureTaskPWG4HighPtTrackQA(prodType, kFALSE, 0, 10, 7, 5, AliVEvent::kEMCEJE);
  }
}

void AddTaskPWG4HighPtTrackQAppRun2(const char *prodType, Bool_t usetriggers) {
  AliPWG4HighPtTrackQA *taskTrackQA05cent10 = ConfigureTaskPWG4HighPtTrackQA(prodType, kFALSE, 0, 10, 0, 5, AliVEvent::kINT7);
  AliPWG4HighPtTrackQA *taskTrackQA75cent10 = ConfigureTaskPWG4HighPtTrackQA(prodType, kFALSE, 0, 10, 7, 5, AliVEvent::kINT7);

  if (usetriggers)
  {
    TString trgnames[4] = {"EG1", "EG2", "EJ1", "EJ2"};
    UInt_t triggerbits[4] = {AliVEvent::kEMCEGA, AliVEvent::kEMCEGA, AliVEvent::kEMCEJE, AliVEvent::kEMCEJE};
    for(int itrg = 0; itrg < 4; itrg++){
      AliPWG4HighPtTrackQA *taskTrackQA05cent10trg = ConfigureTaskPWG4HighPtTrackQA(prodType, kFALSE, 0, 10, 0, 5, triggerbits[itrg], trgnames[itrg]);
      AliPWG4HighPtTrackQA *taskTrackQA75cent10trg = ConfigureTaskPWG4HighPtTrackQA(prodType, kFALSE, 0, 10, 7, 5, triggerbits[itrg], trgnames[itrg]);
    }
  }
}

void AddTaskPWG4HighPtTrackQAAll(const char *prodType, Bool_t isPbPb, Int_t iAODanalysis)
{

  Int_t cent = 10;

  AliPWG4HighPtTrackQA *taskTrackQA00cent10 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 0, 0); //RAA track cuts
  AliPWG4HighPtTrackQA *taskTrackQA01cent10 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 0, 1); //global hybrid unconstrained
  AliPWG4HighPtTrackQA *taskTrackQA70cent10 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 7, 0); //global hybrid constrained category 1
  AliPWG4HighPtTrackQA *taskTrackQA71cent10 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 7, 1); //global hybrid constrained category 2
  AliPWG4HighPtTrackQA *taskTrackQA72cent10 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 7, 2); //global hybrid constrained all

  if (isPbPb)
  {
    for (cent = 0; cent < 4; cent++)
    {
      AliPWG4HighPtTrackQA *taskTrackQA00 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 0, 0);
      AliPWG4HighPtTrackQA *taskTrackQA01 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 0, 1);
      AliPWG4HighPtTrackQA *taskTrackQA70 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 7, 0);
      AliPWG4HighPtTrackQA *taskTrackQA71 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 7, 1);
      AliPWG4HighPtTrackQA *taskTrackQA72 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 7, 2);
    }
  }
}

void AddTaskPWG4HighPtTrackQAAll2011(const char *prodType, Bool_t isPbPb, Int_t iAODanalysis)
{

  Int_t cent = 10;
  UInt_t iPhysicsSelectionFlag = AliVEvent::kMB;
  UInt_t iPhysicsSelectionFlagCentral = AliVEvent::kCentral;
  UInt_t iPhysicsSelectionFlagSemiCentral = AliVEvent::kSemiCentral;

  AliPWG4HighPtTrackQA *taskTrackQA00cent10 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 0, 0, iPhysicsSelectionFlag);
  AliPWG4HighPtTrackQA *taskTrackQA01cent10 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 0, 1, iPhysicsSelectionFlag);
  AliPWG4HighPtTrackQA *taskTrackQA70cent10 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 7, 0, iPhysicsSelectionFlag);
  AliPWG4HighPtTrackQA *taskTrackQA71cent10 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 7, 1, iPhysicsSelectionFlag);
  AliPWG4HighPtTrackQA *taskTrackQA72cent10 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 7, 2, iPhysicsSelectionFlag);

  AliPWG4HighPtTrackQA *taskTrackQA05cent10 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 0, 5, iPhysicsSelectionFlag);
  AliPWG4HighPtTrackQA *taskTrackQA74cent10 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 7, 4, iPhysicsSelectionFlag);
  AliPWG4HighPtTrackQA *taskTrackQA75cent10 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 7, 5, iPhysicsSelectionFlag);

  if (isPbPb)
  {
    for (cent = 0; cent < 4; cent++)
    {
      AliPWG4HighPtTrackQA *taskTrackQA00 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 0, 0, iPhysicsSelectionFlag);
      AliPWG4HighPtTrackQA *taskTrackQA01 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 0, 1, iPhysicsSelectionFlag);
      AliPWG4HighPtTrackQA *taskTrackQA70 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 7, 0, iPhysicsSelectionFlag);
      AliPWG4HighPtTrackQA *taskTrackQA71 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 7, 1, iPhysicsSelectionFlag);
      AliPWG4HighPtTrackQA *taskTrackQA72 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 7, 2, iPhysicsSelectionFlag);
      AliPWG4HighPtTrackQA *taskTrackQA05 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 0, 5, iPhysicsSelectionFlag);
      AliPWG4HighPtTrackQA *taskTrackQA74 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 7, 4, iPhysicsSelectionFlag);
      AliPWG4HighPtTrackQA *taskTrackQA75 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 7, 5, iPhysicsSelectionFlag);

      if (cent == 0)
      {
        AliPWG4HighPtTrackQA *taskTrackQA00 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 0, 0, iPhysicsSelectionFlagCentral);
        AliPWG4HighPtTrackQA *taskTrackQA01 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 0, 1, iPhysicsSelectionFlagCentral);
        AliPWG4HighPtTrackQA *taskTrackQA70 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 7, 0, iPhysicsSelectionFlagCentral);
        AliPWG4HighPtTrackQA *taskTrackQA71 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 7, 1, iPhysicsSelectionFlagCentral);
        AliPWG4HighPtTrackQA *taskTrackQA72 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 7, 2, iPhysicsSelectionFlagCentral);
        AliPWG4HighPtTrackQA *taskTrackQA05 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 0, 5, iPhysicsSelectionFlagCentral);
        AliPWG4HighPtTrackQA *taskTrackQA74 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 7, 4, iPhysicsSelectionFlagCentral);
        AliPWG4HighPtTrackQA *taskTrackQA75 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 7, 5, iPhysicsSelectionFlagCentral);
      }
      else
      {
        AliPWG4HighPtTrackQA *taskTrackQA00 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 0, 0, iPhysicsSelectionFlagSemiCentral);
        AliPWG4HighPtTrackQA *taskTrackQA01 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 0, 1, iPhysicsSelectionFlagSemiCentral);
        AliPWG4HighPtTrackQA *taskTrackQA70 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 7, 0, iPhysicsSelectionFlagSemiCentral);
        AliPWG4HighPtTrackQA *taskTrackQA71 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 7, 1, iPhysicsSelectionFlagSemiCentral);
        AliPWG4HighPtTrackQA *taskTrackQA72 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 7, 2, iPhysicsSelectionFlagSemiCentral);
        AliPWG4HighPtTrackQA *taskTrackQA05 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 0, 5, iPhysicsSelectionFlagSemiCentral);
        AliPWG4HighPtTrackQA *taskTrackQA74 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 7, 4, iPhysicsSelectionFlagSemiCentral);
        AliPWG4HighPtTrackQA *taskTrackQA75 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 7, 5, iPhysicsSelectionFlagSemiCentral);
      }
    }
  }
}

void AddTaskPWG4HighPtTrackQAAllReduced(const char *prodType, Bool_t isPbPb, Int_t iAODanalysis)
{

  int cent = 10;

  if (isPbPb)
  {
    for (cent = 0; cent < 4; cent++)
    {
      AliPWG4HighPtTrackQA *taskTrackQA00 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 0, 0);
      AliPWG4HighPtTrackQA *taskTrackQA01 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 0, 1);
      AliPWG4HighPtTrackQA *taskTrackQA70 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 7, 0);
      AliPWG4HighPtTrackQA *taskTrackQA71 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 7, 1);
      AliPWG4HighPtTrackQA *taskTrackQA72 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 7, 2);
    }
  }
}

void AddTaskPWG4HighPtTrackQALHC11hLTS(const char *prodType, Bool_t isPbPb, Int_t iAODanalysis)
{

  Int_t cent = 10;
  UInt_t iPhysicsSelectionFlag = AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral;
  UInt_t iPhysicsSelectionFlagEMCEJE = AliVEvent::kEMCEJE;

  AliPWG4HighPtTrackQA *taskTrackQA00cent10 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 0, 0, iPhysicsSelectionFlag);
  AliPWG4HighPtTrackQA *taskTrackQA01cent10 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 0, 1, iPhysicsSelectionFlag);
  AliPWG4HighPtTrackQA *taskTrackQA70cent10 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 7, 0, iPhysicsSelectionFlag);
  AliPWG4HighPtTrackQA *taskTrackQA71cent10 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 7, 1, iPhysicsSelectionFlag);
  AliPWG4HighPtTrackQA *taskTrackQA72cent10 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 7, 2, iPhysicsSelectionFlag);
  AliPWG4HighPtTrackQA *taskTrackQA05cent10 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 0, 5, iPhysicsSelectionFlag);
  AliPWG4HighPtTrackQA *taskTrackQA74cent10 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 7, 4, iPhysicsSelectionFlag);
  AliPWG4HighPtTrackQA *taskTrackQA75cent10 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 7, 5, iPhysicsSelectionFlag);

  AliPWG4HighPtTrackQA *taskTrackQAEMCJE00cent10 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 0, 0, iPhysicsSelectionFlagEMCEJE);
  AliPWG4HighPtTrackQA *taskTrackQAEMCJE01cent10 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 0, 1, iPhysicsSelectionFlagEMCEJE);
  AliPWG4HighPtTrackQA *taskTrackQAEMCJE70cent10 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 7, 0, iPhysicsSelectionFlagEMCEJE);
  AliPWG4HighPtTrackQA *taskTrackQAEMCJE71cent10 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 7, 1, iPhysicsSelectionFlagEMCEJE);
  AliPWG4HighPtTrackQA *taskTrackQAEMCJE72cent10 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 7, 2, iPhysicsSelectionFlagEMCEJE);
  AliPWG4HighPtTrackQA *taskTrackQAEMCJE05cent10 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 0, 5, iPhysicsSelectionFlagEMCEJE);
  AliPWG4HighPtTrackQA *taskTrackQAEMCJE74cent10 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 7, 4, iPhysicsSelectionFlagEMCEJE);
  AliPWG4HighPtTrackQA *taskTrackQAEMCJE75cent10 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 7, 5, iPhysicsSelectionFlagEMCEJE);

  if (isPbPb)
  {
    for (cent = 0; cent < 4; cent++)
    {
      AliPWG4HighPtTrackQA *taskTrackQA00 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 0, 0, iPhysicsSelectionFlag);
      AliPWG4HighPtTrackQA *taskTrackQA01 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 0, 1, iPhysicsSelectionFlag);
      AliPWG4HighPtTrackQA *taskTrackQA70 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 7, 0, iPhysicsSelectionFlag);
      AliPWG4HighPtTrackQA *taskTrackQA71 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 7, 1, iPhysicsSelectionFlag);
      AliPWG4HighPtTrackQA *taskTrackQA72 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 7, 2, iPhysicsSelectionFlag);
      AliPWG4HighPtTrackQA *taskTrackQA05 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 0, 5, iPhysicsSelectionFlag);
      AliPWG4HighPtTrackQA *taskTrackQA74 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 7, 4, iPhysicsSelectionFlag);
      AliPWG4HighPtTrackQA *taskTrackQA75 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 7, 5, iPhysicsSelectionFlag);

      AliPWG4HighPtTrackQA *taskTrackQAEMCJE00 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 0, 0, iPhysicsSelectionFlagEMCEJE);
      AliPWG4HighPtTrackQA *taskTrackQAEMCJE01 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 0, 1, iPhysicsSelectionFlagEMCEJE);
      AliPWG4HighPtTrackQA *taskTrackQAEMCJE70 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 7, 0, iPhysicsSelectionFlagEMCEJE);
      AliPWG4HighPtTrackQA *taskTrackQAEMCJE71 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 7, 1, iPhysicsSelectionFlagEMCEJE);
      AliPWG4HighPtTrackQA *taskTrackQAEMCJE72 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 7, 2, iPhysicsSelectionFlagEMCEJE);
      AliPWG4HighPtTrackQA *taskTrackQAEMCJE05 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 0, 5, iPhysicsSelectionFlagEMCEJE);
      AliPWG4HighPtTrackQA *taskTrackQAEMCJE74 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 7, 4, iPhysicsSelectionFlagEMCEJE);
      AliPWG4HighPtTrackQA *taskTrackQAEMCJE75 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 7, 5, iPhysicsSelectionFlagEMCEJE);
    }
  }
}

void AddTaskPWG4HighPtTrackQAAllReduced2011(const char *prodType, Bool_t isPbPb, Int_t iAODanalysis)
{

  int cent = 10;

  UInt_t iPhysicsSelectionFlagCentral = AliVEvent::kCentral;
  UInt_t iPhysicsSelectionFlagSemiCentral = AliVEvent::kSemiCentral;

  AliPWG4HighPtTrackQA *taskTrackQA00C = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 0, 0, iPhysicsSelectionFlagCentral);
  AliPWG4HighPtTrackQA *taskTrackQA01C = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 0, 1, iPhysicsSelectionFlagCentral);
  AliPWG4HighPtTrackQA *taskTrackQA70C = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 7, 0, iPhysicsSelectionFlagCentral);
  AliPWG4HighPtTrackQA *taskTrackQA71C = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 7, 1, iPhysicsSelectionFlagCentral);
  AliPWG4HighPtTrackQA *taskTrackQA72C = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 7, 2, iPhysicsSelectionFlagCentral);
  AliPWG4HighPtTrackQA *taskTrackQA05C = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 0, 5, iPhysicsSelectionFlagCentral);
  AliPWG4HighPtTrackQA *taskTrackQA74C = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 7, 4, iPhysicsSelectionFlagCentral);
  AliPWG4HighPtTrackQA *taskTrackQA75C = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 7, 5, iPhysicsSelectionFlagCentral);

  AliPWG4HighPtTrackQA *taskTrackQA00SC = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 0, 0, iPhysicsSelectionFlagSemiCentral);
  AliPWG4HighPtTrackQA *taskTrackQA01SC = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 0, 1, iPhysicsSelectionFlagSemiCentral);
  AliPWG4HighPtTrackQA *taskTrackQA70SC = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 7, 0, iPhysicsSelectionFlagSemiCentral);
  AliPWG4HighPtTrackQA *taskTrackQA71SC = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 7, 1, iPhysicsSelectionFlagSemiCentral);
  AliPWG4HighPtTrackQA *taskTrackQA72SC = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 7, 2, iPhysicsSelectionFlagSemiCentral);
  AliPWG4HighPtTrackQA *taskTrackQA05SC = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 0, 5, iPhysicsSelectionFlagSemiCentral);
  AliPWG4HighPtTrackQA *taskTrackQA74SC = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 7, 4, iPhysicsSelectionFlagSemiCentral);
  AliPWG4HighPtTrackQA *taskTrackQA75SC = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 7, 5, iPhysicsSelectionFlagSemiCentral);
}

void AddTaskPWG4HighPtTrackQAAOD(const char *prodType, Bool_t isPbPb, Int_t iAODanalysis, Int_t filterBit, Bool_t useSPDTrackletVsClusterBG, Bool_t usetriggers)
{
  UInt_t iPhysicsSelectionFlagMB = AliVEvent::kMB;
  UInt_t iPhysicsSelectionFlagCentral = AliVEvent::kCentral;
  UInt_t iPhysicsSelectionFlagSemiCentral = AliVEvent::kSemiCentral;
  UInt_t iPhysicsSelectionFlagINT7 = AliVEvent::kINT7;
  UInt_t iPhysicsSelectionFlagEMCEJE = AliVEvent::kEMCEJE;

  Int_t cent = 10;

  Int_t filterBit1 = 256; //standard global tracks
  Int_t filterBit2 = 512; //complementary tracks

  Bool_t bIncludeNoITS = kFALSE;

  TString strRunPeriod = TString(prodType);
  strRunPeriod.ToLower();

  if (strRunPeriod == "lhc10h" || strRunPeriod == "lhc11h" ||
      strRunPeriod == "lhc12a" || strRunPeriod == "lhc12b" || strRunPeriod == "lhc12c" || strRunPeriod == "lhc12d" ||
      strRunPeriod == "lhc12e" || strRunPeriod == "lhc12f" || strRunPeriod == "lhc12g" || strRunPeriod == "lhc12g" ||
      strRunPeriod == "lhc12h" || strRunPeriod == "lhc12i" ||
      strRunPeriod == "lhc13b" || strRunPeriod == "lhc13c" || strRunPeriod == "lhc13d" || strRunPeriod == "lhc13e" ||
      strRunPeriod == "lhc13f" || strRunPeriod == "lhc13g" ||
      strRunPeriod == "lhc12a15e" || strRunPeriod == "lhc13b4" || strRunPeriod == "lhc13b4_fix" ||
      strRunPeriod == "lhc13b4_plus" || strRunPeriod == "lhc12a15f" || strRunPeriod.Contains("lhc12a17") || strRunPeriod.Contains("lhc14a1") ||
      strRunPeriod.Contains("lhc16c2") || strRunPeriod.Contains("lhc16e1") || strRunPeriod.Contains("lhc17f8") ||
      (strRunPeriod.Length() == 6 && (strRunPeriod.BeginsWith("lhc15") || strRunPeriod.BeginsWith("lhc16") || strRunPeriod.BeginsWith("lhc17"))) // All run2 data
  )
  {
    filterBit = 768;
    filterBit1 = 256;
    filterBit2 = 512;
    bIncludeNoITS = kFALSE;
    if (strRunPeriod == "lhc10h")
      bIncludeNoITS = kTRUE;
  }
  else if (strRunPeriod == "lhc11a" || strRunPeriod == "lhc10hold" || strRunPeriod == "lhc12a15a" || strRunPeriod.Contains("lhc11a2"))
  {
    filterBit = 272;
    filterBit1 = 16;
    filterBit2 = 256;
    bIncludeNoITS = kTRUE;
  }

  if (isPbPb)
  {

    if (strRunPeriod.Contains("lhc15o"))
    {
      cent = 10;
      AliPWG4HighPtTrackQA *taskTrackQAMB = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 0, 0, iPhysicsSelectionFlagINT7);
      taskTrackQAMB->SetFilterMask(filterBit);
      taskTrackQAMB->SetIncludeNoITS(bIncludeNoITS);

      AliPWG4HighPtTrackQA *taskTrackQAMB1 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 0, 5, iPhysicsSelectionFlagINT7);
      taskTrackQAMB1->SetFilterMask(filterBit1);
      taskTrackQAMB1->SetIncludeNoITS(bIncludeNoITS);

      AliPWG4HighPtTrackQA *taskTrackQAMB2 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 7, 5, iPhysicsSelectionFlagINT7);
      taskTrackQAMB2->SetFilterMask(filterBit2);
      taskTrackQAMB2->SetIncludeNoITS(bIncludeNoITS);
    }
    else
    {

      for (cent = 0; cent < 4; cent++)
      {
        AliPWG4HighPtTrackQA *taskTrackQAMB = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 0, 0, iPhysicsSelectionFlagMB);
        taskTrackQAMB->SetFilterMask(filterBit);
        taskTrackQAMB->SetIncludeNoITS(bIncludeNoITS);

        AliPWG4HighPtTrackQA *taskTrackQAMB1 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 0, 5, iPhysicsSelectionFlagMB);
        taskTrackQAMB1->SetFilterMask(filterBit1);
        taskTrackQAMB1->SetIncludeNoITS(bIncludeNoITS);

        AliPWG4HighPtTrackQA *taskTrackQAMB2 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 7, 5, iPhysicsSelectionFlagMB);
        taskTrackQAMB2->SetFilterMask(filterBit2);
        taskTrackQAMB2->SetIncludeNoITS(bIncludeNoITS);
      }

      cent = 10;

      if (strRunPeriod.Contains("lhc11h"))
      {
        AliPWG4HighPtTrackQA *taskTrackQAC = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 0, 0, iPhysicsSelectionFlagCentral);
        taskTrackQAC->SetFilterMask(filterBit);
        taskTrackQAC->SetIncludeNoITS(bIncludeNoITS);

        AliPWG4HighPtTrackQA *taskTrackQAC1 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 0, 5, iPhysicsSelectionFlagCentral);
        taskTrackQAC1->SetFilterMask(filterBit1);
        taskTrackQAC1->SetIncludeNoITS(bIncludeNoITS);

        AliPWG4HighPtTrackQA *taskTrackQAC2 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 7, 5, iPhysicsSelectionFlagCentral);
        taskTrackQAC2->SetFilterMask(filterBit2);
        taskTrackQAC2->SetIncludeNoITS(bIncludeNoITS);

        AliPWG4HighPtTrackQA *taskTrackQASC = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 0, 0, iPhysicsSelectionFlagSemiCentral);
        taskTrackQASC->SetFilterMask(filterBit);
        taskTrackQASC->SetIncludeNoITS(bIncludeNoITS);

        AliPWG4HighPtTrackQA *taskTrackQASC1 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 0, 5, iPhysicsSelectionFlagSemiCentral);
        taskTrackQASC1->SetFilterMask(filterBit1);
        taskTrackQASC1->SetIncludeNoITS(bIncludeNoITS);

        AliPWG4HighPtTrackQA *taskTrackQASC2 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 7, 5, iPhysicsSelectionFlagSemiCentral);
        taskTrackQASC2->SetFilterMask(filterBit2);
        taskTrackQASC2->SetIncludeNoITS(bIncludeNoITS);
      }
    }
  }
  else
  {
    cent = 10;

    if (strRunPeriod.Contains("lhc13") || strRunPeriod.Contains("lhc12"))
    {
      AliPWG4HighPtTrackQA *taskTrackQAMB = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 0, 0, iPhysicsSelectionFlagINT7);
      taskTrackQAMB->SetFilterMask(filterBit);
      taskTrackQAMB->SetIncludeNoITS(bIncludeNoITS);

      AliPWG4HighPtTrackQA *taskTrackQAMB1 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 0, 5, iPhysicsSelectionFlagINT7);
      taskTrackQAMB1->SetFilterMask(filterBit1);
      taskTrackQAMB1->SetIncludeNoITS(bIncludeNoITS);

      AliPWG4HighPtTrackQA *taskTrackQAMB2 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 7, 5, iPhysicsSelectionFlagINT7);
      taskTrackQAMB2->SetFilterMask(filterBit2);
      taskTrackQAMB2->SetIncludeNoITS(bIncludeNoITS);

      if (strRunPeriod.EqualTo("lhc13d") || strRunPeriod.EqualTo("lhc13e") || strRunPeriod.EqualTo("lhc13f") || strRunPeriod.EqualTo("lhc13g") || strRunPeriod.Contains("lhc12"))
      {
        AliPWG4HighPtTrackQA *taskTrackQAEMCEJE = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 0, 0, iPhysicsSelectionFlagEMCEJE);
        taskTrackQAEMCEJE->SetFilterMask(filterBit);
        taskTrackQAEMCEJE->SetIncludeNoITS(bIncludeNoITS);

        AliPWG4HighPtTrackQA *taskTrackQAEMCEJE1 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 0, 5, iPhysicsSelectionFlagEMCEJE);
        taskTrackQAEMCEJE1->SetFilterMask(filterBit1);
        taskTrackQAEMCEJE1->SetIncludeNoITS(bIncludeNoITS);

        AliPWG4HighPtTrackQA *taskTrackQAEMCEJE2 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 7, 5, iPhysicsSelectionFlagEMCEJE);
        taskTrackQAEMCEJE2->SetFilterMask(filterBit2);
        taskTrackQAEMCEJE2->SetIncludeNoITS(bIncludeNoITS);
      }
    }
    else if(strRunPeriod.Contains("lhc16") || strRunPeriod.Contains("lhc17"))
    {
      AliPWG4HighPtTrackQA *taskTrackQAMB = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 0, 0, iPhysicsSelectionFlagINT7);
      taskTrackQAMB->SetFilterMask(filterBit);
      taskTrackQAMB->SetIncludeNoITS(bIncludeNoITS);

      AliPWG4HighPtTrackQA *taskTrackQAMB1 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 0, 5, iPhysicsSelectionFlagINT7);
      taskTrackQAMB1->SetFilterMask(filterBit1);
      taskTrackQAMB1->SetIncludeNoITS(bIncludeNoITS);

      AliPWG4HighPtTrackQA *taskTrackQAMB2 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 7, 5, iPhysicsSelectionFlagINT7);
      taskTrackQAMB2->SetFilterMask(filterBit2);
      taskTrackQAMB2->SetIncludeNoITS(bIncludeNoITS);

      if(usetriggers) {
        TString trgnames[4] = {"EG1", "EG2", "EJ1", "EJ2"};
        UInt_t triggerbits[4] = {AliVEvent::kEMCEGA, AliVEvent::kEMCEGA, AliVEvent::kEMCEJE, AliVEvent::kEMCEJE};
        for(int itrg = 0; itrg < 4; itrg++){
          AliPWG4HighPtTrackQA *taskTrackQAtrg = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 0, 0, triggerbits[itrg], trgnames[itrg]);
          taskTrackQAtrg->SetFilterMask(filterBit);
          taskTrackQAtrg->SetIncludeNoITS(bIncludeNoITS);
          taskTrackQAtrg->SetUseSPDTrackletVsClusterBG(useSPDTrackletVsClusterBG);

          AliPWG4HighPtTrackQA *taskTrackQAtrg1 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 0, 5, triggerbits[itrg], trgnames[itrg]);
          taskTrackQAtrg1->SetFilterMask(filterBit1);
          taskTrackQAtrg1->SetIncludeNoITS(bIncludeNoITS);
          taskTrackQAtrg1->SetUseSPDTrackletVsClusterBG(useSPDTrackletVsClusterBG);

          AliPWG4HighPtTrackQA *taskTrackQAtrg2 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 7, 5, triggerbits[itrg], trgnames[itrg]);
          taskTrackQAtrg2->SetFilterMask(filterBit2);
          taskTrackQAtrg2->SetIncludeNoITS(bIncludeNoITS);
          taskTrackQAtrg2->SetUseSPDTrackletVsClusterBG(useSPDTrackletVsClusterBG);
        }
      }
    }
    else
    {
      AliPWG4HighPtTrackQA *taskTrackQAMB = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 0, 0, iPhysicsSelectionFlagMB);
      taskTrackQAMB->SetFilterMask(filterBit);
      taskTrackQAMB->SetIncludeNoITS(bIncludeNoITS);
      taskTrackQAMB->SetUseSPDTrackletVsClusterBG(useSPDTrackletVsClusterBG);

      AliPWG4HighPtTrackQA *taskTrackQAMB1 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 0, 5, iPhysicsSelectionFlagMB);
      taskTrackQAMB1->SetFilterMask(filterBit1);
      taskTrackQAMB1->SetIncludeNoITS(bIncludeNoITS);
      taskTrackQAMB1->SetUseSPDTrackletVsClusterBG(useSPDTrackletVsClusterBG);

      AliPWG4HighPtTrackQA *taskTrackQAMB2 = ConfigureTaskPWG4HighPtTrackQA(prodType, isPbPb, iAODanalysis, cent, 7, 5, iPhysicsSelectionFlagMB);
      taskTrackQAMB2->SetFilterMask(filterBit2);
      taskTrackQAMB2->SetIncludeNoITS(bIncludeNoITS);
      taskTrackQAMB2->SetUseSPDTrackletVsClusterBG(useSPDTrackletVsClusterBG);
    }
  }
}

/**
 * @brief Configure the QA task
 *
 * # For ESDs: (AOD see bottom)
 * ------------------------------
 * | trackType: | 0 = global                                                                             |
 * |            | 1 = TPC stand alone                                                                    |
 * |            | 2 = TPC stand alone constrained to SPD vertex                                          |
 * |            | 4 = TPC stand alone constrained to SPD vertex with QA track selection on global tracks |
 * |            | 5 = (Obsolete) TPConly with vtx constraint (cuts on global track)                      |
 * |            | 6 = (Obsolete) Global track with vtx constraint                                        |
 * |            | 7 = Global tracks with vtx constraint (Hybrid tracks, complementary selection)         |
 * | cuts:      | 0 (global) = standard ITSTPC2010 a la RAA analysis                                     |
 * |            | 1 (global) = ITSrefit, no SPD requirements -> standard for jet analysis                |
 * |            | 2 (global) = ITSrefit + no hits in SPD                                                 |
 * |            | 3 (global) = standard ITS tight cuts with nCrossed rows cut for hybrid tracks          |
 * |            | 0 (TPC)    = standard TPC + NClusters>70                                               |
 * |            | 1 (TPC)    = standard TPC + NClusters>0 --> to study new TPC QA recommendations        |
 * |            | 0 (hybrid 5) = constrained TPConly for which no tight ITS is available                 |
 * |            | 0 (hybrid 6) = constrained loose global for which no tight ITS is available            |
 *
 * # Hybrid Tracks in ESDs:
 * --------------------------------
 * 2010 definition:
 * - Type 0 Cuts 1: w/ SPD req. & ITS refit
 * - Type 7 Cuts 1: w/o SPD req. & w/ ITS refit
 * - Type 7 Cuts 2: w/o SPD req. & w/o ITS refit
 * 2011-> definition:
 * - Type 0 Cuts 5: w/ SPD req. & ITS refit
 * - Type 7 Cuts 5: w/o SPD req. & w/ ITS refit
 * For more extensive documentation on the hybrid track definitions, see [Hybrid tracks](https://twiki.cern.ch/twiki/bin/view/ALICE/HybridTracks)
 *
 * For AODs:
 * -----------------
 * The tracktype and cuts numbers are there only for naming purposes. The 'real' work is done by setting the filtermask. The naming for the hybrid tracks is as follows
 * | ESD labelling                               | filterbit (old AODs)  |  filterbit (current)  |
 * -----------------------------------------------------------------------------------------------
 * | tracktype 0 cuts 5 = constrained tracks     |         16            |        256            |
 * | tracktype 7 cuts 5 = complementary tracks   |         256           |        512            |
 * | tracktype 0 cuts 0 = sum of these two       |         272           |        768            |
 *   Filterbit for complementary tracks in 2010 definition contains both the w/o SPD req. & w/o ITS refit. Filterbit are set in AddTaskPWG4HighPtTrackQAAOD based on production name (prodType)
 * 
 * @param prodType                  Name of the production
 * @param isPbPb                    If true the task runs on PbPb
 * @param iAODanalysis              If true the task runs on AODs
 * @param centClass                 Centrality class number
 * @param trackType                 Track type as specified above
 * @param cuts                      Cuts as specified above
 * @param iPhysicsSelectionFlag     Event selection flag as specified above
 * @param emcaltrigger              Name of the EMCAL trigger used in order to separate EMCAL L1 triggers (optional)
 */
AliPWG4HighPtTrackQA *ConfigureTaskPWG4HighPtTrackQA(const char *prodType, Bool_t isPbPb, Int_t iAODanalysis, Int_t centClass, Int_t trackType, Int_t cuts, UInt_t iPhysicsSelectionFlag, TString emcaltrigger)
{

  //Load common track cut class
  TMacro CreateTrackCutsPWGJE(gSystem->ExpandPathName("$ALICE_PHYSICS/PWGJE/macros/CreateTrackCutsPWGJE.C"));

  // Creates HighPtTrackQA analysis task and adds it to the analysis manager.

  // A. Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    Error("AddTaskPWG4HighPtQMC", "No analysis manager to connect to.");
    return NULL;
  }

  // B. Check the analysis type using the event handlers connected to the analysis
  //    manager. The availability of MC handler can also be checked here.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddPWG4TaskHighPtTrackQA", "This task requires an input event handler");
    return NULL;
  }

  // C. Create the task, add it to manager.
  //===========================================================================

  //CREATE THE  CUTS -----------------------------------------------
  //Use AliESDtrackCuts
  AliESDtrackCuts *trackCuts = new AliESDtrackCuts("AliESDtrackCuts", "Standard Cuts");
  AliESDtrackCuts *trackCutsReject = 0x0;
  AliESDtrackCuts *trackCutsTPConly = new AliESDtrackCuts("AliESDtrackCutsTPConly", "TPC only Cuts");

  //Standard Cuts
  //Set track cuts for global tracks
  if (trackType == 0 && cuts == 0)
  {
    // tight global tracks - RAA analysis
    trackCuts = (AliESDtrackCuts *)CreateTrackCutsPWGJE.Exec("1000");
  }
  if (trackType == 0 && cuts == 1)
  {
    //Cuts global tracks with ITSrefit requirement and SPDrequirement for jet analysis
    trackCuts = (AliESDtrackCuts *)CreateTrackCutsPWGJE.Exec("10001006");
  }
  if (trackType == 0 && cuts == 5)
  {
    //Cuts global tracks with ITSrefit requirement and SPDrequirement for jet analysis + NCrossedRowsCut>70 recommended in 2011
    trackCuts = (AliESDtrackCuts *)CreateTrackCutsPWGJE.Exec("10001008");
  }
  if (trackType == 0 && cuts == 6)
  {
    //Cuts global tracks with ITSrefit requirement and no SPDrequirement for jet analysis + NCrossedRowsCut>70 recommended in 2011
    trackCuts = (AliESDtrackCuts *)CreateTrackCutsPWGJE.Exec("10051008");
  }

  if (trackType == 0 && cuts == 2)
  {
    //Cuts global tracks with ITSrefit requirement but without SPD
    trackCuts = (AliESDtrackCuts *)CreateTrackCutsPWGJE.Exec("10011006");
  }
  if (trackType == 7 && cuts == 0)
  {
    // tight global tracks
    trackCuts = (AliESDtrackCuts *)CreateTrackCutsPWGJE.Exec("10041006");
    trackCutsReject = (AliESDtrackCuts *)CreateTrackCutsPWGJE.Exec("1006");
  }
  if (trackType == 7 && cuts == 4)
  {
    // tight global tracks +  NCrossedRowsCut>120 recommended in 2011
    trackCuts = (AliESDtrackCuts *)CreateTrackCutsPWGJE.Exec("10041008");
    trackCutsReject = (AliESDtrackCuts *)CreateTrackCutsPWGJE.Exec("1008");
  }
  if (trackType == 7 && cuts == 1)
  {
    // tight global tracks
    trackCuts = (AliESDtrackCuts *)CreateTrackCutsPWGJE.Exec("10011006");
  }
  if (trackType == 7 && cuts == 5)
  {
    // tight global tracks  + NCrossedRowsCut>70 recommended in 2011
    trackCuts = (AliESDtrackCuts *)CreateTrackCutsPWGJE.Exec("10011008");
  }
  if (trackType == 7 && cuts == 2)
  {
    // no requirements on SPD and ITSrefit failed
    trackCuts = (AliESDtrackCuts *)CreateTrackCutsPWGJE.Exec("10041006");       //no ITSrefit requirement filter 256
    trackCutsReject = (AliESDtrackCuts *)CreateTrackCutsPWGJE.Exec("10001006"); //ITSrefit requirement filter 16
  }
  if (trackType == 7 && cuts == 6)
  {
    // no requirements on SPD and ITSrefit failed
    trackCuts = (AliESDtrackCuts *)CreateTrackCutsPWGJE.Exec("10041008");       //no ITSrefit requirement filter 256
    trackCutsReject = (AliESDtrackCuts *)CreateTrackCutsPWGJE.Exec("10001008"); //ITSrefit requirement filter 16
  }

  if (trackType == 1 && cuts == 0)
  {
    //Set track cuts for TPConly tracks
    trackCuts = (AliESDtrackCuts *)CreateTrackCutsPWGJE.Exec("2001");
  }
  if (trackType == 1 && cuts == 1)
  {
    //Set track cuts for TPConly tracks
    trackCuts = (AliESDtrackCuts *)CreateTrackCutsPWGJE.Exec("10032001");
  }

  if (trackType == 2 && cuts == 0)
  {
    //Set track cuts for TPConly constrained tracks
    trackCuts = (AliESDtrackCuts *)CreateTrackCutsPWGJE.Exec("2001");
  }
  if (trackType == 2 && cuts == 1)
  {
    //Set track cuts for TPConly constrained tracks w/o cut on NClusters or NCrossedRows
    trackCuts = (AliESDtrackCuts *)CreateTrackCutsPWGJE.Exec("10032001");
  }

  if (trackType == 4 && cuts == 0)
  {
    //	      Set track cuts for TPConly constrained tracks
    trackCuts = (AliESDtrackCuts *)CreateTrackCutsPWGJE.Exec("2001");
  }
  if (trackType == 4 && cuts == 1)
  {
    //	      Set track cuts for TPConly constrained tracks
    trackCuts = (AliESDtrackCuts *)CreateTrackCutsPWGJE.Exec("10032001");
  }
  if (trackType == 5 || trackType == 6)
  {
    // tight global tracks
    trackCuts = (AliESDtrackCuts *)CreateTrackCutsPWGJE.Exec("1003");
    trackCutsReject = (AliESDtrackCuts *)CreateTrackCutsPWGJE.Exec("10021003");
    trackCutsTPConly = (AliESDtrackCuts *)CreateTrackCutsPWGJE.Exec("2002");
  }

  trackCuts->SetEtaRange(-0.9, 0.9);
  trackCuts->SetPtRange(0.15, 1e10);
  if (trackCutsReject)
  {
    trackCutsReject->SetEtaRange(-0.9, 0.9);
    trackCutsReject->SetPtRange(0.15, 1e10);
  }
  trackCutsTPConly->SetEtaRange(-0.9, 0.9);
  trackCutsTPConly->SetPtRange(0.15, 1e10);

  TString trigName = "";
  if (iPhysicsSelectionFlag == AliVEvent::kAnyINT)
    trigName += "kAnyINT";
  else if (iPhysicsSelectionFlag == AliVEvent::kAny)
    trigName += "kAny";
  else if (iPhysicsSelectionFlag == AliVEvent::kINT7)
    trigName += "kINT7";
  else if (iPhysicsSelectionFlag == AliVEvent::kINT8)
    trigName += "kINT8";
  else if (iPhysicsSelectionFlag == AliVEvent::kMB)
    trigName += "kMB";
  else if (iPhysicsSelectionFlag == AliVEvent::kCentral)
    trigName += "kCentral";
  else if (iPhysicsSelectionFlag == AliVEvent::kSemiCentral)
    trigName += "kSemiCentral";
  else if (iPhysicsSelectionFlag == AliVEvent::kEMC7)
    trigName += "kEMC7";
  else if (iPhysicsSelectionFlag == AliVEvent::kEMCEJE){
    if(emcaltrigger.Length()) trigName += emcaltrigger;
    else trigName += "kEMCEJE";
  }
  else if (iPhysicsSelectionFlag == AliVEvent::kEMCEGA){
    if(emcaltrigger.Length()) trigName += emcaltrigger;
    else trigName += "kEMCEGA";
  }

  //Create the task
  AliPWG4HighPtTrackQA *taskPWG4TrackQA = new AliPWG4HighPtTrackQA(Form("AliPWG4HighPtTrackQACent%dTrack%dCuts%d%s", centClass, trackType, cuts, trigName.Data()));
  taskPWG4TrackQA->SetTrackType(trackType);
  taskPWG4TrackQA->SetCuts(trackCuts);
  taskPWG4TrackQA->SetCutsITSLoose(trackCutsReject);
  taskPWG4TrackQA->SetCutsTPConly(trackCutsTPConly);
  taskPWG4TrackQA->SetPtMax(200.);

  if (iAODanalysis)
    taskPWG4TrackQA->SetDataType(AliPWG4HighPtTrackQA::kAOD);
  else
    taskPWG4TrackQA->SetDataType(AliPWG4HighPtTrackQA::kESD);

  if (isPbPb)
  {
    taskPWG4TrackQA->SetIsPbPb(kTRUE);
    taskPWG4TrackQA->SetCentralityClass(centClass);
  }
  taskPWG4TrackQA->SelectCollisionCandidates(iPhysicsSelectionFlag);
  if(emcaltrigger.Length()) taskPWG4TrackQA->SetEMCALTrigger(emcaltrigger);

  // E. Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  outputfile += Form(":PWG4_HighPtTrackQACent%dTrackType%dCuts%d%s", centClass, trackType, cuts, trigName.Data());

  AliAnalysisDataContainer *cout_histQAtrack = 0x0;
  TString contName = Form("qa_histsQAtrackCent%dType%dcuts%d%s", centClass, trackType, cuts, trigName.Data());
  cout_histQAtrack = mgr->CreateContainer(contName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, outputfile);

  mgr->AddTask(taskPWG4TrackQA);
  mgr->ConnectInput(taskPWG4TrackQA, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskPWG4TrackQA, 1, cout_histQAtrack);

  // Return task pointer at the end
  return taskPWG4TrackQA;
}
