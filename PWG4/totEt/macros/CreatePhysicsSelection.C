#include "AliPhysicsSelection.h"
#include "AliAnalysisEtSelectionContainer.h"
#include "TFile.h"
#include "TSystem.h"
#include <iostream>

int Init()
{
}

int CreatePhysicsSelection(Int_t runNumber, TString filename);
AliPhysicsSelection* DefaultSelection();


AliPhysicsSelection* GetPhysicsSelection()
{
    //Example of selection CMBS1* triggers (up to 137133), MANUAL SETTINGS REQUIRED
   AliPhysicsSelection * physSel = new AliPhysicsSelection();
   physSel->AddCollisionTriggerClass("+CMBAC-B-NOPF-ALL");
   physSel->AddCollisionTriggerClass("+CMBS1C-B-NOPF-ALL");
   physSel->AddCollisionTriggerClass("+CMBS1A-B-NOPF-ALL");
// This are needed only to fill the statistics tables
   physSel->AddBGTriggerClass("+CMBAC-C-NOPF-ALL");
   physSel->AddBGTriggerClass("+CMBS1C-C-NOPF-ALL");
   physSel->AddBGTriggerClass("+CMBS1A-C-NOPF-ALL");
   physSel->AddBGTriggerClass("+CMBAC-A-NOPF-ALL");
   physSel->AddBGTriggerClass("+CMBS1C-A-NOPF-ALL");
   physSel->AddBGTriggerClass("+CMBS1A-A-NOPF-ALL");
   physSel->AddBGTriggerClass("+CMBAC-E-NOPF-ALL");
   physSel->AddBGTriggerClass("+CMBS1C-E-NOPF-ALL");
   physSel->AddBGTriggerClass("+CMBS1A-E-NOPF-ALL");
   
   return physSel;
   
}

int UpdatePhysicsSelection(Int_t runNumber, TString filename)
{
  
  AliPhysicsSelection *selection = GetPhysicsSelection();
  
  TFile *file = TFile::Open(filename, "UPDATE");
  
  AliAnalysisEtSelectionContainer *container = (dynamic_cast<AliAnalysisEtSelectionContainer*>(file->Get("physicsSelections")));
  if(container)
  {
    container->AddDefaultSelection(DefaultSelection());
    container->AddPhysicsSelection(selection, runNumber);
  
    container->Write();
    file->Close();
  }
  else
  {
    file->Close();
    CreatePhysicsSelection(runNumber, filename);
  }
  return 0;
}

int UpdatePhysicsSelection(AliPhysicsSelection *selection, Int_t runNumber, TString filename)
{
  
  TFile *file = TFile::Open(filename, "UPDATE");
  
  AliAnalysisEtSelectionContainer *container = (dynamic_cast<AliAnalysisEtSelectionContainer*>(file->Get("physicsSelections")));

  if(container)
  {
    container->AddDefaultSelection(DefaultSelection());
    container->AddPhysicsSelection(selection, runNumber);
  
    container->Write();
    file->Close();
  }
  
  else
  {
    file->Close();
    CreatePhysicsSelection(runNumber, filename);
  }
  
  return 0;
}

int CreatePhysicsSelection(Int_t runNumber, TString filename)
{
  TFile *file = TFile::Open(filename, "RECREATE");

  AliAnalysisEtSelectionContainer *container = new AliAnalysisEtSelectionContainer("physicsSelections");
  container->Write();
  file->Close();
  
  UpdatePhysicsSelection(runNumber, filename);
  
  return 0;
}

AliPhysicsSelection* DefaultSelection()
{
  return new AliPhysicsSelection();
}




