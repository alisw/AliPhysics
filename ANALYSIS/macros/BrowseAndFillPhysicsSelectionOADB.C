#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>

// ROOT includes
#include "TFile.h"
#include "TString.h"
#include "TSystem.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TCanvas.h"

// STEER includes
#include "AliVEvent.h"

// ANALYSIS include
#include "AliAnalysisManager.h"

// OADB includes
#include "AliOADBContainer.h"
#include "AliOADBTriggerAnalysis.h"
#include "AliOADBPhysicsSelection.h"
#include "AliOADBFillingScheme.h"
#endif


//FIXME: defaults for 2.76 TeV

void BrowseAndFillPhysicsSelectionOADB(Bool_t fill = kFALSE) {

  // Load stuff
  gSystem->Load("libCore.so");  
  gSystem->Load("libTree.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libMinuit");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");   
  gSystem->Load("libOADB");   



  TString oadbfilename(Form("%s/COMMON/PHYSICSSELECTION/data/physicsSelection.root", AliAnalysisManager::GetOADBPath()));
  if (!fill) {
    TFile * f = new TFile(oadbfilename);
    new TBrowser();
    return;
  }


  AliOADBContainer * oadbContPS = new AliOADBContainer("physSel");
  AliOADBContainer * oadbContFillingScheme = new AliOADBContainer("fillScheme");
  AliOADBContainer * oadbContTriggerAnalysis = new AliOADBContainer("trigAnalysis");

  // Defaults
  AliOADBFillingScheme * fsDefault = new AliOADBFillingScheme("Default");
  fsDefault->SetFillingSchemeName("Default");
  fsDefault->SetBXIDs("B","");
  fsDefault->SetBXIDs("A","" ); 
  fsDefault->SetBXIDs("AC","" ); 
  fsDefault->SetBXIDs("C","");
  fsDefault->SetBXIDs("E",""       );	   
  oadbContFillingScheme->AddDefaultObject(fsDefault);

  // DefaultPP
  AliOADBPhysicsSelection * oadbDefaultPP = new AliOADBPhysicsSelection("oadbDefaultPP");
  oadbDefaultPP->AddCollisionTriggerClass   ( AliVEvent::kINT7,"+CINT7-I-NOPF-ALLNOTRD","B",0);
  oadbDefaultPP->AddBGTriggerClass          ( AliVEvent::kINT7,"+CINT7-AC-NOPF-ALLNOTRD","AC",0);
  oadbDefaultPP->AddBGTriggerClass          ( AliVEvent::kINT7,"+CINT7-E-NOPF-ALLNOTRD","E",0);  
  oadbDefaultPP->SetHardwareTrigger         ( 0,"V0A && V0C");					      
  oadbDefaultPP->SetOfflineTrigger          ( 0,"(V0A && V0C) && !V0ABG && !V0CBG  && !TPCLaserWarmUp");

  oadbDefaultPP->AddCollisionTriggerClass   ( AliVEvent::kMUSH7,"+CMUSH7-B-NOPF-MUON","B",  1);
  oadbDefaultPP->AddBGTriggerClass          ( AliVEvent::kMUSH7,"+CMUSH7-AC-NOPF-MUON","AC",1);
  oadbDefaultPP->AddBGTriggerClass          ( AliVEvent::kMUSH7,"+CMUSH7-E-NOPF-MUON","E",  1);
  oadbDefaultPP->SetHardwareTrigger         ( 1,"V0A && V0C");                                         
  oadbDefaultPP->SetOfflineTrigger          ( 1,"(V0A && V0C) && !V0ABG && !V0CBG  && !TPCLaserWarmUp");
    
  oadbDefaultPP->AddCollisionTriggerClass   ( AliVEvent::kMUL7,"+CMUL7-B-NOPF-MUON","B",  2);
  oadbDefaultPP->AddBGTriggerClass          ( AliVEvent::kMUL7,"+CMUL7-AC-NOPF-MUON","AC",2);
  oadbDefaultPP->AddBGTriggerClass          ( AliVEvent::kMUL7,"+CMUL7-E-NOPF-MUON","E",  2);
  oadbDefaultPP->SetHardwareTrigger         ( 2,"V0A && V0C");                                         
  oadbDefaultPP->SetOfflineTrigger          ( 2,"(V0A && V0C) && !V0ABG && !V0CBG  && !TPCLaserWarmUp");

  oadbDefaultPP->AddCollisionTriggerClass   ( AliVEvent::kMUU7,"+CMUU7-B-NOPF-ALLNOTRD","B",  3);
  oadbDefaultPP->AddBGTriggerClass          ( AliVEvent::kMUU7,"+CMUU7-AC-NOPF-ALLNOTRD","AC",3);
  oadbDefaultPP->AddBGTriggerClass          ( AliVEvent::kMUU7,"+CMUU7-E-NOPF-ALLNOTRD","E",  3);
  oadbDefaultPP->SetHardwareTrigger         ( 3,"V0A && V0C");                                         
  oadbDefaultPP->SetOfflineTrigger          ( 3,"(V0A && V0C) && !V0ABG && !V0CBG  && !TPCLaserWarmUp");
  
  oadbDefaultPP->AddCollisionTriggerClass   ( AliVEvent::kEMC7,"+CEMC7-B-NOPF-ALLNOTRD","B",  4);
  oadbDefaultPP->AddBGTriggerClass          ( AliVEvent::kEMC7,"+CEMC7-AC-NOPF-ALLNOTRD","AC",4);
  oadbDefaultPP->AddBGTriggerClass          ( AliVEvent::kEMC7,"+CEMC7-E-NOPF-ALLNOTRD","E",  4);
  oadbDefaultPP->SetHardwareTrigger         ( 4,"V0A && V0C");                      
  oadbDefaultPP->SetOfflineTrigger          ( 4,"(V0A && V0C) && !V0ABG && !V0CBG  && !TPCLaserWarmUp");
  
  oadbDefaultPP->AddCollisionTriggerClass   ( AliVEvent::kMUS7,"+CMUS7-B-NOPF-MUON","B",  5);
  oadbDefaultPP->AddBGTriggerClass          ( AliVEvent::kMUS7,"+CMUS7-AC-NOPF-MUON","AC",5);
  oadbDefaultPP->AddBGTriggerClass          ( AliVEvent::kMUS7,"+CMUS7-E-NOPF-MUON","E",  5);
  oadbDefaultPP->SetHardwareTrigger         ( 5,"V0A && V0C");                                         
  oadbDefaultPP->SetOfflineTrigger          ( 5,"(V0A && V0C) && !V0ABG && !V0CBG  && !TPCLaserWarmUp");

  oadbDefaultPP->AddCollisionTriggerClass   ( AliVEvent::kPHI7,"+CPHI7-I-NOPF-ALLNOTRD","B",  6);
  oadbDefaultPP->AddCollisionTriggerClass   ( AliVEvent::kPHI7,"+CPHI7-B-NOPF-ALLNOTRD","B",  6);
  oadbDefaultPP->AddBGTriggerClass          ( AliVEvent::kPHI7,"+CPHI7-AC-NOPF-ALLNOTRD","AC",6);
  oadbDefaultPP->AddBGTriggerClass          ( AliVEvent::kPHI7,"+CPHI7-E-NOPF-ALLNOTRD","E",  6);
  oadbDefaultPP->SetHardwareTrigger         ( 6,"V0A && V0C");
  oadbDefaultPP->SetOfflineTrigger          ( 6,"(V0A && V0C) && !V0ABG && !V0CBG  && !TPCLaserWarmUp");
  
  /*
  // PENDING VALIDATION https://savannah.cern.ch/bugs/?86829
  oadbDefaultPP->AddCollisionTriggerClass   ( AliVEvent::kDG5,"+CDG5-I-NOPF-ALLNOTRD","B",  7);
  oadbDefaultPP->AddBGTriggerClass          ( AliVEvent::kDG5,"+CDG5-AC-NOPF-ALLNOTRD","AC",7);
  oadbDefaultPP->AddBGTriggerClass          ( AliVEvent::kDG5,"+CDG5-E-NOPF-ALLNOTRD","E",  7);
  // TODO TOF condition needs to be implemented (Savannah bug #87064)
  oadbDefaultPP->SetHardwareTrigger         ( 7,"!V0A && !V0C && SPDGFOL1 >= 2");
  oadbDefaultPP->SetOfflineTrigger          ( 7,"!V0A && !V0C && !V0ABG && !V0CBG && && SPDGFOL1 >= 2 && !TPCLaserWarmUp");
  */

  oadbContPS->AddDefaultObject(oadbDefaultPP);

  // DefaultPbPb
  AliOADBPhysicsSelection * oadbDefaultPbPb = new AliOADBPhysicsSelection("oadbDefaultPbPb");
  
  Int_t triggerCount = 0;
  oadbDefaultPbPb->AddCollisionTriggerClass   ( AliVEvent::kMB,"+CPBI1-B-NOPF-ALLNOTRD,CPBI2_B1-B-NOPF-ALLNOTRD,CPBI2_B1-B-PF-ALLNOTRD","B",triggerCount);
  oadbDefaultPbPb->AddBGTriggerClass          ( AliVEvent::kMB,"+CPBI1-AC-NOPF-ALLNOTRD,CPBI2_B1-AC-NOPF-ALLNOTRD","AC",triggerCount);
  oadbDefaultPbPb->AddBGTriggerClass          ( AliVEvent::kMB,"+CPBI1-E-NOPF-ALLNOTRD,CPBI2_B1-E-NOPF-ALLNOTRD","E",triggerCount);
  oadbDefaultPbPb->SetHardwareTrigger         ( triggerCount,"V0A && V0C");
  oadbDefaultPbPb->SetOfflineTrigger          ( triggerCount,"V0A && V0C && !V0ABG && !V0CBG && !TPCLaserWarmUp ");

  triggerCount++;
  oadbDefaultPbPb->AddCollisionTriggerClass   ( AliVEvent::kCentral,"+CVHN-B-NOPF-ALLNOTRD,CVHN-B-NOPF-CENTNOTRD,CVHN-B-PF-ALLNOTRD,CVHN-B-PF-CENTNOTRD","B",triggerCount);
  oadbDefaultPbPb->AddBGTriggerClass          ( AliVEvent::kCentral,"+CVHN-AC-NOPF-ALLNOTRD,CVHN-AC-NOPF-CENTNOTRD","AC",triggerCount);
  oadbDefaultPbPb->AddBGTriggerClass          ( AliVEvent::kCentral,"+CVHN-E-NOPF-ALLNOTRD,CVHN-E-NOPF-CENTNOTRD","E",triggerCount);
  oadbDefaultPbPb->SetHardwareTrigger         ( triggerCount,"V0A && V0C && Central");
  oadbDefaultPbPb->SetOfflineTrigger          ( triggerCount,"V0A && V0C && !V0ABG && !V0CBG && !TPCLaserWarmUp ");

  triggerCount++;
  // semicentral includes central ones (and might have different downscaling)
  oadbDefaultPbPb->AddCollisionTriggerClass   ( AliVEvent::kSemiCentral,"+CVHN-B-NOPF-ALLNOTRD,CVHN-B-NOPF-CENTNOTRD,CVHN-B-PF-ALLNOTRD,CVHN-B-PF-CENTNOTRD,CVLN-B-NOPF-ALLNOTRD,CVLN-B-PF-ALLNOTRD","B",triggerCount);
  oadbDefaultPbPb->AddBGTriggerClass          ( AliVEvent::kSemiCentral,"+CVHN-AC-NOPF-ALLNOTRD,CVHN-AC-NOPF-CENTNOTRD,CVLN-AC-NOPF-ALLNOTRD","AC",triggerCount);
  oadbDefaultPbPb->AddBGTriggerClass          ( AliVEvent::kSemiCentral,"+CVHN-E-NOPF-ALLNOTRD,CVHN-E-NOPF-CENTNOTRD,CVLN-E-NOPF-ALLNOTRD","E",triggerCount);
  oadbDefaultPbPb->SetHardwareTrigger         ( triggerCount,"V0A && V0C && SemiCentral && !Central");
  oadbDefaultPbPb->SetOfflineTrigger          ( triggerCount,"V0A && V0C && !V0ABG && !V0CBG && !TPCLaserWarmUp ");

  triggerCount++;
  oadbDefaultPbPb->AddCollisionTriggerClass   ( AliVEvent::kEMCEJE,"+CPBI2EJE-B-NOPF-CENTNOTRD","B",triggerCount);
  oadbDefaultPbPb->AddBGTriggerClass          ( AliVEvent::kEMCEJE,"+CPBI2EJE-AC-NOPF-CENTNOTRD","AC",triggerCount);
  oadbDefaultPbPb->AddBGTriggerClass          ( AliVEvent::kEMCEJE,"+CPBI2EJE-E-NOPF-CENTNOTRD","E",triggerCount);
  oadbDefaultPbPb->SetHardwareTrigger         ( triggerCount,"V0A && V0C");
  // TODO EMC offline check missing, https://savannah.cern.ch/bugs/index.php?87104
  oadbDefaultPbPb->SetOfflineTrigger          ( triggerCount,"V0A && V0C && !V0ABG && !V0CBG && !TPCLaserWarmUp ");
  
  triggerCount++;
  oadbDefaultPbPb->AddCollisionTriggerClass   ( AliVEvent::kEMCEGA,"+CPBI2EGA-B-NOPF-CENTNOTRD","B",triggerCount);
  oadbDefaultPbPb->AddBGTriggerClass          ( AliVEvent::kEMCEGA,"+CPBI2EGA-AC-NOPF-CENTNOTRD","AC",triggerCount);
  oadbDefaultPbPb->AddBGTriggerClass          ( AliVEvent::kEMCEGA,"+CPBI2EGA-E-NOPF-CENTNOTRD","E",triggerCount);
  oadbDefaultPbPb->SetHardwareTrigger         ( triggerCount,"V0A && V0C");
  // TODO EMC offline check missing, https://savannah.cern.ch/bugs/index.php?87104
  oadbDefaultPbPb->SetOfflineTrigger          ( triggerCount,"V0A && V0C && !V0ABG && !V0CBG && !TPCLaserWarmUp ");
  
  triggerCount++;
  oadbDefaultPbPb->AddCollisionTriggerClass   ( AliVEvent::kMUSPB,"+CPBI1MSL-B-NOPF-MUON","B",  triggerCount);
  oadbDefaultPbPb->AddBGTriggerClass          ( AliVEvent::kMUSPB,"+CPBI1MSL-ACE-NOPF-MUON","ACE",triggerCount);
  oadbDefaultPbPb->SetHardwareTrigger         ( triggerCount,"V0A && V0C");                                         
  oadbDefaultPbPb->SetOfflineTrigger          ( triggerCount,"V0A && V0C && !V0ABG && !V0CBG && !TPCLaserWarmUp");

  triggerCount++;
  oadbDefaultPbPb->AddCollisionTriggerClass   ( AliVEvent::kMUSHPB,"+CPBI1MSH-B-NOPF-MUON","B",  triggerCount);
  oadbDefaultPbPb->AddBGTriggerClass          ( AliVEvent::kMUSHPB,"+CPBI1MSH-ACE-NOPF-MUON","ACE",triggerCount);
  oadbDefaultPbPb->SetHardwareTrigger         ( triggerCount,"V0A && V0C");                                         
  oadbDefaultPbPb->SetOfflineTrigger          ( triggerCount,"V0A && V0C && !V0ABG && !V0CBG && !TPCLaserWarmUp");
    
  triggerCount++;
  oadbDefaultPbPb->AddCollisionTriggerClass   ( AliVEvent::kMuonUnlikePB,"+CPBI1MUL-B-NOPF-MUON","B",  triggerCount);
  oadbDefaultPbPb->AddBGTriggerClass          ( AliVEvent::kMuonUnlikePB,"+CPBI1MUL-AC-NOPF-MUON","AC",triggerCount);
  oadbDefaultPbPb->AddBGTriggerClass          ( AliVEvent::kMuonUnlikePB,"+CPBI1MUL-E-NOPF-MUON","E",  triggerCount);
  oadbDefaultPbPb->SetHardwareTrigger         ( triggerCount,"V0A && V0C");                                         
  oadbDefaultPbPb->SetOfflineTrigger          ( triggerCount,"V0A && V0C && !V0ABG && !V0CBG && !TPCLaserWarmUp");

  triggerCount++;
  oadbDefaultPbPb->AddCollisionTriggerClass   ( AliVEvent::kMuonLikePB,"+CPBI1MLL-B-NOPF-MUON","B",  triggerCount);
  oadbDefaultPbPb->AddBGTriggerClass          ( AliVEvent::kMuonLikePB,"+CPBI1MLL-AC-NOPF-MUON","AC",triggerCount);
  oadbDefaultPbPb->AddBGTriggerClass          ( AliVEvent::kMuonLikePB,"+CPBI1MLL-E-NOPF-MUON","E",  triggerCount);
  oadbDefaultPbPb->SetHardwareTrigger         ( triggerCount,"V0A && V0C");                                         
  oadbDefaultPbPb->SetOfflineTrigger          ( triggerCount,"V0A && V0C && !V0ABG && !V0CBG && !TPCLaserWarmUp");
  
  oadbContPS->AddDefaultObject(oadbDefaultPbPb);

  // Trigger analysis defaults
  AliOADBTriggerAnalysis * oadbTrigAnalysis = new AliOADBTriggerAnalysis("Default");
  oadbTrigAnalysis->SetZDCCorrParameters(-0.823, -0.653, 4*0.58, 4*0.5);
  oadbContTriggerAnalysis->AddDefaultObject(oadbTrigAnalysis);  

  // ----- 2009 - 2010 -----
  // ----- proton-proton -----
  // LHC09d+LHC10(abcde)
  AliOADBPhysicsSelection * oadbLHC09d10e = new AliOADBPhysicsSelection("oadbLHC09d10e");
  oadbLHC09d10e->AddCollisionTriggerClass   ( AliVEvent::kMB,"+CINT1B-ABCE-NOPF-ALL","B",0);
  oadbLHC09d10e->AddBGTriggerClass          ( AliVEvent::kMB,"+CINT1A-ABCE-NOPF-ALL","A",0);
  oadbLHC09d10e->AddBGTriggerClass          ( AliVEvent::kMB,"+CINT1C-ABCE-NOPF-ALL","C",0);
  oadbLHC09d10e->AddBGTriggerClass          ( AliVEvent::kMB,"+CINT1-E-NOPF-ALL","E",0);  
  oadbLHC09d10e->SetHardwareTrigger         ( 0,"SPDGFO >= 1 || V0A || V0C");
  oadbLHC09d10e->SetOfflineTrigger          ( 0,"(SPDGFO >= 1 || V0A || V0C) && !V0ABG && !V0CBG  && !TPCLaserWarmUp && !SPDClsVsTrkBG");
  
  oadbLHC09d10e->AddCollisionTriggerClass   ( AliVEvent::kMUON,"+CMUS1B-ABCE-NOPF-MUON","B",1);
  oadbLHC09d10e->AddBGTriggerClass          ( AliVEvent::kMUON,"+CMUS1A-ABCE-NOPF-MUON","A",1);
  oadbLHC09d10e->AddBGTriggerClass          ( AliVEvent::kMUON,"+CMUS1C-ABCE-NOPF-MUON","C",1);
  oadbLHC09d10e->AddBGTriggerClass          ( AliVEvent::kMUON,"+CMUS1-E-NOPF-MUON","E",1);
  oadbLHC09d10e->SetHardwareTrigger         ( 1,"SPDGFO >= 1 || V0A || V0C");                                         
  oadbLHC09d10e->SetOfflineTrigger          ( 1,"(SPDGFO >= 1 || V0A || V0C) && !V0ABG && !V0CBG  && !TPCLaserWarmUp && !SPDClsVsTrkBG");
    
  oadbContPS->AppendObject(oadbLHC09d10e, 104065,118555);
  oadbContPS->AppendObject(oadbLHC09d10e->Clone(), 118562,127711);

  // LHC10c, fill 1069 (problems with the V0 online trigger in ESD)
  AliOADBPhysicsSelection * oadbLHC10cV0Bug = new AliOADBPhysicsSelection("oadbLHC10cV0Bug");
  oadbLHC10cV0Bug->AddCollisionTriggerClass   ( AliVEvent::kMB,"+CINT1B-ABCE-NOPF-ALL","B",0);
  oadbLHC10cV0Bug->AddBGTriggerClass          ( AliVEvent::kMB,"+CINT1A-ABCE-NOPF-ALL","A",0);
  oadbLHC10cV0Bug->AddBGTriggerClass          ( AliVEvent::kMB,"+CINT1C-ABCE-NOPF-ALL","C",0);
  oadbLHC10cV0Bug->AddBGTriggerClass          ( AliVEvent::kMB,"+CINT1-E-NOPF-ALL","E",0);  
  oadbLHC10cV0Bug->SetHardwareTrigger         ( 0,"SPDGFO >= 1 || CTPV0A || CTPV0C");
  oadbLHC10cV0Bug->SetOfflineTrigger          ( 0,"(SPDGFO >= 1 || V0A || V0C) && !V0ABG && !V0CBG  && !TPCLaserWarmUp && !SPDClsVsTrkBG");
  
  oadbLHC10cV0Bug->AddCollisionTriggerClass   ( AliVEvent::kMUON,"+CMUS1B-ABCE-NOPF-MUON","B",1);
  oadbLHC10cV0Bug->AddBGTriggerClass          ( AliVEvent::kMUON,"+CMUS1A-ABCE-NOPF-MUON","A",1);
  oadbLHC10cV0Bug->AddBGTriggerClass          ( AliVEvent::kMUON,"+CMUS1C-ABCE-NOPF-MUON","C",1);
  oadbLHC10cV0Bug->AddBGTriggerClass          ( AliVEvent::kMUON,"+CMUS1-E-NOPF-MUON","E",1);
  oadbLHC10cV0Bug->SetHardwareTrigger         ( 1,"SPDGFO >= 1 || CTPV0A || CTPV0C");                                         
  oadbLHC10cV0Bug->SetOfflineTrigger          ( 1,"(SPDGFO >= 1 || V0A || V0C) && !V0ABG && !V0CBG  && !TPCLaserWarmUp && !SPDClsVsTrkBG");

  oadbContPS->AppendObject(oadbLHC10cV0Bug, 118556,118561);
  

  // LHC10e1
  AliOADBPhysicsSelection * oadbLHC10e1 = new AliOADBPhysicsSelection("oadbLHC10e1");
  oadbLHC10e1->AddCollisionTriggerClass   ( AliVEvent::kMB,"+CINT1-B-NOPF-ALLNOTRD","B",0);
  oadbLHC10e1->AddBGTriggerClass          ( AliVEvent::kMB,"+CINT1-AC-NOPF-ALLNOTRD","AC",0);
  oadbLHC10e1->AddBGTriggerClass          ( AliVEvent::kMB,"+CINT1-E-NOPF-ALLNOTRD","E",0);  
  oadbLHC10e1->SetHardwareTrigger         ( 0,"SPDGFO >= 1 || V0A || V0C");
  oadbLHC10e1->SetOfflineTrigger          ( 0,"(SPDGFO >= 1 || V0A || V0C) && !V0ABG && !V0CBG  && !TPCLaserWarmUp && !SPDClsVsTrkBG");
  
  oadbLHC10e1->AddCollisionTriggerClass   ( AliVEvent::kMUON,"+CMUS1-B-NOPF-ALLNOTRD","B",1);
  oadbLHC10e1->AddBGTriggerClass          ( AliVEvent::kMUON,"+CMUS1-AC-NOPF-ALLNOTRD","AC",1);
  oadbLHC10e1->AddBGTriggerClass          ( AliVEvent::kMUON,"+CMUS1-E-NOPF-ALLNOTRD","E",1);
  oadbLHC10e1->SetHardwareTrigger         ( 1,"SPDGFO >= 1 || V0A || V0C");                                         
  oadbLHC10e1->SetOfflineTrigger          ( 1,"(SPDGFO >= 1 || V0A || V0C) && !V0ABG && !V0CBG  && !TPCLaserWarmUp && !SPDClsVsTrkBG");
    
  oadbLHC10e1->AddCollisionTriggerClass   ( AliVEvent::kHighMult,"+CSH1-B-NOPF-ALLNOTRD","B",2);
  oadbLHC10e1->AddBGTriggerClass          ( AliVEvent::kHighMult,"+CSH1-AC-NOPF-ALLNOTRD","AC",2);
  oadbLHC10e1->AddBGTriggerClass          ( AliVEvent::kHighMult,"+CSH1-E-NOPF-ALLNOTRD","E",2);
  oadbLHC10e1->SetHardwareTrigger         ( 2,"SPDGFO >= 1 || V0A || V0C");                      
  oadbLHC10e1->SetOfflineTrigger          ( 2,"(SPDGFO >= 1 || V0A || V0C) && !V0ABG && !V0CBG  && !TPCLaserWarmUp && !SPDClsVsTrkBG");
    
  oadbContPS->AppendObject(oadbLHC10e1, 127712,127718);

  // lhc10e2
  AliOADBPhysicsSelection * oadbLHC10e2 = new AliOADBPhysicsSelection("oadbLHC10e2");
  oadbLHC10e2->AddCollisionTriggerClass   ( AliVEvent::kMB,"+CINT1B-ABCE-NOPF-ALL","B",0);
  oadbLHC10e2->AddBGTriggerClass          ( AliVEvent::kMB,"+CINT1A-ABCE-NOPF-ALL","A",0);
  oadbLHC10e2->AddBGTriggerClass          ( AliVEvent::kMB,"+CINT1C-ABCE-NOPF-ALL","C",0);
  oadbLHC10e2->AddBGTriggerClass          ( AliVEvent::kMB,"+CINT1-E-NOPF-ALL","E",0);  
  oadbLHC10e2->SetHardwareTrigger         ( 0,"SPDGFO >= 1 || V0A || V0C");
  oadbLHC10e2->SetOfflineTrigger          ( 0,"(SPDGFO >= 1 || V0A || V0C) && !V0ABG && !V0CBG  && !TPCLaserWarmUp && !SPDClsVsTrkBG");
  
  oadbLHC10e2->AddCollisionTriggerClass   ( AliVEvent::kMUON,"+CMUS1B-ABCE-NOPF-MUON","B",1);
  oadbLHC10e2->AddBGTriggerClass          ( AliVEvent::kMUON,"+CMUS1A-ABCE-NOPF-MUON","A",1);
  oadbLHC10e2->AddBGTriggerClass          ( AliVEvent::kMUON,"+CMUS1C-ABCE-NOPF-MUON","C",1);
  oadbLHC10e2->AddBGTriggerClass          ( AliVEvent::kMUON,"+CMUS1-E-NOPF-MUON","E",1);
  oadbLHC10e2->SetHardwareTrigger         ( 1,"SPDGFO >= 1 || V0A || V0C");                                         
  oadbLHC10e2->SetOfflineTrigger          ( 1,"(SPDGFO >= 1 || V0A || V0C) && !V0ABG && !V0CBG  && !TPCLaserWarmUp && !SPDClsVsTrkBG");
    
  oadbContPS->AppendObject(oadbLHC10e2, 127719,127730);

  // LHC10e3
  AliOADBPhysicsSelection * oadbLHC10e3 = new AliOADBPhysicsSelection("oadbLHC10e3");
  oadbLHC10e3->AddCollisionTriggerClass   ( AliVEvent::kMB,"+CINT1-B-NOPF-ALLNOTRD","B",0);
  oadbLHC10e3->AddBGTriggerClass          ( AliVEvent::kMB,"+CINT1-AC-NOPF-ALLNOTRD","AC",0);
  oadbLHC10e3->AddBGTriggerClass          ( AliVEvent::kMB,"+CINT1-E-NOPF-ALLNOTRD","E",0);  
  oadbLHC10e3->SetHardwareTrigger         ( 0,"SPDGFO >= 1 || V0A || V0C");
  oadbLHC10e3->SetOfflineTrigger          ( 0,"(SPDGFO >= 1 || V0A || V0C) && !V0ABG && !V0CBG  && !TPCLaserWarmUp && !SPDClsVsTrkBG");
  
  oadbLHC10e3->AddCollisionTriggerClass   ( AliVEvent::kMUON,"+CMUS1-B-NOPF-ALLNOTRD","B",1);
  oadbLHC10e3->AddBGTriggerClass          ( AliVEvent::kMUON,"+CMUS1-AC-NOPF-ALLNOTRD","AC",1);
  oadbLHC10e3->AddBGTriggerClass          ( AliVEvent::kMUON,"+CMUS1-E-NOPF-ALLNOTRD","E",1);
  oadbLHC10e3->SetHardwareTrigger         ( 1,"SPDGFO >= 1 || V0A || V0C");                                         
  oadbLHC10e3->SetOfflineTrigger          ( 1,"(SPDGFO >= 1 || V0A || V0C) && !V0ABG && !V0CBG  && !TPCLaserWarmUp && !SPDClsVsTrkBG");
    
  oadbLHC10e3->AddCollisionTriggerClass   ( AliVEvent::kHighMult,"+CSH1-B-NOPF-ALLNOTRD","B",2);
  oadbLHC10e3->AddBGTriggerClass          ( AliVEvent::kHighMult,"+CSH1-AC-NOPF-ALLNOTRD","AC",2);
  oadbLHC10e3->AddBGTriggerClass          ( AliVEvent::kHighMult,"+CSH1-E-NOPF-ALLNOTRD","E",2);
  oadbLHC10e3->SetHardwareTrigger         ( 2,"SPDGFO >= 1 || V0A || V0C");                                         
  oadbLHC10e3->SetOfflineTrigger          ( 2,"(SPDGFO >= 1 || V0A || V0C) && !V0ABG && !V0CBG  && !TPCLaserWarmUp && !SPDClsVsTrkBG");
    
  oadbContPS->AppendObject(oadbLHC10e3, 127813,130850);


  // LHC10f1
  AliOADBPhysicsSelection * oadbLHC10f1 = new AliOADBPhysicsSelection("oadbLHC10f1");
  oadbLHC10f1->AddCollisionTriggerClass   ( AliVEvent::kMB,"+CINT1-B-NOPF-ALLNOTRD","B",0);
  oadbLHC10f1->AddBGTriggerClass          ( AliVEvent::kMB,"+CINT1-AC-NOPF-ALLNOTRD","AC",0);
  oadbLHC10f1->AddBGTriggerClass          ( AliVEvent::kMB,"+CINT1-E-NOPF-ALLNOTRD","E",0);  
  oadbLHC10f1->SetHardwareTrigger         ( 0,"SPDGFO >= 1 || V0A || V0C");
  oadbLHC10f1->SetOfflineTrigger          ( 0,"(SPDGFO >= 1 || V0A || V0C) && !V0ABG && !V0CBG  && !TPCLaserWarmUp && !SPDClsVsTrkBG");
  
  oadbLHC10f1->AddCollisionTriggerClass   ( AliVEvent::kMUON,"+CMUS1-B-NOPF-ALLNOTRD","B",1);
  oadbLHC10f1->AddBGTriggerClass          ( AliVEvent::kMUON,"+CMUS1-AC-NOPF-ALLNOTRD","AC",1);
  oadbLHC10f1->AddBGTriggerClass          ( AliVEvent::kMUON,"+CMUS1-E-NOPF-ALLNOTRD","E",1);
  oadbLHC10f1->SetHardwareTrigger         ( 1,"SPDGFO >= 1 || V0A || V0C");                                         
  oadbLHC10f1->SetOfflineTrigger          ( 1,"(SPDGFO >= 1 || V0A || V0C) && !V0ABG && !V0CBG  && !TPCLaserWarmUp && !SPDClsVsTrkBG");
    
  oadbLHC10f1->AddCollisionTriggerClass   ( AliVEvent::kHighMult,"+CSH1-B-NOPF-ALLNOTRD","B",2);
  oadbLHC10f1->AddBGTriggerClass          ( AliVEvent::kHighMult,"+CSH1-AC-NOPF-ALLNOTRD","AC",2);
  oadbLHC10f1->AddBGTriggerClass          ( AliVEvent::kHighMult,"+CSH1-E-NOPF-ALLNOTRD","E",2);
  oadbLHC10f1->SetHardwareTrigger         ( 2,"SPDGFO >= 1 || V0A || V0C");                                         
  oadbLHC10f1->SetOfflineTrigger          ( 2,"(SPDGFO >= 1 || V0A || V0C) && !V0ABG && !V0CBG  && !TPCLaserWarmUp && !SPDClsVsTrkBG");
    
  oadbContPS->AppendObject(oadbLHC10f1, 133004,134690);

  // LHC10f2
  AliOADBPhysicsSelection * oadbLHC10f2 = new AliOADBPhysicsSelection("oadbLHC10f2");
  oadbLHC10f2->AddCollisionTriggerClass   ( AliVEvent::kCINT5,"+CINT5-B-NOPF-ALLNOTRD","B",0);
  oadbLHC10f2->AddBGTriggerClass          ( AliVEvent::kCINT5,"+CINT5-AC-NOPF-ALLNOTRD","AC",0);
  oadbLHC10f2->AddBGTriggerClass          ( AliVEvent::kCINT5,"+CINT5-E-NOPF-ALLNOTRD","E",0);  
  oadbLHC10f2->SetHardwareTrigger         ( 0,"V0A || V0C");
  oadbLHC10f2->SetOfflineTrigger          ( 0,"(V0A || V0C) && !V0ABG && !V0CBG  && !TPCLaserWarmUp && !SPDClsVsTrkBG");
  
  oadbLHC10f2->AddCollisionTriggerClass   ( AliVEvent::kCMUS5,"+CMUS5-B-NOPF-ALLNOTRD","B",1);
  oadbLHC10f2->AddBGTriggerClass          ( AliVEvent::kCMUS5,"+CMUS5-AC-NOPF-ALLNOTRD","AC",1);
  oadbLHC10f2->AddBGTriggerClass          ( AliVEvent::kCMUS5,"+CMUS5-E-NOPF-ALLNOTRD","E",1);
  oadbLHC10f2->SetHardwareTrigger         ( 1,"V0A || V0C");                                         
  oadbLHC10f2->SetOfflineTrigger          ( 1,"(V0A || V0C) && !V0ABG && !V0CBG  && !TPCLaserWarmUp && !SPDClsVsTrkBG");
    
  oadbContPS->AppendObject(oadbLHC10f2, 134776,134931);

  // LHC10f3
  AliOADBPhysicsSelection * oadbLHC10f3 = new AliOADBPhysicsSelection("oadbLHC10f3");
  oadbLHC10f3->AddCollisionTriggerClass   ( AliVEvent::kMB,"+CINT1-B-NOPF-ALLNOTRD","B",0);
  oadbLHC10f3->AddBGTriggerClass          ( AliVEvent::kMB,"+CINT1-AC-NOPF-ALLNOTRD","AC",0);
  oadbLHC10f3->AddBGTriggerClass          ( AliVEvent::kMB,"+CINT1-E-NOPF-ALLNOTRD","E",0);  
  oadbLHC10f3->SetHardwareTrigger         ( 0,"SPDGFO >= 1 || V0A || V0C");
  oadbLHC10f3->SetOfflineTrigger          ( 0,"(SPDGFO >= 1 || V0A || V0C) && !V0ABG && !V0CBG  && !TPCLaserWarmUp && !SPDClsVsTrkBG");
  
  oadbLHC10f3->AddCollisionTriggerClass   ( AliVEvent::kMUON,"+CMUS1-B-NOPF-ALLNOTRD","B",1);
  oadbLHC10f3->AddBGTriggerClass          ( AliVEvent::kMUON,"+CMUS1-AC-NOPF-ALLNOTRD","AC",1);
  oadbLHC10f3->AddBGTriggerClass          ( AliVEvent::kMUON,"+CMUS1-E-NOPF-ALLNOTRD","E",1);
  oadbLHC10f3->SetHardwareTrigger         ( 1,"SPDGFO >= 1 || V0A || V0C");                                         
  oadbLHC10f3->SetOfflineTrigger          ( 1,"(SPDGFO >= 1 || V0A || V0C) && !V0ABG && !V0CBG  && !TPCLaserWarmUp && !SPDClsVsTrkBG");
    
  oadbLHC10f3->AddCollisionTriggerClass   ( AliVEvent::kHighMult,"+CSH1-B-NOPF-ALLNOTRD","B",2);
  oadbLHC10f3->AddBGTriggerClass          ( AliVEvent::kHighMult,"+CSH1-AC-NOPF-ALLNOTRD","AC",2);
  oadbLHC10f3->AddBGTriggerClass          ( AliVEvent::kHighMult,"+CSH1-E-NOPF-ALLNOTRD","E",2);
  oadbLHC10f3->SetHardwareTrigger         ( 1,"SPDGFO >= 1 || V0A || V0C");                                         
  oadbLHC10f3->SetOfflineTrigger          ( 1,"(SPDGFO >= 1 || V0A || V0C) && !V0ABG && !V0CBG  && !TPCLaserWarmUp && !SPDClsVsTrkBG");
    
  oadbContPS->AppendObject(oadbLHC10f3, 135029,135029);


  // LHC10g
  AliOADBPhysicsSelection * oadbLHC10g = new AliOADBPhysicsSelection("oadbLHC10g");
  oadbLHC10g->AddCollisionTriggerClass   ( AliVEvent::kMB,"+CINT1-B-NOPF-ALLNOTRD","B",0);
  oadbLHC10g->AddBGTriggerClass          ( AliVEvent::kMB,"+CINT1-AC-NOPF-ALLNOTRD","AC",0);
  oadbLHC10g->AddBGTriggerClass          ( AliVEvent::kMB,"+CINT1-E-NOPF-ALLNOTRD","E",0);  
  oadbLHC10g->SetHardwareTrigger         ( 0,"SPDGFO >= 1 || V0A || V0C");
  oadbLHC10g->SetOfflineTrigger          ( 0,"(SPDGFO >= 1 || V0A || V0C) && !V0ABG && !V0CBG  && !TPCLaserWarmUp && !SPDClsVsTrkBG");
  
  oadbLHC10g->AddCollisionTriggerClass   ( AliVEvent::kMUON,"+CMUS1-B-NOPF-ALLNOTRD","B",1);
  oadbLHC10g->AddBGTriggerClass          ( AliVEvent::kMUON,"+CMUS1-AC-NOPF-ALLNOTRD","AC",1);
  oadbLHC10g->AddBGTriggerClass          ( AliVEvent::kMUON,"+CMUS1-E-NOPF-ALLNOTRD","E",1);
  oadbLHC10g->SetHardwareTrigger         ( 1,"SPDGFO >= 1 || V0A || V0C");                                         
  oadbLHC10g->SetOfflineTrigger          ( 1,"(SPDGFO >= 1 || V0A || V0C) && !V0ABG && !V0CBG  && !TPCLaserWarmUp && !SPDClsVsTrkBG");
    
  oadbContPS->AppendObject(oadbLHC10g, 135654,136377);


  // filling schemes

  AliOADBFillingScheme * fs4by4a = new AliOADBFillingScheme("4x4a");
  fs4by4a->SetFillingSchemeName("4x4a");
  fs4by4a->SetBXIDs("B"," #2128 #3019");
  fs4by4a->SetBXIDs("A"," #346 #3465" ); 
  fs4by4a->SetBXIDs("AC",""); 
  fs4by4a->SetBXIDs("C"," #1234 #1680");
  fs4by4a->SetBXIDs("E"," #790");	   
  oadbContFillingScheme->AppendObject(fs4by4a, 104065, 104160);

  AliOADBFillingScheme * fs4by4astar = new AliOADBFillingScheme("4x4a*");
  fs4by4astar->SetFillingSchemeName("4x4a*");
  fs4by4astar->SetBXIDs("B"," #2000 #2891");
  fs4by4astar->SetBXIDs("A"," #218 #3337" ); 
  fs4by4astar->SetBXIDs("AC",""); 
  fs4by4astar->SetBXIDs("C"," #1106 #1552");
  fs4by4astar->SetBXIDs("E"," #790");	   
  oadbContFillingScheme->AppendObject(fs4by4astar, 104315, 104321);

  AliOADBFillingScheme * fs4by4b = new AliOADBFillingScheme("4x4b");
  fs4by4b->SetFillingSchemeName("4x4b");
  fs4by4b->SetBXIDs("B"," #2228 #3119");
  fs4by4b->SetBXIDs("A"," #2554 #446" ); 
  fs4by4b->SetBXIDs("AC",""); 
  fs4by4b->SetBXIDs("C"," #1334 #769");
  fs4by4b->SetBXIDs("E"," #790");	   
  oadbContFillingScheme->AppendObject(fs4by4b, 104792, 104803);

  AliOADBFillingScheme * fs4by4c = new AliOADBFillingScheme("4x4c");
  fs4by4c->SetFillingSchemeName("4x4c");
  fs4by4c->SetBXIDs("B"," #3119 #769");
  fs4by4c->SetBXIDs("A"," #2554 #446" ); 
  fs4by4c->SetBXIDs("AC"," "); 
  fs4by4c->SetBXIDs("C"," #1334 #2228");
  fs4by4c->SetBXIDs("E"," #790");	   
  oadbContFillingScheme->AppendObject(fs4by4c, 104824, 104892);

  AliOADBFillingScheme * fs16by16a = new AliOADBFillingScheme("16x16a");
  fs16by16a->SetFillingSchemeName("16x16a");
  fs16by16a->SetBXIDs("B"," #1337 #1418 #2228 #2309 #3119 #3200 #446 #527");
  fs16by16a->SetBXIDs("A"," #1580  #1742  #1904  #2066  #2630  #2792  #2954  #3362" ); 
  fs16by16a->SetBXIDs("AC",""); 
  fs16by16a->SetBXIDs("C"," #845  #1007  #1169   #1577 #3359 #3521 #119  #281");
  fs16by16a->SetBXIDs("E"," #790");	   
  oadbContFillingScheme->AppendObject(fs16by16a, 105143, 105160);

  oadbContFillingScheme->AppendObject(fs4by4c->Clone(), 105256, 105268);

  AliOADBFillingScheme * fsSingle_2b_1_1_1 = new AliOADBFillingScheme("Single_2b_1_1_1");
  fsSingle_2b_1_1_1->SetFillingSchemeName("Single_2b_1_1_1");
  fsSingle_2b_1_1_1->SetBXIDs("B"," #346");
  fsSingle_2b_1_1_1->SetBXIDs("A"," #2131" ); 
  fsSingle_2b_1_1_1->SetBXIDs("AC",""); 
  fsSingle_2b_1_1_1->SetBXIDs("C"," #3019");
  fsSingle_2b_1_1_1->SetBXIDs("E"," #1238");	   
  oadbContFillingScheme->AppendObject(fsSingle_2b_1_1_1, 114786, 116684);

  AliOADBFillingScheme * fsSingle_3b_2_2_2 = new AliOADBFillingScheme("Single_3b_2_2_2");
  fsSingle_3b_2_2_2->SetFillingSchemeName("Single_3b_2_2_2");
  fsSingle_3b_2_2_2->SetBXIDs("B","   #346  #1240");
  fsSingle_3b_2_2_2->SetBXIDs("A","  #2131" ); 
  fsSingle_3b_2_2_2->SetBXIDs("AC",""); 
  fsSingle_3b_2_2_2->SetBXIDs("C"," #3019 ");
  fsSingle_3b_2_2_2->SetBXIDs("E"," #1238");	   
  oadbContFillingScheme->AppendObject(fsSingle_3b_2_2_2, 117048, 117120);


  oadbContFillingScheme->AppendObject(fsSingle_2b_1_1_1, 117220, 118555);

  AliOADBFillingScheme * fsSingle_2b_1_1_1_mis = new AliOADBFillingScheme("Single_2b_1_1_1 - 1 misaligned bx");
  fsSingle_2b_1_1_1_mis->SetFillingSchemeName("Single_2b_1_1_1 - 1 misaligned BX");
  fsSingle_2b_1_1_1_mis->SetBXIDs("B"," #345");
  fsSingle_2b_1_1_1_mis->SetBXIDs("A"," #2130" ); 
  fsSingle_2b_1_1_1_mis->SetBXIDs("AC",""); 
  fsSingle_2b_1_1_1_mis->SetBXIDs("C"," #3018");
  fsSingle_2b_1_1_1_mis->SetBXIDs("E"," #1238");	   
  oadbContFillingScheme->AppendObject(fsSingle_2b_1_1_1_mis, 118556, 118783);

  oadbContFillingScheme->AppendObject(fsSingle_2b_1_1_1->Clone(), 118784, 119163);

  AliOADBFillingScheme * fsSingle_4b_2_2_2 = new AliOADBFillingScheme("Single_4b_2_2_2");
  fsSingle_4b_2_2_2->SetFillingSchemeName("Single_4b_2_2_2");
  fsSingle_4b_2_2_2->SetBXIDs("B","   #669  #3019 ");
  fsSingle_4b_2_2_2->SetBXIDs("A","  #346  #2454 " ); 
  fsSingle_4b_2_2_2->SetBXIDs("AC",""); 
  fsSingle_4b_2_2_2->SetBXIDs("C","  #1234  #2128 ");
  fsSingle_4b_2_2_2->SetBXIDs("E"," #1681 #3463");	   
  oadbContFillingScheme->AppendObject(fsSingle_4b_2_2_2, 119837, 119862);

  AliOADBFillingScheme * fsSingle_6b_3_3_3 = new AliOADBFillingScheme("Single_6b_3_3_3");
  fsSingle_6b_3_3_3->SetFillingSchemeName("Single_6b_3_3_3");
  fsSingle_6b_3_3_3->SetBXIDs("B","   #346  #546  #746 ");
  fsSingle_6b_3_3_3->SetBXIDs("A","  #2131  #2331  #2531 " ); 
  fsSingle_6b_3_3_3->SetBXIDs("AC",""); 
  fsSingle_6b_3_3_3->SetBXIDs("C"," #3019  #3219  #3419");
  fsSingle_6b_3_3_3->SetBXIDs("E"," #1296 #1670");	   
  oadbContFillingScheme->AppendObject(fsSingle_6b_3_3_3, 119902, 120691);

  AliOADBFillingScheme * fsSingle_13b_8_8_8 = new AliOADBFillingScheme("Single_13b_8_8_8");
  fsSingle_13b_8_8_8->SetFillingSchemeName("Single_13b_8_8_8");
  fsSingle_13b_8_8_8->SetBXIDs("B","  #346  #446  #546  #646  #1240  #1340  #1440  #1540");
  fsSingle_13b_8_8_8->SetBXIDs("A","  #946  #2131  #2231  #2331  #2431 " ); 
  fsSingle_13b_8_8_8->SetBXIDs("AC",""); 
  fsSingle_13b_8_8_8->SetBXIDs("C"," #3019  #3119  #3219  #3319  #3519 ");
  fsSingle_13b_8_8_8->SetBXIDs("E"," #1835 #2726");	   
  oadbContFillingScheme->AppendObject(fsSingle_13b_8_8_8, 120741, 122375);

  AliOADBFillingScheme * fs125n_48b_36_16_36 = new AliOADBFillingScheme("125n_48b_36_16_36");
  fs125n_48b_36_16_36->SetFillingSchemeName("125n_48b_36_16_36");
  fs125n_48b_36_16_36->SetBXIDs("B","   #346  #396  #446  #496  #546  #596  #646  #696  #1240  #1290  #1340  #1390  #1440  #1490  #1540  #1590 ");
  fs125n_48b_36_16_36->SetBXIDs("A"," #755  #805  #855  #905  #955  #1005  #1799  #1849  #1899  #2131  #2181  #2231  #2281  #2331  #2381  #2431  #2481  #2531  #2581  #2631  #2846  #3016  #3066  #3116  #3166  #3216  #3266  #3316  #3366  #3425  #3475  #3525 " ); 
  fs125n_48b_36_16_36->SetBXIDs("AC",""); 
  fs125n_48b_36_16_36->SetBXIDs("C","  #3019  #3069  #3119  #3169  #3219  #3269  #3319  #3369  #14  #64  #114  #746  #796  #846  #908  #958  #1008  #1640  #1690  #1740  #2055  #2125  #2175  #2225  #2275  #2325  #2375  #2425  #2475  #2534  #2584  #2634 ");
  fs125n_48b_36_16_36->SetBXIDs("E","");	   
  oadbContFillingScheme->AppendObject(fs125n_48b_36_16_36, 130148, 130375);

  AliOADBFillingScheme * fs1000ns_50b_35_14_35 = new AliOADBFillingScheme("1000ns_50b_35_14_35");
  fs1000ns_50b_35_14_35->SetFillingSchemeName("1000ns_50b_35_14_35");
  fs1000ns_50b_35_14_35->SetBXIDs("B","  #346  #386  #426  #466  #506  #546  #586  #1240  #1280  #1320  #1360  #1400  #1440  #1480 ");
  fs1000ns_50b_35_14_35->SetBXIDs("A","  #626  #666  #706  #746  #786  #826  #866  #1520  #1560  #1600  #1640  #1680  #1720  #1760  #2076  #2131  #2171  #2211  #2251  #2291  #2331  #2371  #2414  #2454  #2494  #2534  #2574  #2614  #2654  #2694  #2734  #2774  #2814 " ); //#2854  #2894  #2934 not present in this run
  fs1000ns_50b_35_14_35->SetBXIDs("AC",""); 
  fs1000ns_50b_35_14_35->SetBXIDs("C"," #3019  #3059  #3099  #3139  #3179  #3219  #3259  #3299  #3339  #3379  #3419  #3459  #3499  #3539  #115  #629  #669  #709  #749  #789  #829  #869  #909  #949  #989  #1029  #1069  #1109  #1149  #1523  #1563  #1603  #1643 "); //#1683  #1723  #1763 not present in this run
  fs1000ns_50b_35_14_35->SetBXIDs("E","");	   
  oadbContFillingScheme->AppendObject(fs1000ns_50b_35_14_35, 130601, 130640);



  // 2010 ----- Heavy Ion -------

  // LHC10h1
  AliOADBPhysicsSelection * oadbLHC10h1 = new AliOADBPhysicsSelection("oadbLHC10h1");
  oadbLHC10h1->AddCollisionTriggerClass   ( AliVEvent::kMB,"+CMBAC-B-NOPF-ALL","B",0);
  oadbLHC10h1->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBAC-A-NOPF-ALL","A",0);
  oadbLHC10h1->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBAC-C-NOPF-ALL","C",0);  
  oadbLHC10h1->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBAC-E-NOPF-ALL","E",0);  
  oadbLHC10h1->SetHardwareTrigger         ( 0,"V0A && V0C");
  oadbLHC10h1->SetOfflineTrigger          ( 0,"(V0A && V0C && SPDGFOL1 > 1) && !V0ABG && !V0CBG  && !TPCLaserWarmUp && ZDCTime");
  
  oadbLHC10h1->AddCollisionTriggerClass   ( AliVEvent::kMB,"+CMBS2A-B-NOPF-ALL","B",1);
  oadbLHC10h1->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBS2A-A-NOPF-ALL","A",1);
  oadbLHC10h1->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBS2A-C-NOPF-ALL","C",1);  
  oadbLHC10h1->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBS2A-E-NOPF-ALL","E",1);  
  oadbLHC10h1->SetHardwareTrigger         ( 1,"SPDGFOL1 > 1 && V0A");
  oadbLHC10h1->SetOfflineTrigger          ( 1,"(V0A && V0C && SPDGFOL1 > 1) && !V0ABG && !V0CBG  && !TPCLaserWarmUp && ZDCTime");

  oadbLHC10h1->AddCollisionTriggerClass   ( AliVEvent::kMB,"+CMBS2C-B-NOPF-ALL","B",2);
  oadbLHC10h1->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBS2C-A-NOPF-ALL","A",2);
  oadbLHC10h1->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBS2C-C-NOPF-ALL","C",2);  
  oadbLHC10h1->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBS2C-E-NOPF-ALL","E",2);  
  oadbLHC10h1->SetHardwareTrigger         ( 2,"SPDGFOL1 > 1 && V0C");
  oadbLHC10h1->SetOfflineTrigger          ( 2,"(V0A && V0C && SPDGFOL1 > 1) && !V0ABG && !V0CBG  && !TPCLaserWarmUp && ZDCTime");
    
  oadbContPS->AppendObject(oadbLHC10h1, 136851,136879);

  // LHC10h2
  AliOADBPhysicsSelection * oadbLHC10h2 = new AliOADBPhysicsSelection("oadbLHC10h2");
  oadbLHC10h2->AddCollisionTriggerClass   ( AliVEvent::kHighMult,"+C0SMH-B-NOPF-ALL","B",0);
  oadbLHC10h2->AddBGTriggerClass          ( AliVEvent::kHighMult,"+C0SMH-A-NOPF-ALL","A",0);
  oadbLHC10h2->AddBGTriggerClass          ( AliVEvent::kHighMult,"+C0SMH-C-NOPF-ALL","C",0);
  oadbLHC10h2->AddBGTriggerClass          ( AliVEvent::kHighMult,"+C0SMH-E-NOPF-ALL","E",0);  
  oadbLHC10h2->SetHardwareTrigger         ( 0,"SPDGFO >= 100");
  oadbLHC10h2->SetOfflineTrigger          ( 0,"SPDGFO >= 100 && !V0ABG && !V0CBG  && !TPCLaserWarmUp && ZDCTime");
  
  oadbLHC10h2->AddCollisionTriggerClass   ( AliVEvent::kMB,"+CMBAC-B-NOPF-ALL","B",1);
  oadbLHC10h2->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBAC-A-NOPF-ALL","A",1);
  oadbLHC10h2->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBAC-C-NOPF-ALL","C",1);
  oadbLHC10h2->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBAC-E-NOPF-ALL","E",1);
  oadbLHC10h2->SetHardwareTrigger         ( 1,"V0A && V0C");
  oadbLHC10h2->SetOfflineTrigger          ( 1,"(V0A && V0C && SPDGFOL1 > 1) && !V0ABG && !V0CBG  && !TPCLaserWarmUp && ZDCTime");
 
  oadbContPS->AppendObject(oadbLHC10h2, 137042,137133);

  // LHC10h3
  AliOADBPhysicsSelection * oadbLHC10h3 = new AliOADBPhysicsSelection("oadbLHC10h3");
  oadbLHC10h3->AddCollisionTriggerClass   ( AliVEvent::kHighMult,"+C0SMH-B-NOPF-ALL","B",0);
  oadbLHC10h3->AddBGTriggerClass          ( AliVEvent::kHighMult,"+C0SMH-A-NOPF-ALL","A",0);
  oadbLHC10h3->AddBGTriggerClass          ( AliVEvent::kHighMult,"+C0SMH-C-NOPF-ALL","C",0);
  oadbLHC10h3->AddBGTriggerClass          ( AliVEvent::kHighMult,"+C0SMH-E-NOPF-ALL","E",0);  
  oadbLHC10h3->SetHardwareTrigger         ( 0,"SPDGFO >= 100");
  oadbLHC10h3->SetOfflineTrigger          ( 0,"SPDGFO >= 100 && !V0ABG && !V0CBG  && !TPCLaserWarmUp && ZDCTime");
  
  oadbLHC10h3->AddCollisionTriggerClass   ( AliVEvent::kMB,"+CMBAC-B-NOPF-ALL","B",1);
  oadbLHC10h3->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBAC-A-NOPF-ALL","A",1);
  oadbLHC10h3->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBAC-C-NOPF-ALL","C",1);
  oadbLHC10h3->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBAC-E-NOPF-ALL","E",1);
  oadbLHC10h3->SetHardwareTrigger         ( 1,"V0A && V0C");
  oadbLHC10h3->SetOfflineTrigger          ( 1,"(V0A && V0C && SPDGFOL1 > 1) && !V0ABG && !V0CBG  && !TPCLaserWarmUp && ZDCTime");
    
  oadbLHC10h3->AddCollisionTriggerClass   ( AliVEvent::kMB,"+CMBS2A-B-NOPF-ALL","B",2);
  oadbLHC10h3->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBS2A-A-NOPF-ALL","A",2);
  oadbLHC10h3->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBS2A-C-NOPF-ALL","C",2);
  oadbLHC10h3->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBS2A-E-NOPF-ALL","E",2);
  oadbLHC10h3->SetHardwareTrigger         ( 2,"SPDGFOL1 > 1 && V0A");
  oadbLHC10h3->SetOfflineTrigger          ( 2,"(V0A && V0C && SPDGFOL1 > 1) && !V0ABG && !V0CBG  && !TPCLaserWarmUp && ZDCTime");

  oadbLHC10h3->AddCollisionTriggerClass   ( AliVEvent::kMB,"+CMBS2C-B-NOPF-ALL","B",3);
  oadbLHC10h3->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBS2C-A-NOPF-ALL","A",3);
  oadbLHC10h3->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBS2C-C-NOPF-ALL","C",3);
  oadbLHC10h3->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBS2C-E-NOPF-ALL","E",3);
  oadbLHC10h3->SetHardwareTrigger         ( 3,"SPDGFOL1 > 1 && V0C");
  oadbLHC10h3->SetOfflineTrigger          ( 3,"(V0A && V0C && SPDGFOL1 > 1) && !V0ABG && !V0CBG  && !TPCLaserWarmUp && ZDCTime");

  oadbContPS->AppendObject(oadbLHC10h3, 137135,137364);

 // LHC10h, run 137365
  AliOADBPhysicsSelection * oadbLHC10h_137365 = new AliOADBPhysicsSelection("oadbLHC10h_x9");
  oadbLHC10h_137365->AddCollisionTriggerClass   ( AliVEvent::kMB,"+CMBAC-B-NOPF-ALL","B",0);
  oadbLHC10h_137365->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBAC-A-NOPF-ALL","A",0);
  oadbLHC10h_137365->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBAC-C-NOPF-ALL","C",0);
  oadbLHC10h_137365->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBAC-E-NOPF-ALL","E",0);
  oadbLHC10h_137365->SetHardwareTrigger         ( 0,"V0A && V0C");
  oadbLHC10h_137365->SetOfflineTrigger          ( 0,"(V0A && V0C && SPDGFOL1 > 1) && !V0ABG && !V0CBG  && !TPCLaserWarmUp && ZDCTime");

  oadbContPS->AppendObject(oadbLHC10h_137365, 137365,137365);

  AliOADBPhysicsSelection * oadbLHC10h4 = new AliOADBPhysicsSelection("oadbLHC10h4");
  oadbLHC10h4->AddCollisionTriggerClass   ( AliVEvent::kHighMult,"+C0SMH-B-NOPF-ALL","B",0);
  oadbLHC10h4->AddBGTriggerClass          ( AliVEvent::kHighMult,"+C0SMH-A-NOPF-ALL","A",0);
  oadbLHC10h4->AddBGTriggerClass          ( AliVEvent::kHighMult,"+C0SMH-C-NOPF-ALL","C",0);
  oadbLHC10h4->AddBGTriggerClass          ( AliVEvent::kHighMult,"+C0SMH-E-NOPF-ALL","E",0);  
  oadbLHC10h4->SetHardwareTrigger         ( 0,"SPDGFO >= 100");
  oadbLHC10h4->SetOfflineTrigger          ( 0,"SPDGFO >= 100 && !V0ABG && !V0CBG  && !TPCLaserWarmUp && ZDCTime");
  
  oadbLHC10h4->AddCollisionTriggerClass   ( AliVEvent::kMB,"+CMBAC-B-NOPF-ALL","B",1);
  oadbLHC10h4->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBAC-A-NOPF-ALL","A",1);
  oadbLHC10h4->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBAC-C-NOPF-ALL","C",1);
  oadbLHC10h4->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBAC-E-NOPF-ALL","E",1);
  oadbLHC10h4->SetHardwareTrigger         ( 1,"V0A && V0C");
  oadbLHC10h4->SetOfflineTrigger          ( 1,"(V0A && V0C && SPDGFOL1 > 1) && !V0ABG && !V0CBG  && !TPCLaserWarmUp && ZDCTime");
    
  oadbLHC10h4->AddCollisionTriggerClass   ( AliVEvent::kMB,"+CMBS2A-B-NOPF-ALL","B",2);
  oadbLHC10h4->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBS2A-A-NOPF-ALL","A",2);
  oadbLHC10h4->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBS2A-C-NOPF-ALL","C",2);
  oadbLHC10h4->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBS2A-E-NOPF-ALL","E",2);
  oadbLHC10h4->SetHardwareTrigger         ( 2,"SPDGFOL1 > 1 && V0A");
  oadbLHC10h4->SetOfflineTrigger          ( 2,"(V0A && V0C && SPDGFOL1 > 1) && !V0ABG && !V0CBG  && !TPCLaserWarmUp && ZDCTime");

  oadbLHC10h4->AddCollisionTriggerClass   ( AliVEvent::kMB,"+CMBS2C-B-NOPF-ALL","B",3);
  oadbLHC10h4->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBS2C-A-NOPF-ALL","A",3);
  oadbLHC10h4->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBS2C-C-NOPF-ALL","C",3);
  oadbLHC10h4->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBS2C-E-NOPF-ALL","E",3);
  oadbLHC10h4->SetHardwareTrigger         ( 3,"SPDGFOL1 > 1 && V0C");
  oadbLHC10h4->SetOfflineTrigger          ( 3,"(V0A && V0C && SPDGFOL1 > 1) && !V0ABG && !V0CBG  && !TPCLaserWarmUp && ZDCTime");

  oadbContPS->AppendObject(oadbLHC10h4, 137366,137595);



  // LHC10h5
  AliOADBPhysicsSelection * oadbLHC10h5 = new AliOADBPhysicsSelection("oadbLHC10h5");
  oadbLHC10h5->AddCollisionTriggerClass   ( AliVEvent::kHighMult,"+C0SMH-B-NOPF-ALL","B",0);
  oadbLHC10h5->AddBGTriggerClass          ( AliVEvent::kHighMult,"+C0SMH-A-NOPF-ALL","A",0);
  oadbLHC10h5->AddBGTriggerClass          ( AliVEvent::kHighMult,"+C0SMH-C-NOPF-ALL","C",0);
  oadbLHC10h5->AddBGTriggerClass          ( AliVEvent::kHighMult,"+C0SMH-E-NOPF-ALL","E",0);  
  oadbLHC10h5->SetHardwareTrigger         ( 0,"SPDGFO >= 100");
  oadbLHC10h5->SetOfflineTrigger          ( 0,"SPDGFO >= 100 && !V0ABG && !V0CBG  && !TPCLaserWarmUp && ZDCTime");
  
  oadbLHC10h5->AddCollisionTriggerClass   ( AliVEvent::kMB,"+CMBAC-B-NOPF-ALL","B",1);
  oadbLHC10h5->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBAC-A-NOPF-ALL","A",1);
  oadbLHC10h5->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBAC-C-NOPF-ALL","C",1);
  oadbLHC10h5->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBAC-E-NOPF-ALL","E",1);
  oadbLHC10h5->SetHardwareTrigger         ( 1,"V0A && V0C");
  oadbLHC10h5->SetOfflineTrigger          ( 1,"(V0A && V0C && SPDGFOL1 > 1) && !V0ABG && !V0CBG  && !TPCLaserWarmUp && ZDCTime");

  oadbContPS->AppendObject(oadbLHC10h5, 137608,137848); // 138982,138983);

  // LHC10h6
  AliOADBPhysicsSelection * oadbLHC10h6 = new AliOADBPhysicsSelection("oadbLHC10h6");
  oadbLHC10h6->AddCollisionTriggerClass   ( AliVEvent::kHighMult,"+C0SMH-B-NOPF-ALL","B",0);
  oadbLHC10h6->AddBGTriggerClass          ( AliVEvent::kHighMult,"+C0SMH-A-NOPF-ALL","A",0);
  oadbLHC10h6->AddBGTriggerClass          ( AliVEvent::kHighMult,"+C0SMH-C-NOPF-ALL","C",0);
  oadbLHC10h6->AddBGTriggerClass          ( AliVEvent::kHighMult,"+C0SMH-E-NOPF-ALL","E",0);  
  oadbLHC10h6->SetHardwareTrigger         ( 0,"SPDGFO >= 100");
  oadbLHC10h6->SetOfflineTrigger          ( 0,"SPDGFO >= 100 && !V0ABG && !V0CBG  && !TPCLaserWarmUp && ZDCTime");
  
  oadbLHC10h6->AddCollisionTriggerClass   ( AliVEvent::kMB,"+CMBACS2-B-NOPF-ALL","B",1);
  oadbLHC10h6->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBACS2-A-NOPF-ALL","A",1);
  oadbLHC10h6->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBACS2-C-NOPF-ALL","C",1);
  oadbLHC10h6->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBACS2-E-NOPF-ALL","E",1);
  oadbLHC10h6->SetHardwareTrigger         ( 1,"(V0A && V0C && SPDGFOL1 > 1)");
  oadbLHC10h6->SetOfflineTrigger          ( 1,"(V0A && V0C && SPDGFOL1 > 1) && !V0ABG && !V0CBG  && !TPCLaserWarmUp && ZDCTime");

  oadbContPS->AppendObject(oadbLHC10h6, 138125,138980);// 139028, 139316);


  AliOADBPhysicsSelection * oadbLHC10h7 = new AliOADBPhysicsSelection("oadbLHC10h7");
  oadbLHC10h7->AddCollisionTriggerClass   ( AliVEvent::kHighMult,"+C0SMH-B-NOPF-ALL","B",0);
  oadbLHC10h7->AddBGTriggerClass          ( AliVEvent::kHighMult,"+C0SMH-A-NOPF-ALL","A",0);
  oadbLHC10h7->AddBGTriggerClass          ( AliVEvent::kHighMult,"+C0SMH-C-NOPF-ALL","C",0);
  oadbLHC10h7->AddBGTriggerClass          ( AliVEvent::kHighMult,"+C0SMH-E-NOPF-ALL","E",0);  
  oadbLHC10h7->SetHardwareTrigger         ( 0,"SPDGFO >= 100");
  oadbLHC10h7->SetOfflineTrigger          ( 0,"SPDGFO >= 100 && !V0ABG && !V0CBG  && !TPCLaserWarmUp && ZDCTime");
  
  oadbLHC10h7->AddCollisionTriggerClass   ( AliVEvent::kMB,"+CMBAC-B-NOPF-ALL","B",1);
  oadbLHC10h7->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBAC-A-NOPF-ALL","A",1);
  oadbLHC10h7->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBAC-C-NOPF-ALL","C",1);
  oadbLHC10h7->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBAC-E-NOPF-ALL","E",1);
  oadbLHC10h7->SetHardwareTrigger         ( 1,"(V0A && V0C)");
  oadbLHC10h7->SetOfflineTrigger          ( 1,"(V0A && V0C && SPDGFOL1 > 1) && !V0ABG && !V0CBG  && !TPCLaserWarmUp && ZDCTime");


  oadbContPS->AppendObject(oadbLHC10h7, 138982,138983);


  // LHC10h8
  AliOADBPhysicsSelection * oadbLHC10h8 = new AliOADBPhysicsSelection("oadbLHC10h8");
  oadbLHC10h8->AddCollisionTriggerClass   ( AliVEvent::kHighMult,"+C0SMH-B-NOPF-ALL","B",0);
  oadbLHC10h8->AddBGTriggerClass          ( AliVEvent::kHighMult,"+C0SMH-A-NOPF-ALL","A",0);
  oadbLHC10h8->AddBGTriggerClass          ( AliVEvent::kHighMult,"+C0SMH-C-NOPF-ALL","C",0);
  oadbLHC10h8->AddBGTriggerClass          ( AliVEvent::kHighMult,"+C0SMH-E-NOPF-ALL","E",0);  
  oadbLHC10h8->SetHardwareTrigger         ( 0,"SPDGFO >= 100");
  oadbLHC10h8->SetOfflineTrigger          ( 0,"SPDGFO >= 100 && !V0ABG && !V0CBG  && !TPCLaserWarmUp && ZDCTime");
  
  oadbLHC10h8->AddCollisionTriggerClass   ( AliVEvent::kMB,"+CMBACS2-B-NOPF-ALL","B",1);
  oadbLHC10h8->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBACS2-A-NOPF-ALL","A",1);
  oadbLHC10h8->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBACS2-C-NOPF-ALL","C",1);
  oadbLHC10h8->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBACS2-E-NOPF-ALL","E",1);
  oadbLHC10h8->SetHardwareTrigger         ( 1,"(V0A && V0C && SPDGFOL1 > 1)");
  oadbLHC10h8->SetOfflineTrigger          ( 1,"(V0A && V0C && SPDGFOL1 > 1) && !V0ABG && !V0CBG  && !TPCLaserWarmUp && ZDCTime");

  oadbContPS->AppendObject(oadbLHC10h8,  139028, 139316);


  // LHC10h9
  AliOADBPhysicsSelection * oadbLHC10h9 = new AliOADBPhysicsSelection("oadbLHC10h9");
  oadbLHC10h9->AddCollisionTriggerClass   ( AliVEvent::kHighMult,"+C0SMH-B-NOPF-ALLNOTRD","B",0);
  oadbLHC10h9->AddBGTriggerClass          ( AliVEvent::kHighMult,"+C0SMH-A-NOPF-ALLNOTRD","A",0);
  oadbLHC10h9->AddBGTriggerClass          ( AliVEvent::kHighMult,"+C0SMH-C-NOPF-ALLNOTRD","C",0);
  oadbLHC10h9->AddBGTriggerClass          ( AliVEvent::kHighMult,"+C0SMH-E-NOPF-ALLNOTRD","E",0);  
  oadbLHC10h9->SetHardwareTrigger         ( 0,"SPDGFO >= 100");
  oadbLHC10h9->SetOfflineTrigger          ( 0,"SPDGFO >= 100 && !V0ABG && !V0CBG  && !TPCLaserWarmUp && ZDCTime");
  
  oadbLHC10h9->AddCollisionTriggerClass   ( AliVEvent::kMB,"+CMBACS2-B-NOPF-ALLNOTRD","B",1);
  oadbLHC10h9->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBACS2-A-NOPF-ALLNOTRD","A",1);
  oadbLHC10h9->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBACS2-C-NOPF-ALLNOTRD","C",1);
  oadbLHC10h9->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBACS2-E-NOPF-ALLNOTRD","E",1);
  oadbLHC10h9->SetHardwareTrigger         ( 1,"(V0A && V0C && SPDGFOL1 > 1)");
  oadbLHC10h9->SetOfflineTrigger          ( 1,"(V0A && V0C && SPDGFOL1 > 1) && !V0ABG && !V0CBG  && !TPCLaserWarmUp && ZDCTime");

  oadbContPS->AppendObject(oadbLHC10h9, 139328,139517);


  // ----- 2011 -----
  // ----- proton-proton -----
  
  // Default LHC11a/b/c, cint1b configuration
  AliOADBPhysicsSelection * oadbLHC11ab = new AliOADBPhysicsSelection("oadbLHC11abc");
  oadbLHC11ab->AddCollisionTriggerClass   ( AliVEvent::kMB,"+CINT1-B-NOPF-ALLNOTRD","B",0);
  oadbLHC11ab->AddBGTriggerClass          ( AliVEvent::kMB,"+CINT1-AC-NOPF-ALLNOTRD","AC",0);
  oadbLHC11ab->AddBGTriggerClass          ( AliVEvent::kMB,"+CINT1-E-NOPF-ALLNOTRD","E",0);  
  oadbLHC11ab->AddCollisionTriggerClass   ( AliVEvent::kMB | AliVEvent::kFastOnly,"+CINT1-B-NOPF-FASTNOTRD -CINT1-B-NOPF-ALLNOTRD","B",0);
  oadbLHC11ab->AddBGTriggerClass          ( AliVEvent::kMB | AliVEvent::kFastOnly,"+CINT1-AC-NOPF-FASTNOTRD -CINT1-AC-NOPF-ALLNOTRD","AC",0);
  oadbLHC11ab->AddBGTriggerClass          ( AliVEvent::kMB | AliVEvent::kFastOnly,"+CINT1-E-NOPF-FASTNOTRD -CINT1-E-NOPF-ALLNOTRD","E",0);  
  oadbLHC11ab->SetHardwareTrigger         ( 0,"SPDGFO >= 1 || V0A || V0C");					      
  oadbLHC11ab->SetOfflineTrigger          ( 0,"(SPDGFO >= 1 || V0A || V0C) && !V0ABG && !V0CBG  && !TPCLaserWarmUp");

  oadbLHC11ab->AddCollisionTriggerClass   ( AliVEvent::kMUON,"+CMUS1-B-NOPF-MUON","B",1);
  oadbLHC11ab->AddBGTriggerClass          ( AliVEvent::kMUON,"+CMUS1-AC-NOPF-MUON","AC",1);
  oadbLHC11ab->AddBGTriggerClass          ( AliVEvent::kMUON,"+CMUS1-E-NOPF-MUON","E",1);
  oadbLHC11ab->SetHardwareTrigger         ( 1,"SPDGFO >= 1 || V0A || V0C");                                         
  oadbLHC11ab->SetOfflineTrigger          ( 1,"(SPDGFO >= 1 || V0A || V0C) && !V0ABG && !V0CBG  && !TPCLaserWarmUp");
    
  oadbLHC11ab->AddCollisionTriggerClass   ( AliVEvent::kHighMult,"+CSH1-B-NOPF-ALLNOTRD","B",2);
  oadbLHC11ab->AddBGTriggerClass          ( AliVEvent::kHighMult,"+CSH1-AC-NOPF-ALLNOTRD","AC",2);
  oadbLHC11ab->AddBGTriggerClass          ( AliVEvent::kHighMult,"+CSH1-E-NOPF-ALLNOTRD","E",2);
  oadbLHC11ab->SetHardwareTrigger         ( 2,"SPDGFO >= 1 || V0A || V0C");                      
  oadbLHC11ab->SetOfflineTrigger          ( 2,"(SPDGFO >= 1 || V0A || V0C) && !V0ABG && !V0CBG  && !TPCLaserWarmUp");

  oadbLHC11ab->AddCollisionTriggerClass   ( AliVEvent::kEMC1,"+CEMC1-B-NOPF-ALLNOTRD","B",2);
  oadbLHC11ab->AddBGTriggerClass          ( AliVEvent::kEMC1,"+CEMC1-AC-NOPF-ALLNOTRD","AC",2);
  oadbLHC11ab->AddBGTriggerClass          ( AliVEvent::kEMC1,"+CEMC1-E-NOPF-ALLNOTRD","E",2);
  // already defined before
  // oadbLHC11ab->SetHardwareTrigger         ( 2,"SPDGFO >= 1 || V0A || V0C");                      
  // oadbLHC11ab->SetOfflineTrigger          ( 2,"(SPDGFO >= 1 || V0A || V0C) && !V0ABG && !V0CBG  && !TPCLaserWarmUp");

  oadbLHC11ab->AddCollisionTriggerClass   ( AliVEvent::kPHI1,"+CPHI1-B-NOPF-ALLNOTRD","B",3);
  oadbLHC11ab->AddBGTriggerClass          ( AliVEvent::kPHI1,"+CPHI1-AC-NOPF-ALLNOTRD","AC",3);
  oadbLHC11ab->AddBGTriggerClass          ( AliVEvent::kPHI1,"+CPHI1-E-NOPF-ALLNOTRD","E",3);
  oadbLHC11ab->SetHardwareTrigger         ( 3,"SPDGFO >= 1 || V0A || V0C");
  oadbLHC11ab->SetOfflineTrigger          ( 3,"(SPDGFO >= 1 || V0A || V0C) && !V0ABG && !V0CBG  && !TPCLaserWarmUp");

  oadbContPS->AppendObject(oadbLHC11ab,  144871, 146856);
  oadbContPS->AppendObject(oadbLHC11ab->Clone(),  146858, 152935);
  oadbContPS->AppendObject(oadbLHC11ab->Clone(),  153360, 153361);
  
  // LHC11a1, "ALL" exceptions, cint1b configuration
  AliOADBPhysicsSelection * oadbLHC11a1 = new AliOADBPhysicsSelection("oadbLHC11a1");
  oadbLHC11a1->AddCollisionTriggerClass   ( AliVEvent::kMB,"+CINT1-B-NOPF-ALL","B",0);
  oadbLHC11a1->AddBGTriggerClass          ( AliVEvent::kMB,"+CINT1-AC-NOPF-ALL","AC",0);
  oadbLHC11a1->AddBGTriggerClass          ( AliVEvent::kMB,"+CINT1-E-NOPF-ALL","E",0);  
  oadbLHC11a1->SetHardwareTrigger         ( 0,"SPDGFO >= 1 || V0A || V0C");					      
  oadbLHC11a1->SetOfflineTrigger          ( 0,"(SPDGFO >= 1 || V0A || V0C) && !V0ABG && !V0CBG  && !TPCLaserWarmUp");

  oadbLHC11a1->AddCollisionTriggerClass   ( AliVEvent::kMUON,"+CMUS1-B-NOPF-ALL","B",1);
  oadbLHC11a1->AddBGTriggerClass          ( AliVEvent::kMUON,"+CMUS1-AC-NOPF-ALL","AC",1);
  oadbLHC11a1->AddBGTriggerClass          ( AliVEvent::kMUON,"+CMUS1-E-NOPF-ALL","E",1);
  oadbLHC11a1->SetHardwareTrigger         ( 1,"SPDGFO >= 1 || V0A || V0C");                                         
  oadbLHC11a1->SetOfflineTrigger          ( 1,"(SPDGFO >= 1 || V0A || V0C) && !V0ABG && !V0CBG  && !TPCLaserWarmUp");
  
  oadbContPS->AppendObject(oadbLHC11a1, 146857,146857);
  
  
  //LHC11c1, MUU7 read out with MUON cluster (default then changed to ALLNOTRD)
  AliOADBPhysicsSelection * oadbLHC11c1 = new AliOADBPhysicsSelection("oadbLHC11c1");
  oadbLHC11c1->AddCollisionTriggerClass   ( AliVEvent::kINT7,"+CINT7-B-NOPF-ALLNOTRD","B",0);
  oadbLHC11c1->AddBGTriggerClass          ( AliVEvent::kINT7,"+CINT7-AC-NOPF-ALLNOTRD","AC",0);
  oadbLHC11c1->AddBGTriggerClass          ( AliVEvent::kINT7,"+CINT7-E-NOPF-ALLNOTRD","E",0);
  oadbLHC11c1->SetHardwareTrigger         ( 0, "V0A && V0C");
  oadbLHC11c1->SetOfflineTrigger          ( 0, "(V0A && V0C) && !V0ABG && !V0CBG && !TPCLaserWarmUp");

  oadbLHC11c1->AddCollisionTriggerClass   ( AliVEvent::kMUSH7,"+CMUSH7-B-NOPF-MUON","B",1);
  oadbLHC11c1->AddBGTriggerClass          ( AliVEvent::kMUSH7,"+CMUSH7-AC-NOPF-MUON","AC",1);
  oadbLHC11c1->AddBGTriggerClass          ( AliVEvent::kMUSH7,"+CMUSH7-E-NOPF-MUON","E",1);
  oadbLHC11c1->SetHardwareTrigger         ( 1, "V0A && V0C");
  oadbLHC11c1->SetOfflineTrigger          ( 1, "(V0A && V0C) && !V0ABG && !V0CBG && !TPCLaserWarmUp");

  oadbLHC11c1->AddCollisionTriggerClass   ( AliVEvent::kMUL7,"+CMUL7-B-NOPF-MUON","B",2);
  oadbLHC11c1->AddBGTriggerClass          ( AliVEvent::kMUL7,"+CMUL7-AC-NOPF-MUON","AC",2);
  oadbLHC11c1->AddBGTriggerClass          ( AliVEvent::kMUL7,"+CMUL7-E-NOPF-MUON","E",2);
  oadbLHC11c1->SetHardwareTrigger         ( 2, "V0A && V0C");
  oadbLHC11c1->SetOfflineTrigger          ( 2, "(V0A && V0C) && !V0ABG && !V0CBG && !TPCLaserWarmUp");

  oadbLHC11c1->AddCollisionTriggerClass   ( AliVEvent::kMUU7,"+CMUU7-B-NOPF-MUON","B",3);
  oadbLHC11c1->AddBGTriggerClass          ( AliVEvent::kMUU7,"+CMUU7-AC-NOPF-MUON","AC",3);
  oadbLHC11c1->AddBGTriggerClass          ( AliVEvent::kMUU7,"+CMUU7-E-NOPF-MUON","E",3);
  oadbLHC11c1->SetHardwareTrigger         ( 3, "V0A && V0C");
  oadbLHC11c1->SetOfflineTrigger          ( 3, "(V0A && V0C) && !V0ABG && !V0CBG && !TPCLaserWarmUp");

  oadbLHC11c1->AddCollisionTriggerClass   ( AliVEvent::kEMC7,"+CEMC7-B-NOPF-ALLNOTRD","B",4);
  oadbLHC11c1->AddBGTriggerClass          ( AliVEvent::kEMC7,"+CEMC7-AC-NOPF-ALLNOTRD","AC",4);
  oadbLHC11c1->AddBGTriggerClass          ( AliVEvent::kEMC7,"+CEMC7-E-NOPF-ALLNOTRD","E",4);
  oadbLHC11c1->SetHardwareTrigger         ( 4, "V0A && V0C");
  oadbLHC11c1->SetOfflineTrigger          ( 4, "(V0A && V0C) && !V0ABG && !V0CBG && !TPCLaserWarmUp");

  oadbLHC11c1->AddCollisionTriggerClass   ( AliVEvent::kPHI7,"+CPHI7-B-NOPF-ALLNOTRD","B",  5);
  oadbLHC11c1->AddBGTriggerClass          ( AliVEvent::kPHI7,"+CPHI7-AC-NOPF-ALLNOTRD","AC",5);
  oadbLHC11c1->AddBGTriggerClass          ( AliVEvent::kPHI7,"+CPHI7-E-NOPF-ALLNOTRD","E",  5);
  oadbLHC11c1->SetHardwareTrigger         ( 5,"V0A && V0C");
  oadbLHC11c1->SetOfflineTrigger          ( 5,"(V0A && V0C) && !V0ABG && !V0CBG  && !TPCLaserWarmUp");

  oadbContPS->AppendObject(oadbLHC11c1, 153056 , 153296);
  oadbContPS->AppendObject(oadbLHC11c1->Clone(), 153362 , 153578);
  oadbContPS->AppendObject(oadbLHC11c1->Clone(), 153583 , 153583);
  oadbContPS->AppendObject(oadbLHC11c1->Clone(), 153587 , 153733);
  
  // Tests with CMUS7
  //LHC11c3
//  AliOADBPhysicsSelection * oadbLHC11c3 = new AliOADBPhysicsSelection("oadbLHC11c3");
//  oadbLHC11c3->AddCollisionTriggerClass   ( AliVEvent::kINT7,"+CINT7-B-NOPF-ALLNOTRD","B",0);
//  oadbLHC11c3->AddBGTriggerClass          ( AliVEvent::kINT7,"+CINT7-AC-NOPF-ALLNOTRD","AC",0);
//  oadbLHC11c3->AddBGTriggerClass          ( AliVEvent::kINT7,"+CINT7-E-NOPF-ALLNOTRD","E",0);
//  oadbLHC11c3->SetHardwareTrigger         ( 0, "V0A && V0C");
//  oadbLHC11c3->SetOfflineTrigger          ( 0, "(V0A && V0C) && !V0ABG && !V0CBG && !TPCLaserWarmUp");
//  
//  oadbLHC11c3->AddCollisionTriggerClass   ( AliVEvent::kMUS7,"+CMUS7-B-NOPF-MUON","B",1);
//  oadbLHC11c3->AddBGTriggerClass          ( AliVEvent::kMUS7,"+CMUS7-AC-NOPF-MUON","AC",1);
//  oadbLHC11c3->AddBGTriggerClass          ( AliVEvent::kMUS7,"+CMUS7-E-NOPF-MUON","E",1);
//  oadbLHC11c3->SetHardwareTrigger         ( 1, "V0A && V0C");
//  oadbLHC11c3->SetOfflineTrigger          ( 1, "(V0A && V0C) && !V0ABG && !V0CBG && !TPCLaserWarmUp");
//  
//  oadbContPS->AppendObject(oadbLHC11c3, 151664 , 151680);
//  
//  //LHC11c4
//  AliOADBPhysicsSelection * oadbLHC11c4 = new AliOADBPhysicsSelection("oadbLHC11c4");
//  oadbLHC11c4->AddCollisionTriggerClass   ( AliVEvent::kINT7,"+CINT7-B-NOPF-ALLNOTRD","B",0);
//  oadbLHC11c4->AddBGTriggerClass          ( AliVEvent::kINT7,"+CINT7-AC-NOPF-ALLNOTRD","AC",0);
//  oadbLHC11c4->AddBGTriggerClass          ( AliVEvent::kINT7,"+CINT7-E-NOPF-ALLNOTRD","E",0);
//  oadbLHC11c4->SetHardwareTrigger         ( 0, "V0A && V0C");
//  oadbLHC11c4->SetOfflineTrigger          ( 0, "(V0A && V0C) && !V0ABG && !V0CBG && !TPCLaserWarmUp");
//  
//  oadbLHC11c4->AddCollisionTriggerClass   ( AliVEvent::kMUS7,"+CMUS7-B-NOPF-MUON","B",1);
//  oadbLHC11c4->AddBGTriggerClass          ( AliVEvent::kMUS7,"+CMUS7-AC-NOPF-MUON","AC",1);
//  oadbLHC11c4->AddBGTriggerClass          ( AliVEvent::kMUS7,"+CMUS7-E-NOPF-MUON","E",1);
//  oadbLHC11c4->SetHardwareTrigger         ( 1, "V0A && V0C");
//  oadbLHC11c4->SetOfflineTrigger          ( 1, "(V0A && V0C) && !V0ABG && !V0CBG && !TPCLaserWarmUp");
//  
//  oadbLHC11c4->AddCollisionTriggerClass   ( AliVEvent::kEMC7,"+CEMC7-B-NOPF-ALLNOTRD","B",2);
//  oadbLHC11c4->AddBGTriggerClass          ( AliVEvent::kEMC7,"+CEMC7-AC-NOPF-ALLNOTRD","AC",2);
//  oadbLHC11c4->AddBGTriggerClass          ( AliVEvent::kEMC7,"+CEMC7-E-NOPF-ALLNOTRD","E",2);
//  oadbLHC11c4->SetHardwareTrigger         ( 2, "V0A && V0C");
//  oadbLHC11c4->SetOfflineTrigger          ( 2, "(V0A && V0C) && !V0ABG && !V0CBG && !TPCLaserWarmUp");
//  
//  oadbContPS->AppendObject(oadbLHC11c4, 152002 , 152007);
//  oadbContPS->AppendObject(oadbLHC11c4->Clone(), 152011 , 152011);
//  oadbContPS->AppendObject(oadbLHC11c4->Clone(), 152082 , 152082);
//  oadbContPS->AppendObject(oadbLHC11c4->Clone(), 152137 , 152137);
  
  

//  //LHC11c5 (Current default)
//  AliOADBPhysicsSelection * oadbLHC11c5 = new AliOADBPhysicsSelection("oadbLHC11c5");
//  oadbLHC11c5->AddCollisionTriggerClass   ( AliVEvent::kINT7,"+CINT7-B-NOPF-ALLNOTRD","B",0);
//  oadbLHC11c5->AddBGTriggerClass          ( AliVEvent::kINT7,"+CINT7-AC-NOPF-ALLNOTRD","AC",0);
//  oadbLHC11c5->AddBGTriggerClass          ( AliVEvent::kINT7,"+CINT7-E-NOPF-ALLNOTRD","E",0);
//  oadbLHC11c5->SetHardwareTrigger         ( 0, "V0A && V0C");
//  oadbLHC11c5->SetOfflineTrigger          ( 0, "(V0A && V0C) && !V0ABG && !V0CBG && !TPCLaserWarmUp");
//
//  oadbLHC11c5->AddCollisionTriggerClass   ( AliVEvent::kMUSH7,"+CMUSH7-B-NOPF-MUON","B",1);
//  oadbLHC11c5->AddBGTriggerClass          ( AliVEvent::kMUSH7,"+CMUSH7-AC-NOPF-MUON","AC",1);
//  oadbLHC11c5->AddBGTriggerClass          ( AliVEvent::kMUSH7,"+CMUSH7-E-NOPF-MUON","E",1);
//  oadbLHC11c5->SetHardwareTrigger         ( 1, "V0A && V0C");
//  oadbLHC11c5->SetOfflineTrigger          ( 1, "(V0A && V0C) && !V0ABG && !V0CBG && !TPCLaserWarmUp");
//
//  oadbLHC11c5->AddCollisionTriggerClass   ( AliVEvent::kMUL7,"+CMUL7-B-NOPF-MUON","B",2);
//  oadbLHC11c5->AddBGTriggerClass          ( AliVEvent::kMUL7,"+CMUL7-AC-NOPF-MUON","AC",2);
//  oadbLHC11c5->AddBGTriggerClass          ( AliVEvent::kMUL7,"+CMUL7-E-NOPF-MUON","E",2);
//  oadbLHC11c5->SetHardwareTrigger         ( 2, "V0A && V0C");
//  oadbLHC11c5->SetOfflineTrigger          ( 2, "(V0A && V0C) && !V0ABG && !V0CBG && !TPCLaserWarmUp");
//
//  oadbLHC11c5->AddCollisionTriggerClass   ( AliVEvent::kMUU7,"+CMUU7-B-NOPF-ALLNOTRD","B",3);
//  oadbLHC11c5->AddBGTriggerClass          ( AliVEvent::kMUU7,"+CMUU7-AC-NOPF-ALLNOTRD","AC",3);
//  oadbLHC11c5->AddBGTriggerClass          ( AliVEvent::kMUU7,"+CMUU7-E-NOPF-ALLNOTRD","E",3);
//  oadbLHC11c5->SetHardwareTrigger         ( 3, "V0A && V0C");
//  oadbLHC11c5->SetOfflineTrigger          ( 3, "(V0A && V0C) && !V0ABG && !V0CBG && !TPCLaserWarmUp");
//
//  oadbLHC11c5->AddCollisionTriggerClass   ( AliVEvent::kEMC7,"+CEMC7-B-NOPF-ALLNOTRD","B",4);
//  oadbLHC11c5->AddBGTriggerClass          ( AliVEvent::kEMC7,"+CEMC7-AC-NOPF-ALLNOTRD","AC",4);
//  oadbLHC11c5->AddBGTriggerClass          ( AliVEvent::kEMC7,"+CEMC7-E-NOPF-ALLNOTRD","E",4);
//  oadbLHC11c5->SetHardwareTrigger         ( 4, "V0A && V0C");
//  oadbLHC11c5->SetOfflineTrigger          ( 4, "(V0A && V0C) && !V0ABG && !V0CBG && !TPCLaserWarmUp");
//
//  oadbContPS->AppendObject(oadbLHC11c5, 153776 , 154570);

  //LHC11d1, CINT7-B and CPHI-7B (before I transition)
  AliOADBPhysicsSelection * oadbLHC11d1 = new AliOADBPhysicsSelection("oadbLHC11d1");

  oadbLHC11d1->AddCollisionTriggerClass   ( AliVEvent::kINT7,"+CINT7-B-NOPF-ALLNOTRD","B",0);
  oadbLHC11d1->AddBGTriggerClass          ( AliVEvent::kINT7,"+CINT7-AC-NOPF-ALLNOTRD","AC",0);
  oadbLHC11d1->AddBGTriggerClass          ( AliVEvent::kINT7,"+CINT7-E-NOPF-ALLNOTRD","E",0);  
  oadbLHC11d1->SetHardwareTrigger         ( 0,"V0A && V0C");					      
  oadbLHC11d1->SetOfflineTrigger          ( 0,"(V0A && V0C) && !V0ABG && !V0CBG  && !TPCLaserWarmUp");

  oadbLHC11d1->AddCollisionTriggerClass   ( AliVEvent::kMUSH7,"+CMUSH7-B-NOPF-MUON","B",  1);
  oadbLHC11d1->AddBGTriggerClass          ( AliVEvent::kMUSH7,"+CMUSH7-AC-NOPF-MUON","AC",1);
  oadbLHC11d1->AddBGTriggerClass          ( AliVEvent::kMUSH7,"+CMUSH7-E-NOPF-MUON","E",  1);
  oadbLHC11d1->SetHardwareTrigger         ( 1,"V0A && V0C");                                         
  oadbLHC11d1->SetOfflineTrigger          ( 1,"(V0A && V0C) && !V0ABG && !V0CBG  && !TPCLaserWarmUp");
    
  oadbLHC11d1->AddCollisionTriggerClass   ( AliVEvent::kMUL7,"+CMUL7-B-NOPF-MUON","B",  2);
  oadbLHC11d1->AddBGTriggerClass          ( AliVEvent::kMUL7,"+CMUL7-AC-NOPF-MUON","AC",2);
  oadbLHC11d1->AddBGTriggerClass          ( AliVEvent::kMUL7,"+CMUL7-E-NOPF-MUON","E",  2);
  oadbLHC11d1->SetHardwareTrigger         ( 2,"V0A && V0C");                                         
  oadbLHC11d1->SetOfflineTrigger          ( 2,"(V0A && V0C) && !V0ABG && !V0CBG  && !TPCLaserWarmUp");

  oadbLHC11d1->AddCollisionTriggerClass   ( AliVEvent::kMUU7,"+CMUU7-B-NOPF-ALLNOTRD","B",  3);
  oadbLHC11d1->AddBGTriggerClass          ( AliVEvent::kMUU7,"+CMUU7-AC-NOPF-ALLNOTRD","AC",3);
  oadbLHC11d1->AddBGTriggerClass          ( AliVEvent::kMUU7,"+CMUU7-E-NOPF-ALLNOTRD","E",  3);
  oadbLHC11d1->SetHardwareTrigger         ( 3,"V0A && V0C");                                         
  oadbLHC11d1->SetOfflineTrigger          ( 3,"(V0A && V0C) && !V0ABG && !V0CBG  && !TPCLaserWarmUp");
  
  oadbLHC11d1->AddCollisionTriggerClass   ( AliVEvent::kEMC7,"+CEMC7-B-NOPF-ALLNOTRD","B",  4);
  oadbLHC11d1->AddBGTriggerClass          ( AliVEvent::kEMC7,"+CEMC7-AC-NOPF-ALLNOTRD","AC",4);
  oadbLHC11d1->AddBGTriggerClass          ( AliVEvent::kEMC7,"+CEMC7-E-NOPF-ALLNOTRD","E",  4);
  oadbLHC11d1->SetHardwareTrigger         ( 4,"V0A && V0C");                      
  oadbLHC11d1->SetOfflineTrigger          ( 4,"(V0A && V0C) && !V0ABG && !V0CBG  && !TPCLaserWarmUp");
  
  oadbLHC11d1->AddCollisionTriggerClass   ( AliVEvent::kMUS7,"+CMUS7-B-NOPF-MUON","B",  5);
  oadbLHC11d1->AddBGTriggerClass          ( AliVEvent::kMUS7,"+CMUS7-AC-NOPF-MUON","AC",5);
  oadbLHC11d1->AddBGTriggerClass          ( AliVEvent::kMUS7,"+CMUS7-E-NOPF-MUON","E",  5);
  oadbLHC11d1->SetHardwareTrigger         ( 5,"V0A && V0C");                                         
  oadbLHC11d1->SetOfflineTrigger          ( 5,"(V0A && V0C) && !V0ABG && !V0CBG  && !TPCLaserWarmUp");  
  
  oadbLHC11d1->AddCollisionTriggerClass   ( AliVEvent::kPHI7,"+CPHI7-B-NOPF-ALLNOTRD","B",  6);
  oadbLHC11d1->AddCollisionTriggerClass   ( AliVEvent::kPHI7,"+CPHI7-I-NOPF-ALLNOTRD","B",  6);
  oadbLHC11d1->AddBGTriggerClass          ( AliVEvent::kPHI7,"+CPHI7-AC-NOPF-ALLNOTRD","AC",6);
  oadbLHC11d1->AddBGTriggerClass          ( AliVEvent::kPHI7,"+CPHI7-E-NOPF-ALLNOTRD","E",  6);
  oadbLHC11d1->SetHardwareTrigger         ( 6,"V0A && V0C");
  oadbLHC11d1->SetOfflineTrigger          ( 6,"(V0A && V0C) && !V0ABG && !V0CBG  && !TPCLaserWarmUp");

  oadbContPS->AppendObject(oadbLHC11d1, 153738, 158793);
  
  // ------------------------------------------------------------------------------------------------------------

  // Trigger Analysis: ZDC timing cuts
 
  AliOADBTriggerAnalysis * oadbTrigAnalysisZDC1 = new AliOADBTriggerAnalysis("ZDCCut1");
  oadbTrigAnalysisZDC1->SetZDCCorrParameters(-66.9, -2.1, 4*0.58, 4*0.5);
  oadbContTriggerAnalysis->AppendObject(oadbTrigAnalysisZDC1, 136851, 137848);  

  AliOADBTriggerAnalysis * oadbTrigAnalysisZDC2 = new AliOADBTriggerAnalysis("ZDCCut2");
  oadbTrigAnalysisZDC2->SetZDCCorrParameters(-66.2, -2.1, 4*0.58, 4*0.5);
  oadbContTriggerAnalysis->AppendObject(oadbTrigAnalysisZDC2, 138125, 138275);  

  AliOADBTriggerAnalysis * oadbTrigAnalysisZDC3 = new AliOADBTriggerAnalysis("ZDCCut3");
  oadbTrigAnalysisZDC3->SetZDCCorrParameters(-65.4, -2.1, 4*0.58, 4*0.5);
  oadbContTriggerAnalysis->AppendObject(oadbTrigAnalysisZDC3, 138359, 138469);  

  AliOADBTriggerAnalysis * oadbTrigAnalysisZDC4 = new AliOADBTriggerAnalysis("ZDCCut4");
  oadbTrigAnalysisZDC4->SetZDCCorrParameters(-67.7, -2.1, 4*0.58, 4*0.5);
  oadbContTriggerAnalysis->AppendObject(oadbTrigAnalysisZDC4, 138533, 138742);  

  AliOADBTriggerAnalysis * oadbTrigAnalysisZDC5 = new AliOADBTriggerAnalysis("ZDCCut5");
  oadbTrigAnalysisZDC5->SetZDCCorrParameters(-67.2, -2.1, 4*0.58, 4*0.5);
  oadbContTriggerAnalysis->AppendObject(oadbTrigAnalysisZDC5, 138795, 138872);  

  AliOADBTriggerAnalysis * oadbTrigAnalysisZDC6 = new AliOADBTriggerAnalysis("ZDCCut6");
  oadbTrigAnalysisZDC6->SetZDCCorrParameters(-65.6, -2.1, 4*0.58, 4*0.5);
  oadbContTriggerAnalysis->AppendObject(oadbTrigAnalysisZDC6, 138924, 139517);  

  oadbTrigAnalysisZDC1->Print();
  oadbTrigAnalysisZDC2->Print();
  oadbTrigAnalysisZDC3->Print();
  oadbTrigAnalysisZDC4->Print();
  oadbTrigAnalysisZDC5->Print();
  oadbTrigAnalysisZDC6->Print();

  
  // ----------------- visualize coverage
  
  //   oadbContPS->List();
  
  // one can get a list of all runs from the RCT with http://alimonitor.cern.ch/configuration/index.jsp?partition=%25&raw_run=&filling_scheme=&filling_config=&fillno=&energy=&intensity_per_bunch=&mu=&interaction_trigger=&rate=&beam_empty_trigger=&empty_empty_trigger=&muon_trigger=&high_multiplicity_trigger=&quality=&muon_quality=&mean_vertex_xyz=&vertex_quality=&comment=&field=&det_aco=&det_emc=&det_fmd=&det_hlt=&det_hmp=&det_mch=&det_mtr=&det_phs=&det_pmd=&det_spd=&det_sdd=&det_ssd=&det_tof=&det_tpc=&det_trd=&det_t00=&det_v00=&det_zdc=&hlt_mode=&changedon=
  Int_t nRuns = 1764;
  Int_t allRuns[] = { 
    159090, 159085, 159076, 159044, 159042, 159040, 158879, 158878, 158877, 158876, 158875, 158868, 158856, 158848, 158844, 158794, 158793, 158792, 158791, 158790, 158788, 158784, 158781, 158780, 158779, 158777, 158776, 158745, 158729, 158722, 158718, 158717, 158714, 158706, 158673, 158626, 158622, 158617, 158615, 158613, 158611, 158608, 158604, 158602, 158598, 158592, 158533, 158531, 158528, 158526, 158522, 158521, 158520, 158518, 158516, 158511, 158509, 158508, 158496, 158495, 158492, 158471, 158469, 158468, 158467, 158341, 158338, 158304, 158303, 158301, 158299, 158293, 158288, 158287, 158285, 158201, 158200, 158196, 158194, 158192, 158191, 158189, 158188, 158185, 158182, 158179, 158177, 158176, 158175, 158173, 158171, 158137, 158136, 158135, 158124, 158120, 158118, 158115, 158114, 158112, 158111, 158110, 158086, 158084, 157976, 157975, 157857, 157848, 157819, 157818, 157770, 157766, 157765, 157734, 157713, 157708, 157707, 157698, 157645, 157569, 157567, 157564, 157562, 157561, 157560, 157496, 157476, 157475, 157277, 157275, 157273, 157272, 157268, 157262, 157261, 157260, 157257, 157227, 157220, 157214, 157212, 157211, 157210, 157209, 157203, 157100, 157098, 157096, 157094, 157092, 157091, 157090, 157087, 157083, 157082, 157079, 157028, 157027, 157025, 157003, 156896, 156895, 156891, 156889, 156860, 156622, 156621, 156481, 156480, 156479, 156477, 155384, 155367, 155365, 155337, 154808, 154796, 154793, 154789, 154787, 154786, 154783, 154780, 154773, 154770, 154763, 154755, 154753, 154750, 154748, 154745, 154742, 154733, 154732, 154726, 154570, 154495, 154485, 154483, 154482, 154480, 154478, 154448, 154385, 154383, 154382, 154296, 154293, 154289, 154286, 154283, 154281, 154273, 154270, 154269, 154266, 154264, 154261, 154257, 154252, 154222, 154221, 154220, 154219, 154211, 154207, 154163, 154158, 154151, 154145, 154143, 154141, 154138, 154136, 154132, 154131, 154130, 154129, 154126, 154125, 154091, 154083, 154081, 154039, 154031, 154030, 154026, 154024, 154018, 153954, 153946, 153944, 153939, 153938, 153935, 153929, 153924, 153916, 153911, 153909, 153906, 153876, 153875, 153873, 153812, 153808, 153807, 153805, 153798, 153796, 153794, 153784, 153781, 153779, 153778, 153777, 153776, 153738, 153733, 153728, 153727, 153726, 153725, 153718, 153709, 153702, 153594, 153591, 153589, 153588, 153587, 153583, 153578, 153571, 153570, 153566, 153560, 153558, 153552, 153548, 153544, 153542, 153541, 153539, 153536, 153533, 153513, 153465, 153415, 153413, 153373, 153371, 153369, 153363, 153362, 153360, 153296, 153234, 153232, 153227, 153223, 153127, 153121, 153120, 153117, 153116, 153115, 153059, 153056, 152935, 152934, 152907, 152823, 152822, 152821, 152820, 152819, 152817, 152780, 152773, 152772, 152751, 152750, 152718, 152717, 152716, 152715, 152708, 152702, 152701, 152698, 152697, 152696, 152695, 152658, 152581, 152513, 152512, 152488, 152455, 152377, 152371, 152369, 152368, 152367, 152334, 152332, 152323, 152322, 152321, 152320, 152319, 152318, 152314, 152313, 152312, 152311, 152310, 152309, 152307, 152306, 152304, 152285, 152284, 152257, 152256, 152214, 152209, 152208, 152207, 152206, 152146, 152138, 152137, 152136, 152091, 152090, 152089, 152087, 152086, 152083, 152082, 152081, 152079, 152078, 152046, 152015, 152011, 152008, 152007, 152005, 152003, 152002, 151852, 151851, 151850, 151849, 151810, 151809, 151752, 151751, 151732, 151724, 151689, 151681, 151680, 151678, 151674, 151672, 151671, 151669, 151666, 151665, 151664, 151661, 151660, 151655, 151638, 151636, 150629, 150518, 150500, 150499, 150440, 150438, 150437, 150434, 150428, 150427, 150423, 150421, 150420, 150375, 150374, 150259, 150256, 150254, 150252, 150250, 150248, 150212, 150211, 150163, 150162, 150160, 150060, 150059, 149975, 149970, 149960, 149931, 149930, 149929, 149927, 149890, 149884, 149883, 149881, 149880, 149790, 149789, 149785, 149781, 149775, 149774, 149773, 149772, 149771, 149769, 149767, 149761, 149757, 149748, 149746, 149740, 149736, 149734, 149732, 149725, 149723, 149722, 149718, 149717, 149688, 149687, 149680, 149678, 149677, 149676, 149675, 149671, 149669, 149665, 149664, 149663, 149656, 149647, 149643, 149608, 149602, 149601, 149597, 149595, 149593, 149591, 149587, 149586, 149581, 149548, 149546, 149542, 149541, 149532, 149531, 149527, 149498, 149497, 149496, 149489, 149487, 149484, 149473, 149469, 149460, 149459, 149451, 149448, 149441, 149436, 149432, 149396, 149395, 149391, 149388, 149387, 149386, 149384, 149380, 149379, 149375, 149367, 149134, 149133, 149130, 149129, 149127, 149113, 149102, 149072, 149071, 149070, 149068, 148857, 148856, 148854, 148853, 148852, 148850, 148847, 148845, 148844, 148843, 148839, 148838, 148835, 148834, 148800, 148797, 148791, 148719, 148711, 148708, 148704, 148697, 148665, 148663, 148659, 148648, 148645, 148630, 148625, 148601, 148592, 148582, 148576, 148569, 148565, 148559, 148556, 148553, 148549, 148547, 148544, 148541, 148538, 148534, 148531, 146860, 146859, 146858, 146857, 146856, 146824, 146817, 146814, 146813, 146812, 146808, 146807, 146806, 146805, 146804, 146803, 146802, 146801, 146748, 146747, 146746, 146689, 146688, 146686, 146459, 146453, 146402, 146401, 146400, 146394, 146393, 146390, 146388, 146386, 146385, 146382, 146381, 146380, 146376, 146374, 146369, 146292, 146287, 146282, 146278, 146277, 146273, 146272, 146223, 146220, 146208, 146205, 146158, 146156, 146153, 146152, 146148, 146147, 146144, 146141, 146105, 146099, 146090, 146089, 146079, 146072, 146071, 146066, 146027, 146026, 146025, 146024, 146023, 146021, 146020, 146018, 146017, 145975, 145974, 145973, 145972, 145967, 145960, 145957, 145956, 145955, 145954, 145681, 145680, 145679, 145674, 145455, 145454, 145448, 145385, 145384, 145383, 145379, 145355, 145354, 145353, 145314, 145309, 145300, 145299, 145294, 145292, 145291, 145290, 145289, 145288, 145182, 145157, 145156, 145145, 145144, 144991, 144871, 139517, 139514, 139513, 139511, 139510, 139507, 139505, 139504, 139503, 139471, 139470, 139467, 139466, 139465, 139441, 139440, 139439, 139438, 139437, 139360, 139329, 139328, 139316, 139314, 139311, 139310, 139309, 139308, 139173, 139172, 139110, 139107, 139105, 139104, 139042, 139038, 139037, 139036, 139034, 139031, 139030, 139029, 139028, 139025, 139024, 138983, 138982, 138980, 138979, 138978, 138977, 138976, 138973, 138972, 138965, 138924, 138872, 138871, 138870, 138837, 138836, 138831, 138830, 138828, 138826, 138796, 138795, 138742, 138740, 138737, 138736, 138732, 138731, 138730, 138666, 138662, 138653, 138652, 138638, 138637, 138624, 138621, 138620, 138583, 138582, 138579, 138578, 138534, 138533, 138469, 138442, 138439, 138438, 138396, 138364, 138359, 138275, 138225, 138201, 138200, 138197, 138192, 138190, 138154, 138153, 138151, 138150, 138126, 138125, 137848, 137847, 137844, 137843, 137752, 137751, 137748, 137724, 137722, 137718, 137704, 137693, 137692, 137691, 137689, 137686, 137685, 137639, 137638, 137609, 137608, 137595, 137549, 137546, 137544, 137541, 137539, 137531, 137530, 137443, 137441, 137440, 137439, 137434, 137432, 137431, 137430, 137370, 137366, 137365, 137243, 137236, 137235, 137232, 137231, 137230, 137165, 137163, 137162, 137161, 137137, 137136, 137135, 137133, 137132, 137125, 137124, 137045, 137042, 136879, 136854, 136851, 136377, 136376, 136375, 136374, 136373, 136372, 136290, 136288, 136199, 136194, 136193, 136192, 136191, 136189, 136181, 136180, 136179, 136177, 136136, 136133, 136130, 136129, 136128, 136125, 136124, 136123, 136118, 136117, 136116, 136115, 136114, 136108, 136107, 136106, 136105, 136104, 136103, 136102, 136101, 135986, 135985, 135983, 135976, 135971, 135969, 135966, 135965, 135958, 135957, 135951, 135949, 135947, 135941, 135900, 135892, 135808, 135795, 135792, 135791, 135788, 135784, 135782, 135780, 135776, 135773, 135771, 135765, 135762, 135761, 135758, 135756, 135751, 135748, 135712, 135711, 135710, 135709, 135707, 135705, 135704, 135701, 135686, 135680, 135674, 135670, 135660, 135658, 135657, 135654, 135029, 134931, 134929, 134925, 134919, 134914, 134905, 134841, 134781, 134780, 134779, 134778, 134691, 134690, 134685, 134681, 134679, 134674, 134672, 134671, 134670, 134667, 134666, 134660, 134657, 134497, 134306, 134305, 134304, 134302, 134301, 134299, 134298, 134297, 134204, 134203, 134202, 134201, 134200, 134199, 134198, 134094, 133985, 133982, 133981, 133979, 133969, 133924, 133920, 133800, 133762, 133761, 133680, 133675, 133674, 133673, 133672, 133671, 133670, 133563, 133419, 133418, 133414, 133330, 133329, 133328, 133327, 133282, 133010, 133007, 133006, 133005, 133004, 130850, 130848, 130847, 130844, 130842, 130840, 130834, 130833, 130831, 130804, 130803, 130802, 130799, 130798, 130795, 130793, 130704, 130696, 130640, 130628, 130627, 130623, 130621, 130620, 130619, 130612, 130609, 130608, 130601, 130526, 130524, 130520, 130519, 130517, 130481, 130480, 130479, 130375, 130369, 130365, 130360, 130358, 130356, 130354, 130353, 130348, 130343, 130342, 130179, 130178, 130172, 130170, 130168, 130158, 130157, 130156, 130151, 130149, 130148, 129983, 129966, 129962, 129961, 129960, 129959, 129763, 129760, 129750, 129748, 129747, 129745, 129744, 129742, 129738, 129736, 129735, 129734, 129731, 129729, 129726, 129725, 129723, 129667, 129666, 129665, 129659, 129655, 129654, 129653, 129652, 129651, 129650, 129649, 129648, 129647, 129641, 129639, 129599, 129598, 129597, 129587, 129586, 129541, 129540, 129536, 129529, 129528, 129527, 129526, 129525, 129524, 129523, 129522, 129521, 129520, 129519, 129516, 129515, 129514, 129513, 129512, 129510, 129508, 129042, 129041, 128913, 128912, 128911, 128910, 128855, 128853, 128850, 128849, 128843, 128836, 128835, 128834, 128833, 128824, 128823, 128820, 128819, 128818, 128817, 128814, 128813, 128778, 128777, 128776, 128678, 128677, 128621, 128615, 128611, 128610, 128609, 128605, 128596, 128594, 128592, 128590, 128589, 128582, 128581, 128507, 128506, 128505, 128504, 128503, 128498, 128495, 128494, 128486, 128483, 128452, 128366, 128359, 128357, 128263, 128262, 128260, 128257, 128256, 128192, 128191, 128190, 128189, 128186, 128185, 128183, 128182, 128180, 128175, 128053, 128050, 127942, 127941, 127940, 127937, 127936, 127935, 127933, 127932, 127931, 127930, 127822, 127820, 127819, 127817, 127815, 127814, 127813, 127730, 127729, 127724, 127723, 127719, 127718, 127715, 127714, 127712, 126437, 126432, 126425, 126424, 126422, 126421, 126420, 126419, 126416, 126411, 126409, 126408, 126407, 126406, 126405, 126404, 126403, 126359, 126352, 126351, 126350, 126285, 126284, 126283, 126177, 126168, 126167, 126162, 126160, 126158, 126097, 126090, 126088, 126087, 126086, 126082, 126081, 126078, 126073, 126008, 126007, 126004, 125855, 125851, 125850, 125849, 125848, 125847, 125844, 125843, 125842, 125841, 125634, 125633, 125632, 125630, 125628, 125296, 125295, 125294, 125292, 125186, 125156, 125140, 125139, 125134, 125133, 125131, 125101, 125100, 125097, 125085, 125083, 125078, 125077, 125076, 125023, 125022, 124886, 124850, 124751, 124750, 124746, 124745, 124702, 124608, 124607, 124606, 124605, 124604, 124603, 124602, 124600, 124400, 124388, 124385, 124383, 124381, 124380, 124378, 124374, 124371, 124367, 124364, 124362, 124360, 124359, 124358, 124355, 124191, 124187, 124186, 124183, 122375, 122374, 122373, 122372, 122195, 121040, 121039, 120829, 120828, 120825, 120824, 120823, 120822, 120821, 120820, 120819, 120818, 120758, 120750, 120746, 120744, 120743, 120742, 120741, 120693, 120692, 120689, 120688, 120687, 120685, 120683, 120682, 120681, 120679, 120678, 120677, 120675, 120672, 120671, 120670, 120625, 120622, 120619, 120617, 120616, 120615, 120614, 120613, 120612, 120611, 120505, 120504, 120503, 120500, 120498, 120496, 120495, 120494, 120493, 120492, 120490, 120489, 120488, 120485, 120483, 120482, 120480, 120479, 120244, 120243, 120242, 120241, 120079, 120078, 120076, 120075, 120073, 120072, 120071, 120069, 120068, 120067, 120066, 120065, 120064, 120004, 120003, 120002, 120001, 120000, 119999, 119998, 119994, 119979, 119978, 119973, 119971, 119970, 119969, 119967, 119965, 119963, 119961, 119959, 119954, 119952, 119948, 119946, 119941, 119935, 119934, 119932, 119926, 119924, 119923, 119917, 119915, 119913, 119909, 119907, 119904, 119903, 119902, 119865, 119862, 119859, 119856, 119854, 119853, 119849, 119846, 119845, 119844, 119843, 119842, 119841, 119840, 119838, 119837, 119164, 119163, 119162, 119161, 119160, 119159, 119158, 119156, 119086, 119085, 119084, 119079, 119077, 119068, 119067, 119061, 119057, 119055, 119048, 119047, 119045, 119043, 119042, 119041, 119040, 119038, 119037, 119035, 119033, 119022, 118903, 118561, 118560, 118558, 118557, 118556, 118518, 118512, 118507, 118506, 118504, 118503, 117223, 117222, 117221, 117220, 117121, 117120, 117119, 117118, 117117, 117116, 117113, 117112, 117110, 117109, 117100, 117099, 117098, 117092, 117086, 117082, 117078, 117077, 117065, 117064, 117063, 117061, 117060, 117059, 117054, 117053, 117052, 117051, 117050, 117049, 117048, 117046, 117045, 117042, 117041, 117039, 117035, 117034, 116948, 116684, 116682, 116681, 116645, 116644, 116643, 116642, 116640, 116611, 116610, 116609, 116574, 116572, 116571, 116562, 116561, 116559, 116403, 116402, 116401, 116288, 116287, 116204, 116203, 116198, 116197, 116134, 116130, 116123, 116118, 116112, 116111, 116102, 116081, 116079, 115892, 115890, 115889, 115888, 115887, 115882, 115881, 115880, 115521, 115516, 115514, 115414, 115413, 115401, 115393, 115345, 115338, 115335, 115328, 115325, 115322, 115318, 115315, 115310, 115237, 115193, 115186, 115173, 115165, 114931, 114930, 114924, 114920, 114919, 114918, 114917, 114916, 114798, 114786, 114785, 114783, 114778, 114757, 114753, 114750, 114747, 114746, 114745, 114744, 114743, 114740, 114737, 114448, 114201, 112980, 105268, 105258, 105257, 105256, 105255, 105254, 105160, 105143, 105108, 104892, 104890, 104876, 104867, 104865, 104861, 104852, 104849, 104845, 104841, 104825, 104824, 104821, 104803, 104802, 104801, 104800, 104799, 104798, 104793, 104792, 104618, 104603, 104439, 104321, 104320, 104316, 104315, 104160, 104159, 104157, 104155, 104065, 101278, 101275, 101268, 101266, 101259, 101235
  };
  
  Int_t min = 1234567890;
  Int_t max = -1;
  for (Int_t i=0; i<oadbContPS->GetNumberOfEntries(); i++)
  {
    Printf("%3d: %d --> %d", i, oadbContPS->LowerLimit(i), oadbContPS->UpperLimit(i));
    min = TMath::Min(min, oadbContPS->LowerLimit(i));
    max = TMath::Max(max, oadbContPS->UpperLimit(i));
  }
  
  for (Int_t i=0; i<nRuns; i++)
  {
    min = TMath::Min(min, allRuns[i]);
    max = TMath::Max(max, allRuns[i]);
  }
  
  
  TH2* coverageHist = new TH2F("coverageHist", "Run coverage;run number", max - min + 1, min - 0.5, max + 0.5, 3, 0, 3);
  coverageHist->SetStats(0);
  coverageHist->GetYaxis()->SetBinLabel(1, "OADB");
  coverageHist->GetYaxis()->SetBinLabel(2, "Physics run");
  coverageHist->GetYaxis()->SetBinLabel(3, "Run w/o OADB");
  
  for (Int_t i=0; i<oadbContPS->GetNumberOfEntries(); i++)
  {
    for (Int_t j=oadbContPS->LowerLimit(i); j <= oadbContPS->UpperLimit(i); j++)
      coverageHist->Fill(j, 0.5, i+1);
  }
  
  for (Int_t i=0; i<nRuns; i++)
  {
    coverageHist->Fill(allRuns[i], 1.5);
    if (coverageHist->GetBinContent(coverageHist->GetXaxis()->FindBin(allRuns[i]), 1) == 0)
      coverageHist->Fill(allRuns[i], 2.5);
  }
  
  new TCanvas;
  coverageHist->Draw("COLZ");
  
  // ----------------- write to file

  oadbContPS->WriteToFile(oadbfilename.Data());
  oadbContFillingScheme->WriteToFile(oadbfilename.Data());
  oadbContTriggerAnalysis->WriteToFile(oadbfilename.Data());

  TFile * fopen = new TFile (oadbfilename); 
  new TBrowser;

}
