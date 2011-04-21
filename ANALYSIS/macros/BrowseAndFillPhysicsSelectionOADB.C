//FIXME: defaults for 2.76 TeV

void BrowseAndFillPhysicsSelectionOADB(Bool_t fill = kFALSE) {

  // Load stuff
  gSystem->Load("libCore.so");  
  gSystem->Load("libTree.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libPhysics.so");
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
  oadbDefaultPP->AddCollisionTriggerClass   ( AliVEvent::kMB,"+CINT1-B-NOPF-ALLNOTRD","B",0);
  oadbDefaultPP->AddBGTriggerClass          ( AliVEvent::kMB,"+CINT1-AC-NOPF-ALLNOTRD","AC",0);
  oadbDefaultPP->AddBGTriggerClass          ( AliVEvent::kMB,"+CINT1-E-NOPF-ALLNOTRD","E",0);  
  oadbDefaultPP->AddCollisionTriggerClass   ( AliVEvent::kMB | AliVEvent::kFastOnly,"+CINT1-B-NOPF-FASTNOTRD -CINT1-B-NOPF-ALLNOTRD","B",0);
  oadbDefaultPP->AddBGTriggerClass          ( AliVEvent::kMB | AliVEvent::kFastOnly,"+CINT1-AC-NOPF-FASTNOTRD -CINT1-AC-NOPF-ALLNOTRD","AC",0);
  oadbDefaultPP->AddBGTriggerClass          ( AliVEvent::kMB | AliVEvent::kFastOnly,"+CINT1-E-NOPF-FASTNOTRD -CINT1-E-NOPF-ALLNOTRD","E",0);  
  oadbDefaultPP->SetHardwareTrigger         ( 0,"SPDGFO >= 1 || V0A || V0C");					      
  oadbDefaultPP->SetOfflineTrigger          ( 0,"(SPDGFO >= 1 || V0A || V0C) && !V0ABG && !V0CBG");

  oadbDefaultPP->AddCollisionTriggerClass   ( AliVEvent::kMUON,"+CMUS1-B-NOPF-MUON","B",1);
  oadbDefaultPP->AddBGTriggerClass          ( AliVEvent::kMUON,"+CMUS1-AC-NOPF-MUON","AC",1);
  oadbDefaultPP->AddBGTriggerClass          ( AliVEvent::kMUON,"+CMUS1-E-NOPF-MUON","E",1);
  oadbDefaultPP->SetHardwareTrigger         ( 1,"SPDGFO >= 1 || V0A || V0C");                                         
  oadbDefaultPP->SetOfflineTrigger          ( 1,"(SPDGFO >= 1 || V0A || V0C) && !V0ABG && !V0CBG");
    
  oadbDefaultPP->AddCollisionTriggerClass   ( AliVEvent::kHighMult,"+CSH1-B-NOPF-ALLNOTRD","B",2);
  oadbDefaultPP->AddBGTriggerClass          ( AliVEvent::kHighMult,"+CSH1-AC-NOPF-ALLNOTRD","AC",2);
  oadbDefaultPP->AddBGTriggerClass          ( AliVEvent::kHighMult,"+CSH1-E-NOPF-ALLNOTRD","E",2);
  oadbDefaultPP->SetHardwareTrigger         ( 2,"SPDGFO >= 1 || V0A || V0C");                      
  oadbDefaultPP->SetOfflineTrigger          ( 2,"(SPDGFO >= 1 || V0A || V0C) && !V0ABG && !V0CBG");

  oadbDefaultPP->AddCollisionTriggerClass   ( AliVEvent::kEMC1,"+CEMC1-B-NOPF-ALLNOTRD","B",2);
  oadbDefaultPP->AddBGTriggerClass          ( AliVEvent::kEMC1,"+CEMC1-AC-NOPF-ALLNOTRD","AC",2);
  oadbDefaultPP->AddBGTriggerClass          ( AliVEvent::kEMC1,"+CEMC1-E-NOPF-ALLNOTRD","E",2);
  oadbDefaultPP->SetHardwareTrigger         ( 2,"SPDGFO >= 1 || V0A || V0C");                      
  oadbDefaultPP->SetOfflineTrigger          ( 2,"(SPDGFO >= 1 || V0A || V0C) && !V0ABG && !V0CBG");


  oadbContPS->AddDefaultObject(oadbDefaultPP);

  // DefaultPbPb
  AliOADBPhysicsSelection * oadbDefaultPbPb = new AliOADBPhysicsSelection("oadbDefaultPbPb");
  oadbDefaultPbPb->AddCollisionTriggerClass   ( AliVEvent::kHighMult,"+C0SMH-B-NOPF-ALLNOTRD","B",0);
  oadbDefaultPbPb->AddBGTriggerClass          ( AliVEvent::kHighMult,"+C0SMH-A-NOPF-ALLNOTRD","A",0);
  oadbDefaultPbPb->AddBGTriggerClass          ( AliVEvent::kHighMult,"+C0SMH-C-NOPF-ALLNOTRD","C",0);
  oadbDefaultPbPb->AddBGTriggerClass          ( AliVEvent::kHighMult,"+C0SMH-E-NOPF-ALLNOTRD","E",0);  
  oadbDefaultPbPb->SetHardwareTrigger         ( 0,"SPDGFO >= 100 && !V0ABG && !V0CBG && ZDCTime");
  oadbDefaultPbPb->SetOfflineTrigger          ( 0,"SPDGFO >= 100");
  
  oadbDefaultPbPb->AddCollisionTriggerClass   ( AliVEvent::kMB,"+CMBACS2-B-NOPF-ALLNOTRD","B",1);
  oadbDefaultPbPb->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBACS2-A-NOPF-ALLNOTRD","A",1);
  oadbDefaultPbPb->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBACS2-C-NOPF-ALLNOTRD","C",1);
  oadbDefaultPbPb->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBACS2-E-NOPF-ALLNOTRD","E",1);
  oadbDefaultPbPb->SetHardwareTrigger         ( 1,"(V0A && V0C && SPDGFO > 1) && !V0ABG && !V0CBG && ZDCTime");
  oadbDefaultPbPb->SetOfflineTrigger          ( 1,"(V0A && V0C && SPDGFO > 1)");

  oadbContPS->AddDefaultObject(oadbDefaultPbPb);

  // Trigger analysis defaults
  AliOADBTriggerAnalysis * oadbTrigAnalysis = new AliOADBTriggerAnalysis("Default");
  oadbTrigAnalysis->SetZDCCorrParameters(-66.5, -2.1, 4*0.58, 4*0.5);
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
  oadbLHC09d10e->SetOfflineTrigger          ( 0,"(SPDGFO >= 1 || V0A || V0C) && !V0ABG && !V0CBG");
  
  oadbLHC09d10e->AddCollisionTriggerClass   ( AliVEvent::kMUON,"+CMUS1B-ABCE-NOPF-MUON","B",1);
  oadbLHC09d10e->AddBGTriggerClass          ( AliVEvent::kMUON,"+CMUS1A-ABCE-NOPF-MUON","A",1);
  oadbLHC09d10e->AddBGTriggerClass          ( AliVEvent::kMUON,"+CMUS1C-ABCE-NOPF-MUON","C",1);
  oadbLHC09d10e->AddBGTriggerClass          ( AliVEvent::kMUON,"+CMUS1-E-NOPF-MUON","E",1);
  oadbLHC09d10e->SetHardwareTrigger         ( 1,"SPDGFO >= 1 || V0A || V0C");                                         
  oadbLHC09d10e->SetOfflineTrigger          ( 1,"(SPDGFO >= 1 || V0A || V0C) && !V0ABG && !V0CBG");
    
  oadbContPS->AppendObject(oadbLHC09d10e, 104065,118555);
  oadbContPS->AppendObject(oadbLHC09d10e->Clone(), 118562,127711);

  // LHC10c, fill 1069 (problems with the V0 online trigger in ESD)
  AliOADBPhysicsSelection * oadbLHC10cV0Bug = new AliOADBPhysicsSelection("oadbLHC10cV0Bug");
  oadbLHC10cV0Bug->AddCollisionTriggerClass   ( AliVEvent::kMB,"+CINT1B-ABCE-NOPF-ALL","B",0);
  oadbLHC10cV0Bug->AddBGTriggerClass          ( AliVEvent::kMB,"+CINT1A-ABCE-NOPF-ALL","A",0);
  oadbLHC10cV0Bug->AddBGTriggerClass          ( AliVEvent::kMB,"+CINT1C-ABCE-NOPF-ALL","C",0);
  oadbLHC10cV0Bug->AddBGTriggerClass          ( AliVEvent::kMB,"+CINT1-E-NOPF-ALL","E",0);  
  oadbLHC10cV0Bug->SetHardwareTrigger         ( 0,"SPDGFO >= 1 || CTPV0A || CTPV0C");
  oadbLHC10cV0Bug->SetOfflineTrigger          ( 0,"(SPDGFO >= 1 || V0A || V0C) && !V0ABG && !V0CBG");
  
  oadbLHC10cV0Bug->AddCollisionTriggerClass   ( AliVEvent::kMUON,"+CMUS1B-ABCE-NOPF-MUON","B",1);
  oadbLHC10cV0Bug->AddBGTriggerClass          ( AliVEvent::kMUON,"+CMUS1A-ABCE-NOPF-MUON","A",1);
  oadbLHC10cV0Bug->AddBGTriggerClass          ( AliVEvent::kMUON,"+CMUS1C-ABCE-NOPF-MUON","C",1);
  oadbLHC10cV0Bug->AddBGTriggerClass          ( AliVEvent::kMUON,"+CMUS1-E-NOPF-MUON","E",1);
  oadbLHC10cV0Bug->SetHardwareTrigger         ( 1,"SPDGFO >= 1 || CTPV0A || CTPV0C");                                         
  oadbLHC10cV0Bug->SetOfflineTrigger          ( 1,"(SPDGFO >= 1 || V0A || V0C) && !V0ABG && !V0CBG");

  oadbContPS->AppendObject(oadbLHC10cV0Bug, 118556,118561);
  

  // LHC10e1
  AliOADBPhysicsSelection * oadbLHC10e1 = new AliOADBPhysicsSelection("oadbLHC10e1");
  oadbLHC10e1->AddCollisionTriggerClass   ( AliVEvent::kMB,"+CINT1-B-NOPF-ALLNOTRD","B",0);
  oadbLHC10e1->AddBGTriggerClass          ( AliVEvent::kMB,"+CINT1-AC-NOPF-ALLNOTRD","AC",0);
  oadbLHC10e1->AddBGTriggerClass          ( AliVEvent::kMB,"+CINT1-E-NOPF-ALLNOTRD","E",0);  
  oadbLHC10e1->SetHardwareTrigger         ( 0,"SPDGFO >= 1 || V0A || V0C");
  oadbLHC10e1->SetOfflineTrigger          ( 0,"(SPDGFO >= 1 || V0A || V0C) && !V0ABG && !V0CBG");
  
  oadbLHC10e1->AddCollisionTriggerClass   ( AliVEvent::kMUON,"+CMUS1-B-NOPF-ALLNOTRD","B",1);
  oadbLHC10e1->AddBGTriggerClass          ( AliVEvent::kMUON,"+CMUS1-AC-NOPF-ALLNOTRD","AC",1);
  oadbLHC10e1->AddBGTriggerClass          ( AliVEvent::kMUON,"+CMUS1-E-NOPF-ALLNOTRD","E",1);
  oadbLHC10e1->SetHardwareTrigger         ( 1,"SPDGFO >= 1 || V0A || V0C");                                         
  oadbLHC10e1->SetOfflineTrigger          ( 1,"(SPDGFO >= 1 || V0A || V0C) && !V0ABG && !V0CBG");
    
  oadbLHC10e1->AddCollisionTriggerClass   ( AliVEvent::kHighMult,"+CSH1-B-NOPF-ALLNOTRD","B",2);
  oadbLHC10e1->AddBGTriggerClass          ( AliVEvent::kHighMult,"+CSH1-AC-NOPF-ALLNOTRD","AC",2);
  oadbLHC10e1->AddBGTriggerClass          ( AliVEvent::kHighMult,"+CSH1-E-NOPF-ALLNOTRD","E",2);
  oadbLHC10e1->SetHardwareTrigger         ( 2,"SPDGFO >= 1 || V0A || V0C");                      
  oadbLHC10e1->SetOfflineTrigger          ( 2,"(SPDGFO >= 1 || V0A || V0C) && !V0ABG && !V0CBG");
    
  oadbContPS->AppendObject(oadbLHC10e1, 127712,127718);

  // lhc10e2
  AliOADBPhysicsSelection * oadbLHC10e2 = new AliOADBPhysicsSelection("oadbLHC10e2");
  oadbLHC10e2->AddCollisionTriggerClass   ( AliVEvent::kMB,"+CINT1B-ABCE-NOPF-ALL","B",0);
  oadbLHC10e2->AddBGTriggerClass          ( AliVEvent::kMB,"+CINT1A-ABCE-NOPF-ALL","A",0);
  oadbLHC10e2->AddBGTriggerClass          ( AliVEvent::kMB,"+CINT1C-ABCE-NOPF-ALL","C",0);
  oadbLHC10e2->AddBGTriggerClass          ( AliVEvent::kMB,"+CINT1-E-NOPF-ALL","E",0);  
  oadbLHC10e2->SetHardwareTrigger         ( 0,"SPDGFO >= 1 || V0A || V0C");
  oadbLHC10e2->SetOfflineTrigger          ( 0,"(SPDGFO >= 1 || V0A || V0C) && !V0ABG && !V0CBG");
  
  oadbLHC10e2->AddCollisionTriggerClass   ( AliVEvent::kMUON,"+CMUS1B-ABCE-NOPF-MUON","B",1);
  oadbLHC10e2->AddBGTriggerClass          ( AliVEvent::kMUON,"+CMUS1A-ABCE-NOPF-MUON","A",1);
  oadbLHC10e2->AddBGTriggerClass          ( AliVEvent::kMUON,"+CMUS1C-ABCE-NOPF-MUON","C",1);
  oadbLHC10e2->AddBGTriggerClass          ( AliVEvent::kMUON,"+CMUS1-E-NOPF-MUON","E",1);
  oadbLHC10e2->SetHardwareTrigger         ( 1,"SPDGFO >= 1 || V0A || V0C");                                         
  oadbLHC10e2->SetOfflineTrigger          ( 1,"(SPDGFO >= 1 || V0A || V0C) && !V0ABG && !V0CBG");
    
  oadbContPS->AppendObject(oadbLHC10e2, 127719,127730);

  // LHC10e3
  AliOADBPhysicsSelection * oadbLHC10e3 = new AliOADBPhysicsSelection("oadbLHC10e3");
  oadbLHC10e3->AddCollisionTriggerClass   ( AliVEvent::kMB,"+CINT1-B-NOPF-ALLNOTRD","B",0);
  oadbLHC10e3->AddBGTriggerClass          ( AliVEvent::kMB,"+CINT1-AC-NOPF-ALLNOTRD","AC",0);
  oadbLHC10e3->AddBGTriggerClass          ( AliVEvent::kMB,"+CINT1-E-NOPF-ALLNOTRD","E",0);  
  oadbLHC10e3->SetHardwareTrigger         ( 0,"SPDGFO >= 1 || V0A || V0C");
  oadbLHC10e3->SetOfflineTrigger          ( 0,"(SPDGFO >= 1 || V0A || V0C) && !V0ABG && !V0CBG");
  
  oadbLHC10e3->AddCollisionTriggerClass   ( AliVEvent::kMUON,"+CMUS1-B-NOPF-ALLNOTRD","B",1);
  oadbLHC10e3->AddBGTriggerClass          ( AliVEvent::kMUON,"+CMUS1-AC-NOPF-ALLNOTRD","AC",1);
  oadbLHC10e3->AddBGTriggerClass          ( AliVEvent::kMUON,"+CMUS1-E-NOPF-ALLNOTRD","E",1);
  oadbLHC10e3->SetHardwareTrigger         ( 1,"SPDGFO >= 1 || V0A || V0C");                                         
  oadbLHC10e3->SetOfflineTrigger          ( 1,"(SPDGFO >= 1 || V0A || V0C) && !V0ABG && !V0CBG");
    
  oadbLHC10e3->AddCollisionTriggerClass   ( AliVEvent::kHighMult,"+CSH1-B-NOPF-ALLNOTRD","B",2);
  oadbLHC10e3->AddBGTriggerClass          ( AliVEvent::kHighMult,"+CSH1-AC-NOPF-ALLNOTRD","AC",2);
  oadbLHC10e3->AddBGTriggerClass          ( AliVEvent::kHighMult,"+CSH1-E-NOPF-ALLNOTRD","E",2);
  oadbLHC10e3->SetHardwareTrigger         ( 2,"SPDGFO >= 1 || V0A || V0C");                                         
  oadbLHC10e3->SetOfflineTrigger          ( 2,"(SPDGFO >= 1 || V0A || V0C) && !V0ABG && !V0CBG");
    
  oadbContPS->AppendObject(oadbLHC10e3, 127813,130850);


  // LHC10f1
  AliOADBPhysicsSelection * oadbLHC10f1 = new AliOADBPhysicsSelection("oadbLHC10f1");
  oadbLHC10f1->AddCollisionTriggerClass   ( AliVEvent::kMB,"+CINT1-B-NOPF-ALLNOTRD","B",0);
  oadbLHC10f1->AddBGTriggerClass          ( AliVEvent::kMB,"+CINT1-AC-NOPF-ALLNOTRD","AC",0);
  oadbLHC10f1->AddBGTriggerClass          ( AliVEvent::kMB,"+CINT1-E-NOPF-ALLNOTRD","E",0);  
  oadbLHC10f1->SetHardwareTrigger         ( 0,"SPDGFO >= 1 || V0A || V0C");
  oadbLHC10f1->SetOfflineTrigger          ( 0,"(SPDGFO >= 1 || V0A || V0C) && !V0ABG && !V0CBG");
  
  oadbLHC10f1->AddCollisionTriggerClass   ( AliVEvent::kMUON,"+CMUS1-B-NOPF-ALLNOTRD","B",1);
  oadbLHC10f1->AddBGTriggerClass          ( AliVEvent::kMUON,"+CMUS1-AC-NOPF-ALLNOTRD","AC",1);
  oadbLHC10f1->AddBGTriggerClass          ( AliVEvent::kMUON,"+CMUS1-E-NOPF-ALLNOTRD","E",1);
  oadbLHC10f1->SetHardwareTrigger         ( 1,"SPDGFO >= 1 || V0A || V0C");                                         
  oadbLHC10f1->SetOfflineTrigger          ( 1,"(SPDGFO >= 1 || V0A || V0C) && !V0ABG && !V0CBG");
    
  oadbLHC10f1->AddCollisionTriggerClass   ( AliVEvent::kHighMult,"+CSH1-B-NOPF-ALLNOTRD","B",2);
  oadbLHC10f1->AddBGTriggerClass          ( AliVEvent::kHighMult,"+CSH1-AC-NOPF-ALLNOTRD","AC",2);
  oadbLHC10f1->AddBGTriggerClass          ( AliVEvent::kHighMult,"+CSH1-E-NOPF-ALLNOTRD","E",2);
  oadbLHC10f1->SetHardwareTrigger         ( 2,"SPDGFO >= 1 || V0A || V0C");                                         
  oadbLHC10f1->SetOfflineTrigger          ( 2,"(SPDGFO >= 1 || V0A || V0C) && !V0ABG && !V0CBG");
    
  oadbContPS->AppendObject(oadbLHC10f1, 133004,134690);

  // LHC10f2
  AliOADBPhysicsSelection * oadbLHC10f2 = new AliOADBPhysicsSelection("oadbLHC10f2");
  oadbLHC10f2->AddCollisionTriggerClass   ( AliVEvent::kCINT5,"+CINT5-B-NOPF-ALLNOTRD","B",0);
  oadbLHC10f2->AddBGTriggerClass          ( AliVEvent::kCINT5,"+CINT5-AC-NOPF-ALLNOTRD","AC",0);
  oadbLHC10f2->AddBGTriggerClass          ( AliVEvent::kCINT5,"+CINT5-E-NOPF-ALLNOTRD","E",0);  
  oadbLHC10f2->SetHardwareTrigger         ( 0,"V0A || V0C");
  oadbLHC10f2->SetOfflineTrigger          ( 0,"(V0A || V0C) && !V0ABG && !V0CBG");
  
  oadbLHC10f2->AddCollisionTriggerClass   ( AliVEvent::kCMUS5,"+CMUS5-B-NOPF-ALLNOTRD","B",1);
  oadbLHC10f2->AddBGTriggerClass          ( AliVEvent::kCMUS5,"+CMUS5-AC-NOPF-ALLNOTRD","AC",1);
  oadbLHC10f2->AddBGTriggerClass          ( AliVEvent::kCMUS5,"+CMUS5-E-NOPF-ALLNOTRD","E",1);
  oadbLHC10f2->SetHardwareTrigger         ( 1,"V0A || V0C");                                         
  oadbLHC10f2->SetOfflineTrigger          ( 1,"(V0A || V0C) && !V0ABG && !V0CBG");
    
  oadbContPS->AppendObject(oadbLHC10f2, 134776,134931);

  // LHC10f3
  AliOADBPhysicsSelection * oadbLHC10f3 = new AliOADBPhysicsSelection("oadbLHC10f3");
  oadbLHC10f3->AddCollisionTriggerClass   ( AliVEvent::kMB,"+CINT1-B-NOPF-ALLNOTRD","B",0);
  oadbLHC10f3->AddBGTriggerClass          ( AliVEvent::kMB,"+CINT1-AC-NOPF-ALLNOTRD","AC",0);
  oadbLHC10f3->AddBGTriggerClass          ( AliVEvent::kMB,"+CINT1-E-NOPF-ALLNOTRD","E",0);  
  oadbLHC10f3->SetHardwareTrigger         ( 0,"SPDGFO >= 1 || V0A || V0C");
  oadbLHC10f3->SetOfflineTrigger          ( 0,"(SPDGFO >= 1 || V0A || V0C) && !V0ABG && !V0CBG");
  
  oadbLHC10f3->AddCollisionTriggerClass   ( AliVEvent::kMUON,"+CMUS1-B-NOPF-ALLNOTRD","B",1);
  oadbLHC10f3->AddBGTriggerClass          ( AliVEvent::kMUON,"+CMUS1-AC-NOPF-ALLNOTRD","AC",1);
  oadbLHC10f3->AddBGTriggerClass          ( AliVEvent::kMUON,"+CMUS1-E-NOPF-ALLNOTRD","E",1);
  oadbLHC10f3->SetHardwareTrigger         ( 1,"SPDGFO >= 1 || V0A || V0C");                                         
  oadbLHC10f3->SetOfflineTrigger          ( 1,"(SPDGFO >= 1 || V0A || V0C) && !V0ABG && !V0CBG");
    
  oadbLHC10f3->AddCollisionTriggerClass   ( AliVEvent::kHighMult,"+CSH1-B-NOPF-ALLNOTRD","B",2);
  oadbLHC10f3->AddBGTriggerClass          ( AliVEvent::kHighMult,"+CSH1-AC-NOPF-ALLNOTRD","AC",2);
  oadbLHC10f3->AddBGTriggerClass          ( AliVEvent::kHighMult,"+CSH1-E-NOPF-ALLNOTRD","E",2);
  oadbLHC10f3->SetHardwareTrigger         ( 1,"SPDGFO >= 1 || V0A || V0C");                                         
  oadbLHC10f3->SetOfflineTrigger          ( 1,"(SPDGFO >= 1 || V0A || V0C) && !V0ABG && !V0CBG");
    
  oadbContPS->AppendObject(oadbLHC10f3, 135029,135029);


  // LHC10g
  AliOADBPhysicsSelection * oadbLHC10g = new AliOADBPhysicsSelection("oadbLHC10g");
  oadbLHC10g->AddCollisionTriggerClass   ( AliVEvent::kMB,"+CINT1-B-NOPF-ALLNOTRD","B",0);
  oadbLHC10g->AddBGTriggerClass          ( AliVEvent::kMB,"+CINT1-AC-NOPF-ALLNOTRD","AC",0);
  oadbLHC10g->AddBGTriggerClass          ( AliVEvent::kMB,"+CINT1-E-NOPF-ALLNOTRD","E",0);  
  oadbLHC10g->SetHardwareTrigger         ( 0,"SPDGFO >= 1 || V0A || V0C");
  oadbLHC10g->SetOfflineTrigger          ( 0,"(SPDGFO >= 1 || V0A || V0C) && !V0ABG && !V0CBG");
  
  oadbLHC10g->AddCollisionTriggerClass   ( AliVEvent::kMUON,"+CMUS1-B-NOPF-ALLNOTRD","B",1);
  oadbLHC10g->AddBGTriggerClass          ( AliVEvent::kMUON,"+CMUS1-AC-NOPF-ALLNOTRD","AC",1);
  oadbLHC10g->AddBGTriggerClass          ( AliVEvent::kMUON,"+CMUS1-E-NOPF-ALLNOTRD","E",1);
  oadbLHC10g->SetHardwareTrigger         ( 1,"SPDGFO >= 1 || V0A || V0C");                                         
  oadbLHC10g->SetOfflineTrigger          ( 1,"(SPDGFO >= 1 || V0A || V0C) && !V0ABG && !V0CBG");
    
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



  // ----- Heavy Ion -------

  // LHC10h1
  AliOADBPhysicsSelection * oadbLHC10h1 = new AliOADBPhysicsSelection("oadbLHC10h1");
  oadbLHC10h1->AddCollisionTriggerClass   ( AliVEvent::kMB,"+CMBAC-B-NOPF-ALL","B",0);
  oadbLHC10h1->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBAC-A-NOPF-ALL","A",0);
  oadbLHC10h1->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBAC-C-NOPF-ALL","C",0);  
  oadbLHC10h1->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBAC-E-NOPF-ALL","E",0);  
  oadbLHC10h1->SetHardwareTrigger         ( 0,"V0A && V0C");
  oadbLHC10h1->SetOfflineTrigger          ( 0,"(V0A && V0C && SPDGFOL1 > 1) && !V0ABG && !V0CBG && ZDCTime");
  
  oadbLHC10h1->AddCollisionTriggerClass   ( AliVEvent::kMB,"+CMBS2A-B-NOPF-ALL","B",1);
  oadbLHC10h1->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBS2A-A-NOPF-ALL","A",1);
  oadbLHC10h1->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBS2A-C-NOPF-ALL","C",1);  
  oadbLHC10h1->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBS2A-E-NOPF-ALL","E",1);  
  oadbLHC10h1->SetHardwareTrigger         ( 1,"SPDGFOL1 > 1 && V0A");
  oadbLHC10h1->SetOfflineTrigger          ( 1,"(V0A && V0C && SPDGFOL1 > 1) && !V0ABG && !V0CBG && ZDCTime");

  oadbLHC10h1->AddCollisionTriggerClass   ( AliVEvent::kMB,"+CMBS2C-B-NOPF-ALL","B",2);
  oadbLHC10h1->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBS2C-A-NOPF-ALL","A",2);
  oadbLHC10h1->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBS2C-C-NOPF-ALL","C",2);  
  oadbLHC10h1->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBS2C-E-NOPF-ALL","E",2);  
  oadbLHC10h1->SetHardwareTrigger         ( 2,"SPDGFOL1 > 1 && V0C");
  oadbLHC10h1->SetOfflineTrigger          ( 2,"(V0A && V0C && SPDGFOL1 > 1) && !V0ABG && !V0CBG && ZDCTime");
    
  oadbContPS->AppendObject(oadbLHC10h1, 136851,136879);

  // LHC10h2
  AliOADBPhysicsSelection * oadbLHC10h2 = new AliOADBPhysicsSelection("oadbLHC10h2");
  oadbLHC10h2->AddCollisionTriggerClass   ( AliVEvent::kHighMult,"+C0SMH-B-NOPF-ALL","B",0);
  oadbLHC10h2->AddBGTriggerClass          ( AliVEvent::kHighMult,"+C0SMH-A-NOPF-ALL","A",0);
  oadbLHC10h2->AddBGTriggerClass          ( AliVEvent::kHighMult,"+C0SMH-C-NOPF-ALL","C",0);
  oadbLHC10h2->AddBGTriggerClass          ( AliVEvent::kHighMult,"+C0SMH-E-NOPF-ALL","E",0);  
  oadbLHC10h2->SetHardwareTrigger         ( 0,"SPDGFO >= 100");
  oadbLHC10h2->SetOfflineTrigger          ( 0,"SPDGFO >= 100 && !V0ABG && !V0CBG && ZDCTime");
  
  oadbLHC10h2->AddCollisionTriggerClass   ( AliVEvent::kMB,"+CMBAC-B-NOPF-ALL","B",1);
  oadbLHC10h2->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBAC-A-NOPF-ALL","A",1);
  oadbLHC10h2->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBAC-C-NOPF-ALL","C",1);
  oadbLHC10h2->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBAC-E-NOPF-ALL","E",1);
  oadbLHC10h2->SetHardwareTrigger         ( 1,"V0A && V0C");
  oadbLHC10h2->SetOfflineTrigger          ( 1,"(V0A && V0C && SPDGFOL1 > 1) && !V0ABG && !V0CBG && ZDCTime");
 
  oadbContPS->AppendObject(oadbLHC10h2, 137042,137133);

  // LHC10h3
  AliOADBPhysicsSelection * oadbLHC10h3 = new AliOADBPhysicsSelection("oadbLHC10h3");
  oadbLHC10h3->AddCollisionTriggerClass   ( AliVEvent::kHighMult,"+C0SMH-B-NOPF-ALL","B",0);
  oadbLHC10h3->AddBGTriggerClass          ( AliVEvent::kHighMult,"+C0SMH-A-NOPF-ALL","A",0);
  oadbLHC10h3->AddBGTriggerClass          ( AliVEvent::kHighMult,"+C0SMH-C-NOPF-ALL","C",0);
  oadbLHC10h3->AddBGTriggerClass          ( AliVEvent::kHighMult,"+C0SMH-E-NOPF-ALL","E",0);  
  oadbLHC10h3->SetHardwareTrigger         ( 0,"SPDGFO >= 100");
  oadbLHC10h3->SetOfflineTrigger          ( 0,"SPDGFO >= 100 && !V0ABG && !V0CBG && ZDCTime");
  
  oadbLHC10h3->AddCollisionTriggerClass   ( AliVEvent::kMB,"+CMBAC-B-NOPF-ALL","B",1);
  oadbLHC10h3->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBAC-A-NOPF-ALL","A",1);
  oadbLHC10h3->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBAC-C-NOPF-ALL","C",1);
  oadbLHC10h3->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBAC-E-NOPF-ALL","E",1);
  oadbLHC10h3->SetHardwareTrigger         ( 1,"V0A && V0C");
  oadbLHC10h3->SetOfflineTrigger          ( 1,"(V0A && V0C && SPDGFOL1 > 1) && !V0ABG && !V0CBG && ZDCTime");
    
  oadbLHC10h3->AddCollisionTriggerClass   ( AliVEvent::kMB,"+CMBS2A-B-NOPF-ALL","B",2);
  oadbLHC10h3->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBS2A-A-NOPF-ALL","A",2);
  oadbLHC10h3->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBS2A-C-NOPF-ALL","C",2);
  oadbLHC10h3->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBS2A-E-NOPF-ALL","E",2);
  oadbLHC10h3->SetHardwareTrigger         ( 2,"SPDGFOL1 > 1 && V0A");
  oadbLHC10h3->SetOfflineTrigger          ( 2,"(V0A && V0C && SPDGFOL1 > 1) && !V0ABG && !V0CBG && ZDCTime");

  oadbLHC10h3->AddCollisionTriggerClass   ( AliVEvent::kMB,"+CMBS2C-B-NOPF-ALL","B",3);
  oadbLHC10h3->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBS2C-A-NOPF-ALL","A",3);
  oadbLHC10h3->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBS2C-C-NOPF-ALL","C",3);
  oadbLHC10h3->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBS2C-E-NOPF-ALL","E",3);
  oadbLHC10h3->SetHardwareTrigger         ( 3,"SPDGFOL1 > 1 && V0C");
  oadbLHC10h3->SetOfflineTrigger          ( 3,"(V0A && V0C && SPDGFOL1 > 1) && !V0ABG && !V0CBG && ZDCTime");

  oadbContPS->AppendObject(oadbLHC10h3, 137135,137364);

 // LHC10h, run 137365
  AliOADBPhysicsSelection * oadbLHC10h_137365 = new AliOADBPhysicsSelection("oadbLHC10h_x9");
  oadbLHC10h_137365->AddCollisionTriggerClass   ( AliVEvent::kMB,"+CMBAC-B-NOPF-ALL","B",0);
  oadbLHC10h_137365->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBAC-A-NOPF-ALL","A",0);
  oadbLHC10h_137365->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBAC-C-NOPF-ALL","C",0);
  oadbLHC10h_137365->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBAC-E-NOPF-ALL","E",0);
  oadbLHC10h_137365->SetHardwareTrigger         ( 0,"V0A && V0C");
  oadbLHC10h_137365->SetOfflineTrigger          ( 0,"(V0A && V0C && SPDGFOL1 > 1) && !V0ABG && !V0CBG && ZDCTime");

  oadbContPS->AppendObject(oadbLHC10h_137365, 137365,137365);

  AliOADBPhysicsSelection * oadbLHC10h4 = new AliOADBPhysicsSelection("oadbLHC10h4");
  oadbLHC10h4->AddCollisionTriggerClass   ( AliVEvent::kHighMult,"+C0SMH-B-NOPF-ALL","B",0);
  oadbLHC10h4->AddBGTriggerClass          ( AliVEvent::kHighMult,"+C0SMH-A-NOPF-ALL","A",0);
  oadbLHC10h4->AddBGTriggerClass          ( AliVEvent::kHighMult,"+C0SMH-C-NOPF-ALL","C",0);
  oadbLHC10h4->AddBGTriggerClass          ( AliVEvent::kHighMult,"+C0SMH-E-NOPF-ALL","E",0);  
  oadbLHC10h4->SetHardwareTrigger         ( 0,"SPDGFO >= 100");
  oadbLHC10h4->SetOfflineTrigger          ( 0,"SPDGFO >= 100 && !V0ABG && !V0CBG && ZDCTime");
  
  oadbLHC10h4->AddCollisionTriggerClass   ( AliVEvent::kMB,"+CMBAC-B-NOPF-ALL","B",1);
  oadbLHC10h4->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBAC-A-NOPF-ALL","A",1);
  oadbLHC10h4->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBAC-C-NOPF-ALL","C",1);
  oadbLHC10h4->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBAC-E-NOPF-ALL","E",1);
  oadbLHC10h4->SetHardwareTrigger         ( 1,"V0A && V0C");
  oadbLHC10h4->SetOfflineTrigger          ( 1,"(V0A && V0C && SPDGFOL1 > 1) && !V0ABG && !V0CBG && ZDCTime");
    
  oadbLHC10h4->AddCollisionTriggerClass   ( AliVEvent::kMB,"+CMBS2A-B-NOPF-ALL","B",2);
  oadbLHC10h4->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBS2A-A-NOPF-ALL","A",2);
  oadbLHC10h4->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBS2A-C-NOPF-ALL","C",2);
  oadbLHC10h4->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBS2A-E-NOPF-ALL","E",2);
  oadbLHC10h4->SetHardwareTrigger         ( 2,"SPDGFOL1 > 1 && V0A");
  oadbLHC10h4->SetOfflineTrigger          ( 2,"(V0A && V0C && SPDGFOL1 > 1) && !V0ABG && !V0CBG && ZDCTime");

  oadbLHC10h4->AddCollisionTriggerClass   ( AliVEvent::kMB,"+CMBS2C-B-NOPF-ALL","B",3);
  oadbLHC10h4->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBS2C-A-NOPF-ALL","A",3);
  oadbLHC10h4->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBS2C-C-NOPF-ALL","C",3);
  oadbLHC10h4->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBS2C-E-NOPF-ALL","E",3);
  oadbLHC10h4->SetHardwareTrigger         ( 3,"SPDGFOL1 > 1 && V0C");
  oadbLHC10h4->SetOfflineTrigger          ( 3,"(V0A && V0C && SPDGFOL1 > 1) && !V0ABG && !V0CBG && ZDCTime");

  oadbContPS->AppendObject(oadbLHC10h4, 137366,137595);



  // LHC10h5
  AliOADBPhysicsSelection * oadbLHC10h5 = new AliOADBPhysicsSelection("oadbLHC10h5");
  oadbLHC10h5->AddCollisionTriggerClass   ( AliVEvent::kHighMult,"+C0SMH-B-NOPF-ALL","B",0);
  oadbLHC10h5->AddBGTriggerClass          ( AliVEvent::kHighMult,"+C0SMH-A-NOPF-ALL","A",0);
  oadbLHC10h5->AddBGTriggerClass          ( AliVEvent::kHighMult,"+C0SMH-C-NOPF-ALL","C",0);
  oadbLHC10h5->AddBGTriggerClass          ( AliVEvent::kHighMult,"+C0SMH-E-NOPF-ALL","E",0);  
  oadbLHC10h5->SetHardwareTrigger         ( 0,"SPDGFO >= 100");
  oadbLHC10h5->SetOfflineTrigger          ( 0,"SPDGFO >= 100 && !V0ABG && !V0CBG && ZDCTime");
  
  oadbLHC10h5->AddCollisionTriggerClass   ( AliVEvent::kMB,"+CMBAC-B-NOPF-ALL","B",1);
  oadbLHC10h5->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBAC-A-NOPF-ALL","A",1);
  oadbLHC10h5->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBAC-C-NOPF-ALL","C",1);
  oadbLHC10h5->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBAC-E-NOPF-ALL","E",1);
  oadbLHC10h5->SetHardwareTrigger         ( 1,"V0A && V0C");
  oadbLHC10h5->SetOfflineTrigger          ( 1,"(V0A && V0C && SPDGFOL1 > 1) && !V0ABG && !V0CBG && ZDCTime");

  oadbContPS->AppendObject(oadbLHC10h5, 137608,137848); // 138982,138983);

  // LHC10h6
  AliOADBPhysicsSelection * oadbLHC10h6 = new AliOADBPhysicsSelection("oadbLHC10h6");
  oadbLHC10h6->AddCollisionTriggerClass   ( AliVEvent::kHighMult,"+C0SMH-B-NOPF-ALL","B",0);
  oadbLHC10h6->AddBGTriggerClass          ( AliVEvent::kHighMult,"+C0SMH-A-NOPF-ALL","A",0);
  oadbLHC10h6->AddBGTriggerClass          ( AliVEvent::kHighMult,"+C0SMH-C-NOPF-ALL","C",0);
  oadbLHC10h6->AddBGTriggerClass          ( AliVEvent::kHighMult,"+C0SMH-E-NOPF-ALL","E",0);  
  oadbLHC10h6->SetHardwareTrigger         ( 0,"SPDGFO >= 100");
  oadbLHC10h6->SetOfflineTrigger          ( 0,"SPDGFO >= 100 && !V0ABG && !V0CBG && ZDCTime");
  
  oadbLHC10h6->AddCollisionTriggerClass   ( AliVEvent::kMB,"+CMBACS2-B-NOPF-ALL","B",1);
  oadbLHC10h6->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBACS2-A-NOPF-ALL","A",1);
  oadbLHC10h6->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBACS2-C-NOPF-ALL","C",1);
  oadbLHC10h6->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBACS2-E-NOPF-ALL","E",1);
  oadbLHC10h6->SetHardwareTrigger         ( 1,"(V0A && V0C && SPDGFOL1 > 1)");
  oadbLHC10h6->SetOfflineTrigger          ( 1,"(V0A && V0C && SPDGFOL1 > 1) && !V0ABG && !V0CBG && ZDCTime");

  oadbContPS->AppendObject(oadbLHC10h6, 138125,138980);// 139028, 139316);


  AliOADBPhysicsSelection * oadbLHC10h7 = new AliOADBPhysicsSelection("oadbLHC10h7");
  oadbLHC10h7->AddCollisionTriggerClass   ( AliVEvent::kHighMult,"+C0SMH-B-NOPF-ALL","B",0);
  oadbLHC10h7->AddBGTriggerClass          ( AliVEvent::kHighMult,"+C0SMH-A-NOPF-ALL","A",0);
  oadbLHC10h7->AddBGTriggerClass          ( AliVEvent::kHighMult,"+C0SMH-C-NOPF-ALL","C",0);
  oadbLHC10h7->AddBGTriggerClass          ( AliVEvent::kHighMult,"+C0SMH-E-NOPF-ALL","E",0);  
  oadbLHC10h7->SetHardwareTrigger         ( 0,"SPDGFO >= 100");
  oadbLHC10h7->SetOfflineTrigger          ( 0,"SPDGFO >= 100 && !V0ABG && !V0CBG && ZDCTime");
  
  oadbLHC10h7->AddCollisionTriggerClass   ( AliVEvent::kMB,"+CMBAC-B-NOPF-ALL","B",1);
  oadbLHC10h7->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBAC-A-NOPF-ALL","A",1);
  oadbLHC10h7->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBAC-C-NOPF-ALL","C",1);
  oadbLHC10h7->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBAC-E-NOPF-ALL","E",1);
  oadbLHC10h7->SetHardwareTrigger         ( 1,"(V0A && V0C)");
  oadbLHC10h7->SetOfflineTrigger          ( 1,"(V0A && V0C && SPDGFOL1 > 1) && !V0ABG && !V0CBG && ZDCTime");


  oadbContPS->AppendObject(oadbLHC10h7, 138982,138983);


  // LHC10h8
  AliOADBPhysicsSelection * oadbLHC10h8 = new AliOADBPhysicsSelection("oadbLHC10h8");
  oadbLHC10h8->AddCollisionTriggerClass   ( AliVEvent::kHighMult,"+C0SMH-B-NOPF-ALL","B",0);
  oadbLHC10h8->AddBGTriggerClass          ( AliVEvent::kHighMult,"+C0SMH-A-NOPF-ALL","A",0);
  oadbLHC10h8->AddBGTriggerClass          ( AliVEvent::kHighMult,"+C0SMH-C-NOPF-ALL","C",0);
  oadbLHC10h8->AddBGTriggerClass          ( AliVEvent::kHighMult,"+C0SMH-E-NOPF-ALL","E",0);  
  oadbLHC10h8->SetHardwareTrigger         ( 0,"SPDGFO >= 100");
  oadbLHC10h8->SetOfflineTrigger          ( 0,"SPDGFO >= 100 && !V0ABG && !V0CBG && ZDCTime");
  
  oadbLHC10h8->AddCollisionTriggerClass   ( AliVEvent::kMB,"+CMBACS2-B-NOPF-ALL","B",1);
  oadbLHC10h8->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBACS2-A-NOPF-ALL","A",1);
  oadbLHC10h8->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBACS2-C-NOPF-ALL","C",1);
  oadbLHC10h8->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBACS2-E-NOPF-ALL","E",1);
  oadbLHC10h8->SetHardwareTrigger         ( 1,"(V0A && V0C && SPDGFOL1 > 1)");
  oadbLHC10h8->SetOfflineTrigger          ( 1,"(V0A && V0C && SPDGFOL1 > 1) && !V0ABG && !V0CBG && ZDCTime");

  oadbContPS->AppendObject(oadbLHC10h8,  139028, 139316);


  // LHC10h9
  AliOADBPhysicsSelection * oadbLHC10h9 = new AliOADBPhysicsSelection("oadbLHC10h9");
  oadbLHC10h9->AddCollisionTriggerClass   ( AliVEvent::kHighMult,"+C0SMH-B-NOPF-ALLNOTRD","B",0);
  oadbLHC10h9->AddBGTriggerClass          ( AliVEvent::kHighMult,"+C0SMH-A-NOPF-ALLNOTRD","A",0);
  oadbLHC10h9->AddBGTriggerClass          ( AliVEvent::kHighMult,"+C0SMH-C-NOPF-ALLNOTRD","C",0);
  oadbLHC10h9->AddBGTriggerClass          ( AliVEvent::kHighMult,"+C0SMH-E-NOPF-ALLNOTRD","E",0);  
  oadbLHC10h9->SetHardwareTrigger         ( 0,"SPDGFO >= 100");
  oadbLHC10h9->SetOfflineTrigger          ( 0,"SPDGFO >= 100 && !V0ABG && !V0CBG && ZDCTime");
  
  oadbLHC10h9->AddCollisionTriggerClass   ( AliVEvent::kMB,"+CMBACS2-B-NOPF-ALLNOTRD","B",1);
  oadbLHC10h9->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBACS2-A-NOPF-ALLNOTRD","A",1);
  oadbLHC10h9->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBACS2-C-NOPF-ALLNOTRD","C",1);
  oadbLHC10h9->AddBGTriggerClass          ( AliVEvent::kMB,"+CMBACS2-E-NOPF-ALLNOTRD","E",1);
  oadbLHC10h9->SetHardwareTrigger         ( 1,"(V0A && V0C && SPDGFOL1 > 1)");
  oadbLHC10h9->SetOfflineTrigger          ( 1,"(V0A && V0C && SPDGFOL1 > 1) && !V0ABG && !V0CBG && ZDCTime");

  oadbContPS->AppendObject(oadbLHC10h9, 139328,139517);


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


  // ----- 2011 -----
  // ----- proton-proton -----


  oadbContPS->WriteToFile(oadbfilename.Data());
  oadbContFillingScheme->WriteToFile(oadbfilename.Data());
  oadbContTriggerAnalysis->WriteToFile(oadbfilename.Data());

  TFile * fopen = new TFile (oadbfilename); 
  new TBrowser;

}
