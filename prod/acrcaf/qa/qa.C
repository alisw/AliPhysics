void qa(Int_t runNumber, Int_t doQASym=1, Int_t doCosmic=0) {
  TStopwatch timer;
  timer.Start();

  gEnv->SetValue("XSec.GSI.DelegProxy","2");
  // Select ROOT version
  TProof::Mgr("aliprod@alicecaf")->SetROOTVersion("v5-26-00");
  // Login to CAF
  TProof::Open("aliprod@alicecaf");

  // Enable AliRoot
 
  gProof->UploadPackage("/afs/cern.ch/alice/caf/sw/ALICE/PARs/v4-18-Release.rec/AF-v4-18-rec.par");
  gProof->EnablePackage("AF-v4-18-rec.par");

  // Enable analysis libs
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gProof->Exec("gSystem->Load(\"libANALYSIS.so\");"     ,kTRUE);
  gProof->Exec("gSystem->Load(\"libANALYSISalice.so\");",kTRUE);

  //____________________________________________//
  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("AnalysisQAManager");
  AliESDInputHandler* esdH = new AliESDInputHandler();
  mgr->SetInputEventHandler(esdH);  
  mgr->SetDebugLevel(3);

  //____________________________________________//
  // 1st Cosmic task
  if (doCosmic) {
    TFile::Cp("$ALICE_ROOT/PWG1/cosmic/AliAnalysisTaskCosmic.h",
              "file:AliAnalysisTaskCosmic.h");
    TFile::Cp("$ALICE_ROOT/PWG1/cosmic/AliAnalysisTaskCosmic.cxx",
              "file:AliAnalysisTaskCosmic.cxx");
    gProof->Load("AliAnalysisTaskCosmic.cxx++g");

    AliAnalysisTaskCosmic *taskcosmic = new AliAnalysisTaskCosmic("TaskCosmic");
    mgr->AddTask(taskcosmic);

    // Create containers for input/output
    AliAnalysisDataContainer *coutputcosmic = mgr->CreateContainer("chist1",TList::Class(),AliAnalysisManager::kOutputContainer,
					    Form("cosmic%d.root",runNumber));
    mgr->ConnectInput  (taskcosmic,  0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput (taskcosmic,  1, coutputcosmic);
  }  

  if (doQASym) {
    //____________________________________________//
    // QA task for central barrel tracking exploiting symmetries
    TFile::Cp("$ALICE_ROOT/PWG1/AliAnalysisTaskQASym.h",
              "file:AliAnalysisTaskQASym.h");
    TFile::Cp("$ALICE_ROOT/PWG1/AliAnalysisTaskQASym.cxx",
              "file:AliAnalysisTaskQASym.cxx");
    gProof->Load("AliAnalysisTaskQASym.cxx++g");
    gROOT->LoadMacro("$ALICE_ROOT/PWG1/PilotTrain/AddTaskQAsym.C");
    AliAnalysisTaskSE* taskqasym = AddTaskQAsym();
    //
    // QA task for VZERO
    TFile::Cp("$ALICE_ROOT/PWG1/AliAnaVZEROQA.h",
              "file:AliAnaVZEROQA.h");
    TFile::Cp("$ALICE_ROOT/PWG1/AliAnaVZEROQA.cxx",
              "file:AliAnaVZEROQA.cxx");
    gProof->Load("AliAnaVZEROQA.cxx++g");
    gROOT->LoadMacro("$ALICE_ROOT/PWG1/PilotTrain/AddTaskVZEROQA.C");
    AliAnalysisTaskSE* taskvzeroqa = AddTaskVZEROQA(runNumber);
    //
    // QA task for vertexing
    TFile::Cp("$ALICE_ROOT/PWG1/global/AliAnalysisTaskVertexESD.h",
              "file:AliAnalysisTaskVertexESD.h");
    TFile::Cp("$ALICE_ROOT/PWG1/global/AliAnalysisTaskVertexESD.cxx",
              "file:AliAnalysisTaskVertexESD.cxx");
    gProof->Load("AliAnalysisTaskVertexESD.cxx++g");
    gROOT->LoadMacro("$ALICE_ROOT/PWG1/macros/AddTaskVertexESD.C");
    AliAnalysisTaskSE* taskvertex = AddTaskVertexESD();
  }
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  mgr->StartAnalysis("proof",
		     Form("/ALIREC/aliprod/run%d",runNumber));

  timer.Stop();
  timer.Print();

  plot(runNumber);
}

void plot(Int_t runNumber, Int_t doQASym=1, Int_t doCosmic=0)
{
  if (doCosmic && gSystem->AccessPathName(Form("cosmic%d.root",runNumber))) {
     printf("Error: doCosmic requested but file cosmic%d.root not found.\n", runNumber);
     doCosmic = 0;
  }   
  if (doQASym && gSystem->AccessPathName("QAsym.root")) {
     printf("Error: doCosmic requested but file QAsym.root not found.\n", runNumber);
     doQASym = 0;
  }   
  if (doCosmic) {
    TFile* f1 = new TFile(Form("cosmic%d.root",runNumber), "read");

    // pt, phi ....
    TCanvas* c1 = new TCanvas("c1", "Cosmic: pt,eta,phi", 10, 10, 1100, 800);
    c1->Divide(2,2);
    c1->cd(1);
    c1->GetPad(1)->SetLogy();
    TList *chist1 = (TList*)f1.Get("chist1");
    fhPtP = (TH1F*) chist1->FindObject("fhPtPC");
    fhPtP->SetLineColor(2);
    fhPtP->Draw();
    fhPtN = (TH1F*) chist1->FindObject("fhPtNC");
    fhPtN->SetLineColor(4);
    fhPtN->Draw("same");

    c1->cd(2);
    fhPhiP = (TH1F*) chist1->FindObject("fhPhiPC");
    fhPhiP->SetLineColor(2);
    fhPhiP->Draw();
    fhPhiN = (TH1F*) chist1->FindObject("fhPhiNC");
    fhPhiN->SetLineColor(4);
    fhPhiN->Draw("same");

    c1->cd(3);
    fhThetaP = (TH1F*) chist1->FindObject("fhThetaPC");
    fhThetaP->SetLineColor(2);
    fhThetaP->Draw();
    fhThetaN = (TH1F*) chist1->FindObject("fhThetaNC");
    fhThetaN->SetLineColor(4);
    fhThetaN->Draw("same");
    
    // delta phi
    TCanvas* c2 = new TCanvas("c2", "Cosmic: #Delta#phi", 10, 10, 1100, 800);
    c2->Divide(2,2);
    // pos/neg charge
    c2->cd(1);
    c2->GetPad(1)->SetLogy();
    fhDPhiP = (TH1F*) chist1->FindObject("fhDPhiPC");
    fhDPhiN = (TH1F*) chist1->FindObject("fhDPhiNC");
    fhDPhiP->SetLineColor(2);
    fhDPhiN->SetLineColor(4);
    fhDPhiP->SetXTitle("#Delta#phi [rad]");
    fhDPhiP->Draw();
    fhDPhiN->Draw("same");

    // pos/neg z
    c2->cd(2);
    c2->GetPad(2)->SetLogy();

    fhDPhiZP = (TH1F*) chist1->FindObject("fhDPhiPZ");
    fhDPhiZN = (TH1F*) chist1->FindObject("fhDPhiNZ");
    fhDPhiBA = (TH1F*) chist1->FindObject("fhDPhi_Bad");
    fhDPhiZP->SetLineColor(2);
    fhDPhiZN->SetLineColor(4);
    fhDPhiBA->SetLineColor(6);

    fhDPhiZP->SetXTitle("#Delta#phi [rad]");
    fhDPhiZP->Draw();
    fhDPhiZN->Draw("same");
    //  fhDPhiBA->Draw("same");

    
    // delta theta
    // pos/neg charge
    c2->cd(3);
    c2->GetPad(3)->SetLogy();
    fhDThetaP = (TH1F*) chist1->FindObject("fhDThetaPC");
    fhDThetaN = (TH1F*) chist1->FindObject("fhDThetaNC");
    fhDThetaP->SetLineColor(2);
    fhDThetaN->SetLineColor(4);
 

    fhDThetaP->SetXTitle("#Delta#theta [rad]");
    fhDThetaP->Draw();
    fhDThetaN->Draw("same");

    // pos/neg z
    c2->cd(4);
    c2->GetPad(4)->SetLogy();

    fhDThetaZP = (TH1F*) chist1->FindObject("fhDThetaPZ");
    fhDThetaZN = (TH1F*) chist1->FindObject("fhDThetaNZ");
    fhDThetaBA = (TH1F*) chist1->FindObject("fhDTheta_Bad");
    fhDThetaZP->SetLineColor(2);
    fhDThetaZN->SetLineColor(4);
    fhDThetaBA->SetLineColor(6);
    fhDThetaZP->SetXTitle("#Delta#theta [rad]");
    fhDThetaZP->Draw();

    fhDThetaZN->Draw("same");
    //    fhDThetaBA->Draw("same");

    
    // delta Pt
    TCanvas* c3 = new TCanvas("c3", "Cosmic: #Delta p_{T}", 10, 10, 1100, 800);
    c3->Divide(2,2);
    // pos/neg charge
    c3->cd(1);
    c3->GetPad(1)->SetLogy();
    fhDPtP = (TH1F*) chist1->FindObject("fhDPtPC");
    fhDPtN = (TH1F*) chist1->FindObject("fhDPtNC");
    fhDPtP->SetLineColor(2);
    fhDPtN->SetLineColor(4);
    
    fhDPtP->Draw();
    fhDPtN->Draw("same");

    // pos/neg z
    c3->cd(2);
    c3->GetPad(2)->SetLogy();

    fhDPtZP = (TH1F*) chist1->FindObject("fhDPtPZ");
    fhDPtZN = (TH1F*) chist1->FindObject("fhDPtNZ");
    fhDPtZP->SetLineColor(2);
    fhDPtZN->SetLineColor(4);
    fhDPtZP->Draw();
    fhDPtZN->Draw("same");

    c3->cd(3);

    fpDPtP = (TH1F*) chist1->FindObject("fpDPtPC");
    fpDPtN = (TH1F*) chist1->FindObject("fpDPtNC");
    fpDPtP->SetLineColor(2);
    fpDPtN->SetLineColor(4);
    fpDPtP->Draw();
    fpDPtN->Draw("same");

    c3->cd(4);

    fpDPtZP = (TH1F*) chist1->FindObject("fpDPtPZ");
    fpDPtZN = (TH1F*) chist1->FindObject("fpDPtNZ");
    fpDPtZP->SetLineColor(2);
    fpDPtZN->SetLineColor(4);
    fpDPtZP->Draw();
    fpDPtZN->Draw("same");
    

    // dZ
    TCanvas* c4 = new TCanvas("c4", "Cosmic: #Delta Z", 10, 10, 1100, 800);
    c4->Divide(2,2);
    // pos/neg charge
    c4->cd(1);
    c4->GetPad(1)->SetLogy();
    fhDZP = (TH1F*) chist1->FindObject("fhDZPC");
    fhDZN = (TH1F*) chist1->FindObject("fhDZNC");
    fhDZP->SetLineColor(2);
    fhDZN->SetLineColor(4);
    
    fhDZP->Draw();
    fhDZN->Draw("same");
  
    // pos/neg z
    c4->cd(2);
    c4->GetPad(2)->SetLogy();

    fhDZZP = (TH1F*) chist1->FindObject("fhDZPZ");
    fhDZZN = (TH1F*) chist1->FindObject("fhDZNZ");
    fhDZBA = (TH1F*) chist1->FindObject("fhDZ_Bad");

    fhDZZP->SetLineColor(2);
    fhDZZN->SetLineColor(4);
    fhDZBA->SetLineColor(6);
  
    fhDZZP->Draw();
    fhDZZN->Draw("same");
    fhDZBA->Draw("same");

    // dX
    // pos/neg charge
    c4->cd(3);
    c4->GetPad(3)->SetLogy();
    fhDXP = (TH1F*) chist1->FindObject("fhDXPC");
    fhDXN = (TH1F*) chist1->FindObject("fhDXNC");
    fhDXP->SetLineColor(2);
    fhDXN->SetLineColor(4);
    
    fhDXP->Draw();
    fhDXN->Draw("same");

    // pos/neg z
    c4->cd(4);
    c4->GetPad(4)->SetLogy();

    fhDXZP = (TH1F*) chist1->FindObject("fhDXPZ");
    fhDXZN = (TH1F*) chist1->FindObject("fhDXNZ");
    fhDXBA = (TH1F*) chist1->FindObject("fhDX_Bad");

    fhDXZP->SetLineColor(2);
    fhDXZN->SetLineColor(4);
    fhDXBA->SetLineColor(6);

    fhDXZP->Draw();
    fhDXZN->Draw("same");
    fhDXBA->Draw("same");

    // dY
    TCanvas* c4a = new TCanvas("c4a", "Cosmic #Delta Y", 10, 10, 1100, 800);
    c4a->Divide(2,2);
    // pos/neg charge
    c4a->cd(1);
    c4a->GetPad(1)->SetLogy();
    fhDYP = (TH1F*) chist1->FindObject("fhDYPC");
    fhDYN = (TH1F*) chist1->FindObject("fhDYNC");
    fhDYP->SetLineColor(2);
    fhDYN->SetLineColor(4);
  
    fhDYP->Draw();
    fhDYN->Draw("same");

    // pos/neg z
    c4a->cd(2);
    c4a->GetPad(2)->SetLogy();

    fhDYZP = (TH1F*) chist1->FindObject("fhDYPZ");
    fhDYZN = (TH1F*) chist1->FindObject("fhDYNZ");
    fhDYBA = (TH1F*) chist1->FindObject("fhDY_Bad");
  
    fhDYZP->SetLineColor(2);
    fhDYZN->SetLineColor(4);
    fhDYBA->SetLineColor(6);

    fhDYZP->Draw();
    fhDYZN->Draw("same");
    fhDYBA->Draw("same");

    // delta Pt
    TCanvas* c5 = new TCanvas("c5", "Cosmic #Delta p_{T} (n-sigma)", 10, 10, 1100, 800);
    c5->Divide(2,2);
    c5->cd(1);

    fpDPtSP = (TH1F*) chist1->FindObject("fpDPtSPC");
    fpDPtSN = (TH1F*) chist1->FindObject("fpDPtSNC");
    fpDPtSP->SetLineColor(2);
    fpDPtSN->SetLineColor(4);
    fpDPtSP->Draw();
    fpDPtSN->Draw("same");


    c5->cd(2);

    fpDPtSPZ = (TH1F*) chist1->FindObject("fpDPtSPZ");
    fpDPtSNZ = (TH1F*) chist1->FindObject("fpDPtSNZ");
    fpDPtSPZ->SetLineColor(2);
    fpDPtSNZ->SetLineColor(4);
    fpDPtSPZ->Draw();
    fpDPtSNZ->Draw("same");

    TCanvas* c6 = new TCanvas("c6", "Cosmic: Dz vs z", 10, 10, 1100, 800);
    fhDZvsZ = (TH2F*) chist1->FindObject("fhDZvsZ");
    fhDZvsZ->SetXTitle("z_{in} * sign(z_{in}) * sign(z_{out}) [cm]");
    fhDZvsZ->SetYTitle("#DeltaZ [cm]");

    gStyle->SetPalette(1);
    
    fhDZvsZ->Draw("colz");
  }
  if (doQASym) {
    TFile* f2 = new TFile("QAsym.root", "read");
    TList *cont = (TList*)f2->Get("QAsymHists");
    if (!cont) {
      printf("Woops... list container changed name inside QAsym.root ?\n");
      return;
    }
    // Global distributions  
    TCanvas* c7 = new TCanvas("c7", "QAsym: #Global quantities", 10, 10, 1100, 800);
    c7->Divide(3,3);
    c7->cd(1);
    gPad->SetLogy();
    fHistRECpt = (TH1F*)cont->FindObject("fHistRECpt");
    fHistRECpt->SetLineColor(2);
    fHistRECpt->Draw();
    c7->cd(2);
    gPad->SetLogy();
    fEta = (TH1F*)cont->FindObject("fEta");
    fEta->SetLineColor(4);
    fEta->Draw();
    c7->cd(3);
    fEtaPhi = (TH2F*)cont->FindObject("fEtaPhi");
    fEtaPhi->SetMarkerColor(2);
    fEtaPhi->SetMarkerStyle(kCross);
    fEtaPhi->Draw();
    c7->cd(4);
    fThetaRec = (TH1F*)cont->FindObject("fThetaRec");
    fThetaRec->SetLineColor(4);
    fThetaRec->Draw();
    c7->cd(5);
    fPhiRec = (TH1F*)cont->FindObject("fPhiRec");
    fPhiRec->SetLineColor(2);
    fPhiRec->Draw();
    c7->cd(6);
    fNumber = (TH1F*)cont->FindObject("fNumber");
    fNumber->SetLineColor(3);
    fNumber->Draw();
    c7->cd(7);
    fVx = (TH1F*)cont->FindObject("fVx");
    fVx->Draw();
    c7->cd(8);
    fVy = (TH1F*)cont->FindObject("fVy");
    fVy->Draw();
    c7->cd(9);
    fVz = (TH1F*)cont->FindObject("fVz");
    fVz->Draw();

    // Other global quantities
    TCanvas* c8 = new TCanvas("c8", "QAsym: #Global quantities #2", 10, 10, 1100, 800);
    c8->Divide(3,2);
    c8->cd(1);
    fEtaPt = (TH1F*)cont->FindObject("fEtaPt");
    fEtaPt->SetLineColor(4);
    fEtaPt->Draw();
    c8->cd(2);
    fQPt = (TH1F*)cont->FindObject("fQPt");
    fQPt->SetLineColor(2);
    fQPt->Draw();
    c8->cd(3);
    fDca = (TH1F*)cont->FindObject("fDca");
    fDca->SetLineColor(3);
    fDca->Draw();
    c8->cd(4);
    gPad->SetLogy();
    fqRec = (TH1F*)cont->FindObject("fqRec");
    fqRec->SetLineColor(2);
    fqRec->Draw();
    c8->cd(5);
    gPad->SetLogy();
    fsigmaPt = (TH1F*)cont->FindObject("fsigmaPt");
    fsigmaPt->SetLineColor(3);
    fsigmaPt->Draw();   
    // ITS DCA
    TCanvas* c9 = new TCanvas("c9", "QAsym: #ITS DCA per layer (red=positive, blue=negative", 10, 10, 1100, 800);
    c9->Divide(4,2);
    for (Int_t i=0; i<7; i++) {
      c9->cd(i+1);
      gPad->SetLogy();
      fSignDcaPos = (TH1F*)cont->FindObject(Form("fSignDcaPos%d", i));
      fSignDcaPos->SetLineColor(4);
      fSignDcaPos->Draw();
      fSignDcaNeg = (TH1F*)cont->FindObject(Form("fSignDcaNeg%d", i));
      fSignDcaNeg->SetLineColor(2);
      fSignDcaNeg->Draw("same");
    }  
    // YIELDs---------- positive and negative particles
    TCanvas* c10 = new TCanvas("c10", "QAsym: #Positive (red) and negative (blue) particles", 10, 10, 1100, 800);
    c10->Divide(3,2);
    c10->cd(1);
    gPad->SetLogy();
    fRecPtPos = (TH1F*)cont->FindObject("fRecPtPos");
    fRecPtPos->SetLineColor(2);
    fRecPtPos->Draw();
    fRecPtNeg = (TH1F*)cont->FindObject("fRecPtNeg");
    fRecPtNeg->SetLineColor(4);
    fRecPtNeg->Draw("same");
    c10->cd(2);
    gPad->SetLogy();
    fRecPhiPos = (TH1F*)cont->FindObject("fRecPhiPos");
    fRecPhiPos->SetLineColor(2);
    fRecPhiPos->Draw();
    fRecPhiNeg = (TH1F*)cont->FindObject("fRecPhiNeg");
    fRecPhiNeg->SetLineColor(4);
    fRecPhiNeg->Draw("same");
    c10->cd(3);
    gPad->SetLogy();
    fRecEtaPos = (TH1F*)cont->FindObject("fRecEtaPos");
    fRecEtaPos->SetLineColor(2);
    fRecEtaPos->Draw();
    fRecEtaNeg = (TH1F*)cont->FindObject("fRecEtaNeg");
    fRecEtaNeg->SetLineColor(4);
    fRecEtaNeg->Draw("same");
    c10->cd(4);
    gPad->SetLogy();
    fRecEtaPtPos = (TH1F*)cont->FindObject("fRecEtaPtPos");
    fRecEtaPtPos->SetLineColor(2);
    fRecEtaPtPos->Draw();
    fRecEtaPtNeg = (TH1F*)cont->FindObject("fRecEtaPtNeg");
    fRecEtaPtNeg->SetLineColor(4);
    fRecEtaPtNeg->Draw("same");
    c10->cd(5);
    gPad->SetLogy();
    fRecDcaPos = (TH1F*)cont->FindObject("fRecDcaPos");
    fRecDcaPos->SetLineColor(2);
    fRecDcaPos->Draw();
    fRecDcaNeg = (TH1F*)cont->FindObject("fRecDcaNeg");
    fRecDcaNeg->SetLineColor(4);
    fRecDcaNeg->Draw("same");
    c10->cd(6);
    gPad->SetLogy();
    fRecDPos = (TH1F*)cont->FindObject("fRecDPos");
    fRecDPos->SetLineColor(2);
    fRecDPos->Draw();
    fRecDNeg = (TH1F*)cont->FindObject("fRecDNeg");
    fRecDNeg->SetLineColor(4);
    fRecDNeg->Draw("same");
 }  
}

