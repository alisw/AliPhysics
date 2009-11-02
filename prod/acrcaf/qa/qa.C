void qa(Int_t runNumber) {
  TStopwatch timer;
  timer.Start();

  gEnv->SetValue("XSec.GSI.DelegProxy","2");
  // Select ROOT version
  TProof::Mgr("aliprod@alicecaf")->SetROOTVersion("v5-24-00b-caf");
  // Login to CAF
  TProof::Open("aliprod@alicecaf");

  // Enable AliRoot
  gProof->UploadPackage("/afs/cern.ch/alice/caf/sw/ALICE/PARs/v4-17-Release.rec/AF-v4-17-rec.par");
  gProof->EnablePackage("AF-v4-17-rec.par");

  // Enable analysis libs
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gProof->Exec("gSystem->Load(\"libANALYSIS.so\");",kTRUE);
  gProof->Exec("gSystem->Load(\"libANALYSISalice.so\");",kTRUE);

  gProof->Load(Form("%s/PWG1/cosmic/AliAnalysisTaskCosmic.cxx++g",
		    gSystem->Getenv("ALICE_ROOT")));

  //____________________________________________//
  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("AnalysisQAManager");
  AliESDInputHandler* esdH = new AliESDInputHandler();
  mgr->SetInputEventHandler(esdH);  
  mgr->SetDebugLevel(10);

  //____________________________________________//
  // 1st Cosmic task
  AliAnalysisTaskCosmic *task1 = new AliAnalysisTaskCosmic("TaskCosmic");
  mgr->AddTask(task1);

  // Create containers for input/output
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("chist1",TList::Class(),AliAnalysisManager::kOutputContainer,
							    Form("run%d.root",runNumber));

  //____________________________________________//
  mgr->ConnectInput  (task1,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput (task1,  1, coutput1);
  
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  mgr->StartAnalysis("proof",
		     Form("/ALIREC/aliprod/run%d",runNumber));

  timer.Stop();
  timer.Print();

  plot(runNumber);
}

void plot(Int_t runNumber)
{
  TFile* f = new TFile(Form("run%d.root",runNumber), "read");

  // pt, phi ....
  TCanvas* c1 = new TCanvas("c1", "pt,eta,phi", 10, 10, 1100, 800);
  c1->Divide(2,2);
  c1->cd(1);
  c1->GetPad(1)->SetLogy();
    
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
  TCanvas* c2 = new TCanvas("c2", "#Delta#phi", 10, 10, 1100, 800);
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
  TCanvas* c3 = new TCanvas("c3", "#Delta p_{T}", 10, 10, 1100, 800);
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
  TCanvas* c4 = new TCanvas("c4", "#Delta Z", 10, 10, 1100, 800);
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
  TCanvas* c4a = new TCanvas("c4a", "#Delta Y", 10, 10, 1100, 800);
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
  TCanvas* c5 = new TCanvas("c5", "#Delta p_{T} (n-sigma)", 10, 10, 1100, 800);
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

  TCanvas* c6 = new TCanvas("c6", "Dz vs z", 10, 10, 1100, 800);
  fhDZvsZ = (TH2F*) chist1->FindObject("fhDZvsZ");
  fhDZvsZ->SetXTitle("z_{in} * sign(z_{in}) * sign(z_{out}) [cm]");
  fhDZvsZ->SetYTitle("#DeltaZ [cm]");

  gStyle->SetPalette(1);
    
  fhDZvsZ->Draw("colz");
}

