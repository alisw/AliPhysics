#define aodV0_cxx
// The class definition in aodV0.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The following methods are defined in this file:
//    Begin():        called everytime a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// Root > T->Process("aodV0.C")
// Root > T->Process("aodV0.C","some options")
// Root > T->Process("aodV0.C+")
//

#include "aodV0.h"
#include <TPDGCode.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TPaveLabel.h>
#include <TPostScript.h>

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliESD.h"
#include "AliPID.h"
#include "AliAODv0.h"
#endif



aodV0::aodV0(TTree *) :
  TSelector(),
  fChain(0),
  fESD(0),
  fHistV0PerEvent(0),
  fHistMassK0(0),
  fHistMassLambda(0),
  fHistMassAntiLambda(0),
  fHistMassLambdaVsProb(0),
  fHistMassLambdaCut(0),
  fHistMassAntiLambdaCut(0),
  fHistPtVsRapK0Short(0),
  fHistPtVsRapLambda(0),
  fHistArmenterosPodolanski(0)
{
  // Constructor. Initialization of pointers
}

aodV0::~aodV0() {
  // Remove all pointers

  delete fESD;

  // histograms are in the output list and deleted when the output
  // list is deleted by the TSelector dtor
}

void aodV0::Begin(TTree *)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();
}

void aodV0::SlaveBegin(TTree * tree)
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).

  Init(tree);
  
  TString option = GetOption();

  // Create histograms on each slave server
  fHistV0PerEvent = new TH1F("h1V0PerEvent", "V^{0} candidates per event", 50, 0, 50);
  fHistV0PerEvent->GetXaxis()->SetTitle("Number of V^{0}");
  fHistV0PerEvent->GetYaxis()->SetTitle("Events");
  fHistV0PerEvent->SetOption("HE");
  fHistV0PerEvent->SetStats(11);

  fHistMassK0 = new TH1F("h1MassK0", "K^{0} candidates", 100, 0.4, 0.6);
  fHistMassK0->GetXaxis()->SetTitle("Invariant Mass_{#pi^{+}#pi^{-}}  (GeV/c^{2})");
  fHistMassK0->SetOption("HE");
  fHistMassK0->SetStats(11);

  fHistMassLambda = new TH1F("h1MassLambda", "#Lambda^{0} candidates", 75, 1.05, 1.2);
  fHistMassLambda->GetXaxis()->SetTitle("Invariant Mass_{p#pi^{-}} (GeV/c^{2})");
  fHistMassLambda->SetOption("HE");
  fHistMassLambda->SetLineStyle(1);
  fHistMassLambda->SetStats(kFALSE);

  fHistMassAntiLambda = new TH1F("h1MassAntiLambda", "#bar{#Lambda}^{0} candidates", 75, 1.05, 1.2);
  fHistMassAntiLambda->GetXaxis()->SetTitle("Invariant Mass_{#bar{p}#pi^{+}} (GeV/c^{2})");
  fHistMassAntiLambda->SetOption("HE");
  fHistMassAntiLambda->SetLineStyle(2);
  fHistMassAntiLambda->SetStats(kFALSE);

  fHistMassLambdaVsProb = new TH2F("h2MassLambdaVsProb", "#Lambda^{0} and #bar{#Lambda}^{0} candidates", 21, -2.5, 102.5, 75, 1.05, 1.2);
  fHistMassLambdaVsProb->GetXaxis()->SetTitle("prob(p) (%)");
  fHistMassLambdaVsProb->GetYaxis()->SetTitle("Invariant Mass_{p#pi} (GeV/c^{2})");
  fHistMassLambdaVsProb->SetOption("BOX");
  fHistMassLambdaVsProb->SetFillStyle(0);
  fHistMassLambdaVsProb->SetStats(kFALSE);

  fHistMassLambdaCut = new TH1F("h1MassLambdaCut", "#Lambda^{0} candidates", 75, 1.05, 1.2);
  fHistMassLambdaCut->GetXaxis()->SetTitle("Invariant Mass_{p#pi^{-}} (GeV/c^{2})");
  fHistMassLambdaCut->SetOption("E");
  fHistMassLambdaCut->SetMarkerStyle(kFullCircle);
  fHistMassLambdaCut->SetMarkerColor(kRed);
  fHistMassLambdaCut->SetStats(kFALSE);

  fHistMassAntiLambdaCut = new TH1F("h1MassAntiLambdaCut", "#bar{#Lambda}^{0} candidates", 75, 1.05, 1.2);
  fHistMassAntiLambdaCut->GetXaxis()->SetTitle("Mass_{#bar{p}#pi^{+}} (GeV/c^{2})");
  fHistMassAntiLambdaCut->SetOption("E");
  fHistMassAntiLambdaCut->SetMarkerStyle(kFullCircle);
  fHistMassAntiLambdaCut->SetMarkerColor(kRed);
  fHistMassAntiLambdaCut->SetStats(kFALSE);

  // AOD Histograms
  fHistPtVsRapK0Short = new TH2F("h2PtVsRapK0Short","K^{0}_{s} phase space",20,-1,1,20,0,10);
  fHistPtVsRapK0Short->SetOption("COL2Z");
  fHistPtVsRapK0Short->SetXTitle("Rapidity");
  fHistPtVsRapK0Short->SetYTitle("p_{t} (GeV/c)");

  fHistPtVsRapLambda = new TH2F("h2PtVsRapLambda","#Lambda phase space",20,-1,1,20,0,10);
  fHistPtVsRapLambda->SetOption("COL2Z");
  fHistPtVsRapLambda->SetXTitle("Rapidity");
  fHistPtVsRapLambda->SetYTitle("p_{t} (GeV/c)");

  fHistArmenterosPodolanski = new TH2F("h2ArmenterosPodolanski","Armenteros-Podolanski phase space",100,-1.0,1.0,50,0,0.5);
  fHistArmenterosPodolanski->SetOption("BOX");
  fHistArmenterosPodolanski->SetXTitle("#alpha");
  fHistArmenterosPodolanski->SetYTitle("p_{t} arm");
}

void aodV0::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses of the tree
  // will be set. It is normaly not necessary to make changes to the
  // generated code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running with PROOF.

  // Set branch addresses
  if (tree == 0) return;
  fChain = tree;

  fChain->SetBranchAddress("ESD",&fESD);
}

Bool_t aodV0::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. Typically here the branch pointers
  // will be retrieved. It is normaly not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed.
  
  return kTRUE;
}


Bool_t aodV0::Process(Long64_t entry)
{
  // The Process() function is called for each entry in the tree (or possibly
  // keyed object in the case of PROOF) to be processed. The entry argument
  // specifies which entry in the currently loaded tree is to be processed.
  // It can be passed to either TTree::GetEntry() or TBranch::GetEntry()
  // to read either all or the required parts of the data. When processing
  // keyed objects with PROOF, the object is already loaded and is available
  // via the fObject pointer.
  //
  // This function should contain the "body" of the analysis. It can contain
  // simple or elaborate selection criteria, run algorithms on the data
  // of the event and typically fill histograms.

  // WARNING when a selector is used with a TChain, you must use
  //  the pointer to the current TTree to call GetEntry(entry).
  //  The entry is always the local entry number in the current tree.
  //  Assuming that fChain is the pointer to the TChain being processed,
  //  use fChain->GetEntry(entry).

  fChain->GetTree()->GetEntry(entry);

  if (!fESD) return kFALSE;

  Int_t mBoDebug = 1;

  // Primary vertex and Position
  const AliESDVertex *lPrimaryVertex = fESD->GetVertex();
  Double_t tPositionPrimaryVertex[3]; lPrimaryVertex->GetXYZ(tPositionPrimaryVertex);
  if (mBoDebug) printf("***BoInfo: Primary Vertex position x=%.3f, y=%.3f, z=%.3f \n",tPositionPrimaryVertex[0],tPositionPrimaryVertex[1],tPositionPrimaryVertex[2]);

  // Following is commented since ntracks is not used (remove warning)
  //  Int_t ntracks = fESD->GetNumberOfTracks();
  Int_t nv0s = fESD->GetNumberOfV0s();
  if (mBoDebug) printf("***BoInfo: Number of V0s =%d \n",nv0s);
  fHistV0PerEvent->Fill(nv0s);

  Double_t mass;
  Double_t pPrior[5] = {0, 0, 0.85, 0.1, 0.05};
  Double_t p[5];
  Double_t pSum, pNorm;

  AliAODv0 *myAODv0 = new AliAODv0();

  // variable for AOD v0
  Double_t dcaPosToPrimVertex = 0, dcaNegToPrimVertex = 0;
  Double_t dcaV0Daughters     = 0, decayLengthV0      = 0;
  Double_t alphaV0            = 0, ptArmV0            = 0;
  Double_t rapK0Short         = 0, rapLambda          = 0;
  Double_t pt                 = 0;

  for (Int_t iV0 = 0; iV0 < nv0s; iV0++) {
    AliESDv0* v0 = fESD->GetV0(iV0);
    if (!v0) continue;
    myAODv0->ResetV0();
    myAODv0->Fill(v0,fESD);

    dcaPosToPrimVertex = myAODv0->DcaPosToPrimVertex();
    dcaNegToPrimVertex = myAODv0->DcaNegToPrimVertex();
    dcaV0Daughters     = myAODv0->DcaV0Daughters();
    decayLengthV0      = myAODv0->DecayLengthV0(tPositionPrimaryVertex);
    alphaV0            = myAODv0->AlphaV0();
    ptArmV0            = myAODv0->PtArmV0();
    rapK0Short         = myAODv0->RapK0Short();
    rapLambda          = myAODv0->RapLambda();
    pt                 = TMath::Sqrt(myAODv0->Ptot2V0());

    fHistPtVsRapK0Short->Fill(rapK0Short,pt);
    fHistPtVsRapLambda->Fill(rapLambda,pt);
    fHistArmenterosPodolanski->Fill(alphaV0,ptArmV0);

    fHistMassK0->Fill(v0->GetEffMass());

    v0->ChangeMassHypothesis(kLambda0);
    mass = v0->GetEffMass();
    fESD->GetTrack(v0->GetPindex())->GetESDpid(p);
    pSum = 0;
    for (Int_t i = 0; i < AliPID::kSPECIES; i++) pSum += p[i] * pPrior[i];
    if (pSum <= 0) pSum = 1.;
    pNorm = p[AliPID::kProton] * pPrior[AliPID::kProton] / pSum;
    fHistMassLambdaVsProb->Fill(100.*pNorm, mass); 
    if (pNorm > 0.1) fHistMassLambda->Fill(mass);
    
    v0->ChangeMassHypothesis(kLambda0Bar);
    mass = v0->GetEffMass();
    fESD->GetTrack(v0->GetNindex())->GetESDpid(p);
    pSum = 0;
    for (Int_t i = 0; i < AliPID::kSPECIES; i++) pSum += p[i] * pPrior[i];
    if (pSum <= 0) pSum = 1.;
    pNorm = p[AliPID::kProton] * pPrior[AliPID::kProton] / pSum;
    fHistMassLambdaVsProb->Fill(100.*pNorm, mass); 
    if (pNorm > 0.1) fHistMassAntiLambda->Fill(mass);
   
  }//loop v0s
  
  return kTRUE;
}



void aodV0::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.

  // Add the histograms to the output on each slave server
  fOutput->Add(fHistV0PerEvent);
  fOutput->Add(fHistMassK0);
  fOutput->Add(fHistMassLambda);
  fOutput->Add(fHistMassAntiLambda);
  fOutput->Add(fHistMassLambdaVsProb);
  fOutput->Add(fHistMassLambdaCut);
  fOutput->Add(fHistMassAntiLambdaCut);
  // Add the AOD histograms to the output on each slave server
  fOutput->Add(fHistPtVsRapK0Short);
  fOutput->Add(fHistPtVsRapLambda);
  fOutput->Add(fHistArmenterosPodolanski);
}

void aodV0::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

  fHistV0PerEvent = dynamic_cast<TH1F*>(fOutput->FindObject("h1V0PerEvent"));
  fHistMassK0 = dynamic_cast<TH1F*>(fOutput->FindObject("h1MassK0"));
  fHistMassLambda = dynamic_cast<TH1F*>(fOutput->FindObject("h1MassLambda"));
  fHistMassAntiLambda = dynamic_cast<TH1F*>(fOutput->FindObject("h1MassAntiLambda"));
  fHistMassLambdaVsProb = dynamic_cast<TH2F*>(fOutput->FindObject("h2MassLambdaVsProb"));
  fHistMassLambdaCut = dynamic_cast<TH1F*>(fOutput->FindObject("h1MassLambdaCut"));
  fHistMassAntiLambdaCut = dynamic_cast<TH1F*>(fOutput->FindObject("h1MassAntiLambdaCut"));

  // AOD Histograms
  fHistPtVsRapK0Short = dynamic_cast<TH2F*>(fOutput->FindObject("h2PtVsRapK0Short"));
  fHistPtVsRapLambda = dynamic_cast<TH2F*>(fOutput->FindObject("h2PtVsRapLambda"));
  fHistArmenterosPodolanski = dynamic_cast<TH2F*>(fOutput->FindObject("h2ArmenterosPodolanski"));

  TFile* file = TFile::Open("aodV0.root", "RECREATE");
  fHistV0PerEvent->Write();
  fHistMassK0->Write();
  fHistMassLambda->Write();
  fHistMassAntiLambda->Write();
  fHistMassLambdaVsProb->Write();
  fHistMassLambdaCut->Write();
  fHistMassAntiLambdaCut->Write();
  // AOD Histograms
  fHistPtVsRapK0Short->Write();      
  fHistPtVsRapLambda->Write();       
  fHistArmenterosPodolanski->Write();

  file->Close();
  delete file;

  // Define postscript


  if (!gROOT->IsBatch()) 
    {
      TCanvas *canOnline = new TCanvas("canOnline","V0s",10,10,610,610);
      canOnline->SetFillColor(10);
      canOnline->SetHighLightColor(10);
      canOnline->Divide(2,1);
 
      canOnline->cd(1);
      fHistMassK0->DrawCopy("E");
      fHistMassK0->Fit("gaus","q","",0.49,0.51);
      TVirtualPad* pad = (TVirtualPad*)canOnline->cd(2);
      pad->Divide(1,2);
      pad->cd(1);
      fHistMassLambda->DrawCopy("E");
      pad->cd(2);
      fHistMassAntiLambda->DrawCopy("E");
    }

  // Following define a postscript output
  int  myAliceRed        = TColor::GetColor(192,0,0);
  TDatime *currentDatime = new TDatime();
  Char_t info[125];
  gStyle->SetPaperSize(20,26);
//   gStyle->SetOptStat(111111);
  gStyle->SetOptStat(11);
  gStyle->SetPalette(1,0);
  TCanvas     *c0 = new TCanvas("c0","AOD V0 Analysis",600,800);
  c0->Range(0,0,20,24);

  TPostScript *ps = new TPostScript("aodV0.ps");

  TPaveLabel  *t0 = new TPaveLabel(4,23.0,16,23.8,"AOD V0 Analysis","br");
  t0->SetFillColor(0);
  t0->SetBorderSize(2);
  t0->SetTextColor(myAliceRed);

  TPaveLabel  *d0 = new TPaveLabel(16.5,23.0,19.5,23.6,currentDatime->AsSQLString(),"br");
  d0->SetFillColor(0);
  d0->SetTextSize(0.4);
  d0->SetBorderSize(0);
  d0->SetTextColor(1);

  TPaveText *ptt1 = new TPaveText(1.5,20.6,8.0,21.4);
  ptt1->SetFillColor(0);
  ptt1->SetBorderSize(2);
  ptt1->SetTextSize(0.03);
  ptt1->SetTextColor(myAliceRed);
  ptt1->SetTextAlign(12);
  ptt1->AddText("Global Information:");

  TPaveText *ptx1 = new TPaveText(2,16,18,21);
  ptx1->SetFillColor(0);
  ptx1->SetBorderSize(1);
  ptx1->SetTextSize(0.02);
  ptx1->SetTextAlign(13);
  ptx1->AddText("");
  TString USER("$USER");
  TString HOST("$HOST");
  TString ROOTSYS("$ROOTSYS");
  TString ALICE_LEVEL("$ALICE_LEVEL");
  gSystem->ExpandPathName(USER);
  gSystem->ExpandPathName(HOST);
  gSystem->ExpandPathName(ROOTSYS);
  gSystem->ExpandPathName(ALICE_LEVEL);
  sprintf(info,"Analysis from %s on %s",USER.Data(),HOST.Data());
  ptx1->AddText(info);
  sprintf(info,"Root Version: %s",ROOTSYS.Data());
  ptx1->AddText(info);
  sprintf(info,"Alice Version: %s",ALICE_LEVEL.Data());
  ptx1->AddText(info);
  sprintf(info,"Number of scanned events : %.0f",fHistV0PerEvent->GetEntries());
  ptx1->AddText(info);

  ps->NewPage();
  t0->Draw();
  d0->Draw();

  ptx1->Draw();
  ptt1->Draw();
  c0->Update();

  delete ptx1;
  delete ptt1;

  delete t0;  
  delete c0;  

  // Rich Primary Vertex 
  ps->NewPage();
  TCanvas* c1 = new TCanvas("c1","Phase Space Page",500,690);
  c1->Update();
  c1->cd();
  c1->Range(0,0,25,18);

  TPaveLabel *ptitle1 = new TPaveLabel(5,17,20,18,"Phase Space Info","br");
  ptitle1->SetFillColor(18);
  ptitle1->SetTextFont(32);
  ptitle1->SetTextSize(0.8);
  ptitle1->SetTextColor(myAliceRed);
  ptitle1->Draw();

  TPad    *c1Pad11 = new TPad("c1Pad11","box",0.01,0.50,0.33,0.90,0);
  TPad    *c1Pad12 = new TPad("c1Pad12","box",0.34,0.50,0.66,0.90,0);
  TPad    *c1Pad13 = new TPad("c1Pad13","box",0.67,0.50,0.99,0.90,0);
  TPad    *c1Pad21 = new TPad("c1Pad21","box",0.01,0.05,0.33,0.45,0);
  TPad    *c1Pad22 = new TPad("c1Pad22","box",0.34,0.05,0.66,0.45,0);
  TPad    *c1Pad23 = new TPad("c1Pad23","box",0.67,0.05,0.99,0.45,0);
  c1Pad11->SetLeftMargin(0.15); c1Pad11->SetBottomMargin(0.15);
  c1Pad11->Draw();
  c1Pad12->SetLeftMargin(0.15); c1Pad12->SetBottomMargin(0.15);
  c1Pad12->Draw();
  c1Pad13->SetLeftMargin(0.15); c1Pad13->SetBottomMargin(0.15);
  c1Pad13->Draw();
  c1Pad21->SetLeftMargin(0.15); c1Pad21->SetBottomMargin(0.15);
  c1Pad21->Draw();
  c1Pad22->SetLeftMargin(0.15); c1Pad22->SetBottomMargin(0.15);
  c1Pad22->Draw();
  c1Pad23->SetLeftMargin(0.15); c1Pad23->SetBottomMargin(0.15);
  c1Pad23->Draw();

  c1Pad11->cd();
  fHistV0PerEvent->Draw("H");
  c1Pad12->cd();
  fHistMassK0->Draw("HE");
  c1Pad13->cd();
  fHistMassLambda->Draw("HE");
  fHistMassAntiLambda->SetLineStyle(2);
  fHistMassAntiLambda->Draw("HESAME");
  c1Pad13->Update();

  gPad->Update();
  c1Pad21->cd();
  fHistPtVsRapK0Short->Draw("COL2Z");
  c1Pad22->cd();
  fHistPtVsRapLambda->Draw("COL2Z");
  c1Pad23->cd();
  fHistArmenterosPodolanski->Draw("BOX");
  c1->Update();
  ps->NewPage();

  delete c1Pad11;
  delete c1Pad12;
  delete c1Pad13;
  delete c1Pad21;
  delete c1Pad22;
  delete c1Pad23;
  delete c1;

  ps->Close();
}
