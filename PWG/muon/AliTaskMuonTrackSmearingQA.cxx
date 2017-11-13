/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
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

// ROOT includes
#include "TROOT.h"
#include "TRootIOCtor.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"

// STEER includes
#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliMCEvent.h"

// ANALYSIS includes
#include "AliAnalysisManager.h"
#include "AliAnalysisMuonUtility.h"

#include "AliTaskMuonTrackSmearingQA.h"

/// \cond CLASSIMP
ClassImp(AliTaskMuonTrackSmearingQA) // Class implementation in ROOT context
/// \endcond


//________________________________________________________________________
AliTaskMuonTrackSmearingQA::AliTaskMuonTrackSmearingQA(TRootIOCtor* ioCtor) :
AliAnalysisTaskSE(),
fGenList(0x0),
fRecList(0x0),
fResList(0x0),
fMuonTrackCuts(0x0),
fcGen(0x0),
fcRec(0x0),
fcRat(0x0),
fcResVsP(0x0)
{
  /// Constructor for streaming
}

//________________________________________________________________________
AliTaskMuonTrackSmearingQA::AliTaskMuonTrackSmearingQA(const char *name, Int_t chosenFunc) :
AliAnalysisTaskSE(name),
fGenList(0x0),
fRecList(0x0),
fResList(0x0),
fMuonTrackCuts(0x0),
fcGen(0x0),
fcRec(0x0),
fcRat(0x0),
fcResVsP(0x0)
{
  /// Constructor
  
  // Output slot #1 writes into a TObjArray container
  DefineOutput(1,TObjArray::Class());
  // Output slot #1 writes into a TObjArray container
  DefineOutput(2,TObjArray::Class());
  // Output slot #1 writes into a TObjArray container
  DefineOutput(3,TObjArray::Class());
  
}


//________________________________________________________________________
AliTaskMuonTrackSmearingQA::~AliTaskMuonTrackSmearingQA()
{
  /// Destructor
  if (!AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
    delete fGenList;
    delete fRecList;
    delete fResList;
  }
  delete fMuonTrackCuts;
  delete fcGen;
  delete fcRec;
  delete fcRat;
  delete fcResVsP;
}

//___________________________________________________________________________
void AliTaskMuonTrackSmearingQA::UserCreateOutputObjects()
{
  /// Create control plots
  
  // pT/eta/phi ranges for generation
  Double_t pTRange[2] = {1.,15.};
  Double_t npTBinPerGeV = 2.;
  Double_t etaRange[2] = {-4.2,-2.3};
  Double_t nEtaBinPerUnit = 50.;
  Double_t phiRange[2] = {0.,360.};
  Double_t nPhiBinPerUnit = 0.25;
  
  // --- histograms for generated and smeared tracks ---
  
  fGenList = new TObjArray(2000);
  fGenList->SetOwner();
  fRecList = new TObjArray(2000);
  fRecList->SetOwner();
  
  Int_t npTBins = (Int_t) (npTBinPerGeV*pTRange[1]);
  TH1F *hpTGen = new TH1F("hpTGen","hpTGen",2*npTBins,0.,2.*pTRange[1]);
  hpTGen->Sumw2();
  fGenList->AddAtAndExpand(hpTGen, kPt);
  TH1F *hpTRec = new TH1F("hpTRec","hpTRec",2*npTBins,0.,2.*pTRange[1]);
  hpTRec->Sumw2();
  fRecList->AddAtAndExpand(hpTRec, kPt);
  
  TH1F *hpGen = new TH1F("hpGen","hpGen",5*npTBins,0.,50.*pTRange[1]);
  hpGen->Sumw2();
  fGenList->AddAtAndExpand(hpGen, kP);
  TH1F *hpRec = new TH1F("hpRec","hpRec",5*npTBins,0.,50.*pTRange[1]);
  hpRec->Sumw2();
  fRecList->AddAtAndExpand(hpRec, kP);
  
  Int_t nEtaBins = (Int_t) (nEtaBinPerUnit*(etaRange[1]-etaRange[0]));
  TH1F *hetaGen = new TH1F("hetaGen","hetaGen",nEtaBins+0.4*nEtaBinPerUnit,etaRange[0]-0.2,etaRange[1]+0.2);
  hetaGen->Sumw2();
  fGenList->AddAtAndExpand(hetaGen, kEta);
  TH1F *hetaRec = new TH1F("hetaRec","hetaRec",nEtaBins+0.4*nEtaBinPerUnit,etaRange[0]-0.2,etaRange[1]+0.2);
  hetaRec->Sumw2();
  fRecList->AddAtAndExpand(hetaRec, kEta);
  
  TH1F *hyGen = new TH1F("hyGen","hyGen",nEtaBins+0.4*nEtaBinPerUnit,etaRange[0]-0.2,etaRange[1]+0.2);
  hyGen->Sumw2();
  fGenList->AddAtAndExpand(hyGen, kY);
  TH1F *hyRec = new TH1F("hyRec","hyRec",nEtaBins+0.4*nEtaBinPerUnit,etaRange[0]-0.2,etaRange[1]+0.2);
  hyRec->Sumw2();
  fRecList->AddAtAndExpand(hyRec, kY);
  
  Int_t nPhiBins = (Int_t) (nPhiBinPerUnit*(phiRange[1]-phiRange[0]));
  TH1F *hphiGen = new TH1F("hphiGen","hphiGen",nPhiBins,phiRange[0],phiRange[1]);
  hphiGen->Sumw2();
  fGenList->AddAtAndExpand(hphiGen, kPhi);
  TH1F *hphiRec = new TH1F("hphiRec","hphiRec",nPhiBins,phiRange[0],phiRange[1]);
  hphiRec->Sumw2();
  fRecList->AddAtAndExpand(hphiRec, kPhi);
  
  // --- histograms for track resolution ---
  
  fResList = new TObjArray(2000);
  fResList->SetOwner();
  
  Double_t deltaPAtVtxEdges[2];
  deltaPAtVtxEdges[0] = -20. - 0.0005*(20.*pTRange[1])*(20.*pTRange[1]);
  deltaPAtVtxEdges[1] = 5. + 0.0005*(20.*pTRange[1])*(20.*pTRange[1]);
  Int_t deltaPAtVtxNBins = (Int_t)(10.*(deltaPAtVtxEdges[1]-deltaPAtVtxEdges[0]));
  if (deltaPAtVtxNBins > 1000) deltaPAtVtxNBins = (deltaPAtVtxNBins+200)/200*200;
  else if (deltaPAtVtxNBins > 100) deltaPAtVtxNBins = (deltaPAtVtxNBins+20)/20*20;
  TH2F *hResPAtVtxVsPIn02degMC = new TH2F("hResPAtVtxVsPIn02deg","#Delta_{p} at vertex versus p in [0,2[ deg MC;p (GeV/c);#Delta_{p} (GeV/c)",2*npTBins,0.,20.*pTRange[1],deltaPAtVtxNBins,deltaPAtVtxEdges[0],deltaPAtVtxEdges[1]);
  hResPAtVtxVsPIn02degMC->SetDirectory(0);
  hResPAtVtxVsPIn02degMC->Sumw2();
  fResList->AddAtAndExpand(hResPAtVtxVsPIn02degMC, kResPAtVtxVsPIn02degMC);
  TH2F *hResPAtVtxVsPIn23deg = new TH2F("hResPAtVtxVsPIn23deg","#Delta_{p} at vertex versus p in [2,3[ deg;p (GeV/c);#Delta_{p} (GeV/c)",2*npTBins,0.,20.*pTRange[1],deltaPAtVtxNBins,deltaPAtVtxEdges[0],deltaPAtVtxEdges[1]);
  hResPAtVtxVsPIn23deg->SetDirectory(0);
  hResPAtVtxVsPIn23deg->Sumw2();
  fResList->AddAtAndExpand(hResPAtVtxVsPIn23deg, kResPAtVtxVsPIn23deg);
  TH2F *hResPAtVtxVsPIn310deg = new TH2F("hResPAtVtxVsPIn310deg","#Delta_{p} at vertex versus p in [3,10[ deg;p (GeV/c);#Delta_{p} (GeV/c)",2*npTBins,0.,20.*pTRange[1],deltaPAtVtxNBins,deltaPAtVtxEdges[0],deltaPAtVtxEdges[1]);
  hResPAtVtxVsPIn310deg->SetDirectory(0);
  hResPAtVtxVsPIn310deg->Sumw2();
  fResList->AddAtAndExpand(hResPAtVtxVsPIn310deg, kResPAtVtxVsPIn310deg);

  TH2F *hResPtAtVtxVsPt = new TH2F("hResPtAtVtxVsPt","#Delta_{p_{t}} at vertex versus p_{t};p_{t} (GeV/c);#Delta_{p_{t}} (GeV/c)",npTBins,0.,pTRange[1],deltaPAtVtxNBins,deltaPAtVtxEdges[0]/10.,deltaPAtVtxEdges[1]/10.);
  hResPtAtVtxVsPt->SetDirectory(0);
  hResPtAtVtxVsPt->Sumw2();
  fResList->AddAtAndExpand(hResPtAtVtxVsPt, kResPtAtVtxVsPt);

  Double_t deltaSlopeAtVtxEdges[2] = {-0.05, 0.05};
  Int_t deltaSlopeAtVtxNBins = 1000;
  TH2F *hResSlopeXAtVtxVsP = new TH2F("hResSlopeXAtVtxVsP","#Delta_{slope_{X}} at vertex versus p;p (GeV/c);#Delta_{slope_{X}}",2*npTBins,0.,20.*pTRange[1], deltaSlopeAtVtxNBins, deltaSlopeAtVtxEdges[0], deltaSlopeAtVtxEdges[1]);
  hResSlopeXAtVtxVsP->SetDirectory(0);
  hResSlopeXAtVtxVsP->Sumw2();
  fResList->AddAtAndExpand(hResSlopeXAtVtxVsP, kResSlopeXAtVtxVsP);
  TH2F *hResSlopeYAtVtxVsP = new TH2F("hResSlopeYAtVtxVsP","#Delta_{slope_{Y}} at vertex versus p;p (GeV/c);#Delta_{slope_{Y}}",2*npTBins,0.,20.*pTRange[1], deltaSlopeAtVtxNBins, deltaSlopeAtVtxEdges[0], deltaSlopeAtVtxEdges[1]);
  hResSlopeYAtVtxVsP->SetDirectory(0);
  hResSlopeYAtVtxVsP->Sumw2();
  fResList->AddAtAndExpand(hResSlopeYAtVtxVsP, kResSlopeYAtVtxVsP);

  Double_t deltaEtaAtVtxEdges[2] = {-0.5, 0.5};
  Int_t deltaEtaAtVtxNBins = 1000;
  TH2F *hResEtaAtVtxVsP = new TH2F("hResEtaAtVtxVsP","#Delta_{eta} at vertex versus p;p (GeV/c);#Delta_{eta}",2*npTBins,0.,20.*pTRange[1], deltaEtaAtVtxNBins, deltaEtaAtVtxEdges[0], deltaEtaAtVtxEdges[1]);
  hResEtaAtVtxVsP->SetDirectory(0);
  hResEtaAtVtxVsP->Sumw2();
  fResList->AddAtAndExpand(hResEtaAtVtxVsP, kResEtaAtVtxVsP);

  TH2F *hResPhiAtVtxVsP = new TH2F("hResPhiAtVtxVsP","#Delta_{phi} at vertex versus p;p (GeV/c);#Delta_{phi}",2*npTBins,0.,20.*pTRange[1], deltaEtaAtVtxNBins, deltaEtaAtVtxEdges[0], deltaEtaAtVtxEdges[1]);
  hResPhiAtVtxVsP->SetDirectory(0);
  hResPhiAtVtxVsP->Sumw2();
  fResList->AddAtAndExpand(hResPhiAtVtxVsP, kResPhiAtVtxVsP);

  // Post data at least once per task to ensure data synchronisation (required for merging)
  PostData(1, fGenList);
  PostData(2, fRecList);
  PostData(3, fResList);
  
}

//________________________________________________________________________
void AliTaskMuonTrackSmearingQA::UserExec(Option_t *)
{
  /// Fill control plots
  
  // get event
  AliVEvent *evt = InputEvent();
  if (!dynamic_cast<AliESDEvent*>(evt) && !dynamic_cast<AliAODEvent*>(evt)) return;
  
  // get MC event
  AliMCEvent *mcEvt = MCEvent();
  if (!mcEvt) return;
  
  Int_t nTracks = AliAnalysisMuonUtility::GetNTracks(evt);
  for (Int_t itrack = 0; itrack < nTracks; itrack++) {
    
    // get smeared track
    AliVParticle *track = AliAnalysisMuonUtility::GetTrack(itrack,evt);
    if (!AliAnalysisMuonUtility::IsMuonTrack(track)) continue;
    if (fMuonTrackCuts && !fMuonTrackCuts->IsSelected(track)) continue;
    Double_t pRec = track->P();
    Double_t pTRec = track->Pt();
    Double_t etaRec = track->Eta();
    Double_t phiRec = track->Phi();
    Double_t thetaAbs = AliAnalysisMuonUtility::GetThetaAbsDeg(track);
    Double_t slopeXRec = track->Px() / track->Pz();
    Double_t slopeYRec = track->Py() / track->Pz();
    
    // get generated track
    if (track->GetLabel() < 0) continue;
    AliVParticle *mctrack = mcEvt->GetTrack(track->GetLabel());
    if (!mctrack || TMath::Abs(mctrack->PdgCode()) != 13 || !AliAnalysisMuonUtility::IsPrimary(mctrack,mcEvt)) continue;
    Double_t pGen = mctrack->P();
    Double_t pTGen = mctrack->Pt();
    Double_t etaGen = mctrack->Eta();
    Double_t phiGen = mctrack->Phi();
    Double_t thetaGen = 180. - mctrack->Theta() * TMath::RadToDeg();
    Double_t slopeXGen = mctrack->Px() / mctrack->Pz();
    Double_t slopeYGen = mctrack->Py() / mctrack->Pz();
    
    // fill histograms for generated and smeared tracks
    ((TH1F*)fGenList->UncheckedAt(kPt))->Fill(pTGen);
    ((TH1F*)fRecList->UncheckedAt(kPt))->Fill(pTRec);
    ((TH1F*)fGenList->UncheckedAt(kP))->Fill(pGen);
    ((TH1F*)fRecList->UncheckedAt(kP))->Fill(pRec);
    ((TH1F*)fGenList->UncheckedAt(kEta))->Fill(etaGen);
    ((TH1F*)fRecList->UncheckedAt(kEta))->Fill(etaRec);
    ((TH1F*)fGenList->UncheckedAt(kY))->Fill(mctrack->Y());
    ((TH1F*)fRecList->UncheckedAt(kY))->Fill(track->Y());
    ((TH1F*)fGenList->UncheckedAt(kPhi))->Fill(phiGen * TMath::RadToDeg());
    ((TH1F*)fRecList->UncheckedAt(kPhi))->Fill(phiRec * TMath::RadToDeg());
    
    // fill histograms for track resolution
    if (thetaGen < 2.) ((TH2F*)fResList->UncheckedAt(kResPAtVtxVsPIn02degMC))->Fill(pGen, pRec-pGen);
    if (thetaAbs > 2. && thetaAbs < 3.) ((TH2F*)fResList->UncheckedAt(kResPAtVtxVsPIn23deg))->Fill(pGen, pRec-pGen);
    else if (thetaAbs >= 3. && thetaAbs < 10.) ((TH2F*)fResList->UncheckedAt(kResPAtVtxVsPIn310deg))->Fill(pGen, pRec-pGen);
    ((TH2F*)fResList->UncheckedAt(kResPtAtVtxVsPt))->Fill(pTGen, pTRec-pTGen);
    ((TH2F*)fResList->UncheckedAt(kResSlopeXAtVtxVsP))->Fill(pGen, slopeXRec-slopeXGen);
    ((TH2F*)fResList->UncheckedAt(kResSlopeYAtVtxVsP))->Fill(pGen, slopeYRec-slopeYGen);
    ((TH2F*)fResList->UncheckedAt(kResEtaAtVtxVsP))->Fill(pGen, etaRec-etaGen);
    Double_t dPhi = phiRec-phiGen;
    if (dPhi < -TMath::Pi()) dPhi += 2.*TMath::Pi();
    else if (dPhi > TMath::Pi()) dPhi -= 2.*TMath::Pi();
    ((TH2F*)fResList->UncheckedAt(kResPhiAtVtxVsP))->Fill(pGen, dPhi);
    
  }
  
  // Post final data. It will be written to a file with option "RECREATE"
  PostData(1, fGenList);
  PostData(2, fRecList);
  PostData(3, fResList);
  
}

//________________________________________________________________________
void AliTaskMuonTrackSmearingQA::Terminate(Option_t *)
{
  /// Draw control plots
  
  fGenList = dynamic_cast<TObjArray*>(GetOutputData(1));
  fRecList = dynamic_cast<TObjArray*>(GetOutputData(2));
  fResList = dynamic_cast<TObjArray*>(GetOutputData(3));
  if (!fGenList || !fRecList || !fResList) return;

  // draw histograms for generated tracks
  fcGen = new TCanvas("cGen","cGen",10,10,600,600);
  fcGen->Divide(2,2);
  gROOT->SetSelectedPad(fcGen->cd(1));
  gPad->SetLogy();
  ((TH1F*)fGenList->UncheckedAt(kPt))->Draw();
  gROOT->SetSelectedPad(fcGen->cd(2));
  gPad->SetLogy();
  ((TH1F*)fGenList->UncheckedAt(kEta))->Draw();
  ((TH1F*)fGenList->UncheckedAt(kY))->SetLineColor(8);
  ((TH1F*)fGenList->UncheckedAt(kY))->Draw("same");
  gROOT->SetSelectedPad(fcGen->cd(3));
  gPad->SetLogy();
  ((TH1F*)fGenList->UncheckedAt(kP))->Draw();
  gROOT->SetSelectedPad(fcGen->cd(4));
  gPad->SetLogy();
  ((TH1F*)fGenList->UncheckedAt(kPhi))->Draw();
  
  // draw histograms for smeared tracks
  fcRec = new TCanvas("cRec","cRec",10,10,600,600);
  fcRec->Divide(2,2);
  gROOT->SetSelectedPad(fcRec->cd(1));
  gPad->SetLogy();
  ((TH1F*)fRecList->UncheckedAt(kPt))->Draw();
  gROOT->SetSelectedPad(fcRec->cd(2));
  gPad->SetLogy();
  ((TH1F*)fRecList->UncheckedAt(kEta))->Draw();
  ((TH1F*)fRecList->UncheckedAt(kY))->SetLineColor(8);
  ((TH1F*)fRecList->UncheckedAt(kY))->Draw("same");
  gROOT->SetSelectedPad(fcRec->cd(3));
  gPad->SetLogy();
  ((TH1F*)fRecList->UncheckedAt(kP))->Draw();
  gROOT->SetSelectedPad(fcRec->cd(4));
  gPad->SetLogy();
  ((TH1F*)fRecList->UncheckedAt(kPhi))->Draw();
  
  // draw smeared over generated ratios
  fcRat = new TCanvas("cRat","cRat",10,10,600,600);
  fcRat->Divide(2,2);
  gROOT->SetSelectedPad(fcRat->cd(1));
  TH1 *hpTRat = static_cast<TH1*>(fRecList->UncheckedAt(kPt)->Clone());
  hpTRat->SetNameTitle("hpTRat","hpTRat");
  hpTRat->Divide(((TH1F*)fGenList->UncheckedAt(kPt)));
  hpTRat->Draw();
  gROOT->SetSelectedPad(fcRat->cd(2));
  TH1 *hetaRat = static_cast<TH1*>(fRecList->UncheckedAt(kEta)->Clone());
  hetaRat->SetNameTitle("hetaRat","hetaRat");
  hetaRat->Divide(((TH1F*)fGenList->UncheckedAt(kEta)));
  hetaRat->Draw();
  TH1 *hyRat = static_cast<TH1*>(fRecList->UncheckedAt(kY)->Clone());
  hyRat->SetNameTitle("hyRat","hyRat");
  hyRat->Divide(((TH1F*)fGenList->UncheckedAt(kY)));
  hyRat->SetLineColor(9);
  hyRat->Draw("same");
  gROOT->SetSelectedPad(fcRat->cd(3));
  TH1 *hpRat = static_cast<TH1*>(fRecList->UncheckedAt(kP)->Clone());
  hpRat->SetNameTitle("hpRat","hpRat");
  hpRat->Divide(((TH1F*)fGenList->UncheckedAt(kP)));
  hpRat->Draw();
  gROOT->SetSelectedPad(fcRat->cd(4));
  TH1 *hphiRat = static_cast<TH1*>(fRecList->UncheckedAt(kPhi)->Clone());
  hphiRat->SetNameTitle("hphiRat","hphiRat");
  hphiRat->Divide(((TH1F*)fGenList->UncheckedAt(kPhi)));
  hphiRat->Draw();
  
  // draw histograms for track resolution
  fcResVsP = new TCanvas("cResVsP","cResVsP",10,10,1200,600);
  fcResVsP->Divide(4,2);
  gROOT->SetSelectedPad(fcResVsP->cd(1));
  gPad->SetLogz();
  ((TH2F*)fResList->UncheckedAt(kResPAtVtxVsPIn02degMC))->Draw("colz");
  gROOT->SetSelectedPad(fcResVsP->cd(2));
  gPad->SetLogz();
  ((TH2F*)fResList->UncheckedAt(kResPAtVtxVsPIn23deg))->Draw("colz");
  gROOT->SetSelectedPad(fcResVsP->cd(3));
  gPad->SetLogz();
  ((TH2F*)fResList->UncheckedAt(kResPAtVtxVsPIn310deg))->Draw("colz");
  gROOT->SetSelectedPad(fcResVsP->cd(4));
  gPad->SetLogz();
  ((TH2F*)fResList->UncheckedAt(kResPtAtVtxVsPt))->Draw("colz");
  gROOT->SetSelectedPad(fcResVsP->cd(5));
  gPad->SetLogz();
  ((TH2F*)fResList->UncheckedAt(kResSlopeXAtVtxVsP))->Draw("colz");
  gROOT->SetSelectedPad(fcResVsP->cd(6));
  gPad->SetLogz();
  ((TH2F*)fResList->UncheckedAt(kResSlopeYAtVtxVsP))->Draw("colz");
  gROOT->SetSelectedPad(fcResVsP->cd(7));
  gPad->SetLogz();
  ((TH2F*)fResList->UncheckedAt(kResEtaAtVtxVsP))->Draw("colz");
  gROOT->SetSelectedPad(fcResVsP->cd(8));
  gPad->SetLogz();
  ((TH2F*)fResList->UncheckedAt(kResPhiAtVtxVsP))->Draw("colz");
  
}

