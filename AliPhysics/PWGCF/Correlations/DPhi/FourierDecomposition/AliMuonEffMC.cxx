// MUON track QA referring AliMuonEffMC.cxx
// Author : Saehanseul Oh

#include "AliMuonEffMC.h"

#include <TList.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THn.h>
#include <TChain.h>
#include <TFile.h>
#include <TParticle.h>

#include "AliStack.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTask.h"
#include "AliESDVertex.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliESDMuonTrack.h"
#include "AliAODTrack.h"
#include "AliESDVertex.h"
#include "AliAODVertex.h"
#include "AliCentrality.h"
#include "AliVParticle.h"
#include "AliMCParticle.h"
#include "AliMCEvent.h"
#include "AliAnalysisHelperJetTasks.h"
#include "AliAODTrack.h"

ClassImp(AliMuonEffMC)

//________________________________________________________________________
AliMuonEffMC::AliMuonEffMC() :
  AliAnalysisTaskSE(), fESD(0), fAOD(0), fMC(0), fStack(0), fCentrality(99), fZVertex(99), 
  fOutputList(0x0), fHEventStat(0), fHEvt(0x0), fIsMc(kTRUE), fPlotMode(0), fMuonCutMask(0), fMuonTrackCuts(0x0),
  fCentralityEstimator("V0M"), fNEtaBins(100), fNpTBins(50), fNCentBins(1), fNZvtxBins(1), fNPhiBins(12), 
  fHFM(0x0), fHPP(0x0), fHMuFM(0x0), fHMuFMrec(0x0), fHMuPP(0x0), fHMuPPrec(0x0), fHMuRec(0x0)
{
  // Constructor
  //DefineInput(0, TChain::Class());
  //DefineOutput(1, TList::Class());
}

//________________________________________________________________________
AliMuonEffMC::AliMuonEffMC(const char *name) :
  AliAnalysisTaskSE(name), fESD(0), fAOD(0), fMC(0), fStack(0), fCentrality(99), fZVertex(99), 
  fOutputList(0x0), fHEventStat(0), fHEvt(0x0), fIsMc(kTRUE), fPlotMode(0), fMuonCutMask(0), fMuonTrackCuts(0x0),
  fCentralityEstimator("V0M"), fNEtaBins(100), fNpTBins(50), fNCentBins(1), fNZvtxBins(1), fNPhiBins(12), 
  fHFM(0x0), fHPP(0x0), fHMuFM(0x0), fHMuFMrec(0x0), fHMuPP(0x0), fHMuPPrec(0x0), fHMuRec(0x0)
{
  // Constructor
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
AliMuonEffMC::~AliMuonEffMC()
{
  //Destructor
  if(fOutputList) delete fOutputList;
}

//________________________________________________________________________
void AliMuonEffMC::UserCreateOutputObjects()
{
  // Create histograms
  // Called once (per slave on PROOF!)
  fOutputList = new TList();
  fOutputList->SetOwner(1);

  fHEventStat = new TH1D("fHEventStat","Event statistics for analysis",18,0,18);
  fHEventStat->GetXaxis()->SetBinLabel(1,"Event");
  fHEventStat->GetXaxis()->SetBinLabel(2,"SelectedEvent");
  fHEventStat->GetXaxis()->SetBinLabel(3,"File");
  fHEventStat->GetXaxis()->SetBinLabel(4,"fSPHighpt");  //!Global Trigger Single plus High p_T
  fHEventStat->GetXaxis()->SetBinLabel(5,"fSPAllpt");   //!Global Trigger Single plus All p_T
  fHEventStat->GetXaxis()->SetBinLabel(6,"fSMLowpt");   //!Global Trigger Single minus Low p_T
  fHEventStat->GetXaxis()->SetBinLabel(7,"fSMHighpt");  //!Global Trigger Single minus High p_T
  fHEventStat->GetXaxis()->SetBinLabel(8,"fSMAllpt");   //!Global Trigger Single minus All p_T
  fHEventStat->GetXaxis()->SetBinLabel(9,"fSULowpt");   //!Global Trigger Single undefined Low p_T
  fHEventStat->GetXaxis()->SetBinLabel(10,"fSUHighpt"); //!Global Trigger Single undefined High p_T
  fHEventStat->GetXaxis()->SetBinLabel(11,"fSUAllpt");  //!Global Trigger Single undefined All p_T
  fHEventStat->GetXaxis()->SetBinLabel(12,"fUSLowpt");  //!Global Trigger UnlikeSign pair Low p_T
  fHEventStat->GetXaxis()->SetBinLabel(13,"fUSHighpt"); //!Global Trigger UnlikeSign pair High p_T
  fHEventStat->GetXaxis()->SetBinLabel(14,"fUSAllpt");  //!Global Trigger UnlikeSign pair All p_T
  fHEventStat->GetXaxis()->SetBinLabel(15,"fLSLowpt");  //!Global Trigger LikeSign pair pair Low p_T
  fHEventStat->GetXaxis()->SetBinLabel(16,"fLSHighpt"); //!Global Trigger LikeSign pair pair High p_T
  fHEventStat->GetXaxis()->SetBinLabel(17,"fLSAllpt");  //!Global Trigger LikeSign pair pair All p_T
  fHEventStat->GetXaxis()->SetBinLabel(18,"fSPLowpt");   //!Global Trigger Single plus Low p_T
  fOutputList->Add(fHEventStat);

  fHEvt = new TH2F("fHEvt", "Event-level variables; Zvtx; Cent", 30, -15, 15, 103, -2, 101);
  fOutputList->Add(fHEvt);

  fMuonTrackCuts = new AliMuonTrackCuts("StdMuonCuts","StdMuonCuts");
  fMuonTrackCuts->SetCustomParamFromRun(197388,"muon_pass2"); // for LHC13b,c,d,e,f 
  fMuonTrackCuts->SetFilterMask(fMuonCutMask);
  AliInfo(Form(" using muon track cuts with bit %u\n", fMuonCutMask));
  
  // Define THn's
  Int_t iTrackBin[6];
  Double_t* trackBins[6];
  const char* trackAxisTitle[6];

  // eta
  Double_t etaBins[fNEtaBins+1];
  for(Int_t i=0; i<=fNEtaBins; i++) { etaBins[i] = (Double_t)(-8.0 + 16.0/fNEtaBins*i); }
  iTrackBin[0] = fNEtaBins;
  trackBins[0] = etaBins;
  trackAxisTitle[0] = "#eta";

   // p_T
  Double_t pTBins[fNpTBins+1];
  for(Int_t i=0; i<=fNpTBins; i++) { pTBins[i] = (Double_t)(5.0/fNpTBins * i); }
  iTrackBin[1] = fNpTBins;
  trackBins[1] = pTBins;
  trackAxisTitle[1] = "p_{T} (GeV/c)";

  // centrality/multiplicity
  Double_t CentBins[fNCentBins+1];
  for (Int_t i=0; i<=fNCentBins; i++) { CentBins[i] = (Double_t)(100.0/fNCentBins * i); }
  iTrackBin[2] = fNCentBins;
  trackBins[2] = CentBins;
  trackAxisTitle[2] = "Mult";

  // phi
  Double_t phiBins[fNPhiBins+1];
  for(Int_t i=0; i<=fNPhiBins; i++) { phiBins[i] = (Double_t)(TMath::TwoPi()/fNPhiBins * i); }
  iTrackBin[3] = fNPhiBins;
  trackBins[3] = phiBins;
  trackAxisTitle[3] = "#phi";

   // charge
  Double_t chargeBins[4] = {-10.0, -0.5, 0.5, 10.0};
  iTrackBin[4] = 3;
  trackBins[4] = chargeBins;
  trackAxisTitle[4] = "charge";

  // species
  Double_t Species[21];
  for(Int_t iSpe=0; iSpe<=20; iSpe++) { Species[iSpe] = (Double_t)iSpe; }
  iTrackBin[5] = 20;
  trackBins[5] = Species;
  trackAxisTitle[5] = "Species";

  if(fIsMc)
  {
    // THn for tracking efficiency
    if(fPlotMode==0)
    {
      fHFM = new THnF("fHFM", "", 6, iTrackBin, 0, 0);
      for (Int_t i=0; i<6; i++)
      {
	fHFM->SetBinEdges(i, trackBins[i]);
	fHFM->GetAxis(i)->SetTitle(trackAxisTitle[i]);
      }
      fOutputList->Add(fHFM); 
    }
    else if(fPlotMode==1)
    {
      fHPP = new THnF("fHPP", "", 6, iTrackBin, 0, 0);
      for (Int_t i=0; i<6; i++)
      {
	fHPP->SetBinEdges(i, trackBins[i]);
	fHPP->GetAxis(i)->SetTitle(trackAxisTitle[i]);
      }
      fOutputList->Add(fHPP); 
    }
 
    else if(fPlotMode==2)
    {
      fHMuFM = new THnF("fHMuFM", "", 6, iTrackBin, 0, 0);
      for (Int_t i=0; i<6; i++)
      {
	fHMuFM->SetBinEdges(i, trackBins[i]);
	fHMuFM->GetAxis(i)->SetTitle(trackAxisTitle[i]);
      }
      fOutputList->Add(fHMuFM);

      fHMuFMrec = (THnF*)fHMuFM->Clone("fHMuFMrec");
      fOutputList->Add(fHMuFMrec);
    }
    else if(fPlotMode==3)
    {
      fHMuPP = new THnF("fHMuPP", "", 6, iTrackBin, 0, 0);
      for (Int_t i=0; i<6; i++)
      {
	fHMuPP->SetBinEdges(i, trackBins[i]);
	fHMuPP->GetAxis(i)->SetTitle(trackAxisTitle[i]);
      }
      fOutputList->Add(fHMuPP);
  
      fHMuPPrec = (THnF*)fHMuPP->Clone("fHMuPPrec");
      fOutputList->Add(fHMuPPrec);
    }
  }
  else
  {
    if(fPlotMode==4)
    {
      fHMuRec = new THnF("fHMuRec", "", 6, iTrackBin, 0, 0);
      for (Int_t i=0; i<6; i++)
      {
	fHMuRec->SetBinEdges(i, trackBins[i]);
	fHMuRec->GetAxis(i)->SetTitle(trackAxisTitle[i]);
      }
      fOutputList->Add(fHMuRec); 
    }
  }
  PostData(1, fOutputList);
}

//________________________________________________________________________
void AliMuonEffMC::UserExec(Option_t *)
{
  // Main loop, Called for each event
  Int_t ntrks = 0; // number of tracks in an event

  if(((TString)InputEvent()->IsA()->GetName())=="AliAODEvent")
  {
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    if (!fAOD) { AliError("AOD event not found. Nothing done!"); return; }
    ntrks = fAOD->GetNumberOfTracks();
  }
  else
  {
    fESD = dynamic_cast<AliESDEvent*>(InputEvent());
    if (!fESD) { AliError("ESD event not found. Nothing done!"); return; }
    ntrks = fESD->GetNumberOfMuonTracks();
  }
       
  fHEventStat->Fill(0.5);
  if(fIsMc)
  {
    fMC = MCEvent();
    if (!fMC) { AliError("MC event not avaliable."); return; }
    fStack = fMC->Stack();
  }

  // Centrality, vertex, other event variables...
  if(fAOD)
  {
    const AliAODVertex* vertex = fAOD->GetPrimaryVertex();
    fZVertex = vertex->GetZ();
    if(fAOD->GetCentrality())  fCentrality = fAOD->GetCentrality()->GetCentralityPercentile(fCentralityEstimator);
  }
  else if(fESD)
  {
    const AliESDVertex* vertex = fESD->GetPrimaryVertex();
    fZVertex = vertex->GetZ();
    if(fESD->GetCentrality()) fCentrality = fESD->GetCentrality()->GetCentralityPercentile(fCentralityEstimator);
  }  

  if ((fESD && !VertexOk(fESD)) || (fAOD && !VertexOk(fAOD))) { 
    //AliInfo(Form("Event REJECTED. z = %.1f", fZVertex));
    return; 
  }
  if (fCentrality > 100. || fCentrality < -1.5) { 
    //AliInfo(Form("Event REJECTED. fCentrality = %.1f", fCentrality));
    return; 
  }
 
  if(fCentrality < 0) fCentrality = 0.5; //ad hoc centrality for pp
  // Fill Event histogram
  fHEvt->Fill(fZVertex, fCentrality);
  fHEventStat->Fill(1.5);
 
  ULong64_t trigword = 0;
  if(fAOD) trigword=fAOD->GetTriggerMask();
  else if(fESD) trigword=fESD->GetTriggerMask();
 
  if (trigword & 0x01) fHEventStat->Fill(17.5);
  if (trigword & 0x02) fHEventStat->Fill(3.5);
  if (trigword & 0x04) fHEventStat->Fill(4.5);
  if (trigword & 0x08) fHEventStat->Fill(5.5);      
  if (trigword & 0x010) fHEventStat->Fill(6.5);
  if (trigword & 0x020) fHEventStat->Fill(7.5);
  if (trigword & 0x040) fHEventStat->Fill(8.5);
  if (trigword & 0x080) fHEventStat->Fill(9.5);
  if (trigword & 0x100) fHEventStat->Fill(10.5);
  if (trigword & 0x200) fHEventStat->Fill(11.5);
  if (trigword & 0x400) fHEventStat->Fill(12.5);
  if (trigword & 0x800) fHEventStat->Fill(13.5);
  if (trigword & 0x1000) fHEventStat->Fill(14.5);
  if (trigword & 0x2000) fHEventStat->Fill(15.5);
  if (trigword & 0x4000) fHEventStat->Fill(16.5);
    
  if(fIsMc && (fPlotMode==0 || fPlotMode==1) && fESD)
  {
    Double_t MEta = 0.0;
    Double_t MPt = 0.0;
    Double_t MPhi = 0.0;
    Double_t MCharge = 0.0;
    Double_t MSpecies = 0.0;
        
    for (Int_t ipart=0; ipart<fMC->GetNumberOfTracks(); ipart++) // generated level loop
    {
      AliMCParticle *McParticle  = (AliMCParticle*)fMC->GetTrack(ipart);
      if(fPlotMode==0) 
      {
	if(ipart < 211 || McParticle->GetMother()!=-1) continue;
	MEta = McParticle->Eta();
	MPt = McParticle->Pt();
	MPhi = McParticle->Phi();
	MCharge = McParticle->Charge();
	MSpecies = GetSpecies(McParticle->PdgCode());
      }
      else if(fPlotMode==1)
      {	  
	if(!fMC->IsPhysicalPrimary(ipart)) continue;
	MEta = McParticle->Eta();
	MPt = McParticle->Pt();
	MPhi = McParticle->Phi();
	MCharge = McParticle->Charge();
	MSpecies = GetSpecies(McParticle->PdgCode());
      }
      Double_t fillArrayM[6] = { MEta, MPt, fCentrality, MPhi, MCharge, MSpecies };
      if(fPlotMode==0) fHFM->Fill(fillArrayM);
      else if(fPlotMode==1) fHPP->Fill(fillArrayM);
    }// end of generated level loop
  } 
    
  if((fPlotMode==2 || fPlotMode==3 || fPlotMode==4))
  {
    for (Int_t iTrack = 0; iTrack<ntrks; iTrack++) // reconstructed level loop
    {
      Int_t label = 0;
      Double_t MuEta = 0.0;
      Double_t MuPt = 0.0;
      Double_t MuPhi = 0.0;
      Double_t MuCharge = 0.0;
      
      Double_t MuMEta = 0.0;
      Double_t MuMPt = 0.0;
      Double_t MuMPhi = 0.0;
      Double_t MuMCharge = 0.0;
      Double_t MuMSpecies = 0.0; 
      
      if(fESD)
      {
	AliESDMuonTrack* muonTrack = fESD->GetMuonTrack(iTrack);
	if(muonTrack)
	{
	  if(!IsGoodMUONtrack(*muonTrack)) continue;
	  MuEta = muonTrack->Eta();
	  MuPt = muonTrack->Pt();
	  MuPhi = muonTrack->Phi();
	  MuCharge = muonTrack->Charge();
	  label =  TMath::Abs(muonTrack->GetLabel());
	  if (label>=fMC->GetNumberOfTracks()) {
	    AliError(Form("Label %d larger than number of particles on stack %d\n",label,fMC->GetNumberOfTracks()));
	    continue;
	  }
	}
      }
      else if(fAOD)
      {
	AliAODTrack* muonTrack = (AliAODTrack*)fAOD->GetTrack(iTrack);
	if(muonTrack)
	{
	  if(!IsGoodMUONtrack(*muonTrack)) continue;
	  MuEta = muonTrack->Eta();
	  MuPt = muonTrack->Pt();
	  MuPhi = muonTrack->Phi();
	  MuCharge = muonTrack->Charge();
	  label =  TMath::Abs(muonTrack->GetLabel());
	  if (label>=fMC->GetNumberOfTracks()) {
	    AliError(Form("Label %d larger than number of particles on stack %d\n",label,fMC->GetNumberOfTracks()));
	    continue;
	  }
        }
      }
          
      if(fIsMc && fESD)
      {
	if(fPlotMode==2) 
	{
	  AliMCParticle* MuMparticle = (AliMCParticle*)fMC->GetTrack(GetFirstMother(label));
	  MuMEta = MuMparticle->Eta();
	  MuMPt = MuMparticle->Pt();
	  MuMPhi = MuMparticle->Phi();
	  MuMCharge = MuMparticle->Charge();
	  MuMSpecies = GetSpecies(MuMparticle->PdgCode());
	}      
	else if(fPlotMode==3)
	{
	  AliMCParticle* MuMparticle = (AliMCParticle*)fMC->GetTrack(GetFirstPPMother(label));
	  MuMEta = MuMparticle->Eta();
	  MuMPt = MuMparticle->Pt();
	  MuMPhi = MuMparticle->Phi();
	  MuMCharge = MuMparticle->Charge();
	  MuMSpecies = GetSpecies(MuMparticle->PdgCode());
	}
      }// end of MC process
      Double_t fillArrayMuRec[6] = { MuEta, MuPt, fCentrality, MuPhi, MuCharge, 0.5 };
      Double_t fillArrayMuRecM[6] = { MuEta, MuPt, fCentrality, MuPhi, MuCharge, MuMSpecies };
      Double_t fillArrayMuM[6] = { MuMEta, MuMPt, fCentrality, MuMPhi, MuMCharge, MuMSpecies };
 
      if(fPlotMode==2) { fHMuFM->Fill(fillArrayMuM); fHMuFMrec->Fill(fillArrayMuRecM); }
      else if(fPlotMode==3) { fHMuPP->Fill(fillArrayMuM); fHMuPPrec->Fill(fillArrayMuRecM); }
      else if(fPlotMode==4) { fHMuRec->Fill(fillArrayMuRec); }  
    }// end of reconstructed loop
  }
  PostData(1, fOutputList);
  return;
}

//________________________________________________________________________
void AliMuonEffMC::Terminate(Option_t *)
{
  // Draw result to the screen
  // Called once at the end of the query
}
//________________________________________________________________________
Bool_t AliMuonEffMC::VertexOk(TObject* obj) const
{
  // Modified from AliAnalyseLeadingTrackUE::VertexSelection()
 
  Int_t nContributors  = 0;
  Double_t zVertex     = 999;
  TString name("");
 
  if (obj->InheritsFrom("AliESDEvent")) {
    AliESDEvent* esdevt = (AliESDEvent*) obj;
    const AliESDVertex* vtx = esdevt->GetPrimaryVertex();
    if (!vtx)
      return 0;
    nContributors = vtx->GetNContributors();
    zVertex       = vtx->GetZ();
    name          = vtx->GetName();
  }
  else if (obj->InheritsFrom("AliAODEvent")) {
    AliAODEvent* aodevt = (AliAODEvent*) obj;
    if (aodevt->GetNumberOfVertices() < 1)
      return 0;
    const AliAODVertex* vtx = aodevt->GetPrimaryVertex();
    nContributors = vtx->GetNContributors();
    zVertex       = vtx->GetZ();
    name          = vtx->GetName();
  }
 
  // Reject if TPC-only vertex
  if (name.CompareTo("TPCVertex")==0)
    return kFALSE;
 
  // Check # contributors and range...
  if( nContributors < 1 || TMath::Abs(zVertex) > 10 ) {
    return kFALSE;
  }
 
  return kTRUE;
}

//________________________________________________________________________
Bool_t AliMuonEffMC::IsGoodMUONtrack(AliESDMuonTrack &track)
{
    return fMuonTrackCuts->IsSelected(&track);
}

//________________________________________________________________________
Bool_t AliMuonEffMC::IsGoodMUONtrack(AliAODTrack &track)
{
    return fMuonTrackCuts->IsSelected(&track);
}

//________________________________________________________________________
Double_t AliMuonEffMC::GetMuonTrackType(AliMCParticle &track)
{
  if(track.GetMother() < fStack->GetNprimary() && track.PdgCode() == 13) return 0.5; 
  else if(track.GetMother() >= fStack->GetNprimary() && track.PdgCode() == 13) return 1.5;
  else return 2.5;
}

//________________________________________________________________________
Double_t AliMuonEffMC::GetSpecies(Int_t PdgCode)
{
  Int_t code = TMath::Abs(PdgCode);
  if(code==13) return 0.5;
  else if(code==211) return 1.5;
  else if(code==321 || code==311) return 2.5;
  else if(code==411 || code==421) return 3.5;
  else if(code==511 || code==521) return 4.5;
  else if(code==213) return 5.5;
  else if(code==313 || code==323) return 6.5;
  else if(code==413 || code==423 || code==431) return 7.5;
  else if(code==513 || code==523 || code==531 || code==533 || code==541 || code==543) return 8.5;
  else if(code==111) return 9.5;
  else if(code==113) return 10.5;
  else if(code==221 || code==331) return 11.5;
  else if(code==223) return 12.5;
  else if(code==2212) return 13.5;
  else if(code==2112) return 14.5;
  else if(code==1114 || code==2114 || code==2214 || code==2224) return 15.5;
  else if(code==3112 || code==3114 || code==3212) return 16.5;
  else if(code>100 && code<1000) return 17.5;
  else if(code>1000 && code<10000) return 18.5;
  else return 19.5;
}

//________________________________________________________________________
Double_t AliMuonEffMC::deltaphi(Double_t phi)
{
  if(phi < -1.0*TMath::Pi()) { return (phi + TMath::TwoPi()); }
  else if(phi > TMath::Pi()) { return (phi - TMath::TwoPi()); }
  else { return phi; }
}

//________________________________________________________________________
Int_t AliMuonEffMC::GetFirstMother(Int_t muonlabel)
{
  if(fAOD) return -1;
  else if(fESD)
  {
    AliMCParticle *McParticle  = (AliMCParticle*)fMC->GetTrack(muonlabel);
    Int_t motherlabel = McParticle->GetMother();
    Int_t nextmother = 0;
    if(motherlabel==-1) return -1;
    else
    {
      while(motherlabel>-1)
      {
	AliMCParticle *MotherParticle  = (AliMCParticle*)fMC->GetTrack(motherlabel);
	nextmother = motherlabel;
	motherlabel = MotherParticle->GetMother();
      }
      return nextmother;
    }
  }
  else return -1;
}

//________________________________________________________________________
Int_t AliMuonEffMC::GetFirstPPMother(Int_t muonlabel)
{
  if(fAOD) return 1;
  else if(fESD)
  {
    AliMCParticle *McParticle  = (AliMCParticle*)fMC->GetTrack(muonlabel);
    if(fMC->IsPhysicalPrimary(muonlabel)) return muonlabel;
    else
    {
      Int_t motherlabel = McParticle->GetMother();
      while(motherlabel > -1)
      {
	AliMCParticle *MotherParticle  = (AliMCParticle*)fMC->GetTrack(motherlabel);
	if(fMC->IsPhysicalPrimary(motherlabel)) break;
	else motherlabel = MotherParticle->GetMother();
      }
      return motherlabel;
    }
  }
  else return -1;
}
