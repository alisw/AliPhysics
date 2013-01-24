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

#include "AliAnalysisManager.h"
#include "AliAnalysisTask.h"
#include "AliESDVertex.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliESDMuonTrack.h"
#include "AliESDVertex.h"
#include "AliAODVertex.h"
#include "AliCentrality.h"
#include "AliVParticle.h"
#include "AliMCParticle.h"
#include "AliAnalysisHelperJetTasks.h"
#include "AliMCEvent.h"

using std::cout;
using std::endl;

ClassImp(AliMuonEffMC)

//________________________________________________________________________
AliMuonEffMC::AliMuonEffMC() :
  AliAnalysisTaskSE(), fESD(0), fAOD(0), fMC(0), fCentrality(99), fZVertex(99), fOutputList(0x0),      
  fHEventStat(0), fHEvt(0x0), fIsMc(kTRUE), fMDProcess(kTRUE), fCentralityEstimator("V0M"),
  fNEtaBins(15), fNpTBins(100), fNCentBins(5), fNZvtxBins(10), fNPhiBins(24),
  fHMuonParGen(0x0), fHMuonDetGen(0x0), fHMuonDetRec(0x0), fHEtcDetRec(0x0)
{
  // Constructor
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());

  for(Int_t i=0; i<4; i++)
  {
    fHMuMotherGenPt[i] = 0x0;
    fHMuMotherRecPt[i] = 0x0;
    fHMuMotherGenPhi[i] = 0x0;
    fHMuMotherRecPhi[i] = 0x0;
    fHMuMotherGenEta[i] = 0x0;
    fHMuMotherRecEta[i] = 0x0;
    fHMuDCA[i] = 0x0;
  }
}

//________________________________________________________________________
AliMuonEffMC::AliMuonEffMC(const char *name) :
  AliAnalysisTaskSE(name), fESD(0), fAOD(0), fMC(0), fCentrality(99), fZVertex(99), fOutputList(0x0),      
  fHEventStat(0), fHEvt(0x0),  fIsMc(kTRUE), fMDProcess(kTRUE), fCentralityEstimator("V0M"),
  fNEtaBins(15), fNpTBins(100), fNCentBins(5), fNZvtxBins(10), fNPhiBins(24),
  fHMuonParGen(0x0), fHMuonDetGen(0x0), fHMuonDetRec(0x0), fHEtcDetRec(0x0)
{
  // Constructor
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());

  for(Int_t i=0; i<4; i++)
  {
    fHMuMotherGenPt[i] = 0x0;
    fHMuMotherRecPt[i] = 0x0;
    fHMuMotherGenPhi[i] = 0x0;
    fHMuMotherRecPhi[i] = 0x0;
    fHMuMotherGenEta[i] = 0x0;
    fHMuMotherRecEta[i] = 0x0;
    fHMuDCA[i] = 0x0;
  }
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


  Int_t iTrackBin[5];
  Double_t* trackBins[5];
  const char* trackAxisTitle[5];

  // eta
  Double_t etaBins[fNEtaBins+1];
  for (Int_t i=0; i<=fNEtaBins; i++) { etaBins[i] = (Double_t)(-4.0 + 1.5/fNEtaBins*i); }
  iTrackBin[0] = fNEtaBins;
  trackBins[0] = etaBins;
  trackAxisTitle[0] = "#eta";

  // p_T
  Double_t pTBins[fNpTBins+1];
  for (Int_t i=0; i<=fNpTBins; i++) { pTBins[i] = (Double_t)(10.0/fNpTBins*i); }
  iTrackBin[1] = fNpTBins;
  trackBins[1] = pTBins;
  trackAxisTitle[1] = "p_{T} (GeV/c)";

  // centrality
  Double_t CentBins[fNCentBins+1];
  for (Int_t i=0; i<=fNCentBins; i++) { CentBins[i] = (Double_t)(100.0/fNCentBins*i); }
  iTrackBin[2] = fNCentBins;
  trackBins[2] = CentBins;
  trackAxisTitle[2] = "Cent";

  // Z-vertex
  Double_t ZvtxBins[fNZvtxBins+1];
  for (Int_t i=0; i<=fNZvtxBins; i++) { ZvtxBins[i] = (Double_t)(-10.0 + 20.0/fNZvtxBins*i); }
  iTrackBin[3] = fNZvtxBins;
  trackBins[3] = ZvtxBins;
  trackAxisTitle[3] = "Zvtx";

  // phi
  Double_t phiBins[fNPhiBins+1];
  for (Int_t i=0; i<=fNPhiBins; i++) { phiBins[i] = (Double_t)(TMath::TwoPi()/fNPhiBins * i); }
  iTrackBin[4] = fNPhiBins;
  trackBins[4] = phiBins;
  trackAxisTitle[4] = "#phi";

  // THn for tracking efficiency
  fHMuonParGen = new THnF("fHMuonParGen", "", 5, iTrackBin, 0, 0);
  for (Int_t i=0; i<5; i++)
  {
    fHMuonParGen->SetBinEdges(i, trackBins[i]);
    fHMuonParGen->GetAxis(i)->SetTitle(trackAxisTitle[i]);
  }
  fHMuonParGen->Sumw2();
  fOutputList->Add(fHMuonParGen);

  fHMuonDetGen = (THnF*) fHMuonParGen->Clone("fHMuonDetGen");
  fHMuonDetGen->Sumw2();
  fOutputList->Add(fHMuonDetGen);

  fHMuonDetRec = (THnF*) fHMuonParGen->Clone("fHMuonDetRec");
  fHMuonDetRec->Sumw2();
  fOutputList->Add(fHMuonDetRec);

  fHEtcDetRec = (THnF*) fHMuonParGen->Clone("fHEtcDetRec");
  fHEtcDetRec->Sumw2();
  fOutputList->Add(fHEtcDetRec);

  if(fMDProcess)
  {
    const char* MotherSpecies[4] = {"Pion","Kaon","D", "Etc"};
    
    for(Int_t i=0; i<4; i++)
    {
      fHMuMotherGenPt[i] = new TH2F(Form("fHMuMotherGenPt_%s",MotherSpecies[i]),";p_{T,muon}^{gen} (GeV/c);p_{T,mother}^{Truth} (GeV/c);",500, 0, 50, 500, 0, 50);
      fOutputList->Add(fHMuMotherGenPt[i]);
      
      fHMuMotherRecPt[i] = new TH2F(Form("fHMuMotherRecPt_%s",MotherSpecies[i]),";p_{T,muon}^{rec} (GeV/c);p_{T,mother}^{Truth} (GeV/c);",500, 0, 50, 500, 0, 50);
      fOutputList->Add(fHMuMotherRecPt[i]);
      
      fHMuMotherGenPhi[i] = new TH2F(Form("fHMuMotherGenPhi_%s",MotherSpecies[i]),";#phi_{gen};mother #phi;",100, 0, TMath::TwoPi(), 100, 0, TMath::TwoPi());
      fOutputList->Add(fHMuMotherGenPhi[i]);
      
      fHMuMotherRecPhi[i] = new TH2F(Form("fHMuMotherRecPhi_%s",MotherSpecies[i]),";#phi_{rec};mother #phi;",100, 0, TMath::TwoPi(), 100, 0, TMath::TwoPi());
      fOutputList->Add(fHMuMotherRecPhi[i]);
      
      fHMuMotherGenEta[i] = new TH2F(Form("fHMuMotherGenEta_%s",MotherSpecies[i]),";#eta_{gen};mother #eta;",100, -5., -1., 100, -5., -1.);
      fOutputList->Add(fHMuMotherGenEta[i]);
      
      fHMuMotherRecEta[i] = new TH2F(Form("fHMuMotherRecEta_%s",MotherSpecies[i]),";#eta_{rec};mother #eta;",100, -5., -1., 100, -5., -1.);
      fOutputList->Add(fHMuMotherRecEta[i]);
      
      fHMuDCA[i] =  new TH1F(Form("fHMuDCA_%s",MotherSpecies[i]), ";DCA", 100, 0, 50);
      fOutputList->Add(fHMuDCA[i]);
    }
  }
  PostData(1, fOutputList);
}

//________________________________________________________________________
void AliMuonEffMC::UserExec(Option_t *) 
{
  // Main loop, Called for each event
  fESD = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!fESD) 
  {
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  }

  if (!fESD && !fAOD) 
  {
    AliError("Neither fESD nor fAOD available");
    return;
  }
	
  fHEventStat->Fill(0.5);
  if(fIsMc) 
  {
    fMC = MCEvent();
    if (!fMC) {
      printf("ERROR: fMC not available\n");
      return;
    }
  }

  const AliESDVertex* vertex = fESD->GetPrimaryVertex();
  fZVertex = vertex->GetZ();
  if(fESD->GetCentrality()) fCentrality = fESD->GetCentrality()->GetCentralityPercentile("fCentralityEstimator");

  // Fill Event histogram
  fHEvt->Fill(fZVertex, fCentrality);

  Int_t iVerb = 0;

   // Centrality, vertex, other event variables...
  if (!VertexOk(fESD)) { 
    if (iVerb>1)
      AliInfo(Form("Event REJECTED. z = %.1f", fZVertex));  
    return; 
  }

  if (fCentrality > 100. || fCentrality < -1.5) { 
    if (iVerb>1)
      AliInfo(Form("Event REJECTED. fCentrality = %.1f", fCentrality)); 
    return; 
  }
 
  if(fCentrality < 0) fCentrality = 1.0;
  fHEventStat->Fill(1.5);
 
  ULong64_t trigword=fESD->GetTriggerMask();
 
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

  if(fIsMc)
  {
    for (Int_t iTrack = 0; iTrack<fESD->GetNumberOfMuonTracks(); iTrack++) 
    { 
      // loop on muon tracks
      AliESDMuonTrack* muonTrack = fESD->GetMuonTrack(iTrack);
      if(muonTrack)
      {
	if(!IsGoodMUONtrack(*muonTrack)) continue;
	Double_t fillArrayDetRec[5] = { muonTrack->Eta(), muonTrack->Pt(), fCentrality, fZVertex, muonTrack->Phi() };
	fHMuonDetRec->Fill(fillArrayDetRec); 
	
	Int_t label = TMath::Abs(muonTrack->GetLabel());
	AliMCParticle *McParticle  = (AliMCParticle*)fMC->GetTrack(label);
	if(TMath::Abs(McParticle->PdgCode()) != 13) 
	{
	  fHEtcDetRec->Fill(fillArrayDetRec); 
	  continue;
	}
	Double_t fillArrayDetGen[5] = { McParticle->Eta(), McParticle->Pt(), fCentrality, fZVertex, McParticle->Phi() };
	fHMuonDetGen->Fill(fillArrayDetGen);
	
	if(fMDProcess)
	{
	  if(McParticle->GetMother() > 0)
	  {
	    AliMCParticle *MotherParticle  = (AliMCParticle*)fMC->GetTrack(McParticle->GetMother());
	    Int_t motherlabel = TMath::Abs(MotherParticle->PdgCode());
	    
	    if(motherlabel==411 || motherlabel==413 || motherlabel==421 || motherlabel==423 || motherlabel==431 || motherlabel==433 || motherlabel==10413 || motherlabel==10411 || motherlabel==10423 || motherlabel==10421 || motherlabel==10433 || motherlabel==10431 || motherlabel==20413 || motherlabel==415 || motherlabel==20423 || motherlabel==425 || motherlabel==20433 || motherlabel==435)
	    {  
	      fHMuMotherGenPt[2]->Fill(McParticle->Pt(), MotherParticle->Pt());
	      fHMuMotherRecPt[2]->Fill(muonTrack->Pt(), MotherParticle->Pt());
	      fHMuMotherGenPhi[2]->Fill(McParticle->Phi(), MotherParticle->Phi());
	      fHMuMotherRecPhi[2]->Fill(muonTrack->Phi(), MotherParticle->Phi());
	      fHMuMotherGenEta[2]->Fill(McParticle->Eta(), MotherParticle->Eta());
	      fHMuMotherRecEta[2]->Fill(muonTrack->Eta(), MotherParticle->Eta());
	      fHMuDCA[2]->Fill(muonTrack->GetDCA());
	    }
	    else if(motherlabel==211) 
	    {
	      fHMuMotherGenPt[0]->Fill(McParticle->Pt(), MotherParticle->Pt());
	      fHMuMotherRecPt[0]->Fill(muonTrack->Pt(), MotherParticle->Pt());
	      fHMuMotherGenPhi[0]->Fill(McParticle->Phi(), MotherParticle->Phi());
	      fHMuMotherRecPhi[0]->Fill(muonTrack->Phi(), MotherParticle->Phi());
	      fHMuMotherGenEta[0]->Fill(McParticle->Eta(), MotherParticle->Eta());
	      fHMuMotherRecEta[0]->Fill(muonTrack->Eta(), MotherParticle->Eta());
	      fHMuDCA[0]->Fill(muonTrack->GetDCA());
	    }
	    else if(motherlabel==321)
	    {
	      fHMuMotherGenPt[1]->Fill(McParticle->Pt(), MotherParticle->Pt());
	      fHMuMotherRecPt[1]->Fill(muonTrack->Pt(), MotherParticle->Pt());
	      fHMuMotherGenPhi[1]->Fill(McParticle->Phi(), MotherParticle->Phi());
	      fHMuMotherRecPhi[1]->Fill(muonTrack->Phi(), MotherParticle->Phi());
	      fHMuMotherGenEta[1]->Fill(McParticle->Eta(), MotherParticle->Eta());
	      fHMuMotherRecEta[1]->Fill(muonTrack->Eta(), MotherParticle->Eta());
	      fHMuDCA[1]->Fill(muonTrack->GetDCA());
	    }	
	    else 
	    {
	      fHMuMotherGenPt[3]->Fill(McParticle->Pt(), MotherParticle->Pt());
	      fHMuMotherRecPt[3]->Fill(muonTrack->Pt(), MotherParticle->Pt());
	      fHMuMotherGenPhi[3]->Fill(McParticle->Phi(), MotherParticle->Phi());
	      fHMuMotherRecPhi[3]->Fill(muonTrack->Phi(), MotherParticle->Phi());
	      fHMuMotherGenEta[3]->Fill(McParticle->Eta(), MotherParticle->Eta());
	      fHMuMotherRecEta[3]->Fill(muonTrack->Eta(), MotherParticle->Eta());
	      fHMuDCA[3]->Fill(muonTrack->GetDCA());
	    }
	  }    
	} // (mother hadron) : (daughter muon) QA 
      } 
    }
    
    for (Int_t ipart=0; ipart<fMC->GetNumberOfTracks(); ipart++)
    {
      AliMCParticle *McParticle  = (AliMCParticle*)fMC->GetTrack(ipart);
      if(TMath::Abs(McParticle->PdgCode())==13) 
      {
	Double_t fillArrayParGen[5] = { McParticle->Eta(), McParticle->Pt(), fCentrality, fZVertex, McParticle->Phi() };
	fHMuonParGen->Fill(fillArrayParGen);
      }
    }  
  }
  PostData(1, fOutputList);
  return;
}

//________________________________________________________________________
void AliMuonEffMC::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query

  fOutputList = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutputList) {
    AliError("Output list not available");
    return;
  }
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
  // Applying track cuts for MUON tracks
  if(!track.ContainTrackerData()) return kFALSE;
  if(!track.ContainTriggerData()) return kFALSE;

  Double_t thetaTrackAbsEnd = TMath::ATan(track.GetRAtAbsorberEnd()/505.) * TMath::RadToDeg();
  Double_t eta = track.Eta();

  // Theta cut at absorber end
  if(thetaTrackAbsEnd <= 2. || thetaTrackAbsEnd >= 10.) return kFALSE;
  // Eta cut
  if(eta <= -4. || eta >= -2.5) return kFALSE;
  if(track.GetMatchTrigger() <= 0) return kFALSE;

  return kTRUE;
}
