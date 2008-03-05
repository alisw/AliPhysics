#include "AliMuonAnalysis.h"
//________________________________
///////////////////////////////////////////////////////////
//
// class AliMuonAnalysis
//
// MUON Analysis
//
//
// 
// finck@subatech.in2p3.fr
//
///////////////////////////////////////////////////////////
/*********************************************************/

#include <TString.h>
#include <TParticle.h>

#include <AliStack.h>
#include <AliAOD.h>
#include <AliAODParticle.h>
#include <AliAODParticleCut.h>

#include <AliESDtrack.h>
#include <AliESD.h>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"

ClassImp(AliMuonAnalysis)

AliMuonAnalysis::AliMuonAnalysis():
  fHistoFile(0x0),
  fHPtMuon(0x0),
  fHPtMuonPlus(0x0),
  fHPtMuonMinus(0x0),
  fHPMuon(0x0),
  fHInvMassAll(0x0),
  fHRapMuon(0x0),
  fHRapResonance(0x0),
  fHPtResonance(0x0),
  fHInvMassAllvsPt(0x0),
  fPartCut(0x0)
{
 //ctor
}

/*********************************************************/
AliMuonAnalysis::~AliMuonAnalysis()
{
 //dtor
  delete fPartCut;
  delete fHistoFile;
  delete fHPtMuon;
  delete fHPtMuonPlus;
  delete fHPtMuonMinus;
  delete fHPMuon;
  delete fHInvMassAll;
  delete fHRapMuon;
  delete fHRapResonance;
  delete fHPtResonance;
  delete fHInvMassAllvsPt;
}
/*********************************************************/

Int_t AliMuonAnalysis::Init()
{
  //Initilizes analysis
  Info("Init","Histo initialized for MUON Analysis");

  fHistoFile         = new TFile("MUONmassPlot.root", "RECREATE");
  fHPtMuon           = new TH1F("hPtMuon", "Muon Pt (GeV/c)", 100, 0., 20.);
  fHPMuon            = new TH1F("hPMuon", "Muon P (GeV/c)", 100, 0., 200.);
  fHPtMuonPlus       = new TH1F("hPtMuonPlus", "Muon+ Pt (GeV/c)", 100, 0., 20.);
  fHPtMuonMinus      = new TH1F("hPtMuonMinus", "Muon- Pt (GeV/c)", 100, 0., 20.);
  fHInvMassAll       = new TH1F("hInvMassAll", "Mu+Mu- invariant mass (GeV/c2)", 480, 0., 12.);
  fHRapMuon          = new TH1F("hRapMuon"," Muon Rapidity",50,-4.5,-2);
  fHRapResonance     = new TH1F("hRapResonance"," Resonance Rapidity",50,-4.5,-2);
  fHPtResonance      = new TH1F("hPtResonance", "Resonance Pt (GeV/c)", 100, 0., 20.);
  fHInvMassAllvsPt = new TH2F("hInvMassAll_vs_Pt","hInvMassAll_vs_Pt",480,0.,12.,80,0.,20.);

  return 0; 
}
/*********************************************************/

Int_t AliMuonAnalysis::ProcessEvent(AliAOD* aodrec, AliAOD* aodsim)
{
  //
  // process the event
  // 
  if (aodrec) {
    GetInvMass(aodrec);
    //    Info("ProcessEvent","Inv Mass Rec");
  }  

  if (aodsim) {
    //     Info("ProcessEvent","aodsim not implemented");
  }  
  
  return 0;
  
}

/*********************************************************/

Int_t AliMuonAnalysis::Finish()
{
  //Finish analysis and writes results
  Info("Finish","Histo writing for MUON Analysis");

  fHistoFile->Write();
  fHistoFile->Close();

  return 0;
}
/*********************************************************/

void AliMuonAnalysis::GetInvMass(AliAOD* aod)
{
  // get the invariant mass distribution
  // from the oad events
 
  TLorentzVector lorV1, lorV2, lorVtot;
  Float_t massMin = 9.17;
  Float_t massMax = 9.77;
  Int_t charge1, charge2;

//returns flow parameters: v2 and event plane
  if (aod == 0x0) {
     Error("AliMuonAnalysis::GetInvMass","Pointer to AOD is NULL");
     return;
  }
   
  Int_t nPart = aod->GetNumberOfParticles();
  
  for (Int_t iPart1 = 0; iPart1 < nPart; iPart1++)  {
    AliAODParticle* aodPart1 =  (AliAODParticle*)aod->GetParticle(iPart1);

    if (aodPart1 == 0x0) {
      Error("AliMuonAnalysis::GetInvMass","Cannot get particle %d", iPart1);
      continue;
    }

    lorV1.SetPxPyPzE(aodPart1->Px(),
			 aodPart1->Py(),
			 aodPart1->Pz(),
			 aodPart1->E());

    fHPtMuon->Fill(lorV1.Pt());
    fHPMuon->Fill(lorV1.P());

    charge1 = TMath::Sign(1,aodPart1->GetPdgCode());

    if (charge1 > 0) {
      fHPtMuonPlus->Fill(lorV1.Pt());
    } else {
      fHPtMuonMinus->Fill(lorV1.Pt());
    }
    fHRapMuon->Fill(lorV1.Rapidity());
    for (Int_t iPart2 = iPart1 + 1; iPart2 < nPart; iPart2++)  {

      AliAODParticle* aodPart2 =  (AliAODParticle*)aod->GetParticle(iPart2);

      lorV2.SetPxPyPzE(aodPart2->Px(),
			 aodPart2->Py(),
			 aodPart2->Pz(),
			 aodPart2->E());

      charge2 = TMath::Sign(1,aodPart2->GetPdgCode());

      if ((charge1 * charge2) == -1) {

	lorVtot = lorV1;
	lorVtot += lorV2;
	Float_t invMass = lorVtot.M();

	fHInvMassAll->Fill(invMass);
	fHInvMassAllvsPt->Fill(invMass,lorVtot.Pt());

	if (invMass > massMin && invMass < massMax) {
	  fHRapResonance->Fill(lorVtot.Rapidity());
	  fHPtResonance->Fill(lorVtot.Pt());
	}
      }

    }
  }
}
