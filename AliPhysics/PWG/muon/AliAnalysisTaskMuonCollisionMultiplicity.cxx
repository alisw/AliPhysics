/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

/// Compute the number of Muons tracks as a function of the SPD tracklets multiplicity
/// Compare with Monte Carlo tracks
/// Author Matthieu LENHARDT - SUBATECH, Nantes


//PWG3/muon includes
#include "AliAnalysisTaskMuonCollisionMultiplicity.h"

//STEER includes
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliESDEvent.h"
#include "AliESDVertex.h"
#include "AliAODTrack.h"
#include "AliESDMuonTrack.h"
#include "AliAODTracklets.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisManager.h"
#include "AliAODDimuon.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliMultiplicity.h"

//ROOT includes
#include <TH2D.h>
#include <THnSparse.h>
#include <TChain.h>
#include <TList.h>
#include <TArrayD.h>
#include <Riostream.h>
#include <TParticle.h>
#include <TLorentzVector.h>

//___________________________________________________
AliAnalysisTaskMuonCollisionMultiplicity::AliAnalysisTaskMuonCollisionMultiplicity()
  :
  AliAnalysisTaskSE(),
  fIsInit(kFALSE),
  fAOD(0x0),
  fESD(0x0),
  fEtaCut(0),
  fTrackletMultiplicity(0),
  fTriggerList(0),
  fSingleMuonList(0),
  fDimuonList(0),
  fMonteCarloList(0)
{
  ///Default Constructor
}


//___________________________________________________
AliAnalysisTaskMuonCollisionMultiplicity::AliAnalysisTaskMuonCollisionMultiplicity(const AliAnalysisTaskMuonCollisionMultiplicity& src)
  :
  AliAnalysisTaskSE(),
  fIsInit(kFALSE),
  fAOD(0x0),
  fESD(0x0),
  fEtaCut(0),
  fTrackletMultiplicity(0),
  fTriggerList(0),
  fSingleMuonList(0),
  fDimuonList(0),
  fMonteCarloList(0)
{
  /// copy ctor
  src.Copy(*this);

}

//___________________________________________________
AliAnalysisTaskMuonCollisionMultiplicity& AliAnalysisTaskMuonCollisionMultiplicity::operator=(const AliAnalysisTaskMuonCollisionMultiplicity& src)
{
  /// assignement operator
  if ( this != &src ) 
  {
    src.Copy(*this);
  }
  return *this;
}




//___________________________________________________
AliAnalysisTaskMuonCollisionMultiplicity::AliAnalysisTaskMuonCollisionMultiplicity(const Char_t *name)
  :
  AliAnalysisTaskSE(name),
  fIsInit(kFALSE),
  fAOD(0x0),
  fESD(0x0),
  fEtaCut(0),
  fTrackletMultiplicity(0),
  fTriggerList(0),
  fSingleMuonList(0),
  fDimuonList(0),
  fMonteCarloList(0)
{
  // Define Inputs and outputs
  //DefineInput(0, TChain::Class());
  //DefineInput(1, TChain::Class());
  DefineOutput(0, TList::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
  DefineOutput(3, TList::Class());
  DefineOutput(4, TList::Class());
}


//______________________________________________________________________________
AliAnalysisTaskMuonCollisionMultiplicity::~AliAnalysisTaskMuonCollisionMultiplicity()
{
// Destructor.
  delete fAOD;
  delete fESD;
  delete fTriggerList;
  delete fSingleMuonList;
  delete fDimuonList;
  delete fMonteCarloList;
}



//________________________________________________________________________
void AliAnalysisTaskMuonCollisionMultiplicity::UserCreateOutputObjects()
{
  // Initialise the object, and open the file.
  if (!fIsInit)
    Init();
  OpenFile(0);
}




//________________________________________________________________________
void AliAnalysisTaskMuonCollisionMultiplicity::UserExec(Option_t */*option*/)
{
  // Execute the analysis task
  fAOD = 0x0;
  fESD = 0x0;

  if (!fIsInit)
    Init();

  fAOD = dynamic_cast<AliAODEvent *> (InputEvent());
  if (!fAOD)
    fESD = dynamic_cast<AliESDEvent *> (InputEvent());
  
 
  if (fAOD)
    CheckEventAOD();

  if (fESD)
    CheckEventESD();

  PostData(1, fTriggerList);
  PostData(2, fSingleMuonList);
  PostData(3, fDimuonList);
  PostData(4, fMonteCarloList);
}



//________________________________________________________________________
void AliAnalysisTaskMuonCollisionMultiplicity::NotifyRun()
{
  // Notify run
}


//________________________________________________________________________
void AliAnalysisTaskMuonCollisionMultiplicity::FinishTaskOutput()
{
  // Finish the task
}


//__________________________________________________________________________
Bool_t AliAnalysisTaskMuonCollisionMultiplicity::CheckEventAOD()
{
  // Check if the AOD event pass the cuts
  AliAODVertex *vertex = fAOD->GetPrimaryVertex();
  
  if (!vertex)
     return kFALSE;

  if (vertex->GetNContributors() < 1)
    return kFALSE;
  
  ComputeMultiplicity();
  
  // Variables use to determine the type of trigger :
  // 0 for minimum bias : CINT1B, CINT1-B or MB1
  // 1 for muon events : CMUS1B, CMUS1-B or MULow
  // -1 for everything else
  TString trigger = fAOD->GetFiredTriggerClasses();
  Int_t triggerClass = -1; 

  if (trigger.Contains("CINT1B") || trigger.Contains("CINT1-B") || trigger.Contains("MB1"))
    triggerClass = 0;

  if (trigger.Contains("CMUS1B") || trigger.Contains("CMUS1-B") || trigger.Contains("MULow"))
    triggerClass = 1;

  if (triggerClass >= 0)
    FillHistosAOD(triggerClass);
  
  return kTRUE;
}



//__________________________________________________________________________
Bool_t AliAnalysisTaskMuonCollisionMultiplicity::CheckEventESD()
{
  // Check if the ESD event pass the cuts
  const AliESDVertex *vertex = fESD->GetPrimaryVertex();
  
  if (!vertex)
    return kFALSE;

  if (vertex->GetNContributors() < 1)
    return kFALSE;

  ComputeMultiplicity();

  // Variables use to determine the type of trigger :
  // 0 for minimum bias : CINT1B, CINT1-B or MB1
  // 1 for muon events : CMUS1B, CMUS1-B or MULow
  // -1 for everything else
  TString trigger = fESD->GetFiredTriggerClasses();
  Int_t triggerClass = -1; 

  if (trigger.Contains("CINT1B") || trigger.Contains("CINT1-B") || trigger.Contains("MB1") || trigger.Contains("CMBAC-B") || trigger.Contains("CMBACS2-B"))
    triggerClass = 0;

  if (trigger.Contains("CMUS1B") || trigger.Contains("CMUS1-B") || trigger.Contains("MULow"))
    triggerClass = 1;

  if (triggerClass >= 0)
    FillHistosESD(triggerClass);

  if (fMCEvent)
    FillHistosMC();

  return kTRUE;
}




//________________________________________________________________________
void AliAnalysisTaskMuonCollisionMultiplicity::FillHistosAOD(Int_t triggerClass)
{
  // Fill histos for AOD events
  Int_t nTracks = fAOD->GetNumberOfTracks();
  Int_t nDimuons = fAOD->GetNDimuons();
  
  // Fill histos
  Double_t vertexPosition = fAOD->GetPrimaryVertex()->GetZ();
  Double_t pileUp = !(fAOD->IsPileupFromSPD());
  
  Double_t valuesTrigger[3] = {static_cast<Double_t>(fTrackletMultiplicity), vertexPosition, pileUp};
  ((THnSparseD *)fTriggerList->At(triggerClass))->Fill(valuesTrigger);
  

  // Loop on the muons tracks
  for (Int_t ii = 0; ii < nTracks; ii++){
    AliAODTrack * aodtrack = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(ii));
    if(!aodtrack) AliFatal("Not a standard AOD");
    if (IsUsableMuon(aodtrack))
      {
	Double_t matchTrigger = aodtrack->GetMatchTrigger();
	if (matchTrigger > 1.0)
	  matchTrigger = 1.0;                                            // We don't care what type of trigger it is (low or high pT)
	
	Double_t thetaAbs = (180.0 / TMath::Pi()) * TMath::ATan(aodtrack->GetRAtAbsorberEnd()/505.0);
	Double_t eta = aodtrack->Eta();
	Double_t p = aodtrack->P();
	Double_t pT = aodtrack->Pt();
	// For the p used in the pDCA, we want the mean between the p before the muon go through the absorber (p corrected) and after (p uncorrected)
	// However p uncorrected is not saved in the AODs
	// Instead we define p uncorrected as p corrected minus the mean p lost when a muon go through the absorber
	// p lost for muons (it depends of theta_abs) :
	// 2.0 < theta_abs < 3.0 ---> 2.98 GeV
	// 3.0 < theta_abs < 10.0 --->	2.4 GeV
	// No correction applied otherwise
	Double_t pDCA = p*aodtrack->DCA();
	if (2.0 < thetaAbs && thetaAbs < 3.0)
	  pDCA = (p-2.98/2.0) * aodtrack->DCA();
	if (3.0 < thetaAbs && thetaAbs < 10.0)
	  pDCA = (p-2.4/2.0) * aodtrack->DCA();
	
	Double_t valuesMuon[9] = {static_cast<Double_t>(fTrackletMultiplicity), vertexPosition, pileUp, matchTrigger, thetaAbs, eta, pDCA, static_cast<Double_t>(p), static_cast<Double_t>(pT)};
	((THnSparseD *)fSingleMuonList->At(triggerClass))->Fill(valuesMuon);
      }
  }
  // Loop on Dimuons
  for (Int_t ii = 0; ii < nDimuons; ii++)
    if (fAOD->GetDimuon(ii)->Charge() == 0.0)
      if (IsUsableMuon(fAOD->GetDimuon(ii)->GetMu(0)))
        if (IsUsableMuon(fAOD->GetDimuon(ii)->GetMu(1)))
	  {
	    Double_t matchTrigger1 = fAOD->GetDimuon(ii)->GetMu(0)->GetMatchTrigger();
	    if (matchTrigger1 > 0.0)
	      matchTrigger1 = 1.0;
	    Double_t matchTrigger2 = fAOD->GetDimuon(ii)->GetMu(1)->GetMatchTrigger();
	    if (matchTrigger2 > 0.0)
	      matchTrigger2 = 1.0;
	    
	    Double_t nMatchTrigger = matchTrigger1 + matchTrigger2;
	    
	    Double_t thetaAbs1 = (180.0 / TMath::Pi()) * TMath::ATan(fAOD->GetDimuon(ii)->GetMu(0)->GetRAtAbsorberEnd()/505.0);
	    Double_t thetaAbs2 = (180.0 / TMath::Pi()) * TMath::ATan(fAOD->GetDimuon(ii)->GetMu(1)->GetRAtAbsorberEnd()/505.0);
	    Double_t eta1 = fAOD->GetDimuon(ii)->GetMu(0)->Eta();
	    Double_t eta2 = fAOD->GetDimuon(ii)->GetMu(1)->Eta();
	    
	    Double_t p1 = fAOD->GetDimuon(ii)->GetMu(0)->P();
	    Double_t p2 = fAOD->GetDimuon(ii)->GetMu(1)->P();
	    // See the explanation on how the pDCA is computed in the single muon loop
	    Double_t pDCA1 = p1*fAOD->GetDimuon(ii)->GetMu(0)->DCA();
	    if (2.0 < thetaAbs1 && thetaAbs1 < 3.0)
	      pDCA1 = (p1-2.98/2.0) * fAOD->GetDimuon(ii)->GetMu(0)->DCA();
	    if (3.0 < thetaAbs1 && thetaAbs1 < 10.0)
	      pDCA1 = (p1-2.4/2.0) * fAOD->GetDimuon(ii)->GetMu(0)->DCA();
	    Double_t pDCA2 = p2*fAOD->GetDimuon(ii)->GetMu(1)->DCA();
	    if (2.0 < thetaAbs2 && thetaAbs2 < 3.0)
	      pDCA2 = (p2-2.98/2.0) * fAOD->GetDimuon(ii)->GetMu(1)->DCA();
	    if (3.0 < thetaAbs2 && thetaAbs2 < 10.0)
	      pDCA2 = (p2-2.4/2.0) * fAOD->GetDimuon(ii)->GetMu(1)->DCA();
	    
	    Double_t y = fAOD->GetDimuon(ii)->Y();
	    Double_t p = fAOD->GetDimuon(ii)->P();
	    Double_t pT = fAOD->GetDimuon(ii)->Pt();
	    Double_t M = fAOD->GetDimuon(ii)->M();
	    
	    Double_t valuesDimuon[18] = {static_cast<Double_t>(fTrackletMultiplicity), vertexPosition, pileUp, matchTrigger1, matchTrigger2, static_cast<Double_t>(nMatchTrigger), thetaAbs1, thetaAbs2,
					 eta1, eta2, pDCA1, pDCA2, static_cast<Double_t>(p1), p2, y, p, pT, M};
	    ((THnSparseD *)fDimuonList->At(triggerClass))->Fill(valuesDimuon);
	  }
}



void AliAnalysisTaskMuonCollisionMultiplicity::FillHistosESD(Int_t triggerClass)
{
  // Fill the histos for ESD events 
  Int_t nTracks = fESD->GetNumberOfMuonTracks();

  Double_t vertexPosition = fESD->GetPrimaryVertex()->GetZ();
  Double_t pileUp = !(fESD->IsPileupFromSPD());
  
  Double_t valuesTrigger[3] = {static_cast<Double_t>(fTrackletMultiplicity), vertexPosition, pileUp};
  ((THnSparseD *)fTriggerList->At(triggerClass))->Fill(valuesTrigger);
  

  // Loop on the muons tracks
  for (Int_t ii = 0; ii < nTracks; ii++)
    if (IsUsableMuon(fESD->GetMuonTrack(ii)))
      {
	Double_t matchTrigger1 = fESD->GetMuonTrack(ii)->GetMatchTrigger();
	if (matchTrigger1 > 1.0)
	  matchTrigger1 = 1.0;                                            // We don't care what type of trigger it is (low or high pT)
	
	Double_t thetaAbs1 = (180.0 / TMath::Pi()) * TMath::ATan(fESD->GetMuonTrack(ii)->GetRAtAbsorberEnd()/505.0);
	Double_t eta1 = fESD->GetMuonTrack(ii)->Eta();
	Double_t p1 = fESD->GetMuonTrack(ii)->P();
	Double_t pT1 = fESD->GetMuonTrack(ii)->Pt();

	// For the p used in the pDCA, we want the mean between the p before the muon go through the absorber (p corrected) and after (p uncorrected)
	// These are respectively AliESDMuonTrack::P() and AliESDMuonTrack::PUncorrected()
	Double_t pUncor1 = fESD->GetMuonTrack(ii)->PUncorrected();
	Double_t pDCA1 = (p1+pUncor1) * fESD->GetMuonTrack(ii)->GetDCA() / 2.0;
	
	Double_t valuesMuon[9] = {static_cast<Double_t>(fTrackletMultiplicity), vertexPosition, pileUp, matchTrigger1, thetaAbs1, eta1, pDCA1, p1, pT1};
	((THnSparseD *)fSingleMuonList->At(triggerClass))->Fill(valuesMuon);
	
	// Second loop on muons, to fill the dimuons histos
	for (Int_t jj = ii+1; jj < nTracks; jj++)
	  if (IsUsableMuon(fESD->GetMuonTrack(jj)))
	    if (fESD->GetMuonTrack(ii)->Charge() + fESD->GetMuonTrack(jj)->Charge() == 0.0)
	      {
		Double_t matchTrigger2 = fESD->GetMuonTrack(jj)->GetMatchTrigger();
		if (matchTrigger2 > 0.0)
		  matchTrigger2 = 1.0;
		
		Double_t nMatchTrigger = matchTrigger1 + matchTrigger2;

		Double_t thetaAbs2 = (180.0 / TMath::Pi()) * TMath::ATan(fESD->GetMuonTrack(jj)->GetRAtAbsorberEnd()/505.0);
		Double_t eta2 = fESD->GetMuonTrack(jj)->Eta();
		Double_t p2 = fESD->GetMuonTrack(jj)->P();

		// For the p used in the pDCA, we want the mean between the p before the muon go through the absorber (p corrected) and after (p uncorrected)
		// These are respectively AliESDMuonTrack::P() and AliESDMuonTrack::PUncorrected()
		Double_t pUncor2 = fESD->GetMuonTrack(jj)->PUncorrected();
		Double_t pDCA2 = (p2+pUncor2) * fESD->GetMuonTrack(jj)->GetDCA() / 2.0;
		
		// To compute the p, pT and M of the dimuon, we need a TLorentz vector of the dimuon
		Double_t E = fESD->GetMuonTrack(ii)->E() + fESD->GetMuonTrack(jj)->E();
		Double_t pX = fESD->GetMuonTrack(ii)->Px() + fESD->GetMuonTrack(jj)->Px();
		Double_t pY = fESD->GetMuonTrack(ii)->Py() + fESD->GetMuonTrack(jj)->Py();
		Double_t pZ = fESD->GetMuonTrack(ii)->Pz() + fESD->GetMuonTrack(jj)->Pz();
		TLorentzVector *dimuonVector = new TLorentzVector(pX, pY, pZ, E);
		dimuonVector->SetPxPyPzE(pX, pY, pZ, E);
		
		Double_t y = dimuonVector->Rapidity();
		Double_t p = dimuonVector->P();
		Double_t pT = TMath::Sqrt(pX*pX + pY*pY);
		Double_t M = dimuonVector->M();
		
		Double_t valuesDimuon[18] = {static_cast<Double_t>(fTrackletMultiplicity), vertexPosition, pileUp, matchTrigger1, matchTrigger2, static_cast<Double_t>(nMatchTrigger), thetaAbs1, thetaAbs2,
					     eta1, eta2, pDCA1, pDCA2, p1, p2, y, p, pT, M};
		((THnSparseD *)fDimuonList->At(triggerClass))->Fill(valuesDimuon);
		delete dimuonVector;
	      }
      }
  

  // Since this is the ESD, we are also going to fill a the V0amp vs multiplicity histos
  if (triggerClass == 0)
    {
      AliESDVZERO *v0 = fESD->GetVZEROData();
      Float_t multV0 = 0;
      Float_t multV0A = 0;
      Float_t multV0C = 0;

      for (Int_t ii = 0; ii < 64; ii++)
	{
	  multV0 += v0->GetMultiplicity(ii);
	  if (ii < 32)
	    {
	      multV0C += v0->GetMultiplicityV0C(ii);
	      multV0A += v0->GetMultiplicityV0A(ii);
	    }
	}

      ((TH2D *)fTriggerList->At(2))->Fill(fTrackletMultiplicity, multV0);
      ((TH2D *)fTriggerList->At(3))->Fill(fTrackletMultiplicity, multV0A);
      ((TH2D *)fTriggerList->At(4))->Fill(fTrackletMultiplicity, multV0C);
    }


}


//________________________________________________________________________
void AliAnalysisTaskMuonCollisionMultiplicity::FillHistosMC()
{
  // Fill the histo of the correlation between MC tracks and ESD tracks

  Int_t multiplicityGenerated = 0;
  const AliESDVertex *vertex = fESD->GetPrimaryVertexTracks();

  for (Int_t nn = 0; nn < fMCEvent->GetNumberOfTracks(); nn++)
    {
      AliMCParticle *particle = (AliMCParticle *) fMCEvent->GetTrack(nn);
      Bool_t isGoodMult = kTRUE;
      
      if (particle->Particle()->GetStatusCode() != 1)
	isGoodMult = kFALSE;
      
      if (particle->Charge() == 0)
	isGoodMult = kFALSE;

      if (TMath::Abs(particle->Eta()) > 1.6)
	isGoodMult = kFALSE;

      // Check if the particle is a pion, kaon, proton, electron or muon
      if (TMath::Abs(particle->PdgCode()) != 211 && TMath::Abs(particle->PdgCode()) != 321 && TMath::Abs(particle->PdgCode()) != 2212 &&
	  TMath::Abs(particle->PdgCode()) != 11 && TMath::Abs(particle->PdgCode()) != 13)
	isGoodMult = kFALSE;

      // Check if the distance to vertex is inferior to 1 cm
      Double_t distanceToVertex = TMath::Sqrt((particle->Xv() - vertex->GetX())*(particle->Xv() - vertex->GetX()) + 
					      (particle->Yv() - vertex->GetY())*(particle->Yv() - vertex->GetY()) + 
					      (particle->Zv() - vertex->GetZ())*(particle->Zv() - vertex->GetZ()));
      if (distanceToVertex > 1.0)
	isGoodMult = kFALSE;
	
      if (isGoodMult)
	multiplicityGenerated += 1;
    }

  ((TH2D *)fMonteCarloList->At(0))->Fill(multiplicityGenerated, fTrackletMultiplicity);
}



//________________________________________________________________________
void AliAnalysisTaskMuonCollisionMultiplicity::ComputeMultiplicity()
{
  // Compute the collision multiplicity based on AOD or ESD tracklets

  Int_t multiplicity = 0;
  
  if (fAOD)
    {
      AliAODTracklets *tracklets = fAOD->GetTracklets();
      Int_t nTracklets = tracklets->GetNumberOfTracklets();
      for (Int_t nn = 0; nn < nTracklets; nn++)
	{
	  Double_t theta = tracklets->GetTheta(nn);
	  Double_t eta = -TMath::Log(TMath::Tan(theta/2.0));
	  
	  if (TMath::Abs(eta) < fEtaCut)
	    multiplicity += 1;
	}
    }


  if (fESD)
    {
      // In ESDs, we use the EstimateMultiplicity function
      Int_t multTracklets = 0;
      Int_t multTPC = 0;
      Int_t multITSSA = 0;
      fESD->EstimateMultiplicity(multTracklets, multTPC, multITSSA, fEtaCut);
      multiplicity = multTracklets;
    }

  fTrackletMultiplicity =  multiplicity;
}



//________________________________________________________________________
Bool_t AliAnalysisTaskMuonCollisionMultiplicity::IsUsableMuon(AliAODTrack *track)
{
  // Check if the track is a usable muon track
  // Cuts applied :
  // - is it a muon track?
  // - does it have a pT > 0.0?

  Bool_t isGood = kFALSE;

  if (!track->IsMuonTrack())
    return isGood;

  if (!(track->Pt() > 0.0))
    return isGood;

  isGood = kTRUE;
  return isGood;
}



//________________________________________________________________________
Bool_t AliAnalysisTaskMuonCollisionMultiplicity::IsUsableMuon(AliESDMuonTrack *track)
{
  // Check if the track is a usable muon track
  // Cuts applied :
  // - is it a muon track?
  // - does it have a pT > 0.0?

  Bool_t isGood = kFALSE;

  if (!track->ContainTrackerData())
    return isGood;

  if (!(track->Pt() > 0.0))
    return isGood;

  isGood = kTRUE;
  return isGood;
}



//________________________________________________________________________
void AliAnalysisTaskMuonCollisionMultiplicity::Init()
{
  // Initialize the object

  fTriggerList = new TList();
  fSingleMuonList = new TList();
  fDimuonList = new TList();
  fMonteCarloList = new TList();

  fTriggerList->SetOwner();
  fSingleMuonList->SetOwner();
  fDimuonList->SetOwner();
  fMonteCarloList->SetOwner();



  // Trigger histos
  // dimension 0 : multiplicity of the event
  // dimension 1 : z vertex of the event
  // dimension 2 : is it an event without pile up (0 for no, 1 for yes)?
  Int_t nBinsTrigger[3] =       {  150,    60,   2};
  Double_t minRangeTrigger[3] = {  0.0, -30.0, 0.0};
  Double_t maxRangeTrigger[3] = {150.0,  30.0, 2.0};
  THnSparseD *CINT1B = new THnSparseD ("CINT1B", "CINT1B", 3, nBinsTrigger, minRangeTrigger, maxRangeTrigger);
  THnSparseD *CMUS1B = new THnSparseD ("CMUS1B", "CMUS1B", 3, nBinsTrigger, minRangeTrigger, maxRangeTrigger);
  CINT1B->Sumw2();
  CMUS1B->Sumw2();

  TH2D *CompSPDV0 = new TH2D ("CompSPDV0", "CompSPDV0", 150, 0.0, 150.0, 2000, 0.0, 2000.0);
  CompSPDV0->Sumw2();
  TH2D *CompSPDV0A = new TH2D ("CompSPDV0A", "CompSPDV0A", 150, 0.0, 150.0, 1000, 0.0, 1000.0);
  CompSPDV0A->Sumw2();
  TH2D *CompSPDV0C = new TH2D ("CompSPDV0C", "CompSPDV0C", 150, 0.0, 150.0, 1000, 0.0, 1000.0);
  CompSPDV0C->Sumw2();

  fTriggerList->AddAt(CINT1B, 0);
  fTriggerList->AddAt(CMUS1B, 1);
  fTriggerList->AddAt(CompSPDV0, 2);
  fTriggerList->AddAt(CompSPDV0A, 3);
  fTriggerList->AddAt(CompSPDV0C, 4);




  // Muons histos
  // dimension 0 : multiplicity of the event
  // dimension 1 : z vertex of the event
  // dimension 2 : is it an event without pile up (0 for no, 1 for yes)?
  // dimension 3 : does the muon match the trigger (0 for no, 1 for yes)?
  // dimension 4 : theta_abs of the muon
  // dimension 5 : eta of the muon
  // dimension 6 : p DCA of the muon
  // dimension 7 : p of the muon
  // dimension 8 : pT of the muon

  Int_t nBinsMuon[9] =       {  150,    60,   2,   2,  110,   35,   300,   300,  300};
  Double_t minRangeMuon[9] = {  0.0, -30.0, 0.0, 0.0,  0.0, -5.0,   0.0,   0.0,  0.0};
  Double_t maxRangeMuon[9] = {150.0,  30.0, 2.0, 2.0, 11.0, -1.5, 300.0, 150.0, 30.0};

  THnSparseD *muonCINT1B = new THnSparseD("muonCINT1B", "muonCINT1B", 9, nBinsMuon, minRangeMuon, maxRangeMuon);
  THnSparseD *muonCMUS1B = new THnSparseD("muonCMUS1B", "muonCMUS1B", 9, nBinsMuon, minRangeMuon, maxRangeMuon);
  muonCINT1B->Sumw2();
  muonCMUS1B->Sumw2();

  fSingleMuonList->AddAt(muonCINT1B, 0);
  fSingleMuonList->AddAt(muonCMUS1B, 1);


  // Dimuons histos
  // dimension 0  : multiplicity of the event
  // dimension 1  : z range (0 for no, 1 for yes)?
  // dimension 2  : is it an event without pile up (0 for no, 1 for yes)?
  // dimension 3  : does the first muon match the trigger (0 for no, 1 for yes)?
  // dimension 4  : does the second muon match the trigger (0 for no, 1 for yes)?
  // dimension 5  : number of muons matching the trigger in the dimuon
  // dimension 6  : theta_abs of the first muon
  // dimension 7  : theta_abs of the second muon
  // dimension 8  : eta of the first muon
  // dimension 9  : eta of the second muon
  // dimension 10 : p DCA of the first muon
  // dimension 11 : p DCA of the second muon
  // dimension 12 : p of the first muon
  // dimension 13 : p of the second muon
  // dimension 14 : y of the dimuon
  // dimension 15 : p of the dimuon
  // dimension 16 : pT of the dimuon
  // dimension 17 : invariant mass of the dimuon

  Int_t nBinsDimuon[18] =       {  150,    60,   2,   2,   2,   3,  110,  110,   35,   35,   200,   200,   300,   300,   35,   300,  300,  375};
  Double_t minRangeDimuon[18] = {  0.0, -30.0, 0.0, 0.0, 0.0, 0.0,  0.0,  0.0, -5.0, -5.0,   0.0,   0.0,   0.0,   0.0, -5.0,   0.0,  0.0,  0.0};
  Double_t maxRangeDimuon[18] = {150.0,  30.0, 2.0, 2.0, 2.0, 3.0, 11.0, 11.0, -1.5, -1.5, 600.0, 600.0, 150.0, 150.0, -1.5, 300.0, 30.0, 15.0};

  THnSparseD *dimuonCINT1B = new THnSparseD("dimuonCINT1B", "dimuonCINT1B", 18, nBinsDimuon, minRangeDimuon, maxRangeDimuon);
  THnSparseD *dimuonCMUS1B = new THnSparseD("dimuonCMUS1B", "dimuonCMUS1B", 18, nBinsDimuon, minRangeDimuon, maxRangeDimuon);

  dimuonCINT1B->Sumw2();
  dimuonCMUS1B->Sumw2();

  fDimuonList->AddAt(dimuonCINT1B, 0);
  fDimuonList->AddAt(dimuonCMUS1B, 1);

  // MonteCarlo Histo
  TH2D *correlGenerReco = new TH2D("correlGenerReco", "correlGenerReco", 250, 0.0, 250.0, 250, 0.0, 250.0);
  correlGenerReco->GetXaxis()->SetTitle("N ch gener");
  correlGenerReco->GetYaxis()->SetTitle("N reco tracklets");

  correlGenerReco->Sumw2();

  fMonteCarloList->AddAt(correlGenerReco, 0);

  fIsInit = kTRUE;
}





//________________________________________________________________________
void AliAnalysisTaskMuonCollisionMultiplicity::Terminate(Option_t */*option*/)
{
//Terminate analysis
  
  fTriggerList = (TList *) GetOutputData(0);
  fSingleMuonList = (TList *) GetOutputData(1);
  fDimuonList = (TList *) GetOutputData(2);
}
