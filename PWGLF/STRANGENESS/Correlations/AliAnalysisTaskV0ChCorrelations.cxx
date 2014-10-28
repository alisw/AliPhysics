/**************************************************************************
 * Copyright(c) 1998-2013, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: Marek Bombara                                                  *
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

/* The task selects candidates for K0s, Lambdas and AntiLambdas (trigger particles)
 * and calculates correlations with charged unidentified particles (associated particles) in phi and eta. 
 * The task works with AOD (with or without MC info) events only and containes also mixing for acceptance corrections.
 * Last update edited by Marek Bombara, August 2013, Marek.Bombara@cern.ch
 */

#include <TCanvas.h>
#include <TList.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <THnSparse.h>
#include <TObjArray.h>
#include <TObject.h>

#include "AliLog.h"
#include "AliAnalysisManager.h"
#include "AliAODTrack.h"
#include "AliAODEvent.h"
#include "AliAODv0.h"
#include "AliAODcascade.h"
#include "AliAODVertex.h"

#include "AliMCEvent.h"
#include "AliMCVertex.h" 
#include "AliGenEventHeader.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliAnalyseLeadingTrackUE.h"

#include "AliAODPid.h"
#include "AliPIDResponse.h"
#include "AliEventPoolManager.h"
#include "AliCentrality.h"

#include "AliAnalysisTaskV0ChCorrelations.h"
#include "AliPhysicsSelectionTask.h"

#include <AliMultiInputEventHandler.h>
#include <AliMixInputEventHandler.h>

ClassImp(AliAnalysisTaskV0ChCorrelations)
ClassImp(AliV0ChBasicParticle)

//________________________________________________________________________
AliAnalysisTaskV0ChCorrelations::AliAnalysisTaskV0ChCorrelations(const char *name) // All data members should be initialised here
   : AliAnalysisTaskSE(name),
     fAnalysisMC(kFALSE),
	 fFillMixed(kTRUE),
	 fMixingTracks(50000),
	 fPoolMgr(0x0),
     
     fDcaDToPV(0.5),
	 fDcaV0D(0.1),
	 fWithChCh(kFALSE),
	 fOStatus(1),

	 fOutput(0),
	 fPIDResponse(0),
     fHistCentVtx(0),
     fHistMultiMain(0),
	 fHistMassK0(0),
	 fHistMassLambda(0),
	 fHistMassAntiLambda(0),
 
	 fHistdPhidEtaSib(0),
	 fHistdPhidEtaMix(0),
	 fHistNTrigSib(0),

   	 fHistMCPtCentTrig(0),
     fHistRCPtCentTrig(0),
     fHistMCPtCentAs(0),
     fHistRCPtCentAs(0),
     fHistRCPtCentAll(0),
	 
	 fHistTemp(0),
	 fHistTemp2(0)// The last in the above list should not have a comma after it
{
   // Constructor
   // Define input and output slots here (never in the dummy constructor)
   // Input slot #0 works with a TChain - it is connected to the default input container
   // Output slot #1 writes into a TH1 container
   DefineOutput(1, TList::Class());                                            // for output list
}

//________________________________________________________________________
AliAnalysisTaskV0ChCorrelations::~AliAnalysisTaskV0ChCorrelations()
{
   // Destructor. Clean-up the output list, but not the histograms that are put inside
   // (the list is owner and will clean-up these histograms). Protect in PROOF case.
   if (fOutput && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
      delete fOutput;
   }
}

//________________________________________________________________________
void AliAnalysisTaskV0ChCorrelations::UserCreateOutputObjects()
{
   // Create histograms

   fOutput = new TList();
   fOutput->SetOwner();  // IMPORTANT!
   // defining bins for centrality
   Int_t nCentralityBins  = 9;
   Double_t centBins[] = {0.,10.,20.,30.,40.,50.,60.,70.,80.,90.};
   const Double_t* centralityBins = centBins;
   // defining bins for Z vertex
   Int_t nZvtxBins  = 7;
   //Double_t vertexBins[] = {-10.,-8.,-6.,-4.,-2.,0.,2.,4.,6.,8.,10.};
   Double_t vertexBins[] = {-7.,-5.,-3.,-1.,1.,3.,5.,7.};
   const Double_t* zvtxBins = vertexBins;
   // pt bins of associated particles for the analysis
   Int_t nPtBins = 7;
   const Double_t PtBins[8] = {2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0}; 
   //Int_t nPtBins = 1;
   //const Double_t PtBins[2] = {3.0,15.0}; 
   // pt bins of trigger particles for the analysis
   Int_t nPtBinsV0 = 11;
   //const Double_t PtBinsV0[2] = {6.0,15.0}; 
   const Double_t PtBinsV0[12] = {4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0}; 
   // V0 candidate: 1 - K0sig, 2 - Lamsig, 3 - Alamsig, 4 - K0bg, 5 - Lambg, 6 - Alambg
   Int_t nTrigC = 7;
   const Double_t TrigC[8] = {0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5};
   
   // Create general histograms 
   fHistCentVtx = new TH2F("fHistCentVtx", "Centrality vs. Z vertex", nCentralityBins, centralityBins, nZvtxBins, zvtxBins);
   fHistMultiMain = new TH1F("fHistMultiMain", "Multiplicity of main events", 2000, 0, 2000);

   const Int_t mrBins[3] = {nPtBinsV0, nCentralityBins, nTrigC};
   const Double_t mrMin[3] = {PtBinsV0[0], centralityBins[0], TrigC[0]};
   const Double_t mrMax[3] = {PtBinsV0[11], centralityBins[9], TrigC[6]};

   // Create histograms for reconstruction track and V0 efficiency
   fHistMCPtCentTrig = new THnSparseF("fHistMCPtCentTrig", "MC Pt vs. Cent. Trig", 3, mrBins, mrMin, mrMax);
   fHistRCPtCentTrig = new THnSparseF("fHistRCPtCentTrig", "Rec Pt vs. Cent. Trig", 3, mrBins, mrMin, mrMax);

   // pt bins of associated particles for the efficiency
   Int_t nPtBinsAs = 13;
   const Double_t PtBinsAs[14] = {2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0};

   fHistMCPtCentAs = new TH2D("fHistMCPtCentAs", "MC Pt vs. Cent. Assoc", nPtBinsAs, PtBinsAs[0], PtBinsAs[13], nCentralityBins, centBins[0], centBins[9]);
   fHistRCPtCentAs = new TH2D("fHistRCPtCentAs", "Rec Pt vs. Cent. Assoc", nPtBinsAs, PtBinsAs[0], PtBinsAs[13], nCentralityBins, centBins[0], centBins[9]);
   fHistRCPtCentAll = new TH2D("fHistRCPtCentAll", "Rec Pt vs. Cent. All Prim+Sec", nPtBinsAs, PtBinsAs[0], PtBinsAs[13], nCentralityBins, centBins[0], centBins[9]);

   // defining bins for mass distributions
   Int_t nBins = 200;
   Double_t mk0Min = 0.40;
   Double_t mk0Max = 0.60;
   Double_t mlaMin = 1.07;
   Double_t mlaMax = 1.17;
   Double_t malMin = 1.07;
   Double_t malMax = 1.17;
   
   const Int_t spBins[3] = {nBins, nPtBinsV0, nCentralityBins};
   const Double_t spMinK0[3] = {mk0Min, PtBinsV0[0], centralityBins[0]};
   const Double_t spMaxK0[3] = {mk0Max, PtBinsV0[11], centralityBins[9]};
   const Double_t spMinLa[3] = {mlaMin, PtBinsV0[0], centralityBins[0]};
   const Double_t spMaxLa[3] = {mlaMax, PtBinsV0[11], centralityBins[9]};
   const Double_t spMinAl[3] = {malMin, PtBinsV0[0], centralityBins[0]};
   const Double_t spMaxAl[3] = {malMax, PtBinsV0[11], centralityBins[9]};
   // Create mass histograms
   fHistMassK0 = new THnSparseF("fHistMassK0","V0 mass for K0 hypothesis", 3, spBins, spMinK0, spMaxK0);
   fHistMassLambda = new THnSparseF("fHistMassLambda","V0 mass for Lambda hypothesis", 3, spBins, spMinLa, spMaxLa);
   fHistMassAntiLambda = new THnSparseF("fHistMassAntiLambda","V0 mass for AntiLambda hypothesis", 3, spBins, spMinAl, spMaxAl);
   
   // defining bins for dPhi distributions
   const Int_t nbPhiBins = 72;
   //const Double_t kPi = TMath::Pi();
   //Double_t PhiMin = -kPi/2.;
   //Double_t PhiMax = -kPi/2. + 2*kPi;
   Double_t PhiMin = -1.57;
   Double_t PhiMax = 4.71;
   Double_t PhiBins[nbPhiBins+1] = {0.};
   PhiBins[0] = PhiMin;
   for (Int_t i=0; i<nbPhiBins; i++) { PhiBins[i+1] = PhiBins[i] + (PhiMax - PhiMin)/nbPhiBins; }

   // defining bins for dEta distributions
   const Int_t nbEtaBins = 40;
   Double_t EtaMin = -2.0;
   Double_t EtaMax = 2.0;
   Double_t EtaBins[nbEtaBins+1] = {0.};
   EtaBins[0] = EtaMin;
   for (Int_t i=0; i<nbEtaBins; i++) { EtaBins[i+1] = EtaBins[i] + (EtaMax - EtaMin)/nbEtaBins; }

   const Int_t corBins[7] = {nbPhiBins, nbEtaBins, nPtBinsV0, nPtBins, nCentralityBins, nZvtxBins, nTrigC};
   const Double_t corMin[7] = {PhiBins[0], EtaBins[0], PtBinsV0[0], PtBins[0], centralityBins[0], zvtxBins[0], TrigC[0]};
   const Double_t corMax[7] = {PhiBins[72], EtaBins[40], PtBinsV0[11], PtBins[7], centralityBins[9], zvtxBins[7], TrigC[7]};
   // Create correlation histograms
   fHistdPhidEtaSib = new THnSparseF("fHistdPhidEtaSib","dPhi vs. dEta siblings", 7, corBins, corMin, corMax); 
   fHistdPhidEtaMix = new THnSparseF("fHistdPhidEtaMix","dPhi vs. dEta mixed", 7, corBins, corMin, corMax);

   // Create histograms for counting the numbers of trigger particles
   const Int_t corNTrigBins[5] = {nPtBinsV0, nCentralityBins, nZvtxBins, nTrigC};
   const Double_t corNTrigMin[5] = {PtBinsV0[0], centBins[0], vertexBins[0], TrigC[0]};
   const Double_t corNTrigMax[5] = {PtBinsV0[11], centBins[9], vertexBins[7], TrigC[7]};
   fHistNTrigSib = new THnSparseF("fHistNTrigSib","Number trigger sib", 4, corNTrigBins, corNTrigMin, corNTrigMax);

   // Histograms for debugging
   fHistTemp = new TH1D("fHistTemp", "Temporary", 100, -10., 10.);
   fHistTemp2 = new TH1D("fHistTemp2", "Temporary2", 100, -10., 10.);

   fOutput->Add(fHistCentVtx);

   fOutput->Add(fHistMultiMain);
   fOutput->Add(fHistMassK0);
   fOutput->Add(fHistMassLambda);
   fOutput->Add(fHistMassAntiLambda);
   
   fOutput->Add(fHistdPhidEtaSib);
   fOutput->Add(fHistdPhidEtaMix);
   fOutput->Add(fHistNTrigSib);

   fOutput->Add(fHistMCPtCentTrig);
   fOutput->Add(fHistRCPtCentTrig);
   fOutput->Add(fHistMCPtCentAs);
   fOutput->Add(fHistRCPtCentAs);
   fOutput->Add(fHistRCPtCentAll);
   
   fOutput->Add(fHistTemp);
   fOutput->Add(fHistTemp2);

   PostData(1, fOutput); // Post data for ALL output slots >0 here, to get at least an empty histogram

   // Settings for event mixing -------------------------------------
   Int_t trackDepth = fMixingTracks;
   Int_t poolSize   = 200;  // Maximum number of events, ignored in the present implemented of AliEventPoolManager

   fPoolMgr = new AliEventPoolManager(poolSize, trackDepth, nCentralityBins, centBins, nZvtxBins, vertexBins);
   //----------------------------------------------
}

//________________________________________________________________________
void AliAnalysisTaskV0ChCorrelations::Terminate(Option_t *)
{
   // Draw result to screen, or perform fitting, normalizations
   // Called once at the end of the query

   fOutput = dynamic_cast<TList*>(GetOutputData(1));
   if (!fOutput) { AliError("Could not retrieve TList fOutput"); return; }

   // NEW HISTO should be retrieved from the TList container in the above way,
   // so it is available to draw on a canvas such as below
}

//_________________________________________________________________________
void AliAnalysisTaskV0ChCorrelations::UserExec(Option_t *)
{
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler *inEvMain = (AliInputEventHandler*)mgr->GetInputEventHandler();

    // physics selection
    UInt_t maskIsSelected = inEvMain->IsEventSelected();
   // 2010 data trigger selection
    //Bool_t isSelected = (maskIsSelected & AliVEvent::kMB);
   // 2011 data trigger selection
    Bool_t isSelected = ((maskIsSelected & AliVEvent::kMB) || (maskIsSelected & AliVEvent::kCentral) || (maskIsSelected & AliVEvent::kSemiCentral));
    if (!isSelected) return;
	// calculating the event types
    if (maskIsSelected & AliVEvent::kMB) fHistTemp->Fill(2);
    if (maskIsSelected & AliVEvent::kCentral) fHistTemp->Fill(4);
    if (maskIsSelected & AliVEvent::kSemiCentral) fHistTemp->Fill(6);

    AliAODEvent* aod = (AliAODEvent*)inEvMain->GetEvent();
    fPIDResponse = inEvMain->GetPIDResponse();

  // pt intervals for trigger particles  
  const Double_t kPi = TMath::Pi();
  Double_t PtTrigMin = 4.0;
  Double_t PtTrigMax = 15.0;
  // pt interval for associated particles
  Double_t PtAssocMin = 2.0;
  fHistMultiMain->Fill(aod->GetNumberOfTracks());

  // Vertex cut
  Double_t cutPrimVertex = 7.0;
  AliAODVertex *myPrimVertex = aod->GetPrimaryVertex();
  if (!myPrimVertex) return;
  if ( ( TMath::Abs(myPrimVertex->GetZ()) ) >= cutPrimVertex) return ;
  Double_t lPVx = myPrimVertex->GetX();
  Double_t lPVy = myPrimVertex->GetY();
  Double_t lPVz = myPrimVertex->GetZ();

  if (TMath::Abs(lPVx)<10e-5 && TMath::Abs(lPVy)<10e-5 && TMath::Abs(lPVz)<10e-5) return;

  // Centrality definition
  Double_t lCent = 0.0;
  AliCentrality *centralityObj = 0;
  centralityObj = ((AliVAODHeader*)aod->GetHeader())->GetCentralityP();
  lCent = centralityObj->GetCentralityPercentile("V0M");
  if ((lCent < 0.)||(lCent > 90.)) return;
  fHistCentVtx->Fill(lCent,lPVz);

//=========== MC loop ===============================
  if (fAnalysisMC)
  {
    AliAODMCHeader *aodMCheader = (AliAODMCHeader*)aod->FindListObject(AliAODMCHeader::StdBranchName());
    Float_t vzMC = aodMCheader->GetVtxZ();
    if (TMath::Abs(vzMC) >= 7.) return;
    //retrieve MC particles from event
    TClonesArray *mcArray = (TClonesArray*)aod->FindListObject(AliAODMCParticle::StdBranchName());
    if(!mcArray){
        Printf("No MC particle branch found");
        return;
    }
    
    Int_t nMCAllTracks = mcArray->GetEntriesFast();
    // new tracks array - without injected signal
    TObjArray * mcTracks = new TObjArray;
    //selectedMCTracks->SetOwner(kTRUE);
  
    for (Int_t i = 0; i < nMCAllTracks; i++)
    { 
        AliAODMCParticle *mcTrack = (AliAODMCParticle*)mcArray->At(i);
        if (!mcTrack) {
            Error("ReadEventAODMC", "Could not receive particle %d", i);
            continue;
        }
        mcTracks->Add(mcTrack);
    }

    // get labels for injected particles -------
    TObject* mc = mcArray;
    //TObjArray * mcTracks = mcArray;
    Int_t skipParticlesAbove = 0;
    //AliGenEventHeader* eventHeader = 0;
    //Int_t headers = 0;
    //headers = aodMCheader->GetNCocktailHeaders();
    //eventHeader = aodMCheader->GetCocktailHeader(0);
    AliGenEventHeader* eventHeader = aodMCheader->GetCocktailHeader(0);
    Int_t headers = aodMCheader->GetNCocktailHeaders();
    if (!eventHeader)
    {
      // We avoid AliFatal here, because the AOD productions sometimes have events where the MC header is missing
      // (due to unreadable Kinematics) and we don't want to loose the whole job because of a few events
      AliError("First event header not found. Skipping this event.");
      return;
    }
    skipParticlesAbove = eventHeader->NProduced();
    //cout << "ski label = " << skipParticlesAbove << endl;
    AliInfo(Form("Injected signals in this event (%d headers). Keeping events of %s.Will skip particles/tracks above %d.", headers, eventHeader->ClassName(), skipParticlesAbove));
    //cout << "before MC compressing = " << mcTracks->GetEntriesFast() << endl;
    RemovingInjectedSignal(mcTracks,mc,skipParticlesAbove);
    //cout << "after MC compressing = " << mcTracks->GetEntriesFast() << endl;
    // -------------

    Int_t nMCTracks = mcTracks->GetEntriesFast();
    //cout << "number of MC tracks = " << nMCTracks << endl;
    for (Int_t iMC = 0; iMC<nMCTracks; iMC++)
    {
      AliAODMCParticle *mcTrack = (AliAODMCParticle*)mcTracks->At(iMC);
        if (!mcTrack) {
            Error("ReadEventAODMC", "Could not receive particle %d", iMC);
            continue;
      }
      // track part
      Double_t mcTrackEta = mcTrack->Eta();
      Double_t mcTrackPt = mcTrack->Pt();
      Bool_t TrIsPrim = mcTrack->IsPhysicalPrimary();
      Bool_t TrEtaMax = TMath::Abs(mcTrackEta)<0.9;
      Bool_t TrPtMin = mcTrackPt>PtAssocMin;
      Bool_t TrCharge = (mcTrack->Charge())!=0;
      if (TrIsPrim && TrEtaMax && TrPtMin && TrCharge) fHistMCPtCentAs->Fill(mcTrackPt,lCent);
      // V0 part
      Int_t mcMotherPdg = 0;
	  Int_t mcPartPdg = mcTrack->GetPdgCode();

      // Keep only K0s, Lambda and AntiLambda
	  if ((mcPartPdg != 310) && (mcPartPdg != 3122) && (mcPartPdg != (-3122))) continue;

      //cout << " mc pdg is " << mcPartPdg << endl;

      Bool_t IsK0 = mcPartPdg==310;
      Bool_t IsLambda = mcPartPdg==3122;
      Bool_t IsAntiLambda = mcPartPdg==-3122;
	  Bool_t IsSigma = kFALSE;
	  Int_t mcMotherLabel = mcTrack->GetMother();
      AliAODMCParticle *mcMother = (AliAODMCParticle*)mcArray->At(mcMotherLabel);
	  if (mcMotherLabel < 0) {mcMotherPdg = 0;} else {mcMotherPdg = mcMother->GetPdgCode();}
	  //if ((mcMotherLabel >= 0) && mcMother) 
	  //if ((mcMotherLabel >= 0) && mcMother) 
	  //{
      	//  Bool_t IsFromCascade = ((mcMotherPdg==3312)||(mcMotherPdg==-3312)||(mcMotherPdg==3334)||(mcMotherPdg==-3334));
          Bool_t IsFromSigma = ((mcMotherPdg==3212)||(mcMotherPdg==-3212));
          IsFromSigma = IsFromSigma || ((mcMotherPdg==3224)||(mcMotherPdg==-3224));
          IsFromSigma = IsFromSigma || ((mcMotherPdg==3214)||(mcMotherPdg==-3214));
          IsFromSigma = IsFromSigma || ((mcMotherPdg==3114)||(mcMotherPdg==-3114));
		  if ((IsFromSigma) && (mcMother->IsPhysicalPrimary()) && (IsLambda || IsAntiLambda)) IsSigma = kTRUE;
      	  Double_t mcRapidity = mcTrack->Y();
      	  Bool_t V0RapMax = TMath::Abs(mcRapidity)<0.75;
          Bool_t PtInterval = ((mcTrackPt>PtTrigMin)&&(mcTrackPt<PtTrigMax));
          //Bool_t PtInterval = kTRUE;
          IsK0 = IsK0 && (mcTrack->IsPhysicalPrimary());
          IsLambda = IsLambda && (mcTrack->IsPhysicalPrimary() || IsSigma);
          IsAntiLambda = IsAntiLambda && (mcTrack->IsPhysicalPrimary() || IsSigma);
      	  Double_t mcK0[3] = {mcTrackPt, lCent, 1};
      	  Double_t mcLa[3] = {mcTrackPt, lCent, 2};
      	  Double_t mcAl[3] = {mcTrackPt, lCent, 3};
      	  if (IsK0 && V0RapMax && PtInterval) fHistMCPtCentTrig->Fill(mcK0); 
      	  if (IsLambda && V0RapMax && PtInterval) fHistMCPtCentTrig->Fill(mcLa);
      	  if (IsAntiLambda && V0RapMax && PtInterval) fHistMCPtCentTrig->Fill(mcAl);
	   //}
    }
    // ------- access the real data -----------
    Int_t nTracks = aod->GetNumberOfTracks();
    // new tracks array - without injected signal
    TObjArray * selectedMCTracks = new TObjArray;
    //selectedMCTracks->SetOwner(kTRUE);
  
    for (Int_t i = 0; i < nTracks; i++)
    {
        AliAODTrack* tr = dynamic_cast<AliAODTrack*>(aod->GetTrack(i));
        if(!tr) AliFatal("Not a standard AOD"); 
        selectedMCTracks->Add(tr);
    }
   // cout << "before compressing = " << selectedMCTracks->GetEntriesFast() << endl;
    RemovingInjectedSignal(selectedMCTracks,mc,skipParticlesAbove);
    //AliAnalyseLeadingTrackUE::RemoveInjectedSignals(selectedMCTracks,mc,skipParticlesAbove);
    //RemoveInjectedSignals(selectedMCTracks,mc,skipParticlesAbove);
   // cout << "after compressing = " << selectedMCTracks->GetEntriesFast() << endl;
    //-----------------------
    //cout << "before getting label -4" << endl;
    Int_t nRecTracks = selectedMCTracks->GetEntriesFast();
    //cout << "before getting label -3" << endl;
    for (Int_t i = 0; i < nRecTracks; i++)
    {
      //AliAODTrack* tras = dynamic_cast<AliAODTrack*>(aod->GetTrack(i));
      AliAODTrack* tras = (AliAODTrack*)selectedMCTracks->At(i);
      //cout << "before getting label -2" << endl;
      //cout << " and i = " << i << endl;
      //cout << "pt = " << tras->Pt() << " and i = " << i << endl;
      if ((tras->Pt())<PtAssocMin) continue; 
      //cout << "before getting label -1" << endl;
      if (!(IsMyGoodPrimaryTrack(tras))) continue;
      //cout << "before getting label 0" << endl;
      Int_t AssocLabel = tras->GetLabel();
      if (AssocLabel<=0) continue;
      //cout << "before getting label 1" << endl;
      Bool_t isPhyPrim = static_cast<AliAODMCParticle*>(mcArray->At(tras->GetLabel()))->IsPhysicalPrimary();
      //Bool_t isPhyPrim = static_cast<AliAODMCParticle*>(mcArray->At(tras->GetLabel()))->IsPhysicalPrimary();
      //if(!(static_cast<AliAODMCParticle*>(mcArray->At(tras->GetLabel()))->IsPhysicalPrimary())) continue;
	     //cout << "before getting label 2" << endl;
    Double_t mcPt = static_cast<AliAODMCParticle*>(mcArray->At(tras->GetLabel()))->Pt();
	     //cout << "before getting label 3" << endl;
    //Double_t mcEta = static_cast<AliAODMCParticle*>(mcArray->At(tras->GetLabel()))->Eta();
	  //if (mcPt<PtAssocMin) continue;
	  //if (TMath::Abs(mcEta)>0.8) continue;
      //Double_t trPt = tras->Pt();
      fHistRCPtCentAll->Fill(mcPt,lCent);
      if (isPhyPrim) fHistRCPtCentAs->Fill(mcPt,lCent);
    }
    // ------- end of real data access, for V0s see the main V0 loop -------- 
  }
//============= End of MC loop ======================
	
	// Track selection loop
	//--------------------------------
	Int_t nTracks = aod->GetNumberOfTracks();
	// new tracks array
	TObjArray * selectedTracks = new TObjArray;
	selectedTracks->SetOwner(kTRUE);
	
	Bool_t ChChWith = GetWithChCh();
    TObjArray * selectedChargedTriggers = new TObjArray;
    selectedChargedTriggers->SetOwner(kTRUE);
    for (Int_t i = 0; i < nTracks; i++)
    {
        AliAODTrack* tr = dynamic_cast<AliAODTrack*>(aod->GetTrack(i));
        if(!tr) AliFatal("Not a standard AOD");
        if ((tr->Pt())<PtAssocMin) continue;
        if (!(IsMyGoodPrimaryTrack(tr))) continue;
        selectedTracks->Add(tr);
        // saving the Charged trigger particles
        if ((tr->Pt()>4.)&&(tr->Pt()<15.))
        {
           selectedChargedTriggers->Add(new AliV0ChBasicParticle(tr->Eta(), tr->Phi(), tr->Pt(), 7));
        }
    }

	//---------------------------------
	
	// V0 loop for reconstructed event
    TObjArray * selectedV0s = new TObjArray;
	selectedV0s->SetOwner(kTRUE);

	Int_t nV0s = aod->GetNumberOfV0s();

	for (Int_t i = 0; i < nV0s; i++)
	{ // start of V0 slection loop
		AliAODv0* aodV0 = dynamic_cast<AliAODv0 *>(aod->GetV0(i));
		if (!aodV0) {
         AliError(Form("ERROR: Could not retrieve aodv0 %d", i));
         continue;
        }
		if (aodV0->GetOnFlyStatus()) fHistTemp->Fill(6.);
		if (!aodV0->GetOnFlyStatus()) fHistTemp->Fill(8.);
		//cout << "pt of v0: " << aodV0->Pt() << endl;
		if (((aodV0->Pt())<PtTrigMin)||((aodV0->Pt())>PtTrigMax)) continue;
        // get daughters
        const AliAODTrack *myTrackPos;
        const AliAODTrack *myTrackNeg;
   	    AliAODTrack *myTrackNegTest=dynamic_cast<AliAODTrack *>(aodV0->GetDaughter(1));
        AliAODTrack *myTrackPosTest=dynamic_cast<AliAODTrack *>(aodV0->GetDaughter(0));

        if (!myTrackPosTest || !myTrackNegTest) {
           Printf("strange analysis::UserExec:: Error:Could not retreive one of the daughter track\n");
           continue;
        }

        if( myTrackPosTest->Charge() ==1){
           myTrackPos = myTrackPosTest;
           myTrackNeg = myTrackNegTest;
        }

        if( myTrackPosTest->Charge() ==-1){
            myTrackPos = myTrackNegTest;
            myTrackNeg = myTrackPosTest;
        }

        // effective mass calculations for each hypothesis
		Double_t lInvMassK0 = aodV0->MassK0Short();
        Double_t lInvMassAntiLambda = aodV0->MassAntiLambda();
        Double_t lInvMassLambda = aodV0->MassLambda();

		// calculations for c*tau cut--------------------------------------
	//	Double_t lDVx = aodV0->GetSecVtxX();
    //    Double_t lDVy = aodV0->GetSecVtxY();
    //    Double_t lDVz = aodV0->GetSecVtxZ();
    //    const Double_t kLambdaMass = 1.115683;
    //    const Double_t kK0Mass = 0.497648;
    //    Double_t cutcTauLam = 3*7.89;
    //    Double_t cutcTauK0 = 3*2.68;
     //   Double_t lV0DecayLength = TMath::Sqrt(TMath::Power(lDVx - lPVx,2) + TMath::Power(lDVy- lPVy,2) + TMath::Power(lDVz - lPVz,2 ));
      //  Double_t lPV0 = TMath::Sqrt((aodV0->Pt())*(aodV0->Pt())+(aodV0->Pz())*(aodV0->Pz()));
     //   Double_t lcTauLam = (lV0DecayLength*kLambdaMass)/lPV0;
     //   Double_t lcTauK0 = (lV0DecayLength*kK0Mass)/lPV0;
        // sc - standard cuts
		//Bool_t cutK0sc = (lcTauK0<cutcTauK0);
        //Bool_t cutLambdasc = (lcTauLam<cutcTauLam);
        //Bool_t cutAntiLambdasc = (lcTauLam<cutcTauLam);
		Bool_t cutK0sc = kTRUE;
        Bool_t cutLambdasc = kTRUE;
        Bool_t cutAntiLambdasc = kTRUE;

		//------------------------------------------------

		// preparations for daughter's PID cut------------------------------
        Float_t nSigmaPosPion   = 0.;
        Float_t nSigmaNegPion   = 0.;
        Float_t nSigmaPosProton = 0.;
        Float_t nSigmaNegProton = 0.;

        const AliAODPid *pPid = myTrackPos->GetDetPid();
        const AliAODPid *nPid = myTrackNeg->GetDetPid();

        if (pPid)
        {
            Double_t pdMom = pPid->GetTPCmomentum();
            if (pdMom<1.)
            {
                nSigmaPosPion = fPIDResponse->NumberOfSigmasTPC(myTrackPos, AliPID::kPion);
                nSigmaPosProton = fPIDResponse->NumberOfSigmasTPC(myTrackPos, AliPID::kProton);
            }
        }

        if (nPid)
        {
            Double_t ndMom = nPid->GetTPCmomentum();
            if (ndMom<1.)
            {
                nSigmaNegPion = fPIDResponse->NumberOfSigmasTPC(myTrackPos, AliPID::kPion);
                nSigmaNegProton = fPIDResponse->NumberOfSigmasTPC(myTrackPos, AliPID::kProton);
            }
        }
		Bool_t bpPion = TMath::Abs(nSigmaPosPion) <= 3.;
        Bool_t bpProton = TMath::Abs(nSigmaPosProton) <= 3.;
        Bool_t bnPion = TMath::Abs(nSigmaNegPion) <= 3.;
        Bool_t bnProton = TMath::Abs(nSigmaNegProton) <= 3.;

        Bool_t cutK0Pid = (bpPion && bnPion);
        Bool_t cutLambdaPid = (bpProton && bnPion);
        Bool_t cutAntiLambdaPid = (bpPion && bnProton);
        //--------------------------------------------------
        // mass cuts
        Bool_t cutMassLambda = ((lInvMassLambda>1.10) && (lInvMassLambda<1.13));
        Bool_t cutMassAntiLambda = ((lInvMassAntiLambda>1.10) && (lInvMassAntiLambda<1.13));
        Bool_t cutMassK0 = (lInvMassK0>0.47) && (lInvMassK0<0.53);

        cutK0sc = cutK0Pid && (!cutMassLambda) && (!cutMassAntiLambda);
		cutLambdasc = cutLambdaPid && (!cutMassK0); 
		cutAntiLambdasc = cutAntiLambdaPid && (!cutMassK0);
        // special cut related to AP diagram for K0s
		Bool_t k0APcut = (aodV0->PtArmV0()>(TMath::Abs(0.2*aodV0->AlphaV0())));
		cutK0sc = cutK0sc && k0APcut;
        // fill the mass histograms

        Int_t oStatus = GetOStatus();

        if (!IsMyGoodV0(aod,aodV0,myTrackPos,myTrackNeg,oStatus)) continue;
		Double_t spK0[3] = {lInvMassK0, aodV0->Pt(), lCent};
		Double_t spLa[3] = {lInvMassLambda, aodV0->Pt(), lCent};
		Double_t spAl[3] = {lInvMassAntiLambda, aodV0->Pt(), lCent};
        if (cutK0sc) fHistMassK0->Fill(spK0);
		if (cutLambdasc) fHistMassLambda->Fill(spLa);
		if (cutAntiLambdasc) fHistMassAntiLambda->Fill(spAl);
        // select final V0s for correlation, selected background to study its contribution to correlation function
        // the values for signal might change in the future ...
        Bool_t K0Signal = (lInvMassK0>0.48)&&(lInvMassK0<0.52);
        Bool_t K0Bckg = ((lInvMassK0>0.40)&&(lInvMassK0<0.44)) || ((lInvMassK0>0.56)&&(lInvMassK0<0.60));

        Bool_t LamSignal = (lInvMassLambda>1.108)&&(lInvMassLambda<1.125);
        Bool_t LamBckg = ((lInvMassLambda>1.090)&&(lInvMassLambda<1.100)) || ((lInvMassLambda>1.135)&&(lInvMassLambda<1.145));

        Bool_t ALamSignal = (lInvMassAntiLambda>1.108)&&(lInvMassAntiLambda<1.125);
        Bool_t ALamBckg = ((lInvMassAntiLambda>1.090)&&(lInvMassAntiLambda<1.100)) || ((lInvMassAntiLambda>1.135)&&(lInvMassAntiLambda<1.145));
 
		// Fill selected V0 particle array 
        if ((cutK0sc)&&(K0Signal))
        {
            selectedV0s->Add(new AliV0ChBasicParticle(aodV0->Eta(), aodV0->Phi(), aodV0->Pt(), 1.));
            Double_t nTrigK0Sig[4] = {aodV0->Pt(), lCent, lPVz, 1.};
            fHistNTrigSib->Fill(nTrigK0Sig);
        }
        if ((cutK0sc)&&(K0Bckg))
        {
            selectedV0s->Add(new AliV0ChBasicParticle(aodV0->Eta(), aodV0->Phi(), aodV0->Pt(), 4.));
            Double_t nTrigK0Bkg[4] = {aodV0->Pt(), lCent, lPVz, 4.};
            fHistNTrigSib->Fill(nTrigK0Bkg);
        }
        if ((cutLambdasc)&&(LamSignal))
        {
            selectedV0s->Add(new AliV0ChBasicParticle(aodV0->Eta(), aodV0->Phi(), aodV0->Pt(), 2.));
            Double_t nTrigLaSig[4] = {aodV0->Pt(), lCent, lPVz, 2.};
            fHistNTrigSib->Fill(nTrigLaSig);
        }
        if ((cutLambdasc)&&(LamBckg))
        {
            selectedV0s->Add(new AliV0ChBasicParticle(aodV0->Eta(), aodV0->Phi(), aodV0->Pt(), 5.));
            Double_t nTrigLaBkg[4] = {aodV0->Pt(), lCent, lPVz, 5.};
            fHistNTrigSib->Fill(nTrigLaBkg);
        }
        if ((cutAntiLambdasc)&&(ALamSignal))
        {
            selectedV0s->Add(new AliV0ChBasicParticle(aodV0->Eta(), aodV0->Phi(), aodV0->Pt(), 3.));
            Double_t nTrigAlSig[4] = {aodV0->Pt(), lCent, lPVz, 3.};
            fHistNTrigSib->Fill(nTrigAlSig);
        }
        if ((cutAntiLambdasc)&&(ALamBckg))
        {
            selectedV0s->Add(new AliV0ChBasicParticle(aodV0->Eta(), aodV0->Phi(), aodV0->Pt(), 6.));
            Double_t nTrigAlBkg[4] = {aodV0->Pt(), lCent, lPVz, 6.};
            fHistNTrigSib->Fill(nTrigAlBkg);
        }


    	//===== MC part for V0 reconstruction efficiency ==============
    	if (fAnalysisMC)
    	{
      		TClonesArray *mcArray = (TClonesArray*)aod->FindListObject(AliAODMCParticle::StdBranchName());
      		if(!mcArray){
        	Printf("No MC particle branch found");
        	return;
      		}

			Int_t MotherOfMotherPdg = 0;

      		Int_t myTrackPosLabel = TMath::Abs(myTrackPos->GetLabel());
      		Int_t myTrackNegLabel = TMath::Abs(myTrackNeg->GetLabel());
      		AliAODMCParticle *mcPosTrack = (AliAODMCParticle*)mcArray->At(myTrackPosLabel);
			if (!mcPosTrack) continue;
      		AliAODMCParticle *mcNegTrack = (AliAODMCParticle*)mcArray->At(myTrackNegLabel);
			if (!mcNegTrack) continue;

      		Int_t PosTrackPdg = mcPosTrack->GetPdgCode();
      		Int_t NegTrackPdg = mcNegTrack->GetPdgCode();
        	//if (!mcPosTrack->IsPrimary()) continue;
        	//if (!mcNegTrack->IsPrimary()) continue;

      		Int_t myTrackPosMotherLabel = mcPosTrack->GetMother();
      		Int_t myTrackNegMotherLabel = mcNegTrack->GetMother();

      		if ((myTrackPosMotherLabel==-1)||(myTrackNegMotherLabel==-1)) continue;

      		AliAODMCParticle *mcPosMother = (AliAODMCParticle*)mcArray->At(myTrackPosMotherLabel);
			if (!mcPosMother) continue;	
      		Int_t MotherPdg = mcPosMother->GetPdgCode();
			Int_t MotherOfMother = mcPosMother->GetMother();
			//if (MotherOfMother == -1) MotherOfMotherPdg = 0;

			if (myTrackPosMotherLabel!=myTrackNegMotherLabel) continue;
        	//if (!mcPosMother->IsPrimary()) continue;
			Bool_t IsK0FromMC = (mcPosMother->IsPhysicalPrimary())&&(MotherPdg==310)&&(PosTrackPdg==211)&&(NegTrackPdg==-211);
			Bool_t IsLambdaFromMC = (mcPosMother->IsPhysicalPrimary())&&(MotherPdg==3122)&&(PosTrackPdg==2212)&&(NegTrackPdg==-211);
			Bool_t IsAntiLambdaFromMC = (mcPosMother->IsPhysicalPrimary())&&(MotherPdg==-3122)&&(PosTrackPdg==211)&&(NegTrackPdg==-2212);

			Bool_t ComeFromSigma = kFALSE;
			Bool_t ComeFromSigmaLa = kFALSE;
			Bool_t ComeFromSigmaAl = kFALSE;
		
		    if (MotherOfMother != -1)
			{
      			AliAODMCParticle *mcPosMotherOfMother = (AliAODMCParticle*)mcArray->At(MotherOfMother);
				MotherOfMotherPdg = mcPosMotherOfMother->GetPdgCode();
            	Int_t MoMPdg = TMath::Abs(MotherOfMotherPdg); 
				ComeFromSigma = (mcPosMotherOfMother->IsPhysicalPrimary())&&((MoMPdg==3212)||(MoMPdg==3224)||(MoMPdg==3214)||(MoMPdg==3114));
				ComeFromSigmaLa = ComeFromSigma && (MotherPdg==3122)&&(PosTrackPdg==2212)&&(NegTrackPdg==-211);
				ComeFromSigmaAl = ComeFromSigma && (MotherPdg==-3122)&&(PosTrackPdg==211)&&(NegTrackPdg==-2212);
			}

            IsLambdaFromMC = IsLambdaFromMC || ComeFromSigmaLa;
            IsAntiLambdaFromMC = IsAntiLambdaFromMC || ComeFromSigmaAl;
			
      		Double_t RecMotherPt = aodV0->Pt();
			//cout << "Pt of rec v0 = " << RecMotherPt << " nMC = " << mcArray->GetEntries() << endl;
			//cout << "Pos Label = " << myTrackPosLabel << " Neg Label = " << myTrackNegLabel << endl;
      		Double_t rcK0[3] = {RecMotherPt, lCent, 1};
      		Double_t rcLa[3] = {RecMotherPt, lCent, 2};
      		Double_t rcAl[3] = {RecMotherPt, lCent, 3};
      		if ((cutK0sc)&&(K0Signal)&&(IsK0FromMC)) fHistRCPtCentTrig->Fill(rcK0);
      		if ((cutLambdasc)&&(LamSignal)&&(IsLambdaFromMC)) fHistRCPtCentTrig->Fill(rcLa);
      		if ((cutAntiLambdasc)&&(ALamSignal)&&(IsAntiLambdaFromMC)) fHistRCPtCentTrig->Fill(rcAl);

    	}

    	//===== End of the MC part for V0 reconstruction efficiency ===

		Int_t nSelectedTracks = selectedTracks->GetEntries();
		// Correlation part
		//===================================
		for (Int_t j = 0; j < nSelectedTracks; j++)
		{
			AliAODTrack* atr = (AliAODTrack*) selectedTracks->At(j);
			if ((atr->Pt())>=(aodV0->Pt())) continue;
			Double_t dEta = atr->Eta() - aodV0->Eta();
			Double_t dPhi = atr->Phi() - aodV0->Phi();
			if (dPhi > (1.5*kPi)) dPhi -= 2.0*kPi;
			if (dPhi < (-0.5*kPi)) dPhi += 2.0*kPi;
			// removing autocorrelations 
			//----------------------------------
			Int_t negID = myTrackNeg->GetID();
			Int_t posID = myTrackPos->GetID();
			Int_t atrID = atr->GetID();

			if ((TMath::Abs(negID)+1)==(TMath::Abs(atrID))) continue;
			if ((TMath::Abs(posID)+1)==(TMath::Abs(atrID))) continue;
			//----------------------------------

			fHistNTrigSib->Sumw2();
			fHistdPhidEtaSib->Sumw2();
			// Filling correlation histograms and histograms for triggers counting
			//----------------- K0 ---------------------
			if ((cutK0sc)&&(K0Signal)) 
			{
				Double_t spK0Sig[7] = {dPhi, dEta, aodV0->Pt(), atr->Pt(), lCent, lPVz, 1.};
				fHistdPhidEtaSib->Fill(spK0Sig);
			}

			if ((cutK0sc)&&(K0Bckg)) 
			{
				Double_t spK0Bkg[7] = {dPhi, dEta, aodV0->Pt(), atr->Pt(), lCent, lPVz, 4.};
				fHistdPhidEtaSib->Fill(spK0Bkg);
			}
			//---------------- Lambda -------------------
			if ((cutLambdasc)&&(LamSignal)) 
			{
				Double_t spLaSig[7] = {dPhi, dEta, aodV0->Pt(), atr->Pt(), lCent, lPVz, 2.};
				fHistdPhidEtaSib->Fill(spLaSig);
			}

			if ((cutLambdasc)&&(LamBckg)) 
			{
				Double_t spLaBkg[7] = {dPhi, dEta, aodV0->Pt(), atr->Pt(), lCent, lPVz, 5.};
				fHistdPhidEtaSib->Fill(spLaBkg);
			}
			//------------- AntiLambda -------------------
			if ((cutAntiLambdasc)&&(ALamSignal)) 
			{
				Double_t spAlSig[7] = {dPhi, dEta, aodV0->Pt(), atr->Pt(), lCent, lPVz, 3.};
				fHistdPhidEtaSib->Fill(spAlSig);
			}

			if ((cutAntiLambdasc)&&(ALamBckg)) 
			{
				Double_t spAlBkg[7] = {dPhi, dEta, aodV0->Pt(), atr->Pt(), lCent, lPVz, 6.};
				fHistdPhidEtaSib->Fill(spAlBkg);
			}
	
		} // end of correlation loop
		//===================================


	} // end of V0 selection loop

	// ===========================================
    // Ch-Ch correlation part

    if (ChChWith)
    {
        Int_t nSelectedChargedTriggers = selectedChargedTriggers->GetEntries();
        for (Int_t i = 0; i < nSelectedChargedTriggers; i++)
        {
            AliV0ChBasicParticle* chTrig = (AliV0ChBasicParticle*) selectedChargedTriggers->At(i);
            Double_t chTrigPhi = chTrig->Phi();
            Double_t chTrigEta = chTrig->Eta();
            Double_t chTrigPt = chTrig->Pt();

            Double_t nTrigChSig[4] = {chTrigPt, lCent, lPVz, 7.};
            fHistNTrigSib->Fill(nTrigChSig);

            Int_t nSelectedTracks = selectedTracks->GetEntries();
            for (Int_t j = 0; j < nSelectedTracks; j++)
            {
                AliAODTrack* atr = (AliAODTrack*) selectedTracks->At(j);
                if ((atr->Pt())>=chTrigPt) continue;
                Double_t dEta = atr->Eta() - chTrigEta;
                Double_t dPhi = atr->Phi() - chTrigPhi;
                if (dPhi > (1.5*kPi)) dPhi -= 2.0*kPi;
                if (dPhi < (-0.5*kPi)) dPhi += 2.0*kPi;

                Double_t spCh[7] = {dPhi, dEta, chTrigPt, atr->Pt(), lCent, lPVz, 7.};
                fHistdPhidEtaSib->Fill(spCh);
            }
        }
    }
    // end of Ch-Ch correlation part

	// Mixing ==============================================
	
	fHistdPhidEtaMix->Sumw2();
    AliEventPool* pool = fPoolMgr->GetEventPool(lCent, lPVz);
	if (!pool) AliFatal(Form("No pool found for centrality = %f, zVtx = %f", lCent, lPVz));
	//pool->SetDebug(1);
    if (pool->IsReady() || pool->NTracksInPool() > fMixingTracks / 10 || pool->GetCurrentNEvents() >= 5)
	{
		Int_t nMix = pool->GetCurrentNEvents();
		for (Int_t jMix=0; jMix<nMix; jMix++)
		{// loop through mixing events
			TObjArray* bgTracks = pool->GetEvent(jMix);
			for (Int_t i=0; i<selectedV0s->GetEntriesFast(); i++)
			{// loop through selected V0 particles
				AliV0ChBasicParticle* trig = (AliV0ChBasicParticle*) selectedV0s->At(i);
				Double_t trigPhi = trig->Phi();
				Double_t trigEta = trig->Eta();
				Double_t trigPt = trig->Pt();
				Short_t trigC = trig->WhichCandidate();
				for (Int_t j=0; j<bgTracks->GetEntriesFast(); j++)
				{ // mixing tracks loop 
					AliVParticle* assoc = (AliVParticle*) bgTracks->At(j);
					// be careful tracks may have bigger pt than v0s.
					if ( ( (assoc->Pt())>=trigPt )||( (assoc->Pt())<PtAssocMin ) ) continue;
					Double_t dEtaMix = assoc->Eta() - trigEta;
					Double_t dPhiMix = assoc->Phi() - trigPhi;
					if (dPhiMix > (1.5*kPi)) dPhiMix -= 2.0*kPi;
					if (dPhiMix < (-0.5*kPi)) dPhiMix += 2.0*kPi;
				
					Double_t spMix[7] = {dPhiMix, dEtaMix, trigPt, assoc->Pt(), lCent, lPVz, static_cast<Double_t>(trigC)};
				    fHistdPhidEtaMix->Fill(spMix);
					
				} // end of mixing track loop
			}// end of loop through selected V0 particles
			if (ChChWith)
            {
                for (Int_t i=0; i<selectedChargedTriggers->GetEntriesFast(); i++)
                {// loop through selected charged trigger particles
                    AliV0ChBasicParticle* trig = (AliV0ChBasicParticle*) selectedChargedTriggers->At(i);
                    Double_t trigPhi = trig->Phi();
                    Double_t trigEta = trig->Eta();
                    Double_t trigPt = trig->Pt();
                    Short_t trigC = trig->WhichCandidate();
                    for (Int_t j=0; j<bgTracks->GetEntriesFast(); j++)
                    { // mixing tracks loop 
                        AliVParticle* assoc = (AliVParticle*) bgTracks->At(j);
                        // be careful tracks may have bigger pt than v0s.
                        if ( ( (assoc->Pt())>=trigPt )||( (assoc->Pt())<PtAssocMin ) ) continue;

                        Double_t dEtaMix = assoc->Eta() - trigEta;
                        Double_t dPhiMix = assoc->Phi() - trigPhi;
                        if (dPhiMix > (1.5*kPi)) dPhiMix -= 2.0*kPi;
                        if (dPhiMix < (-0.5*kPi)) dPhiMix += 2.0*kPi;

                        Double_t spMix[7] = {dPhiMix, dEtaMix, trigPt, assoc->Pt(), lCent, lPVz, static_cast<Double_t>(trigC)};
                        fHistdPhidEtaMix->Fill(spMix);

                    } // end of mixing track loop
                }// end of loop through selected Ch particles
            } // end of ChCh mixing

		}// end of loop of mixing events
	}

	TObjArray* tracksClone = (TObjArray*) selectedTracks->Clone();
	tracksClone->SetOwner(kTRUE);
	pool->UpdatePool(tracksClone);
	//delete selectedtracks;

    PostData(1, fOutput);

}
//___________________________________________
Bool_t AliAnalysisTaskV0ChCorrelations::IsMyGoodPrimaryTrack(const AliAODTrack *t)
{
	// Pseudorapidity cut
	if (TMath::Abs(t->Eta())>0.9) return kFALSE;
	// Should correspond to set of cuts suitable for correlation analysis, 768 - hybrid tracks for 2011 data
    if (!t->TestFilterBit(768)) return kFALSE;

	return kTRUE;
}
//_____________________________________________
Bool_t AliAnalysisTaskV0ChCorrelations::IsMyGoodDaughterTrack(const AliAODTrack *t)
{
	// TPC refit
	if (!t->IsOn(AliAODTrack::kTPCrefit)) return kFALSE;
	// Minimum number of clusters
	Float_t nCrossedRowsTPC = t->GetTPCClusterInfo(2,1);
	if (nCrossedRowsTPC < 70) return kFALSE;
	Int_t findable=t->GetTPCNclsF();
	if (findable <= 0) return kFALSE;
	if (nCrossedRowsTPC/findable < 0.8) return kFALSE;

	return kTRUE;
}
//______________________________________________
Bool_t AliAnalysisTaskV0ChCorrelations::IsMyGoodV0(const AliAODEvent* aod, const AliAODv0* aodV0, const AliAODTrack* myTrackPos, const AliAODTrack* myTrackNeg, Int_t oSta)
{
	if (!aodV0) {
       AliError(Form("ERROR: Could not retrieve aodV0"));
       return kFALSE;
    }

	// Offline reconstructed V0 only
    if (oSta==1) {if (aodV0->GetOnFlyStatus()) return kFALSE;}
    if (oSta==3) {if (!aodV0->GetOnFlyStatus()) return kFALSE;}
   
	if (oSta==2) 
	{
		if (aodV0->GetOnFlyStatus()) 
		{
			return kTRUE;
		} else {
			return kFALSE;
		}
	}
    if (oSta==4) 
	{
		if (!aodV0->GetOnFlyStatus()) 
		{
			return kTRUE;
		} else {
			return kFALSE;
		}
	}
    
    // Get daughters and check them
	AliAODTrack *myTrackNegTest=dynamic_cast<AliAODTrack *>(aodV0->GetDaughter(1));
	AliAODTrack *myTrackPosTest=dynamic_cast<AliAODTrack *>(aodV0->GetDaughter(0));
	
	if (!myTrackPosTest || !myTrackNegTest) {
		Printf("strange analysis::UserExec:: Error:Could not retreive one of the daughter track\n");
		return kFALSE;
	}

    if( myTrackPosTest->Charge() ==1){
            myTrackPos = myTrackPosTest;
            myTrackNeg = myTrackNegTest;
    }

    if( myTrackPosTest->Charge() ==-1){
            myTrackPos = myTrackNegTest;
            myTrackNeg = myTrackPosTest;
    }

	// Track cuts for daughter tracks
    if ( !(IsMyGoodDaughterTrack(myTrackPos)) || !(IsMyGoodDaughterTrack(myTrackNeg)) ) return kFALSE;

	// Unlike signs of daughters
    if (myTrackNegTest->Charge() == myTrackPosTest->Charge()) return kFALSE;

	// Rapidity cut
	//Double_t lCutRap = 0.75;
	//Double_t lRapK0s = aodV0->Y(310);
	//Double_t lRapLambda = aodV0->Y(3122);
	//Double_t lRapAntiLambda = aodV0->Y(-3122);

	// Pseudorapidity cut - there are high pt V0
	Double_t lCutEta = 0.75;
	Double_t lEtaV0 = aodV0->Eta();
	if (TMath::Abs(lEtaV0)>=lCutEta) return kFALSE;
	//if (TMath::Abs(lRapK0s)>=lCutRap) return kFALSE;
	//if (TMath::Abs(lRapLambda)>=lCutRap) return kFALSE;
	//if (TMath::Abs(lRapAntiLambda)>=lCutRap) return kFALSE;
    
    // getting global variables
	Float_t dcaDaughtersToPrimVtx = GetDcaDToPV();
	Float_t dcaBetweenDaughters = GetDcaV0D();

    // DCA of daughter track to Primary Vertex
    Float_t xyn=aodV0->DcaNegToPrimVertex();
    if (TMath::Abs(xyn)<dcaDaughtersToPrimVtx) return kFALSE;
    Float_t xyp=aodV0->DcaPosToPrimVertex();
    if (TMath::Abs(xyp)<dcaDaughtersToPrimVtx) return kFALSE;

	// DCA of daughter tracks 
    Double_t dca=aodV0->DcaV0Daughters();
    if (dca>dcaBetweenDaughters) return kFALSE;

	// Cosinus of pointing angle
    Double_t cpa=aodV0->CosPointingAngle(aod->GetPrimaryVertex());
    if (cpa<0.998) return kFALSE;

	// Fiducial volume cut
    Double_t xyz[3]; aodV0->GetSecondaryVtx(xyz);
    Double_t r2=xyz[0]*xyz[0] + xyz[1]*xyz[1];
    if (r2<0.9*0.9) return kFALSE;
    if (r2>100*100) return kFALSE;

	// c*tau cut - in main V0 loop - depends on particle hypothesis

	// Minimum pt of daughters
    Double_t  lMomPos[3] = {999,999,999};
	Double_t  lMomNeg[3] = {999,999,999};

	lMomPos[0] = aodV0->MomPosX();
    lMomPos[1] = aodV0->MomPosY();
    lMomPos[2] = aodV0->MomPosZ();

    lMomNeg[0] = aodV0->MomNegX();
    lMomNeg[1] = aodV0->MomNegY();
    lMomNeg[2] = aodV0->MomNegZ();

    Double_t lPtPos = TMath::Sqrt(lMomPos[0]*lMomPos[0] + lMomPos[1]*lMomPos[1]);
    Double_t lPtNeg = TMath::Sqrt(lMomNeg[0]*lMomNeg[0] + lMomNeg[1]*lMomNeg[1]);
	
	Double_t cutMinPtDaughter = 0.160;
	if (lPtPos<cutMinPtDaughter || lPtNeg<cutMinPtDaughter) return kFALSE;

	// Daughter PID cut - in main V0 loop - depends on particle hypothesis

	return kTRUE;
}

void AliAnalysisTaskV0ChCorrelations::RemovingInjectedSignal(TObjArray* tracks, TObject* mcObj, Int_t maxLabel)
{
  // remove injected signals (primaries above <maxLabel>)
  // <tracks> can be the following cases:
  // a. tracks: in this case the label is taken and then case b.
  // b. particles: the first stable mother is searched and checked if it is <= <maxLabel>
  // <mcObj> can be AOD (TClonesArray) or ESD (AliMCEvent)

  TClonesArray* arrayMC = 0;
  AliMCEvent* mcEvent = 0;
  if (mcObj->InheritsFrom("AliMCEvent"))
    mcEvent = static_cast<AliMCEvent*>(mcObj);
  else if (mcObj->InheritsFrom("TClonesArray"))
    arrayMC = static_cast<TClonesArray*>(mcObj);
  else
  {
    mcObj->Dump();
    AliFatal("Invalid object passed");
  }

  Int_t before = tracks->GetEntriesFast();

  for (Int_t i=0; i<before; ++i)
  {
    AliVParticle* part = (AliVParticle*) tracks->At(i);

    if (part->InheritsFrom("AliESDtrack")) part = mcEvent->GetTrack(TMath::Abs(part->GetLabel()));
    if (part->InheritsFrom("AliAODTrack"))
    {
      part =(AliVParticle*)arrayMC->At(TMath::Abs(part->GetLabel()));
      //cout << "toto musi byt len pri reco trackoch" << endl;
    }

    AliVParticle* mother = part;
    if (mcEvent)
    {
      while (!mcEvent->IsPhysicalPrimary(mother->GetLabel()))
      {
    if (((AliMCParticle*)mother)->GetMother() < 0)
    {
      mother = 0;
      break;
    }

    mother = (AliMCParticle*) mcEvent->GetTrack(((AliMCParticle*)mother)->GetMother());
    if (!mother)
      break;
      }
    }
    else
    {
      // find the primary mother
     while (!((AliAODMCParticle*)mother)->IsPhysicalPrimary())
      {
    if (((AliAODMCParticle*)mother)->GetMother() < 0)
    {
      mother = 0;
      break;
    }

    mother = (AliVParticle*) arrayMC->At(((AliAODMCParticle*)mother)->GetMother());
    if (!mother)
      break;
      }
    }

    if (!mother)
    {
      //Printf("WARNING: No mother found for particle %d:", part->GetLabel());
      continue;
    }

    if (mother->GetLabel() > maxLabel)
    {
//      Printf("Removing %d with label %d", i, part->GetLabel()); part->Dump();
      TObject* object = tracks->RemoveAt(i);
      if (tracks->IsOwner())
    delete object;
    }
  }

  tracks->Compress();

  AliInfo(Form("Reduced from %d to %d", before, tracks->GetEntriesFast()));

}

