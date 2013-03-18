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
 * The task works with AOD events only and containes also mixing for acceptance corrections.
 * Last update edited by Marek Bombara, January 2013, Marek.Bombara@cern.ch
 */

#include <TCanvas.h>
#include <TList.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <THnSparse.h>

#include "AliLog.h"
#include "AliAnalysisManager.h"
#include "AliAODTrack.h"
#include "AliAODEvent.h"
#include "AliAODv0.h"
#include "AliAODcascade.h"
#include "AliAODVertex.h"

#include "AliMCEvent.h"
#include "AliMCVertex.h" 
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"

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
   fAnalysisMC(kTRUE),
	 fFillMixed(kTRUE),
	 fMixingTracks(500),
	 fPoolMgr(0x0),
     
     fDcaDToPV(0.5),
	 fDcaV0D(0.1),

	 fOutput(0),
	 fPIDResponse(0),
     fHistCentVtx(0),
     fHistMultiMain(0),
	 fHistMassK0(0),
	 fHistMassLambda(0),
	 fHistMassAntiLambda(0),
 
	 fHistdPhidEtaSib(0),
	 fHistdPhidEtaMix(0),
	 fHistTrigSib(0),
	 fHistTrigMix(0),

   fHistMCPtCent(0),
   fHistRCPtCent(0),
	 
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
   Int_t nZvtxBins  = 20;
   Double_t vertexBins[] = {-10.,-9.,-8.,-7.,-6.,-5.,-4.,-3.,-2.,-1.,0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.};
   const Double_t* zvtxBins = vertexBins;
   // pt bins of associated particles for the analysis
   Int_t nPtBins = 6;
   const Double_t PtBins[7] = {3.0,4.0,5.0,6.0,7.0,8.0,9.0}; 
   //Int_t nPtBins = 1;
   //const Double_t PtBins[2] = {3.0,15.0}; 
   // pt bins of trigger particles for the analysis
   Int_t nPtBinsV0 = 9;
   //const Double_t PtBinsV0[2] = {6.0,15.0}; 
   const Double_t PtBinsV0[10] = {6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0}; 
   // V0 candidate: 1 - K0sig, 2 - Lamsig, 3 - Alamsig, 4 - K0bg, 5 - Lambg, 6 - Alambg
   Int_t nTrigC = 6;
   const Double_t TrigC[7] = {0.5,1.5,2.5,3.5,4.5,5.5,6.5};
   
   // Create general histograms 
   fHistCentVtx = new TH2F("fHistCentVtx", "Centrality vs. Z vertex", nCentralityBins, centralityBins, nZvtxBins, zvtxBins);
   fHistMultiMain = new TH1F("fHistMultiMain", "Multiplicity of main events", 2000, 0, 2000);

   // Create histograms for reconstruction track efficiency
   fHistMCPtCent = new TH2D("fHistMCPtCent", "MC Pt vs. Cent.", nPtBins, PtBins[0], PtBinsV0[9], nCentralityBins, centBins[0], centBins[9]);
   fHistRCPtCent = new TH2D("fHistRCPtCent", "Rec Pt vs. Cent.", nPtBins, PtBins[0], PtBinsV0[9], nCentralityBins, centBins[0], centBins[9]);

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
   const Double_t spMaxK0[3] = {mk0Max, PtBinsV0[9], centralityBins[9]};
   const Double_t spMinLa[3] = {mlaMin, PtBinsV0[0], centralityBins[0]};
   const Double_t spMaxLa[3] = {mlaMax, PtBinsV0[9], centralityBins[9]};
   const Double_t spMinAl[3] = {malMin, PtBinsV0[0], centralityBins[0]};
   const Double_t spMaxAl[3] = {malMax, PtBinsV0[9], centralityBins[9]};
   // Create mass histograms
   fHistMassK0 = new THnSparseF("fHistMassK0","V0 mass for K0 hypothesis", 3, spBins, spMinK0, spMaxK0);
   fHistMassLambda = new THnSparseF("fHistMassLambda","V0 mass for Lambda hypothesis", 3, spBins, spMinLa, spMaxLa);
   fHistMassAntiLambda = new THnSparseF("fHistMassAntiLambda","V0 mass for AntiLambda hypothesis", 3, spBins, spMinAl, spMaxAl);
   
   // defining bins for dPhi distributions
   const Int_t nbPhiBins = 144;
   const Double_t kPi = TMath::Pi();
   Double_t PhiMin = -kPi/2.;
   Double_t PhiMax = -kPi/2. + 2*kPi;
   Double_t PhiBins[nbPhiBins+1] = {0.};
   PhiBins[0] = PhiMin;
   for (Int_t i=0; i<nbPhiBins; i++) { PhiBins[i+1] = PhiBins[i] + (PhiMax - PhiMin)/nbPhiBins; }

   // defining bins for dEta distributions
   const Int_t nbEtaBins = 72;
   Double_t EtaMin = -1.8;
   Double_t EtaMax = 1.8;
   Double_t EtaBins[nbEtaBins+1] = {0.};
   EtaBins[0] = EtaMin;
   for (Int_t i=0; i<nbEtaBins; i++) { EtaBins[i+1] = EtaBins[i] + (EtaMax - EtaMin)/nbEtaBins; }

   const Int_t corBins[7] = {nbPhiBins, nbEtaBins, nPtBinsV0, nPtBins, nCentralityBins, nZvtxBins, nTrigC};
   const Double_t corMin[7] = {PhiBins[0], EtaBins[0], PtBinsV0[0], PtBins[0], centralityBins[0], zvtxBins[0], TrigC[0]};
   const Double_t corMax[7] = {PhiBins[144], EtaBins[72], PtBinsV0[9], PtBins[6], centralityBins[9], zvtxBins[20], TrigC[6]};
   // Create correlation histograms
   fHistdPhidEtaSib = new THnSparseF("fHistdPhidEtaSib","dPhi vs. dEta siblings", 7, corBins, corMin, corMax); 
   fHistdPhidEtaMix = new THnSparseF("fHistdPhidEtaMix","dPhi vs. dEta mixed", 7, corBins, corMin, corMax);

   const Int_t corTrigBins[5] = {nPtBinsV0, nPtBins, nCentralityBins, nZvtxBins, nTrigC};
   const Double_t corTrigMin[5] = {PtBinsV0[0], PtBins[0], centBins[0], vertexBins[0], TrigC[0]};
   const Double_t corTrigMax[5] = {PtBinsV0[9], PtBins[6], centBins[9], vertexBins[20], TrigC[6]};
   // Create histograms for trigger particles
   fHistTrigSib = new THnSparseF("fHistTrigSib","pt trigger sib", 5, corTrigBins, corTrigMin, corTrigMax); 
   fHistTrigMix = new THnSparseF("fHistTrigMix","pt trigger mix", 5, corTrigBins, corTrigMin, corTrigMax); 
   
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
   fOutput->Add(fHistTrigSib);
   fOutput->Add(fHistTrigMix);

   fOutput->Add(fHistMCPtCent);
   fOutput->Add(fHistRCPtCent);
   
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

    AliAODEvent* aod = (AliAODEvent*)inEvMain->GetEvent();
    fPIDResponse = inEvMain->GetPIDResponse();

  // pt intervals for trigger particles  
  const Double_t kPi = TMath::Pi();
  Double_t PtTrigMin = 6.0;
  Double_t PtTrigMax = 15.0;
  // pt interval for associated particles
  Double_t PtAssocMin = 3.0;
  fHistMultiMain->Fill(aod->GetNumberOfTracks());

  // Vertex cut
  Double_t cutPrimVertex = 10.0;
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
  centralityObj = aod->GetHeader()->GetCentralityP();
  lCent = centralityObj->GetCentralityPercentile("CL1");

  if ((lCent < 0.)||(lCent > 90.)) return;
  fHistCentVtx->Fill(lCent,lPVz);

//=========== MC loop ===============================
  if (fAnalysisMC)
  {
    AliAODMCHeader *aodMCheader = (AliAODMCHeader*)aod->FindListObject(AliAODMCHeader::StdBranchName());
    Float_t vzMC = aodMCheader->GetVtxZ();
    if (TMath::Abs(vzMC) >= 10.) return;
    //retreive MC particles from event
    TClonesArray *mcArray = (TClonesArray*)aod->FindListObject(AliAODMCParticle::StdBranchName());
    if(!mcArray){
        Printf("No MC particle branch found");
        return;
    }
   
    Int_t nMCTracks = mcArray->GetEntriesFast();
    //cout << "number of MC tracks = " << nMCTracks << endl;
    for (Int_t iMC = 0; iMC<nMCTracks; iMC++)
    {
      AliAODMCParticle *mcTrack = (AliAODMCParticle*)mcArray->At(iMC);
        if (!mcTrack) {
            Error("ReadEventAODMC", "Could not receive particle %d", iMC);
            continue;
      }
      
      Double_t mcTrackEta = mcTrack->Eta();
      Double_t mcTrackPt = mcTrack->Pt();
      if (!(mcTrack->IsPhysicalPrimary())) continue;
      if (TMath::Abs(mcTrackEta)>0.8) continue;
      if (mcTrackPt<PtAssocMin) continue;
      if ((mcTrack->Charge())==0) continue;
      fHistMCPtCent->Fill(mcTrackPt,lCent);
    }
    // ------- access the real data -----------
    Int_t nRecTracks = aod->GetNumberOfTracks();
    for (Int_t i = 0; i < nRecTracks; i++)
    {
      AliAODTrack* tras = aod->GetTrack(i);
      if ((tras->Pt())<PtAssocMin) continue; 
      if (!(IsMyGoodPrimaryTrack(tras))) continue;
      Int_t AssocLabel = tras->GetLabel();
      if (AssocLabel<=0) continue;
      if(!(static_cast<AliAODMCParticle*>(mcArray->At(tras->GetLabel()))->IsPhysicalPrimary())) continue;
      Double_t trPt = tras->Pt();
      fHistRCPtCent->Fill(trPt,lCent);
    }
    // ------- end of real data access -------- 
  }
//============= End of MC loop ======================
	
	// Track selection loop
	//--------------------------------
	Int_t nTracks = aod->GetNumberOfTracks();
	// new tracks array
	TObjArray * selectedTracks = new TObjArray;
	selectedTracks->SetOwner(kTRUE);
	for (Int_t i = 0; i < nTracks; i++)
	{
		AliAODTrack* tr = aod->GetTrack(i);
		if ((tr->Pt())<PtAssocMin) continue; 
		if (!(IsMyGoodPrimaryTrack(tr))) continue;
		selectedTracks->Add(tr);
		//if (tr->Pt()>6.) fHistTemp->Fill(tr->Pt());
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

//        if (!IsMyGoodV0CutSet0(aod,aodV0,myTrackPos,myTrackNeg)) continue;

        // effective mass calculations for each hypothesis
		Double_t lInvMassK0 = aodV0->MassK0Short();
        Double_t lInvMassAntiLambda = aodV0->MassAntiLambda();
        Double_t lInvMassLambda = aodV0->MassLambda();

		// calculations for c*tau cut--------------------------------------
		Double_t lDVx = aodV0->GetSecVtxX();
        Double_t lDVy = aodV0->GetSecVtxY();
        Double_t lDVz = aodV0->GetSecVtxZ();
        const Double_t kLambdaMass = 1.115683;
        const Double_t kK0Mass = 0.497648;
        Double_t cutcTauLam = 3*7.89;
        Double_t cutcTauK0 = 3*2.68;
        Double_t lV0DecayLength = TMath::Sqrt(TMath::Power(lDVx - lPVx,2) + TMath::Power(lDVy- lPVy,2) + TMath::Power(lDVz - lPVz,2 ));
        Double_t lPV0 = TMath::Sqrt((aodV0->Pt())*(aodV0->Pt())+(aodV0->Pz())*(aodV0->Pz()));
        Double_t lcTauLam = (lV0DecayLength*kLambdaMass)/lPV0;
        Double_t lcTauK0 = (lV0DecayLength*kK0Mass)/lPV0;
        // sc - standard cuts
		Bool_t cutK0sc = (lcTauK0<cutcTauK0);
        Bool_t cutLambdasc = (lcTauLam<cutcTauLam);
        Bool_t cutAntiLambdasc = (lcTauLam<cutcTauLam);

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
		Bool_t k0APcut = (aodV0->Pt()>(TMath::Abs(0.2*aodV0->AlphaV0())));
		cutK0sc = cutK0sc && k0APcut;
        // fill the mass histograms
	
        if (!IsMyGoodV0(aod,aodV0,myTrackPos,myTrackNeg)) continue;
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
		if ((cutK0sc)&&(K0Signal)) selectedV0s->Add(new AliV0ChBasicParticle(aodV0->Eta(), aodV0->Phi(), aodV0->Pt(), 1));
		if ((cutK0sc)&&(K0Bckg)) selectedV0s->Add(new AliV0ChBasicParticle(aodV0->Eta(), aodV0->Phi(), aodV0->Pt(), 4));
		if ((cutLambdasc)&&(LamSignal)) selectedV0s->Add(new AliV0ChBasicParticle(aodV0->Eta(), aodV0->Phi(), aodV0->Pt(), 2));
		if ((cutLambdasc)&&(LamBckg)) selectedV0s->Add(new AliV0ChBasicParticle(aodV0->Eta(), aodV0->Phi(), aodV0->Pt(), 5));
		if ((cutAntiLambdasc)&&(ALamSignal)) selectedV0s->Add(new AliV0ChBasicParticle(aodV0->Eta(), aodV0->Phi(), aodV0->Pt(), 3));
		if ((cutAntiLambdasc)&&(ALamBckg)) selectedV0s->Add(new AliV0ChBasicParticle(aodV0->Eta(), aodV0->Phi(), aodV0->Pt(), 6));

		// preparation for calculation of the number of triggers for each pt_assoc bin
		// in the future probably better to make the number of bins global
        Int_t NumberOfTriggersK0Sig[6] = {0.};
        Int_t NumberOfTriggersK0Bkg[6] = {0.};
        Int_t NumberOfTriggersLaSig[6] = {0.};
        Int_t NumberOfTriggersLaBkg[6] = {0.};
        Int_t NumberOfTriggersAlSig[6] = {0.};
        Int_t NumberOfTriggersAlBkg[6] = {0.};

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

			fHistTrigSib->Sumw2();
			fHistdPhidEtaSib->Sumw2();
			// Filling correlation histograms and histograms for triggers counting
			//----------------- K0 ---------------------
			if ((cutK0sc)&&(K0Signal)) 
			{
				for (Int_t k=0; k<6; k++)
				{
					if ((NumberOfTriggersK0Sig[k]==0)&&(atr->Pt()>(k+3)))
					{
						Double_t spK0TrigSig[5] = {aodV0->Pt(), k+3.5, lCent, lPVz, 1.};
						fHistTrigSib->Fill(spK0TrigSig);
						NumberOfTriggersK0Sig[k]=1;
					}
				}
				Double_t spK0Sig[7] = {dPhi, dEta, aodV0->Pt(), atr->Pt(), lCent, lPVz, 1.};
				fHistdPhidEtaSib->Fill(spK0Sig);
			}

			if ((cutK0sc)&&(K0Bckg)) 
			{
				for (Int_t k=0; k<6; k++)
				{
					if ((NumberOfTriggersK0Bkg[k]==0)&&(atr->Pt()>(k+3)))
					{
						Double_t spK0TrigBkg[5] = {aodV0->Pt(), k+3.5, lCent, lPVz, 4.};
						fHistTrigSib->Fill(spK0TrigBkg);
						NumberOfTriggersK0Bkg[k]=1;
					}
				}
				Double_t spK0Bkg[7] = {dPhi, dEta, aodV0->Pt(), atr->Pt(), lCent, lPVz, 4.};
				fHistdPhidEtaSib->Fill(spK0Bkg);
			}
			//---------------- Lambda -------------------
			if ((cutLambdasc)&&(LamSignal)) 
			{
				for (Int_t k=0; k<6; k++)
				{
					if ((NumberOfTriggersLaSig[k]==0)&&(atr->Pt()>(k+3)))
					{
						Double_t spLaTrigSig[5] = {aodV0->Pt(), k+3.5, lCent, lPVz, 2.};
						fHistTrigSib->Fill(spLaTrigSig);
						NumberOfTriggersLaSig[k]=1;
					}
				}
				Double_t spLaSig[7] = {dPhi, dEta, aodV0->Pt(), atr->Pt(), lCent, lPVz, 2.};
				fHistdPhidEtaSib->Fill(spLaSig);
			}

			if ((cutLambdasc)&&(LamBckg)) 
			{
				for (Int_t k=0; k<6; k++)
				{
					if ((NumberOfTriggersLaBkg[k]==0)&&(atr->Pt()>(k+3)))
					{
						Double_t spLaTrigBkg[5] = {aodV0->Pt(), k+3.5, lCent, lPVz, 5.};
						fHistTrigSib->Fill(spLaTrigBkg);
						NumberOfTriggersLaBkg[k]=1;
					}
				}
				Double_t spLaBkg[7] = {dPhi, dEta, aodV0->Pt(), atr->Pt(), lCent, lPVz, 5.};
				fHistdPhidEtaSib->Fill(spLaBkg);
			}
			//------------- AntiLambda -------------------
			if ((cutAntiLambdasc)&&(ALamSignal)) 
			{
				for (Int_t k=0; k<6; k++)
				{
					if ((NumberOfTriggersAlSig[k]==0)&&(atr->Pt()>(k+3)))
					{
						Double_t spAlTrigSig[5] = {aodV0->Pt(), k+3.5, lCent, lPVz, 3.};
						fHistTrigSib->Fill(spAlTrigSig);
						NumberOfTriggersAlSig[k]=1;
					}
				}
				Double_t spAlSig[7] = {dPhi, dEta, aodV0->Pt(), atr->Pt(), lCent, lPVz, 3.};
				fHistdPhidEtaSib->Fill(spAlSig);
			}

			if ((cutAntiLambdasc)&&(ALamBckg)) 
			{
				for (Int_t k=0; k<6; k++)
				{
					if ((NumberOfTriggersAlBkg[k]==0)&&(atr->Pt()>(k+3)))
					{
						Double_t spAlTrigBkg[5] = {aodV0->Pt(), k+3.5, lCent, lPVz, 6.};
						fHistTrigSib->Fill(spAlTrigBkg);
						NumberOfTriggersAlBkg[k]=1;
					}
				}
				Double_t spAlBkg[7] = {dPhi, dEta, aodV0->Pt(), atr->Pt(), lCent, lPVz, 6.};
				fHistdPhidEtaSib->Fill(spAlBkg);
			}
	
		} // end of correlation loop
		//===================================


	} // end of V0 selection loop

	// Mixing ==============================================
	
	fHistTrigMix->Sumw2();
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
				
				    if (j==0) 
				    { 
					    Double_t spTrigMix[4] = {trigPt, lCent, lPVz, trigC};
					    fHistTrigMix->Fill(spTrigMix);	
				    }
					Double_t spMix[7] = {dPhiMix, dEtaMix, trigPt, assoc->Pt(), lCent, lPVz, trigC};
				    fHistdPhidEtaMix->Fill(spMix);
					
				} // end of mixing track loop
			}// end of loop through selected V0 particles
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
	if (TMath::Abs(t->Eta())>0.8) return kFALSE;
	// Should correspond to set of cuts suitable for correlation analysis
    if (!t->TestFilterBit(1<<7)) return kFALSE;

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
Bool_t AliAnalysisTaskV0ChCorrelations::IsMyGoodV0(const AliAODEvent* aod, const AliAODv0* aodV0, const AliAODTrack* myTrackPos, const AliAODTrack* myTrackNeg)
{
	if (!aodV0) {
       AliError(Form("ERROR: Could not retrieve aodV0"));
       return kFALSE;
    }

    // Rapidity cut
	Double_t lCutRap = 0.75;
	Double_t lRapK0s = aodV0->Y(310);
	Double_t lRapLambda = aodV0->Y(3122);
	Double_t lRapAntiLambda = aodV0->Y(-3122);

	if (TMath::Abs(lRapK0s)>=lCutRap) return kFALSE;
	if (TMath::Abs(lRapLambda)>=lCutRap) return kFALSE;
	if (TMath::Abs(lRapAntiLambda)>=lCutRap) return kFALSE;
    
	// Offline reconstructed V0 only
    if (aodV0->GetOnFlyStatus()) return kFALSE;

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

	// Track cuts for daugher tracks
    if ( !(IsMyGoodDaughterTrack(myTrackPos)) || !(IsMyGoodDaughterTrack(myTrackNeg)) ) return kFALSE;

	// Unlike signs of daughters
    if (myTrackNegTest->Charge() == myTrackPosTest->Charge()) return kFALSE;

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

