// AliAnalysisTaskMultPbTracks

// Author: Michele Floris, CERN
// TODO:
// - Add chi2/cluster plot for primary, secondaries and fakes

#include "AliAnalysisTaskMultPbTracks.h"
#include "AliESDInputHandler.h"
#include "AliAnalysisMultPbTrackHistoManager.h"
#include "AliAnalysisManager.h"
#include "AliESDtrackCuts.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "TH1I.h"
#include "TH3D.h"
#include "TH2D.h"
#include "AliMCParticle.h"
#include "AliGenEventHeader.h"
#include "AliCentrality.h"
#include "AliMultiplicity.h"
#include "TFile.h"
#include <iostream>
#include "AliAnalysisMultPbCentralitySelector.h"
#include "AliTriggerAnalysis.h"
#include "AliPIDResponse.h"

using namespace std;

ClassImp(AliAnalysisTaskMultPbTracks)

AliAnalysisTaskMultPbTracks::AliAnalysisTaskMultPbTracks()
: AliAnalysisTaskSE("TaskMultPbTracks"),
  fESD(0),fHistoManager(0),fCentrSelector(0),fTrackCuts(0),fTrackCutsNoDCA(0),fOfflineTrigger(0), fIsMC(0),fIsTPCOnly(0), fTriggerAnalysis(0),fPIDResponse(0)
{
  // constructor

  DefineOutput(1, AliAnalysisMultPbTrackHistoManager::Class());
  DefineOutput(2, AliESDtrackCuts::Class());
  DefineOutput(3, AliAnalysisMultPbCentralitySelector::Class());

  fHistoManager = new AliAnalysisMultPbTrackHistoManager("histoManager","Hitograms, Multiplicity, Track analysis");
  if(fIsMC) fHistoManager->SetSuffix("MC");

}
AliAnalysisTaskMultPbTracks::AliAnalysisTaskMultPbTracks(const char * name)
  : AliAnalysisTaskSE(name),
    fESD(0),fHistoManager(0),fCentrSelector(0),fTrackCuts(0),fTrackCutsNoDCA(0),fOfflineTrigger(0),fIsMC(0),fIsTPCOnly(0), fTriggerAnalysis(0),fPIDResponse(0)
{
  //
  // Standard constructur which should be used
  //

  DefineOutput(1, AliAnalysisMultPbTrackHistoManager::Class());
  DefineOutput(2, AliESDtrackCuts::Class());
  DefineOutput(3, AliAnalysisMultPbCentralitySelector::Class());

  fHistoManager = new AliAnalysisMultPbTrackHistoManager("histoManager","Hitograms, Multiplicity, Track analysis");
  if(fIsMC) fHistoManager->SetSuffix("MC");

}

AliAnalysisTaskMultPbTracks::AliAnalysisTaskMultPbTracks(const AliAnalysisTaskMultPbTracks& obj) : 
  AliAnalysisTaskSE(obj) ,fESD (0), fHistoManager(0), fCentrSelector(0), fTrackCuts(0),fTrackCutsNoDCA(0),fOfflineTrigger(0),fIsMC(0),fIsTPCOnly(0), fTriggerAnalysis(0),fPIDResponse(0)
{
  //copy ctor
  fESD = obj.fESD ;
  fHistoManager= obj.fHistoManager;
  fTrackCuts  = obj.fTrackCuts;
  fTrackCutsNoDCA  = obj.fTrackCutsNoDCA;
  fCentrSelector = obj.fCentrSelector;
  fOfflineTrigger = obj.fOfflineTrigger;
  fIsMC = obj.fIsMC;
  fIsTPCOnly = obj.fIsTPCOnly;
  fTriggerAnalysis = obj.fTriggerAnalysis;

}

AliAnalysisTaskMultPbTracks::~AliAnalysisTaskMultPbTracks(){
  // destructor

  if(!AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
    if(fHistoManager) {
      delete fHistoManager;
      fHistoManager = 0;
      delete fTriggerAnalysis;
      fTriggerAnalysis=0;
    }
  }
  // Histo list should not be destroyed: fListWrapper is owner!

}
void AliAnalysisTaskMultPbTracks::UserCreateOutputObjects()
{
  // Called once

  // For the DCA distributions, we want to use exactly the same cuts
  // as for the other distributions, with the exception of the DCA cut
  fTrackCutsNoDCA = new AliESDtrackCuts(*fTrackCuts); // clone cuts
  // Reset all DCA cuts; FIXME: is this all?
  fTrackCutsNoDCA->SetMaxDCAToVertexXY();
  fTrackCutsNoDCA->SetMaxDCAToVertexZ ();
  fTrackCutsNoDCA->SetMaxDCAToVertexXYPtDep();
  fTrackCutsNoDCA->SetMaxDCAToVertexZPtDep();

  fTriggerAnalysis = new AliTriggerAnalysis();
  fTriggerAnalysis->SetAnalyzeMC(fIsMC);

  //The common PID object can then be retrieved from the input handler. This can naturally be done in the UserCreateOutputObjects:
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  fPIDResponse = inputHandler->GetPIDResponse();

}


void AliAnalysisTaskMultPbTracks::UserExec(Option_t *)
{
  // User code


  /* PostData(0) is taken care of by AliAnalysisTaskSE */
  PostData(1,fHistoManager);
  PostData(2,fTrackCuts);
  PostData(3,fCentrSelector);

  // Cache (some) histogram pointers to avoid loop in the histo manager
  static TH3D * hTracks  [AliAnalysisMultPbTrackHistoManager::kNHistos];
  static TH2D * hDCA     [AliAnalysisMultPbTrackHistoManager::kNHistos];
  static TH1D * hNTracks [AliAnalysisMultPbTrackHistoManager::kNHistos];
  hTracks[AliAnalysisMultPbTrackHistoManager::kHistoGen]        = fHistoManager->GetHistoPtEtaVz(AliAnalysisMultPbTrackHistoManager::kHistoGen       );
  hTracks[AliAnalysisMultPbTrackHistoManager::kHistoRec]        = fHistoManager->GetHistoPtEtaVz(AliAnalysisMultPbTrackHistoManager::kHistoRec       );
  hTracks[AliAnalysisMultPbTrackHistoManager::kHistoRecPrim]    = fHistoManager->GetHistoPtEtaVz(AliAnalysisMultPbTrackHistoManager::kHistoRecPrim   );
  hTracks[AliAnalysisMultPbTrackHistoManager::kHistoRecFake]    = fHistoManager->GetHistoPtEtaVz(AliAnalysisMultPbTrackHistoManager::kHistoRecFake   );
  hTracks[AliAnalysisMultPbTrackHistoManager::kHistoRecSecMat]  = fHistoManager->GetHistoPtEtaVz(AliAnalysisMultPbTrackHistoManager::kHistoRecSecMat );
  hTracks[AliAnalysisMultPbTrackHistoManager::kHistoRecSecWeak] = fHistoManager->GetHistoPtEtaVz(AliAnalysisMultPbTrackHistoManager::kHistoRecSecWeak);

  hDCA[AliAnalysisMultPbTrackHistoManager::kHistoGen]        = fHistoManager->GetHistoDCA(AliAnalysisMultPbTrackHistoManager::kHistoGen       );
  hDCA[AliAnalysisMultPbTrackHistoManager::kHistoRec]        = fHistoManager->GetHistoDCA(AliAnalysisMultPbTrackHistoManager::kHistoRec       );
  hDCA[AliAnalysisMultPbTrackHistoManager::kHistoRecPrim]    = fHistoManager->GetHistoDCA(AliAnalysisMultPbTrackHistoManager::kHistoRecPrim   );
  hDCA[AliAnalysisMultPbTrackHistoManager::kHistoRecFake]    = fHistoManager->GetHistoDCA(AliAnalysisMultPbTrackHistoManager::kHistoRecFake   );
  hDCA[AliAnalysisMultPbTrackHistoManager::kHistoRecSecMat]  = fHistoManager->GetHistoDCA(AliAnalysisMultPbTrackHistoManager::kHistoRecSecMat );
  hDCA[AliAnalysisMultPbTrackHistoManager::kHistoRecSecWeak] = fHistoManager->GetHistoDCA(AliAnalysisMultPbTrackHistoManager::kHistoRecSecWeak);

  hNTracks[AliAnalysisMultPbTrackHistoManager::kHistoGen]        = fHistoManager->GetHistoMult(AliAnalysisMultPbTrackHistoManager::kHistoGen );
  hNTracks[AliAnalysisMultPbTrackHistoManager::kHistoRec]        = fHistoManager->GetHistoMult(AliAnalysisMultPbTrackHistoManager::kHistoRec );


  fHistoManager->GetHistoPtEvent(AliAnalysisMultPbTrackHistoManager::kHistoRec)->Reset();

  fESD = dynamic_cast<AliESDEvent*>(fInputEvent);
  if (strcmp(fESD->ClassName(),"AliESDEvent")) {
    AliFatal("Not processing ESDs");
  }
  
  fHistoManager->GetHistoStats()->Fill(AliAnalysisMultPbTrackHistoManager::kStatAll);


  // Centrality selection
  Bool_t isCentralitySelected = fCentrSelector->IsCentralityBinSelected(fESD,fTrackCuts);  
  if(!isCentralitySelected) return;

  // AliESDCentrality *centrality = fESD->GetCentrality();
  // if(!centrality) {
  //   AliFatal("Centrality object not available"); 
  // }
  // else {
  //   Int_t centrBin = centrality->GetCentralityClass5(fCentralityEstimator.Data()) ;
    
  //   if (centrBin != fCentrBin && fCentrBin != -1) return;
  // }

  fHistoManager->GetHistoStats()->Fill(AliAnalysisMultPbTrackHistoManager::kStatCentr);

  Bool_t isSelected = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & fOfflineTrigger);

  if(!isSelected) return;
  fHistoManager->GetHistoStats()->Fill(AliAnalysisMultPbTrackHistoManager::kStatPhysSel);


  if (fIsMC) {
    // Get MC vertex
    //FIXME: which vertex do I take for MC?
    TArrayF   vertex;
    fMCEvent->GenEventHeader()->PrimaryVertex(vertex);
    Float_t zvGen = vertex[2];    

    if (!fMCEvent) {
      AliError("No MC info found");
    } else {
      
      //loop on the MC event, only over primaries, which are always
      //      the first in stack.
      // WRONG: D&B decays are produced in the transport... Need To Loop over all particles
      //      Int_t nMCTracks = fMCEvent->GetNumberOfPrimaries();
      Int_t nMCTracks = fMCEvent->GetNumberOfTracks();
      Int_t nPhysicalPrimaries = 0;
      for (Int_t ipart=0; ipart<nMCTracks; ipart++) { 
	
	AliMCParticle *mcPart  = (AliMCParticle*)fMCEvent->GetTrack(ipart);
	
	// We don't care about neutrals and non-physical primaries
	if(mcPart->Charge() == 0) continue;

	//check if current particle is a physical primary
	if(!IsPhysicalPrimaryAndTransportBit(ipart)) continue;
	if(TMath::Abs(mcPart->Zv()-zvGen)>1e-6) {
	  // This cures a bug in Hijing
	  // A little hack here: I put those in the underflow bin of the process type to keep track of them
	  fHistoManager->GetHistoProcess(AliAnalysisMultPbTrackHistoManager::kHistoGen)->Fill(-1);
	  continue;
	}

	nPhysicalPrimaries++;
	// Fill species histo and particle species
	fHistoManager->GetHistoProcess(AliAnalysisMultPbTrackHistoManager::kHistoGen)->Fill(mcPart->Particle()->GetUniqueID());
	fHistoManager->FillParticleID(AliAnalysisMultPbTrackHistoManager::kHistoGen, mcPart);
	
	//	Float_t zv = vtxESD->GetZ();
	// Fill generated histo
	hTracks[AliAnalysisMultPbTrackHistoManager::kHistoGen]->Fill(mcPart->Pt(),mcPart->Eta(),zvGen);
	Int_t partCode = fHistoManager->GetLocalParticleID(mcPart);
	//	cout << "Part " << partCode << endl;
	if (partCode == AliAnalysisMultPbTrackHistoManager::kPartPiPlus  || 
	    partCode == AliAnalysisMultPbTrackHistoManager::kPartPiMinus ||
	    partCode == AliAnalysisMultPbTrackHistoManager::kPartKPlus   || 
	    partCode == AliAnalysisMultPbTrackHistoManager::kPartKMinus  ||
	    partCode == AliAnalysisMultPbTrackHistoManager::kPartP       || 
	    partCode == AliAnalysisMultPbTrackHistoManager::kPartPBar  
	    ){
	  //	  cout << "Found " << partCode << endl;
	  
	  fHistoManager->GetHistoPtEtaVz(AliAnalysisMultPbTrackHistoManager::kHistoGen, partCode);
	}
      }
      hNTracks[AliAnalysisMultPbTrackHistoManager::kHistoGen]->Fill(nPhysicalPrimaries);
      fHistoManager->GetHistoVzEvent(AliAnalysisMultPbTrackHistoManager::kHistoGen)->Fill(zvGen);
  
    
    }
  }
  

  const AliESDVertex* vtxESD = fESD->GetPrimaryVertex();
  if(!vtxESD) return;
  // TODO: leave the cuts here or find a better place?)
  // Quality cut on vertexer Z, as suggested by Francesco Prino
  if(vtxESD->IsFromVertexerZ()) {
    if (vtxESD->GetNContributors() <= 0) return;
    if (vtxESD->GetDispersion() > 0.04) return;
    if (vtxESD->GetZRes() > 0.25) return;
  }
  // "Beam gas" vertex cut
  const AliESDVertex * vtxESDTPC= fESD->GetPrimaryVertexTPC(); 
  if(vtxESDTPC->GetNContributors()<1) return;
  if (vtxESDTPC->GetNContributors()<(-10.+0.25*fESD->GetMultiplicity()->GetNumberOfITSClusters(0)))     return;  

  // Fill  statistics
  fHistoManager->GetHistoStats()->Fill(AliAnalysisMultPbTrackHistoManager::kStatVtx);
  
  // ZDC cut, only ZNs
  Bool_t zdcA   = fTriggerAnalysis->ZDCTDCTrigger(fESD, AliTriggerAnalysis::kASide, kTRUE, kFALSE) ; 
  Bool_t zdcC   = fTriggerAnalysis->ZDCTDCTrigger(fESD, AliTriggerAnalysis::kCSide, kTRUE, kFALSE) ;			

  if (!(zdcA && zdcC) && (!fIsMC)) return;
  fHistoManager->GetHistoStats()->Fill(AliAnalysisMultPbTrackHistoManager::kStatZDCCut);


  if(TMath::Abs(vtxESD->GetZ()) > 10) return;
  fHistoManager->GetHistoStats()->Fill(AliAnalysisMultPbTrackHistoManager::kStatVtxRangeCut);

  // FIXME
  Float_t ntracksOk = fTrackCuts->CountAcceptedTracks(fESD);
  const AliMultiplicity * mult = fESD->GetMultiplicity();
  if (ntracksOk < 1000 &&  mult->GetNumberOfITSClusters(1) > 4500) {
    Float_t dummy;
    Printf("IEV: %d, Orbit: %x, Period: %d, BC: %d\n",fESD->GetEventNumberInFile(), fESD->GetOrbitNumber(),fESD->GetPeriodNumber(),fESD->GetBunchCrossNumber());
    printf("%s ----> Processing event # %lld\n", CurrentFileName(), Entry());
    cout << "File " << AliAnalysisManager::GetAnalysisManager()->GetTree()->GetCurrentFile()->GetName() << endl;   
    cout << "Nt " << ntracksOk << ", "<< fESD->GetNumberOfTracks() << " V0 " <<  fCentrSelector->GetCorrV0(fESD,dummy)  << " SPD1 " << mult->GetNumberOfITSClusters(1) << endl;       
  }


  
  // Fill Vertex

  fHistoManager->GetHistoVzEvent(AliAnalysisMultPbTrackHistoManager::kHistoRec)->Fill(vtxESD->GetZ());

  // loop on tracks
  Int_t ntracks = fESD->GetNumberOfTracks();
  Int_t acceptedTracks = 0;

  for (Int_t itrack = 0; itrack<ntracks; itrack++) {    
    AliESDtrack * esdTrack = fIsTPCOnly ? fTrackCuts->GetTPCOnlyTrack(fESD,itrack) : fESD->GetTrack(itrack);// FIXMEL TrackCuts or TrackCutsNoDCA for the TPC?
    if (!esdTrack) continue;
          
    // Fill DCA distibutions, without the DCA cut, Fill the other stuff, with all cuts!
    Bool_t accepted = fTrackCuts->AcceptTrack(esdTrack);
    Bool_t acceptedNoDCA = fTrackCutsNoDCA->AcceptTrack(esdTrack);

    // accepted = accepted && ((fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kElectron) > 2) ||
    // 			    (fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kPion    ) < 1) || 
    // 			    (fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kProton  ) < 1) || 
    // 			    (fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kKaon    ) < 1) 
    // 			    //			    (esdTrack->P() > 1.2)
    // 			    );//FIXME SKIP ELECTRONS below p = 1.2 gev, make configurable. Keep particles if they are in the crossing

    // acceptedNoDCA = acceptedNoDCA && ((fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kElectron) > 2) ||
    // 				      (fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kPion    ) < 1) || 
    // 				      (fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kProton  ) < 1) || 
    // 				      (fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kKaon    ) < 1) 
    // 				      //				      (esdTrack->P() > 1.2)
    // 				      );//FIXME SKIP ELECTRONS below p = 1.2 gev, make configurable. Keep particles if they are in the crossing

    if(accepted) acceptedTracks++;

    // Compute weighted offset
    // TODO: other choiches for the DCA?
    Double_t weightedDCA = 10;
    

#if defined WEIGHTED_DCA
    Double_t b = fESD->GetMagneticField();
    Double_t dca[2];
    Double_t cov[3];
    if (esdTrack->PropagateToDCA(vtxESD, b,3., dca, cov)) {
      Double_t det = cov[0]*cov[2]-cov[1]*cov[1]; 
      if (det<=0) {
	AliError("DCA Covariance matrix is not positive definite!");
      }
      else {
	weightedDCA = (dca[0]*dca[0]*cov[2]+dca[1]*dca[1]*cov[0]-2*dca[0]*dca[1]*cov[1])/det; 
	weightedDCA = weightedDCA>0 ? TMath::Sqrt(weightedDCA) : 10;
      }
      //      printf("dR:%e dZ%e  sigR2:%e sigRZ:%e sigZ2:%e\n",dca[0],dca[1],cov[0],cov[1],cov[2]);
    }
    
#elif defined TRANSVERSE_DCA
    Float_t xz[2]; 
    esdTrack->GetDZ(vtxESD->GetX(),vtxESD->GetY(),vtxESD->GetZ(), fESD->GetMagneticField(), xz); 
    weightedDCA = xz[0];
#endif
    

    // Alternative: p*DCA
    // Float_t xz[2]; 
    // esdTrack->GetDZ(vtxESD->GetX(),vtxESD->GetY(),vtxESD->GetZ(), fESD->GetMagneticField(), xz); 
    // Float_t dca = xz[0]*esdTrack->P();

  // This cures a bug in Hijing (displaced primaries)
  if (fIsMC) {
    Int_t label = TMath::Abs(esdTrack->GetLabel()); // no fakes!!!    
    if (IsPhysicalPrimaryAndTransportBit(label)) {
      AliMCParticle *mcPart  =  (AliMCParticle*)fMCEvent->GetTrack(label);
      if(!mcPart) AliFatal("No particle");
      TArrayF   vertex;
      fMCEvent->GenEventHeader()->PrimaryVertex(vertex);
      Float_t zvGen = vertex[2];    
      if(TMath::Abs(mcPart->Zv()-zvGen)>1e-6) continue;    
    }
  }
  // for each track
  if(accepted) {
    hTracks[AliAnalysisMultPbTrackHistoManager::kHistoRec]->Fill(esdTrack->Pt(),esdTrack->Eta(),vtxESD->GetZ());
    if (TMath::Abs(esdTrack->Eta())<0.5 && TMath::Abs(vtxESD->GetZ())<7)
      fHistoManager->GetHistoPtEvent(AliAnalysisMultPbTrackHistoManager::kHistoRec)->Fill(esdTrack->Pt());
  }
  if(acceptedNoDCA)
    hDCA[AliAnalysisMultPbTrackHistoManager::kHistoRec]->Fill(weightedDCA,esdTrack->Pt());

    // Fill separately primaries and secondaries
    // FIXME: fakes? Is this the correct way to account for them?
    // Get label and corresponding mcPart;
    if (fIsMC) {
      Int_t label = TMath::Abs(esdTrack->GetLabel()); // no fakes!!!
      AliMCParticle *mcPart  = label < 0 ? 0 : (AliMCParticle*)fMCEvent->GetTrack(label);
      if (!mcPart)  {
	if(accepted)
	  hTracks[AliAnalysisMultPbTrackHistoManager::kHistoRecFake]->Fill(esdTrack->Pt(),esdTrack->Eta(),vtxESD->GetZ());
	if(acceptedNoDCA)
	  hDCA[AliAnalysisMultPbTrackHistoManager::kHistoRecFake]->Fill(weightedDCA,esdTrack->Pt());
      }
      else {
	if(IsPhysicalPrimaryAndTransportBit(label)) {
	  // Fill species histo
	  fHistoManager->GetHistoProcess(AliAnalysisMultPbTrackHistoManager::kHistoRecPrim)->Fill(mcPart->Particle()->GetUniqueID());
	  if(accepted) {
	    hTracks[AliAnalysisMultPbTrackHistoManager::kHistoRecPrim]->Fill(esdTrack->Pt(),esdTrack->Eta(),vtxESD->GetZ());
	    Int_t partCode = fHistoManager->GetLocalParticleID(mcPart);
	    if (partCode == AliAnalysisMultPbTrackHistoManager::kPartPiPlus  || 
		partCode == AliAnalysisMultPbTrackHistoManager::kPartPiMinus ||
		partCode == AliAnalysisMultPbTrackHistoManager::kPartKPlus   || 
		partCode == AliAnalysisMultPbTrackHistoManager::kPartKMinus  ||
		partCode == AliAnalysisMultPbTrackHistoManager::kPartP       || 
		partCode == AliAnalysisMultPbTrackHistoManager::kPartPBar  
		)
	      fHistoManager->GetHistoPtEtaVz(AliAnalysisMultPbTrackHistoManager::kHistoRecPrim, partCode);
	  }
	  if(acceptedNoDCA)
	    hDCA[AliAnalysisMultPbTrackHistoManager::kHistoRecPrim]->Fill(weightedDCA,esdTrack->Pt());
	} 
	else {
	  Int_t mfl=0;
	  Int_t indexMoth=mcPart->Particle()->GetFirstMother();
	  if(indexMoth>=0){
	    TParticle* moth = fMCEvent->Stack()->Particle(indexMoth);
	    Float_t codemoth = TMath::Abs(moth->GetPdgCode());
	    mfl = Int_t (codemoth/ TMath::Power(10, Int_t(TMath::Log10(codemoth))));
	  }
	  if(mfl==3){ // strangeness
	    fHistoManager->GetHistoProcess(AliAnalysisMultPbTrackHistoManager::kHistoRecSecWeak)->Fill(mcPart->Particle()->GetUniqueID());
	    if(accepted)
	      hTracks[AliAnalysisMultPbTrackHistoManager::kHistoRecSecWeak]->Fill(esdTrack->Pt(),esdTrack->Eta(),vtxESD->GetZ());
	    if(acceptedNoDCA)
	      hDCA[AliAnalysisMultPbTrackHistoManager::kHistoRecSecWeak]->Fill(weightedDCA,esdTrack->Pt());	  
	  }else{ // material
	    fHistoManager->GetHistoProcess(AliAnalysisMultPbTrackHistoManager::kHistoRecSecMat)->Fill(mcPart->Particle()->GetUniqueID());
	    if(accepted)
	      hTracks[AliAnalysisMultPbTrackHistoManager::kHistoRecSecMat]->Fill(esdTrack->Pt(),esdTrack->Eta(),vtxESD->GetZ());
	    if(acceptedNoDCA)
	      hDCA[AliAnalysisMultPbTrackHistoManager::kHistoRecSecMat]->Fill(weightedDCA,esdTrack->Pt());	  

	  }


	}
      }
    }


  }
  //  cout << acceptedTracks << endl;
  
  hNTracks[AliAnalysisMultPbTrackHistoManager::kHistoRec]  ->Fill(acceptedTracks);

  // FIXME: this used to be filled with rescaled V0. I'm using raw V0 to test some newer production here.
  AliESDVZERO* esdV0 = fESD->GetVZEROData();
  Float_t multV0A=esdV0->GetMTotV0A();
  Float_t multV0C=esdV0->GetMTotV0C();
  Float_t multV0 = multV0A+multV0C;
  fHistoManager->GetHistoV0vsNtracks(AliAnalysisMultPbTrackHistoManager::kHistoRec)->Fill(acceptedTracks, multV0);
  // FIXME: the old version:
  // Float_t v0;
  // fHistoManager->GetHistoV0vsNtracks(AliAnalysisMultPbTrackHistoManager::kHistoRec)->Fill(acceptedTracks, fCentrSelector->GetCorrV0(fESD,v0));

  // Fill <pt> analysis histos
  fHistoManager->GetHistoPtEvent(AliAnalysisMultPbTrackHistoManager::kHistoRec)->Scale(1,"width");// correct bin width
  if (fHistoManager->GetHistoPtEvent(AliAnalysisMultPbTrackHistoManager::kHistoRec)->GetEntries()>0)
    fHistoManager->GetHistoMeanPt(AliAnalysisMultPbTrackHistoManager::kHistoRec)->Fill(fHistoManager->GetHistoPtEvent(AliAnalysisMultPbTrackHistoManager::kHistoRec)->GetMean());
  if (fHistoManager->GetHistoPtEvent(AliAnalysisMultPbTrackHistoManager::kHistoRec)->GetMean() >  
      fHistoManager->GetHistoPtEvent(AliAnalysisMultPbTrackHistoManager::kHistoRecHighestMeanPt)->GetMean()) {
    // Found a new highest pt: first reset than add it to the container for the highest
    fHistoManager->GetHistoPtEvent(AliAnalysisMultPbTrackHistoManager::kHistoRecHighestMeanPt)->Reset();
    fHistoManager->GetHistoPtEvent(AliAnalysisMultPbTrackHistoManager::kHistoRecHighestMeanPt)
      ->Add(fHistoManager->GetHistoPtEvent(AliAnalysisMultPbTrackHistoManager::kHistoRec));
  }
  if ((fHistoManager->GetHistoPtEvent(AliAnalysisMultPbTrackHistoManager::kHistoRec)->GetMean() <  
       fHistoManager->GetHistoPtEvent(AliAnalysisMultPbTrackHistoManager::kHistoRecLowestMeanPt)->GetMean() &&
       fHistoManager->GetHistoPtEvent(AliAnalysisMultPbTrackHistoManager::kHistoRec)->GetEntries()>0) ||
      !(fHistoManager->GetHistoPtEvent(AliAnalysisMultPbTrackHistoManager::kHistoRecLowestMeanPt)->GetEntries()>0)) {
    // Found a new lowest pt: first reset than add it to the container for the lowest
    fHistoManager->GetHistoPtEvent(AliAnalysisMultPbTrackHistoManager::kHistoRecLowestMeanPt)->Reset();
    fHistoManager->GetHistoPtEvent(AliAnalysisMultPbTrackHistoManager::kHistoRecLowestMeanPt)
      ->Add(fHistoManager->GetHistoPtEvent(AliAnalysisMultPbTrackHistoManager::kHistoRec));
  }
}

void   AliAnalysisTaskMultPbTracks::Terminate(Option_t *){
  // terminate

}


Bool_t AliAnalysisTaskMultPbTracks::IsPhysicalPrimaryAndTransportBit(Int_t ipart) {
  // Check if it is a primary and if it was transported
  Bool_t physprim=fMCEvent->IsPhysicalPrimary(ipart);
  if (!physprim) return kFALSE;
  Bool_t transported = ((AliMCParticle*)fMCEvent->GetTrack(ipart))->Particle()->TestBit(kTransportBit);
  if(!transported) return kFALSE;

  return kTRUE;

}

// void AliAnalysisTaskEvil::PrintProcInfo()
// {
//   ProcInfo_t info;
//   gSystem->GetProcInfo(&info);
//   AliInfo(Form("fMemResident=%10ld kB  fMemVirtual=%10ld kB",info.fMemResident,info.fMemVirtual));
// }
