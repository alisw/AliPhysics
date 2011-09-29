#include <exception>
#include "TChain.h"
#include "TTree.h"
#include "TVector3.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TH2.h"
#include "TH1.h"
#include "TH3.h"

#include "AliLog.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskPi0Reconstruction.h"

#include "AliESDEvent.h"
#include "AliESDInputHandler.h"

#include "AliAODCentrality.h"
#include "AliTPCPIDResponse.h"
#include "AliAODFmdCluster.h"
#include "AliKFParticle.h"
#include "AliTracker.h"
#include "AliV0Reader.h"
#include "AliAODv0.h"
#include "AliAODPid.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliAODHandler.h"
#include "AliConversionAODBGHandlerRP.h"
#include "AliV0Reader.h"
#include "AliV0ReaderV1.h"
#include "AliCDBEntry.h"
#include "TObjArray.h"
#include "AliEventplane.h"
#include "AliCentrality.h"
#include "AliVCluster.h"
#include "AliESDCaloCluster.h"
#include <iostream>

#include <exception>

// Author Daniel Lohner (Daniel.Lohner@cern.ch)

using namespace std;

ClassImp(AliAnalysisTaskPi0Reconstruction)


//________________________________________________________________________
AliAnalysisTaskPi0Reconstruction::AliAnalysisTaskPi0Reconstruction(const char *name) : AliV0ReaderV1(name),
    fDeltaAODBranchName("GammaConv_gamma"),
    fBGHandler(NULL),
    Pi0MassRange(NULL),
    kEventMixing(kTRUE),
    fIsHeavyIon(kTRUE),
    fRandomizer(NULL),
    kUseSatelliteAODs(kTRUE),
    fPHOSGammas(NULL),
    fEMCALGammas(NULL),
    fPi0Candidates(NULL),
    fBGPi0s(NULL),
    fMCTruePi0s(NULL),
    fRapidityCut(0.9),
    fAlphaCut(1),
    kUseOnlyTaggedPhotons(kFALSE),
    fNRandomEventsForBGRotation(15)
{

    fRandomizer= new TRandom3();

    // If V0 Reader is available use it
    SetUseAODConversionPhoton(kTRUE);

    // Default Values
    Pi0MassRange=new Double_t[2];
    Pi0MassRange[0]=0.;
    Pi0MassRange[1]=0.3;

    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
}
//________________________________________________________________________
void AliAnalysisTaskPi0Reconstruction::SetBGHandler(AliConversionAODBGHandlerRP *bgHandler)
{
    fBGHandler=bgHandler;
    fNCentralityBins=fBGHandler->GetNCentralityBins();
}
//________________________________________________________________________
void AliAnalysisTaskPi0Reconstruction::SetDefaultBGHandler()
{
    AliInfo("Setting up default BGHandler");
    fBGHandler=new AliConversionAODBGHandlerRP();
    fNCentralityBins=fBGHandler->GetNCentralityBins();
}

//________________________________________________________________________
AliAnalysisTaskPi0Reconstruction::~AliAnalysisTaskPi0Reconstruction(){

    if(fBGHandler){delete fBGHandler;
	fBGHandler=0x0;}
    if(fRandomizer){delete fRandomizer;fRandomizer=0x0;}

    if(Pi0MassRange)delete[] Pi0MassRange;

}

//________________________________________________________________________
void AliAnalysisTaskPi0Reconstruction::UserCreateOutputObjects()
{
  // Create histograms
  // Called once
    // Create the output container
 /* if(fOutputList != NULL){
    delete fOutputList;
    fOutputList = NULL;
  }
  if(fOutputList == NULL){
    fOutputList = new TList();
    fOutputList->SetOwner(kTRUE);
    }*/
    AliV0ReaderV1::UserCreateOutputObjects();

    TList *fPi0List=new TList();
    fPi0List->SetName("Pi0Reconstruction");
    fPi0List->SetOwner(kTRUE);
    fOutputList->Add(fPi0List);


    hPi0Cuts=new TH1F("Pi0Cand_Cuts","Pi0Cand_Cuts" ,10,-0.5,9.5);
    fPi0List->Add(hPi0Cuts);
    hPi0BGCuts=new TH1F("Pi0Cand_BG_Cuts","Pi0Cand_BG_Cuts" ,10,-0.5,9.5);
    fPi0List->Add(hPi0BGCuts);
    hPi0Alpha=new TH1F("Pi0Cand_Alpha" ,"Pi0Cand Alpha" ,200,0,1);
    fPi0List->Add(hPi0Alpha);
    hPi0OpeningAngle=new TH1F("Pi0Cand_OpeningAngle" ,"Pi0Cand OpeningAngle" ,400,0,TMath::Pi());
    fPi0List->Add(hPi0OpeningAngle);
    hPi0Rapidity=new TH1F("Pi0Cand_Rapidity" ,"Pi0Cand_Rapidity" ,100,-1,1);
    fPi0List->Add(hPi0Rapidity);
    fPi0List->Add(hPi0BGCuts);
    hPi0AlphaRejected=new TH1F("Pi0Cand_Alpha_Rejected" ,"Pi0Cand Alpha Rejected" ,200,0,1);
    fPi0List->Add(hPi0AlphaRejected);
    hPi0OpeningAngleRejected=new TH1F("Pi0Cand_OpeningAngle_Rejected" ,"Pi0Cand OpeningAngle Rejected" ,400,0,TMath::Pi());
    fPi0List->Add(hPi0OpeningAngleRejected);
    hPi0RapidityRejected=new TH1F("Pi0Cand_Rapidity_Rejected" ,"Pi0Cand_Rapidity_Rejected" ,100,-1,1);
    fPi0List->Add(hPi0RapidityRejected);

    hPool=new TH3F("BGPool","BGPool",fNCentralityBins,-0.5,fNCentralityBins-0.5,7,-0.5,6.5,fBGHandler->GetNRPBins(),-0.5,fBGHandler->GetNRPBins()-0.5);
    fPi0List->Add(hPool);

    // MC

    if(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()){

	Int_t kGCnXBinsSpectra = Int_t((Pi0MassRange[1]-Pi0MassRange[0])*500);  //500 for range 0 - 1
	Double_t kGCfirstXBinSpectra = Pi0MassRange[0];
	Double_t kGClastXBinSpectra = Pi0MassRange[1];

	Int_t kGCnYBinsSpectra = 250;
	Double_t kGCfirstYBinSpectra = 0.;
	Double_t kGClastYBinSpectra = 25.;

	hPi0TRUE=new TH1F*[fNCentralityBins];
	hPi0RECOTRUE=new TH2F*[fNCentralityBins];

	for(Int_t m=0;m<fNCentralityBins;m++){

	    hPi0TRUE[m]=new TH1F(Form("%d_TRUE_Pi0_Pt",m),Form("%d_TRUE_Pi0_Pt",m),kGCnYBinsSpectra, kGCfirstYBinSpectra, kGClastYBinSpectra);
	    hPi0TRUE[m]->Sumw2();
	    fPi0List->Add(hPi0TRUE[m]);

	    hPi0RECOTRUE[m]=new TH2F(Form("%d_RECOTRUE_Pi0_Pt",m),Form("%d_RECOTRUE_Pi0_Pt",m),kGCnXBinsSpectra, kGCfirstXBinSpectra, kGClastXBinSpectra,kGCnYBinsSpectra, kGCfirstYBinSpectra, kGClastYBinSpectra);
	    hPi0RECOTRUE[m]->Sumw2();
	    fPi0List->Add(hPi0RECOTRUE[m]);

	}
    }

  PostData(1, fOutputList);
}

//________________________________________________________________________
void AliAnalysisTaskPi0Reconstruction::UserExec(Option_t *) 
{


    AliV0ReaderV1::UserExec("");

    // Check if BG Handler exists

    if(!fBGHandler){SetDefaultBGHandler();}
    
    // Process Event Info

    if(fEventIsSelected){

	GetConversionGammas();

	CalculatePi0Candidates();

	CalculateBackground();
      
	ProcessMCMesons();
    }

    PostData(1, fOutputList);
}

//________________________________________________________________________

void AliAnalysisTaskPi0Reconstruction::ProcessMCMesons(){

    if(fMCStack){

	if(fMCTruePi0s == NULL){
	fMCTruePi0s = new TClonesArray("TParticle",100);
	}
	fMCTruePi0s->Delete();//Reset the TClonesArray

	for (Int_t part=0; part<fMCStack->GetNprimary(); part++) {

	    TParticle *fMCMother = fMCStack->Particle(part);

	    if(fMCMother){

		// Process only Pions

		if(fMCMother->GetPdgCode()==111){

		    // Histogram for pi0->2y efficiency (take acceptance)
		    if(IsMCMesonInReconstructionAcceptance(fMCMother)){
			hPi0TRUE[fCentralityBin]->Fill(fMCMother->Pt());
			fMCTruePi0s->AddAt(fMCMother,fMCTruePi0s->GetEntriesFast());
		    }
		}
	    }
	}
    }
}

//________________________________________________________________________
Bool_t AliAnalysisTaskPi0Reconstruction::IsMCMesonInReconstructionAcceptance(TParticle *fMCMother){

    // Select only pi0->2y
    if(fMCMother->GetNDaughters()!=2){return kFALSE;}

    // PseudoRapidity Cut
    if(TMath::Abs(fMCMother->Eta())>fRapidityCut){return kFALSE;}

    Int_t NGammaAcceptance=0;

    for(Int_t i=0;i<2;i++){
	TParticle *MDaughter=fMCStack->Particle(fMCMother->GetDaughter(i));

	// Is Daughter a Photon?
	if(MDaughter->GetPdgCode()==22){

	    // Gamma Acceptance Cut
	    if(IsMCConversionGammaInAcceptance(MDaughter)){
		NGammaAcceptance++;}
	}
    }

    if(NGammaAcceptance!=2){return kFALSE;}

    return kTRUE;
}

//________________________________________________________________________

Bool_t AliAnalysisTaskPi0Reconstruction::IsTruePi0(AliAODConversionMother *pi0){

    if(fMCStack){

	AliAODConversionPhoton *gamma0=dynamic_cast<AliAODConversionPhoton*>(fConversionGammas->At(pi0->GetLabel(0)));
	AliAODConversionPhoton *gamma1=dynamic_cast<AliAODConversionPhoton*>(fConversionGammas->At(pi0->GetLabel(1)));

        TParticle *mcgam0=gamma0->GetMCParticle(fMCStack);
	TParticle *mcgam1=gamma1->GetMCParticle(fMCStack);
        TParticle *mcpi0=NULL;

	    if(mcgam0&&mcgam1){
		if(mcgam0->GetMother(0)==mcgam1->GetMother(0)){
		    mcpi0=fMCStack->Particle(mcgam0->GetMother(0));

		    if(mcpi0){
			if(mcpi0->GetPdgCode()==111){return kTRUE;}
		    }
		}
	    }

    }
return kFALSE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskPi0Reconstruction::IsMCPhotonReconstructed(TParticle *MDaughter){

    for(Int_t ii=0;ii<fConversionGammas->GetEntriesFast();ii++){

	AliAODConversionPhoton *gamma=(AliAODConversionPhoton *)fConversionGammas->At(ii);

	TParticle *fMCParticle=gamma->GetMCParticle(fMCStack);

	if(fMCParticle){

	    if(fMCParticle->GetDaughter(0)==MDaughter->GetDaughter(0)&&fMCParticle->GetDaughter(1)==MDaughter->GetDaughter(1)){
		return kTRUE;
	    }
	}
    }
    return kFALSE;
}

//________________________________________________________________________
void AliAnalysisTaskPi0Reconstruction::CalculatePi0Candidates(){

    if(fPi0Candidates == NULL){
	fPi0Candidates = new TClonesArray("AliAODConversionMother",100);
    }
    fPi0Candidates->Delete();//Reset the TClonesArray

    if(fConversionGammas->GetEntriesFast()>1){

	for(Int_t firstGammaIndex=0;firstGammaIndex<fConversionGammas->GetEntriesFast()-1;firstGammaIndex++){

	    AliAODConversionPhoton *gamma0=dynamic_cast<AliAODConversionPhoton*>(fConversionGammas->At(firstGammaIndex));

	    // Combine Photons

	    for(Int_t secondGammaIndex=firstGammaIndex+1;secondGammaIndex<fConversionGammas->GetEntriesFast();secondGammaIndex++){

		AliAODConversionPhoton *gamma1=dynamic_cast<AliAODConversionPhoton*>(fConversionGammas->At(secondGammaIndex));

		//Check for same Electron ID
		if(gamma0->GetTrackLabelPositive()==gamma1->GetTrackLabelPositive()||gamma0->GetTrackLabelNegative()==gamma1->GetTrackLabelNegative()
		   ||gamma0->GetTrackLabelNegative()==gamma1->GetTrackLabelPositive()||gamma0->GetTrackLabelPositive()==gamma1->GetTrackLabelNegative())continue;

		AliAODConversionMother *pi0cand=new AliAODConversionMother(gamma0,gamma1);
		pi0cand->SetLabels(firstGammaIndex,secondGammaIndex);


		if((CheckPi0(pi0cand))){
		    if (pi0cand->M() > Pi0MassRange[0] && pi0cand->M() < Pi0MassRange[1] ){

                        // Add Pi0 to Stack
			new((*fPi0Candidates)[fPi0Candidates->GetEntriesFast()]) AliAODConversionMother(*pi0cand);

			// Test MC truth
			if(IsTruePi0(pi0cand)){
			    hPi0RECOTRUE[fCentralityBin]->Fill(pi0cand->M(),pi0cand->Pt());
			}
                         
		    }
		}

		delete pi0cand;
		pi0cand=0x0;

	    }

	}

    }
}


//________________________________________________________________________
void AliAnalysisTaskPi0Reconstruction::Terminate(Option_t *) 
{
 
    // Draw result to the screen
    fOutputList->Print();
    // Called once at the end of the query
}


//________________________________________________________________________
void AliAnalysisTaskPi0Reconstruction::RotateParticle(AliAODConversionPhoton *gamma,Double_t angle){
    Int_t fNDegreesPMBackground=15;
    Double_t nRadiansPM = fNDegreesPMBackground*TMath::Pi()/180;
    Double_t rotationValue = fRandomizer->Rndm()*2*nRadiansPM + TMath::Pi()-nRadiansPM;
    gamma->RotateZ(rotationValue);

}

//________________________________________________________________________

void AliAnalysisTaskPi0Reconstruction::CalculateBackground(){

    if(fBGPi0s == NULL){
	fBGPi0s = new TClonesArray("AliAODConversionMother",100);
    }
    fBGPi0s->Delete();//Reset the TClonesArray

    AliAODConversionMother *BGcandidate=NULL;
    TClonesArray *currentEventGammas=fConversionGammas;

    //Rotation Method
    if(!kEventMixing){

	// Correct for the number of rotations
	// BG is for rotation the same, except for factor NRotations
	Double_t weight=1./Double_t(fNRandomEventsForBGRotation);

	for(Int_t firstGammaIndex=0;firstGammaIndex<currentEventGammas->GetEntriesFast();firstGammaIndex++){

	    AliAODConversionPhoton *gamma1=dynamic_cast<AliAODConversionPhoton*>(currentEventGammas->At(firstGammaIndex));

	    for(Int_t secondGammaIndex=firstGammaIndex+1;secondGammaIndex<currentEventGammas->GetEntriesFast();secondGammaIndex++){
		for(Int_t nRandom=0;nRandom<fNRandomEventsForBGRotation;nRandom++){
	
		    AliAODConversionPhoton *gamma2=dynamic_cast<AliAODConversionPhoton*>(currentEventGammas->At(secondGammaIndex));

		    RotateParticle(gamma2);

		    BGcandidate=new AliAODConversionMother(gamma1,gamma2);

		    if(CheckPi0(BGcandidate,kFALSE)){

			if (BGcandidate->M() > Pi0MassRange[0] && BGcandidate->M() < Pi0MassRange[1] ){

			    new((*fBGPi0s)[fBGPi0s->GetEntriesFast()]) AliAODConversionMother(*BGcandidate);

                            dynamic_cast<AliAODConversionMother*>(fBGPi0s->At(fBGPi0s->GetEntriesFast()-1))->SetWeight(weight);
			}
		    }
		    delete BGcandidate;
                    BGcandidate=NULL;

		}

	    }
	}
    }

    else{
    // Do Event Mixing

  	Int_t psibin=fBGHandler->GetRPBinIndex(fEPAngle);

	for(Int_t nEventsInBG=0;nEventsInBG <fBGHandler->GetNBGEvents(psibin,fBGHandler->GetZBinIndex(fVertexZ),fCentralityBin);nEventsInBG++){
         
	    AliGammaConversionPhotonVector *previousEventGammas = fBGHandler->GetBGGoodGammas(psibin,fBGHandler->GetZBinIndex(fVertexZ),fCentralityBin,nEventsInBG);

	    if(previousEventGammas){

		// test weighted background
		Double_t weight=1.0;///Double_t(fBGHandler->GetNBGEvents(psibin,zbin,fCentralityBin));
		// Correct for the number of eventmixing:
		// N gammas -> (N-1) + (N-2) +(N-3) ...+ (N-(N-1))  using sum formula sum(i)=N*(N-1)/2  -> N*(N-1)/2
		// real combinations (since you cannot combine a photon with its own)
		// but BG leads to N_{a}*N_{b} combinations
		weight*=0.5*(Double_t(currentEventGammas->GetEntriesFast()-1))/Double_t(previousEventGammas->size());

		for(Int_t iCurrent=0;iCurrent<currentEventGammas->GetEntriesFast();iCurrent++){
		    AliAODConversionPhoton *currentEventGamma = (AliAODConversionPhoton*)(currentEventGammas->At(iCurrent));
		    for(UInt_t iPrevious=0;iPrevious<previousEventGammas->size();iPrevious++){
          
			AliAODConversionPhoton *previousEventGamma = (AliAODConversionPhoton*)(previousEventGammas->at(iPrevious));

			BGcandidate=new AliAODConversionMother(previousEventGamma,currentEventGamma);

			if(CheckPi0(BGcandidate,kFALSE)){

			    if (BGcandidate->M() > Pi0MassRange[0] && BGcandidate->M() < Pi0MassRange[1] ){

				new((*fBGPi0s)[fBGPi0s->GetEntriesFast()]) AliAODConversionMother(*BGcandidate);
				dynamic_cast<AliAODConversionMother*>(fBGPi0s->At(fBGPi0s->GetEntriesFast()-1))->SetWeight(weight);
			    }
			}
			delete BGcandidate;
			BGcandidate=NULL;
		    }
		}
	    }
	}

        // Add Current Event to BG POOL
	fBGHandler->AddEvent(fConversionGammas,fEPAngle,fVertexZ,Int_t(fCentrality));
	hPool->Fill(fCentralityBin,fBGHandler->GetZBinIndex(fVertexZ),fBGHandler->GetRPBinIndex(fEPAngle));

    }
}

//________________________________________________________________________
Bool_t AliAnalysisTaskPi0Reconstruction::GetPHOSGammas(){

    // TClones Array for local working copy

   /*
    if(fPHOSGammas == NULL){
	fPHOSGammas = new TClonesArray("AliAODConversionPhoton",100);
    }
    fPHOSGammas->Delete();//Reset the TClonesArray

    //AliESDCaloCluster *clu=NULL;
    AliVCluster *clu=NULL;
    TLorentzVector pPHOS;
    AliAODConversionPhoton *gammaPHOS=NULL;

    Double_t vtx[3];
    vtx[0] = fESDEvent->GetPrimaryVertex()->GetX();
    vtx[1] = fESDEvent->GetPrimaryVertex()->GetY();
    vtx[2] = fESDEvent->GetPrimaryVertex()->GetZ();


  for (Int_t i=0; i<fESDEvent->GetNumberOfCaloClusters(); i++){
      clu = fESDEvent->GetCaloCluster(i);

      // Cuts from AliAnalysisTaskCaloConv
      if ( !clu->IsPHOS() || clu->E()<0.1 ) continue;
      clu ->GetMomentum(pPHOS ,vtx);
      if(pPHOS.Energy()<0.25){
	  continue ;
      }
      if(clu->GetNCells()<=2){
	  continue ;
      }

    Bool_t isNeutral = kTRUE ;
    Bool_t isDispOK = kTRUE ;
    Bool_t isTOFOK = kTRUE ;
    Int_t iMod,iX,iZ ;

    isNeutral = clu->GetEmcCpvDistance()>5. ;  //To be improved
    isDispOK = kFALSE ;
    Double_t l0=clu->GetM02(),l1=clu->GetM20() ;
    if(l1>= 0   && l0>= 0   && l1 < 0.1 && l0 < 0.1) isDispOK=kFALSE ;
    if(l1>= 0   && l0 > 0.5 && l1 < 0.1 && l0 < 1.5) isDispOK=kTRUE ;
    if(l1>= 0   && l0 > 2.0 && l1 < 0.1 && l0 < 2.7) isDispOK=kFALSE ;
    if(l1>= 0   && l0 > 2.7 && l1 < 0.1 && l0 < 4.0) isDispOK=kFALSE ;
    if(l1 > 0.1 && l1 < 0.7 && l0 > 0.7 && l0 < 2.1) isDispOK=kTRUE ;
    if(l1 > 0.1 && l1 < 0.3 && l0 > 3.0 && l0 < 5.0) isDispOK=kFALSE  ;
    if(l1 > 0.3 && l1 < 0.7 && l0 > 2.5 && l0 < 4.0) isDispOK=kFALSE ;
    if(l1 > 0.7 && l1 < 1.3 && l0 > 1.0 && l0 < 1.6) isDispOK=kTRUE ;
    if(l1 > 0.7 && l1 < 1.3 && l0 > 1.6 && l0 < 3.5) isDispOK=kTRUE ;
    if(l1 > 1.3 && l1 < 3.5 && l0 > 1.3 && l0 < 3.5) isDispOK=kTRUE ;

    Float_t xyz[3] = {0,0,0};
    clu->GetPosition(xyz);   //Global position in ALICE system
    TVector3 global(xyz) ;
    Int_t relid[4] ;
    if(!fPHOSgeom->GlobalPos2RelId(global,relid)){
	printf("PHOS_beyond: x=%f, y=%f, z=%f \n",xyz[0],xyz[1],xyz[2]) ;
	continue ;
    }
    iMod=relid[0] ;
    iX=relid[2];
    iZ=relid[3] ;
    if(!IsGoodChannel("PHOS",iMod,iX,iZ))
	continue ;

    // Add selected cluster to Photonqueue
    gammaPHOS=new AliAODConversionPhoton();
    gammaPHOS->SetXYZM(pPHOS.Px(),pPHOS.Py(),pPHOS.Pz(),0.);
    new((*fPHOSGammas)[fPHOSGammas->GetEntriesFast()]) AliAODConversionPhoton(*gammaPHOS);
  }

  if(fPHOSGammas->GetEntries()){return kTRUE;}
*/
  return kFALSE;


}

//________________________________________________________________________
Bool_t AliAnalysisTaskPi0Reconstruction::GetConversionGammas(){

    // ESD mode

    if(fESDEvent){

	fConversionGammas=GetReconstructedGammas();

	if(!fConversionGammas)AliError("No Gammas from V0 Reader");
    }


    // AOD mode

    if(fAODEvent&&kUseSatelliteAODs){

	if(fConversionGammas == NULL){
	    fConversionGammas = new TClonesArray("AliAODConversionPhoton",100);
	}
	fConversionGammas->Delete();//Reset the TClonesArray
   
	//Get Gammas from satellite AOD gamma branch

	AliAODConversionPhoton *gamma=0x0;
    
	TClonesArray *fInputGammas=dynamic_cast<TClonesArray*>(fAODEvent->FindListObject(fDeltaAODBranchName.Data()));
	if(!fInputGammas){FindDeltaAODBranchName();
	    fInputGammas=dynamic_cast<TClonesArray*>(fAODEvent->FindListObject(fDeltaAODBranchName.Data()));}
	if(!fInputGammas){AliError("No Gamma Satellites found");return kFALSE;}

	// Apply Selection Cuts to Gammas and create local working copy
	if(fInputGammas){
	    for(Int_t i=0;i<fInputGammas->GetEntriesFast();i++){
		gamma=dynamic_cast<AliAODConversionPhoton*>(fInputGammas->At(i));
		if(gamma){
		    if(IsGammaCandidate(gamma)){
			new((*fConversionGammas)[fConversionGammas->GetEntriesFast()]) AliAODConversionPhoton(*gamma);}
		}
	    }
	}

    }

    if(fConversionGammas->GetEntries()){return kTRUE;}

    return kFALSE;
}

//________________________________________________________________________
void AliAnalysisTaskPi0Reconstruction::FindDeltaAODBranchName(){
    TList *list=fAODEvent->GetList();
    for(Int_t ii=0;ii<list->GetEntries();ii++){
	TString name((list->At(ii))->GetName());
	if(name.BeginsWith("GammaConv")&&name.EndsWith("gamma")){
	    fDeltaAODBranchName=name;
            AliInfo(Form("Set DeltaAOD BranchName to: %s",fDeltaAODBranchName.Data()));
	    return;
	}
    }
}


///________________________________________________________________________
Bool_t AliAnalysisTaskPi0Reconstruction::CheckPi0(AliAODConversionMother *pi0,Bool_t IsSignal)
{

    Bool_t passcuts=kTRUE;

    TH1 *hist=0x0;

    if(IsSignal){hist=hPi0Cuts;
    }
    else{hist=hPi0BGCuts;}

    // Undefined Rapidity -> Floating Point exception
    if((pi0->E()+pi0->Pz())/(pi0->E()-pi0->Pz())<=0){
        hist->Fill(0);
	passcuts=kFALSE;
    }
    else{
	// PseudoRapidity Cut
	if(TMath::Abs(pi0->PseudoRapidity())>fRapidityCut){
	    hist->Fill(1);
	    passcuts=kFALSE;
	}
    }

    // Opening Angle Cut
    if(pi0->GetOpeningAngle()<0.005){
        hist->Fill(2);
	passcuts=kFALSE;
    }

    // Alpha (Energy Asymmetry) Cut
    if(pi0->GetAlpha()>fAlphaCut){
	hist->Fill(3);
	passcuts=kFALSE;
    }

    // Fill Histograms

    if(passcuts){
	hist->Fill(9);
        if(IsSignal){
	hPi0Rapidity->Fill(pi0->PseudoRapidity());
	hPi0Alpha->Fill(pi0->GetAlpha());
	hPi0OpeningAngle->Fill(pi0->GetOpeningAngle());

	}
    }
    else{
	hist->Fill(8);
	if(IsSignal){
	    hPi0RapidityRejected->Fill(pi0->PseudoRapidity());
	    hPi0AlphaRejected->Fill(pi0->GetAlpha());
            hPi0OpeningAngleRejected->Fill(pi0->GetOpeningAngle());
	}
    }

    return passcuts;
}
