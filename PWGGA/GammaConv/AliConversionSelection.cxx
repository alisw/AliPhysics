#include "AliConversionSelection.h"
#include "AliAODHeader.h"
#include "AliVVZERO.h"
#include "AliMultiplicity.h"
#include <iostream>


// Author Daniel Lohner (Daniel.Lohner@cern.ch)

using namespace std;

ClassImp(AliConversionSelection)


//________________________________________________________________________
AliConversionSelection::AliConversionSelection(AliConversionCuts *convCut, AliConversionMesonCuts *mesonCut) : TObject(),
    fInputEvent(NULL),
    fMCEvent(NULL),
    fConversionCut(convCut),
    fMesonCut(mesonCut),
    fESDTrackCuts(NULL),
    fGoodGammas(NULL),
    fPi0Candidates(NULL),
    fBGPi0s(NULL),
    fRandomizer(NULL),
    fBGHandler(NULL),
    fCurrentEventNumber(-1),
    fIsOwner(kFALSE)
{
    // Default Values
    fInvMassRange[0]=0.05;
    fInvMassRange[1]=0.3;
}

//________________________________________________________________________
AliConversionSelection::AliConversionSelection(TString convCut, TString mesonCut) : TObject(),
    fInputEvent(NULL),
    fMCEvent(NULL),
    fConversionCut(NULL),
    fMesonCut(NULL),
    fESDTrackCuts(NULL),
    fGoodGammas(NULL),
    fPi0Candidates(NULL),
    fBGPi0s(NULL),
    fRandomizer(NULL),
    fBGHandler(NULL),
    fCurrentEventNumber(-1),
    fIsOwner(kTRUE)
{
    // Default Values
    fInvMassRange[0]=0.05;
    fInvMassRange[1]=0.3;

    fConversionCut = new AliConversionCuts();
    fConversionCut -> InitializeCutsFromCutString(convCut.Data());
    fMesonCut = new AliConversionMesonCuts();
    fMesonCut -> InitializeCutsFromCutString(mesonCut.Data());

}


//________________________________________________________________________
AliConversionSelection::AliConversionSelection(const AliConversionSelection &ref) : TObject(ref),
    fInputEvent(NULL),
    fMCEvent(NULL),
    fConversionCut(NULL),
    fMesonCut(NULL),
    fESDTrackCuts(NULL),
    fGoodGammas(NULL),
    fPi0Candidates(NULL),
    fBGPi0s(NULL),
    fRandomizer(NULL),
    fBGHandler(NULL),
    fCurrentEventNumber(-1),
    fIsOwner(kTRUE)
{
    // Copy Constructor

    fConversionCut=new AliConversionCuts(*ref.fConversionCut);
    fMesonCut=new AliConversionMesonCuts(*ref.fMesonCut);

    fInvMassRange[0]=ref.fInvMassRange[0];
    fInvMassRange[1]=ref.fInvMassRange[1];
}

//________________________________________________________________________
AliConversionSelection::~AliConversionSelection(){

    if(fBGHandler){
	delete fBGHandler;
	fBGHandler=NULL;
    }
    if(fRandomizer){
	delete fRandomizer;
	fRandomizer=NULL;
    }
    if(fPi0Candidates){
	delete fPi0Candidates;
	fPi0Candidates=NULL;
    }
    if(fBGPi0s){
	delete fBGPi0s;
	fBGPi0s=NULL;
    }
    if(fESDTrackCuts){
	delete fESDTrackCuts;
	fESDTrackCuts=NULL;
    }
    if(fIsOwner){
	if(fConversionCut){
	    delete fConversionCut;
	    fConversionCut=NULL;
	}
	if(fMesonCut){
	    delete fMesonCut;
	    fMesonCut=NULL;
	}
    }
}

//________________________________________________________________________
Bool_t AliConversionSelection::ProcessEvent(TClonesArray *photons,AliVEvent *inputEvent,AliMCEvent *mcEvent){
    fInputEvent=inputEvent;
    fMCEvent=mcEvent;

    // Protection
    Int_t eventnumber=GetEventNumber(inputEvent);
    if(eventnumber==fCurrentEventNumber){
	AliWarning("Event already analyzed! Return.");
	return kFALSE;
    }
    else{
	fCurrentEventNumber=eventnumber;
    }

    // Initialize and Reset Arrays
    if(fGoodGammas == NULL){
	fGoodGammas=new TObjArray(30);
    }
    fGoodGammas->Clear();

    if(fPi0Candidates == NULL){
	fPi0Candidates = new TClonesArray("AliAODConversionMother",100);
    }
    fPi0Candidates->Delete();

    if(fBGPi0s == NULL){
	fBGPi0s = new TClonesArray("AliAODConversionMother",100);
    }
    fBGPi0s->Delete();


    if(!photons||!fInputEvent)return kFALSE;
    
    if(!fConversionCut->EventIsSelected(fInputEvent,fMCEvent))return kFALSE;

    // Select photons
    for(Int_t i = 0; i < photons->GetEntriesFast(); i++) {
	AliAODConversionPhoton* gamma =dynamic_cast<AliAODConversionPhoton*>(photons->At(i));
	if(!gamma) continue;
	if(!fConversionCut->PhotonIsSelected(gamma,fInputEvent))continue;
	fGoodGammas->Add(gamma);
    }

    // Do MC Smearing
    Double_t *fUnsmearedPx=NULL;
    Double_t *fUnsmearedPy=NULL;
    Double_t *fUnsmearedPz=NULL;
    Double_t *fUnsmearedE=NULL;
   
    if(fMesonCut->UseMCPSmearing() && fMCEvent){
	fUnsmearedPx = new Double_t[fGoodGammas->GetEntries()]; // Store unsmeared Momenta
	fUnsmearedPy = new Double_t[fGoodGammas->GetEntries()];
	fUnsmearedPz = new Double_t[fGoodGammas->GetEntries()];
	fUnsmearedE =  new Double_t[fGoodGammas->GetEntries()];

	for(Int_t gamma=0;gamma<fGoodGammas->GetEntries();gamma++){ // Smear the AODPhotons in MC
	    fUnsmearedPx[gamma] = ((AliAODConversionPhoton*)fGoodGammas->At(gamma))->Px();
	    fUnsmearedPy[gamma] = ((AliAODConversionPhoton*)fGoodGammas->At(gamma))->Py();
	    fUnsmearedPz[gamma] = ((AliAODConversionPhoton*)fGoodGammas->At(gamma))->Pz();
	    fUnsmearedE[gamma] =  ((AliAODConversionPhoton*)fGoodGammas->At(gamma))->E();
	    fMesonCut->SmearParticle(dynamic_cast<AliAODConversionPhoton*>(fGoodGammas->At(gamma)));
	}
    }

    // Reconstruct Pi0 and BG
    CalculatePi0Candidates();
    CalculateBackground();
    if(fBGHandler)fBGHandler->AddEvent(fGoodGammas,fInputEvent);

    // Undo MC Smearing
    if(fMesonCut->UseMCPSmearing() && fMCEvent){
	for(Int_t gamma=0;gamma<fGoodGammas->GetEntries();gamma++){ // Smear the AODPhotons in MC
	    ((AliAODConversionPhoton*)fGoodGammas->At(gamma))->SetPx(fUnsmearedPx[gamma]); // Reset Unsmeared Momenta
	    ((AliAODConversionPhoton*)fGoodGammas->At(gamma))->SetPy(fUnsmearedPy[gamma]);
	    ((AliAODConversionPhoton*)fGoodGammas->At(gamma))->SetPz(fUnsmearedPz[gamma]);
	    ((AliAODConversionPhoton*)fGoodGammas->At(gamma))->SetE(fUnsmearedE[gamma]);
	}
	delete[] fUnsmearedPx; fUnsmearedPx = 0x0;
	delete[] fUnsmearedPy; fUnsmearedPy = 0x0;
	delete[] fUnsmearedPz; fUnsmearedPz = 0x0;
	delete[] fUnsmearedE;  fUnsmearedE  = 0x0;
    }

    return kTRUE;
}

//________________________________________________________________________
AliAODConversionMother* AliConversionSelection::GetPi0(Int_t index){

    if(index>=0&&index<GetNumberOfPi0s()){
	return dynamic_cast<AliAODConversionMother*>(fPi0Candidates->At(index));
    }
    return NULL;
}

//________________________________________________________________________
AliAODConversionMother* AliConversionSelection::GetBG(Int_t index){

    if(index>=0&&index<GetNumberOfBGs()){
	return dynamic_cast<AliAODConversionMother*>(fBGPi0s->At(index));
    }
    return NULL;
}

//________________________________________________________________________
AliAODConversionPhoton* AliConversionSelection::GetPhoton(Int_t index){

    if(index>=0&&index<GetNumberOfPhotons()){
	return dynamic_cast<AliAODConversionPhoton*>(fGoodGammas->At(index));
    }
    return NULL;
}

//________________________________________________________________________
void AliConversionSelection::CalculatePi0Candidates(){

    // Conversion Gammas
    if(fGoodGammas->GetEntriesFast()>1){

	for(Int_t firstGammaIndex=0;firstGammaIndex<fGoodGammas->GetEntriesFast()-1;firstGammaIndex++){

	    AliAODConversionPhoton *gamma0=dynamic_cast<AliAODConversionPhoton*>(fGoodGammas->At(firstGammaIndex));
      if (gamma0==NULL) continue;
	    // Combine Photons

	    for(Int_t secondGammaIndex=firstGammaIndex+1;secondGammaIndex<fGoodGammas->GetEntriesFast();secondGammaIndex++){

		AliAODConversionPhoton *gamma1=dynamic_cast<AliAODConversionPhoton*>(fGoodGammas->At(secondGammaIndex));
      if (gamma1==NULL) continue;
		//Check for same Electron ID
		if(gamma0->GetTrackLabelPositive()==gamma1->GetTrackLabelPositive()||gamma0->GetTrackLabelNegative()==gamma1->GetTrackLabelNegative()
		   ||gamma0->GetTrackLabelNegative()==gamma1->GetTrackLabelPositive()||gamma0->GetTrackLabelPositive()==gamma1->GetTrackLabelNegative())continue;

		AliAODConversionMother pi0cand(gamma0,gamma1);
		pi0cand.SetLabels(firstGammaIndex,secondGammaIndex);

		// Set MC Label

		if(fMCEvent){

		    TParticle *mcgam0=gamma0->GetMCParticle(fMCEvent->Stack());
		    TParticle *mcgam1=gamma1->GetMCParticle(fMCEvent->Stack());

		    if(mcgam0&&mcgam1){
			// Have same Mother?

			if(mcgam0->GetMother(0)==mcgam1->GetMother(0)){

			    pi0cand.SetMCLabel(mcgam0->GetMother(0));
			}
		    }
		}

		if((fMesonCut->MesonIsSelected(&pi0cand))){
		    if(MesonInMassWindow(&pi0cand)){

			// Add Pi0 to Stack
			new((*fPi0Candidates)[fPi0Candidates->GetEntriesFast()]) AliAODConversionMother(pi0cand);
		    }
		}
	    }
	}
    }
}

//________________________________________________________________________
Bool_t AliConversionSelection::MesonInMassWindow(AliAODConversionMother *pi0cand)
{
    if (pi0cand->M() > fInvMassRange[0] && pi0cand->M() < fInvMassRange[1] ){
	return kTRUE;
    }
    return kFALSE;
}

//________________________________________________________________________
void AliConversionSelection::RotateParticle(AliAODConversionPhoton *gamma,Int_t nDegreesPMBackground){

    if(!fRandomizer){
	fRandomizer=new TRandom3();
	fRandomizer->SetSeed(0);
    }
    Double_t nRadiansPM = nDegreesPMBackground*TMath::Pi()/180;
    Double_t rotationValue = fRandomizer->Rndm()*2*nRadiansPM + TMath::Pi()-nRadiansPM;
    gamma->RotateZ(rotationValue);
}

//________________________________________________________________________

void AliConversionSelection::CalculateBackground(){

    //Rotation Method
    if(fMesonCut->UseRotationMethod()){

	// Correct for the number of rotations
	// BG is for rotation the same, except for factor NRotations
	Double_t weight=1./Double_t(fMesonCut->GetNumberOfBGEvents());

	for(Int_t firstGammaIndex=0;firstGammaIndex<fGoodGammas->GetEntriesFast();firstGammaIndex++){

	    AliAODConversionPhoton *gamma0=dynamic_cast<AliAODConversionPhoton*>(fGoodGammas->At(firstGammaIndex));
       if (gamma0 ==NULL) continue;
	    for(Int_t secondGammaIndex=firstGammaIndex+1;secondGammaIndex<fGoodGammas->GetEntriesFast();secondGammaIndex++){
		AliAODConversionPhoton *gamma1=dynamic_cast<AliAODConversionPhoton*>(fGoodGammas->At(secondGammaIndex));
      if (gamma1==NULL) continue;
		if(!fConversionCut->PhotonIsSelected(gamma1,fInputEvent))continue;
		for(Int_t nRandom=0;nRandom<fMesonCut->GetNumberOfBGEvents();nRandom++){

		    RotateParticle(gamma1,fMesonCut->NDegreesRotation());

		    AliAODConversionMother BGcandidate(gamma0,gamma1);

		    if(fMesonCut->MesonIsSelected(&BGcandidate,kFALSE)){
			if(MesonInMassWindow(&BGcandidate)){

			    new((*fBGPi0s)[fBGPi0s->GetEntriesFast()]) AliAODConversionMother(BGcandidate);
             
			    dynamic_cast<AliAODConversionMother*>(fBGPi0s->At(fBGPi0s->GetEntriesFast()-1))->SetWeight(weight);
			}
		    }
		}
	    }
	}
    }

    else{
	// Do Event Mixing

	if(fBGHandler==NULL){
	    fBGHandler=new AliConversionAODBGHandlerRP(fConversionCut->IsHeavyIon(),fMesonCut->UseTrackMultiplicity());
	}

	for(Int_t nEventsInBG=0;nEventsInBG <fBGHandler->GetNBGEvents(fGoodGammas,fInputEvent);nEventsInBG++){
         
	    AliGammaConversionPhotonVector *previousEventGammas = fBGHandler->GetBGGoodGammas(fGoodGammas,fInputEvent,nEventsInBG);

	    if(previousEventGammas){

		// test weighted background
		Double_t weight=1.0;
		// Correct for the number of eventmixing:
		// N gammas -> (N-1) + (N-2) +(N-3) ...+ (N-(N-1))  using sum formula sum(i)=N*(N-1)/2  -> N*(N-1)/2
		// real combinations (since you cannot combine a photon with its own)
		// but BG leads to N_{a}*N_{b} combinations
		weight*=0.5*(Double_t(fGoodGammas->GetEntriesFast()-1))/Double_t(previousEventGammas->size());

		for(Int_t iCurrent=0;iCurrent<fGoodGammas->GetEntriesFast();iCurrent++){

		    AliAODConversionPhoton *gamma0 = (AliAODConversionPhoton*)(fGoodGammas->At(iCurrent));

		    for(UInt_t iPrevious=0;iPrevious<previousEventGammas->size();iPrevious++){

			AliAODConversionPhoton *gamma1 = (AliAODConversionPhoton*)(previousEventGammas->at(iPrevious));

			AliAODConversionMother BGcandidate(gamma0,gamma1);

			if(fMesonCut->MesonIsSelected(&BGcandidate,kFALSE)){
			    if(MesonInMassWindow(&BGcandidate)){

				new((*fBGPi0s)[fBGPi0s->GetEntriesFast()]) AliAODConversionMother(BGcandidate);
				dynamic_cast<AliAODConversionMother*>(fBGPi0s->At(fBGPi0s->GetEntriesFast()-1))->SetWeight(weight);
			    }
			}
		    }
		}
	    }
	}
    }
}
//________________________________________________________________________
Double_t AliConversionSelection::GetMultiplicity(AliVEvent *inputEvent){

    switch(fConversionCut->GetMultiplicityMethod())
    {
    case 0:
	return Double_t(GetNumberOfPhotons());
    case 1:
	return Double_t(GetNumberOfChargedTracks(inputEvent));
    case 2:
	return GetVZEROMult(inputEvent);
    case 3:
	return GetSPDMult(inputEvent);
    case 9:
        return 1; // if mult is used as a weight, this number can be used to switch off weighting
    default:
	return 0;
    }
}

//________________________________________________________________________
Int_t AliConversionSelection::GetNumberOfChargedTracks(AliVEvent *inputEvent){

    Int_t ntracks = 0;
       
    AliESDEvent * esdEvent = dynamic_cast<AliESDEvent*>(inputEvent);
    if(esdEvent) {
	if(!fESDTrackCuts){
	    fESDTrackCuts= AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE);
	    fESDTrackCuts->SetMaxDCAToVertexZ(2);
            Double_t etamax=fConversionCut->GetEtaCut();
	    fESDTrackCuts->SetEtaRange(-etamax, etamax);
	    fESDTrackCuts->SetPtRange(0.15);
	}
	for(Int_t iTracks = 0; iTracks < inputEvent->GetNumberOfTracks(); iTracks++){
	    AliESDtrack* currentTrack = esdEvent->GetTrack(iTracks);
	    if(!currentTrack) continue;
	    if(fESDTrackCuts->AcceptTrack(currentTrack))ntracks++;
	}
    } else {
	for(Int_t ii=0; ii<inputEvent->GetNumberOfTracks(); ii++) {
	    AliVTrack * track = dynamic_cast<AliVTrack*>(inputEvent->GetTrack(ii));
       if (track==NULL) continue;
	    if(TMath::Abs(track->Eta())>fConversionCut->GetEtaCut())continue;
	    if(track)ntracks++;
	}
    }

    return ntracks;
}

//________________________________________________________________________
Double_t AliConversionSelection::GetVZEROMult(AliVEvent *inputEvent){

    AliVVZERO *vzero=inputEvent->GetVZEROData();
    Double_t multV0A=vzero->GetMTotV0A();
    Double_t multV0C=vzero->GetMTotV0C();
    Double_t mult=multV0A+multV0C;

    return mult;
}

//________________________________________________________________________
Double_t AliConversionSelection::GetSPDMult(AliVEvent *inputEvent){

    AliESDEvent * esdEvent = dynamic_cast<AliESDEvent*>(inputEvent);
    if(esdEvent) {
	const AliMultiplicity *esdmult=esdEvent->GetMultiplicity();
	return esdmult->GetNumberOfITSClusters(1);
    } else {
	// AOD implementation
        AliAODHeader *header=(AliAODHeader*)inputEvent->GetHeader();
	return header->GetNumberOfITSClusters(1);
    }
    return 0;
}

//________________________________________________________________________
Int_t AliConversionSelection::GetEventNumber(AliVEvent *inputEvent){

    AliESDEvent * esdEvent = dynamic_cast<AliESDEvent*>(inputEvent);
    if(esdEvent) {
	return esdEvent->GetEventNumberInFile();
    }
    else{
	AliAODHeader *header=(AliAODHeader*)inputEvent->GetHeader();
	return header->GetEventNumberESDFile();
    }
    return 0;
}

//________________________________________________________________________
TString AliConversionSelection::GetCutString(){
    TString a= Form("%s%s",fConversionCut->GetCutNumber().Data(),fMesonCut->GetCutNumber().Data());
    return a;
}

