#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliAODConversionPhoton.h"
#include "AliKFConversionPhoton.h"

using namespace std;

ClassImp(AliAODConversionPhoton)

AliAODConversionPhoton::AliAODConversionPhoton() :
AliAODConversionParticle(),
AliConversionPhotonBase(),
fDCArPrimVtx(0),
fDCAzPrimVtx(0),
fCaloPhoton(0),
fNCaloPhotonMCLabels(0),
fCaloPhotonMCFlags(0)
{
	// initialize calo photon MC labels
	for (Int_t i =0; i<50; i++){
		fCaloPhotonMCLabels[i]=-1;		
	}
  //Standard constructor
}

AliAODConversionPhoton::AliAODConversionPhoton(AliKFConversionPhoton *kfphoton) :
AliAODConversionParticle(kfphoton),
AliConversionPhotonBase(*((AliConversionPhotonBase*)kfphoton)),
fDCArPrimVtx(0),
fDCAzPrimVtx(0),
fCaloPhoton(0),
fNCaloPhotonMCLabels(0),
fCaloPhotonMCFlags(0)
{
    //Constructor from kfphoton

    // puts the mass to zero and store dilepton mass
	    SetMass(kfphoton->M());

	//SetE(P());
		
	// initialize calo photon MC labels
	for (Int_t i =0; i<50; i++){
		fCaloPhotonMCLabels[i]=-1;		
	}
	
}

AliAODConversionPhoton::AliAODConversionPhoton(TLorentzVector *vec) :
AliAODConversionParticle(vec),
AliConversionPhotonBase(),
fDCArPrimVtx(0),
fDCAzPrimVtx(0),
fCaloPhoton(0),
fNCaloPhotonMCLabels(0),
fCaloPhotonMCFlags(0)
{
    //Constructor from TLorentzVector

	// initialize calo photon MC labels
	for (Int_t i =0; i<50; i++){
		fCaloPhotonMCLabels[i]=-1;		
	}
}



AliAODConversionPhoton::AliAODConversionPhoton(const AliAODConversionPhoton & original) :
AliAODConversionParticle(original),
AliConversionPhotonBase(original),
fDCArPrimVtx(original.fDCArPrimVtx),
fDCAzPrimVtx(original.fDCAzPrimVtx),
fCaloPhoton(original.fCaloPhoton),
fNCaloPhotonMCLabels(original.fNCaloPhotonMCLabels),
fCaloPhotonMCFlags(original.fCaloPhotonMCFlags)
{
	//Copy constructor
	
	// initialize calo photon MC labels
	for (Int_t i =0; i<50; i++){
		fCaloPhotonMCLabels[i]=original.fCaloPhotonMCLabels[i];		
	}
}

AliAODConversionPhoton::~AliAODConversionPhoton()
{
  // empty standard destructor
}

AliAODConversionPhoton & AliAODConversionPhoton::operator = (const AliAODConversionPhoton & /*source*/)
{
  // assignment operator
  return *this;
}

///________________________________________________________________________
void AliAODConversionPhoton::CalculateDistanceOfClossetApproachToPrimVtx(const AliVVertex* primVertex ){

   Double_t primCo[3] = {primVertex->GetX(),primVertex->GetY(),primVertex->GetZ()};

   Double_t absoluteP = TMath::Sqrt(TMath::Power(GetPx(),2) + TMath::Power(GetPy(),2) + TMath::Power(GetPz(),2));   
   Double_t p[3] = {GetPx()/absoluteP,GetPy()/absoluteP,GetPz()/absoluteP};
   Double_t CP[3];
   
   CP[0] =  fConversionPoint[0] - primCo[0];
   CP[1] =  fConversionPoint[1] - primCo[1];
   CP[2] =  fConversionPoint[2] - primCo[2];
   
   Double_t Lambda = - (CP[0]*p[0]+CP[1]*p[1]+CP[2]*p[2])/(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
   
   Double_t S[3];
   S[0] = fConversionPoint[0] + p[0]*Lambda;
   S[1] = fConversionPoint[1] + p[1]*Lambda;
   S[2] = fConversionPoint[2] + p[2]*Lambda;
  
   fDCArPrimVtx = TMath::Sqrt( TMath::Power(primCo[0]-S[0],2) + TMath::Power(primCo[1]-S[1],2));
   fDCAzPrimVtx = primCo[2]-S[2];
    
   return;
}


void AliAODConversionPhoton::SetCaloPhotonMCFlags(AliStack *MCStack){
	
	Bool_t isPhoton = kFALSE;						// largest contribution to cluster is photon 
	Bool_t isElectron = kFALSE;						// largest contribution to cluster is electron
	Bool_t isConversion = kFALSE;					// largest contribution to cluster is converted electron
	Bool_t isConversionFullyContained = kFALSE;		// largest contribution to cluster is converted electron, second electron has been found in same cluster
	Bool_t isMerged = kFALSE;						// largest contribution to cluster is photon, second photon or electron from dalitz has been found in same cluster 
	Bool_t isMergedPartConv = kFALSE;				// cluster contains more than one particle from the same decay and at least one of the particles came from a conversion
	Bool_t isDalitz = kFALSE;						// this cluster was created by a particle stemming from a dality decay
	Bool_t isDalitzMerged = kFALSE;					// this cluster was created by a particle stemming from a dality decay and more than one particle of the dalitz decay is contained in the cluster
	Bool_t isPhotonWithElecMother = kFALSE;			// this cluster is from a photon with an electron as mother
	Bool_t isShower = kFALSE;						// this cluster contains as a largest contribution a particle from a shower or radiative process
	
	
	TParticle* Photon;
	if (fNCaloPhotonMCLabels==0) return;
	Photon = MCStack->Particle(GetCaloPhotonMCLabel(0));
	
	if(Photon == NULL){
		return;
	}


	Int_t particleMotherLabel = Photon->GetMother(0);
	Int_t particleGrandMotherLabel = -1; 
	Int_t particleMotherPDG = -1; 
	Int_t particleGrandMotherPDG = -1; 
	Int_t particleMotherNDaugthers = 0;
	Int_t particleGrandMotherNDaugthers = 0;
	if (particleMotherLabel > -1){
		particleMotherNDaugthers = MCStack->Particle(Photon->GetMother(0))->GetNDaughters();
		particleGrandMotherLabel = MCStack->Particle(Photon->GetMother(0))->GetMother(0);
		particleMotherPDG = MCStack->Particle(Photon->GetMother(0))->GetPdgCode();
		if (particleGrandMotherLabel > -1){
			particleGrandMotherPDG = MCStack->Particle(MCStack->Particle(Photon->GetMother(0))->GetMother(0))->GetPdgCode();
			particleGrandMotherNDaugthers = MCStack->Particle(MCStack->Particle(Photon->GetMother(0))->GetMother(0))->GetNDaughters();
		}	
	}
		
	// Check wether the first contribution was photon	
	if(MCStack->Particle(GetCaloPhotonMCLabel(0))->GetPdgCode() == 22){
		isPhoton=kTRUE;
		// did it decay via the dalitz channel
		if (particleMotherLabel > -1 && particleMotherNDaugthers == 3) isDalitz = kTRUE;
		// Test wether particle stems from a shower or radiation
		if (abs(particleMotherPDG) == 11){									// check whether photon stems from electron
			isPhotonWithElecMother = kTRUE;
			if (particleGrandMotherLabel > -1){	 								// test whether first particle has a grandmother
				if (abs(particleGrandMotherPDG) == 22 ) isShower = kTRUE;	// check whether grandmother is a photon (meaning this is most likely a shower)
			}	
		}			
	}
	// Check wether the first contribution was electron
	if( abs(MCStack->Particle(GetCaloPhotonMCLabel(0))->GetPdgCode()) == 11 ){
		isElectron=kTRUE;	
		if (particleMotherLabel > -1) {
			// was it a conversion
			if (abs(particleMotherPDG) == 22) isConversion = kTRUE;
			// did it decay via the dalitz channel
			if (particleGrandMotherLabel > -1 && particleGrandMotherNDaugthers == 3 ) isDalitz = kTRUE;
		}
		if (particleGrandMotherLabel > -1){										// check whether electron has a grandmother
			if (abs(particleGrandMotherPDG) == 11 ||  abs(particleGrandMotherPDG) == 22){	// test wether electron has photon or electron as grandmother (meaning will most likely be a shower)
				isShower = kTRUE; 
			}	
		}
	}
	
	// check wether there were other contributions to the cluster
	if (fNCaloPhotonMCLabels>1){
		// largest contribution was from photon and is not from shower or electron mother
		if (isPhoton && (!isShower || !isPhotonWithElecMother )){
// 			cout << "largest contribution is photon, with mother: " <<  particleMotherLabel <<"(" << particleMotherPDG<< ")" << "\t with " << particleMotherNDaugthers << " daugthers and grand mother: " << particleGrandMotherLabel  <<"(" << particleGrandMotherPDG << ")" << endl;
			
			TParticle* dummyPart =NULL;
			for (Int_t i = 1; i< fNCaloPhotonMCLabels; i++){
				if (i > 49) continue;													// abort if more than 20 entries to the cluster have been checked (more are not stored in these objects)
				dummyPart = MCStack->Particle(GetCaloPhotonMCLabel(i));
				Int_t dummyPartMotherLabel = dummyPart->GetMother(0);
				Int_t dummyPartGrandMotherLabel = -1;
				Int_t dummyPartMotherPDG = -1; 
				Int_t dummyPartGrandMotherPDG = -1; 
					
				// check wether this particle has a mother & obtain the pdg code
				if (dummyPartMotherLabel > -1){
					dummyPartGrandMotherLabel = MCStack->Particle(dummyPart->GetMother(0))->GetMother(0);
					dummyPartMotherPDG = MCStack->Particle(dummyPart->GetMother(0))->GetPdgCode(); 
					// check wether this particle has a grandmother & obtain its pdg code
					if (dummyPartGrandMotherLabel > -1){
						dummyPartGrandMotherPDG = MCStack->Particle(MCStack->Particle(dummyPart->GetMother(0))->GetMother(0))->GetPdgCode();
					}
				}
				
				if (particleMotherLabel > -1){  												// test whether first particle has a mother
					if (dummyPartMotherLabel == particleMotherLabel) isMerged = kTRUE;			// test whether current and first particle have the same mother => i.e. other gamma from decay or dalitz electron
					if (dummyPartGrandMotherLabel > -1){	 									// test whether first particle has a grandmother
						// check wether particle is an electron from a conversion of a photon from the original mother
						if (abs(dummyPart->GetPdgCode()) == 11 && dummyPartGrandMotherLabel == particleMotherLabel ) isMergedPartConv = kTRUE; 	
						// check wether particle is an electron from a dalitz decay from the original mother
						if (abs(dummyPart->GetPdgCode()) == 11 && dummyPartMotherLabel == particleMotherLabel ) isDalitzMerged = kTRUE;
					}	
				}	
// 				cout << i << "\t" <<GetCaloPhotonMCLabel(i) << "\t" << dummyPart->GetPdgCode() << "\t mother is: " <<  dummyPartMotherLabel <<"(" << dummyPartMotherPDG<< ")"<< "\t grandmother is: " <<  dummyPartGrandMotherLabel<<"(" << dummyPartGrandMotherPDG<< ")"<<  endl;
				
			}	
		}
		// largest contribution was from electron & not a from a shower
		if (isElectron && !isShower){
// 			cout << "largest contribution is electron, with mother: " <<  particleMotherLabel <<"(" << particleMotherPDG<< ")" << "\t and grand mother: " << particleGrandMotherLabel  <<"(" << particleGrandMotherPDG << ")" << endl;
			TParticle* dummyPart =NULL;
			for (Int_t i = 1; i< fNCaloPhotonMCLabels; i++){
				if (i > 49) continue;
				
				dummyPart = MCStack->Particle(GetCaloPhotonMCLabel(i));
				Int_t dummyPartMotherLabel = dummyPart->GetMother(0);
				Int_t dummyPartGrandMotherLabel = -1;
				Int_t dummyPartMotherPDG = -1; 
				Int_t dummyPartGrandMotherPDG = -1; 
				
				// check wether this particle has a mother & obtain the pdg code
				if (dummyPartMotherLabel > -1){
					dummyPartGrandMotherLabel = MCStack->Particle(dummyPart->GetMother(0))->GetMother(0);
					dummyPartMotherPDG = MCStack->Particle(dummyPart->GetMother(0))->GetPdgCode();
					// check wether this particle has a grandmother & obtain its pdg code
					if (dummyPartGrandMotherLabel > -1){
						dummyPartGrandMotherPDG = MCStack->Particle(MCStack->Particle(dummyPart->GetMother(0))->GetMother(0))->GetPdgCode();
					}
				}
				
				if (particleMotherLabel > -1){																			// test whether first particle has a mother
					if (isConversion && dummyPartMotherLabel == particleMotherLabel) isConversionFullyContained = kTRUE;		// test whether conversion is fully contained in cluster 

					if (particleGrandMotherLabel > -1){																	// test whether first particle has a grandmother															
						if (abs(dummyPart->GetPdgCode()) == 22){													// test whether this particle is a photon
							// check wether orginal electron and this photon stem from the same particle and electron stems from conversion
							if( dummyPartMotherLabel == particleGrandMotherLabel && (abs(dummyPartMotherPDG) != 11 ||  abs(dummyPartMotherPDG) != 22 )   ) isMergedPartConv = kTRUE;
							// check wether orginal electron and this photon stem from the same particle and electron originated in dalitz
							if( dummyPartMotherLabel == particleMotherLabel && (abs(dummyPartMotherPDG) != 11 ||  abs(dummyPartMotherPDG) != 22 )   ) isDalitzMerged = kTRUE;
						}
						if (abs(dummyPart->GetPdgCode()) == 11) {
							// check wether orginal electron and this electron stem from the same particle and electron stems from conversion
							if( dummyPartGrandMotherLabel == particleGrandMotherLabel &&  (abs(dummyPartGrandMotherPDG) != 11 ||  abs(dummyPartGrandMotherPDG) != 22 )   ) isMergedPartConv = kTRUE;
							// check wether orginal electron and this electron stem from the same particle and electron originated in dalitz decay
							if( dummyPartMotherLabel == particleMotherLabel && abs(particleMotherPDG) != 22 &&  (abs(dummyPartGrandMotherPDG) != 11 ||  abs(dummyPartGrandMotherPDG) != 22 )   ) isDalitzMerged = kTRUE;
						}
						
					}	
				}	
// 				cout << i << "\t" << GetCaloPhotonMCLabel(i) << "\t" << dummyPart->GetPdgCode() << "\t mother is: " <<  dummyPartMotherLabel <<"(" << dummyPartMotherPDG<< ")"<< "\t grandmother is: " <<  dummyPartGrandMotherLabel<<"(" << dummyPartGrandMotherPDG<< ")"<<  endl;
				
			}	
		}
	}		
	fCaloPhotonMCFlags = isPhoton *1 + isElectron *2 + isConversion*4+ isConversionFullyContained *8 + isMerged *16 + isMergedPartConv*32 + isDalitz *64 + isDalitzMerged *128 + isPhotonWithElecMother *256 + isShower * 512;		
// 	cout << "isPhoton: \t" << isPhoton << "\t isElectron: \t" << isElectron << "\t isConversion: \t" << isConversion <<  "\t isConversionFullyContained: \t" << isConversionFullyContained << "\t isMerged: \t" << isMerged << "\t isMergedPartConv: \t" << isMergedPartConv << "\t isPhotonWithElecMother: \t" << isPhotonWithElecMother<< "\t isDalitz: \t" << isDalitz << "\t isDalitzMerged: \t" << isDalitzMerged << "\t isShower: \t" << isShower << "\t"<< fCaloPhotonMCFlags<< endl;		
	
}

void AliAODConversionPhoton::SetCaloPhotonMCFlagsAOD(AliVEvent* event){
	
	TClonesArray *AODMCTrackArray = dynamic_cast<TClonesArray*>(event->FindListObject(AliAODMCParticle::StdBranchName()));
	if (!AODMCTrackArray) return;
	
	
	Bool_t isPhoton = kFALSE;						// largest contribution to cluster is photon 
	Bool_t isElectron = kFALSE;						// largest contribution to cluster is electron
	Bool_t isConversion = kFALSE;					// largest contribution to cluster is converted electron
	Bool_t isConversionFullyContained = kFALSE;		// largest contribution to cluster is converted electron, second electron has been found in same cluster
	Bool_t isMerged = kFALSE;						// largest contribution to cluster is photon, second photon or electron from dalitz has been found in same cluster 
	Bool_t isMergedPartConv = kFALSE;				// cluster contains more than one particle from the same decay and at least one of the particles came from a conversion
	Bool_t isDalitz = kFALSE;						// this cluster was created by a particle stemming from a dality decay
	Bool_t isDalitzMerged = kFALSE;					// this cluster was created by a particle stemming from a dality decay and more than one particle of the dalitz decay is contained in the cluster
	Bool_t isPhotonWithElecMother = kFALSE;			// this cluster is from a photon with an electron as mother
	Bool_t isShower = kFALSE;						// this cluster contains as a largest contribution a particle from a shower or radiative process
	
	
	AliAODMCParticle* Photon;
	if (fNCaloPhotonMCLabels==0) return;
	Photon = (AliAODMCParticle*) AODMCTrackArray->At(GetCaloPhotonMCLabel(0));
	
	if(Photon == NULL){
		return;
	}

	AliAODMCParticle* PhotonMother;
	AliAODMCParticle* PhotonGrandMother;
	Int_t particleMotherLabel = Photon->GetMother();
	Int_t particleGrandMotherLabel = -1; 
	Int_t particleMotherPDG = -1; 
	Int_t particleGrandMotherPDG = -1; 
	Int_t particleMotherNDaugthers = 0;
	Int_t particleGrandMotherNDaugthers = 0;
	if (particleMotherLabel > -1){
		PhotonMother = (AliAODMCParticle*) AODMCTrackArray->At(Photon->GetMother());
		particleMotherNDaugthers = PhotonMother->GetNDaughters();
		particleGrandMotherLabel = PhotonMother->GetMother();
		particleMotherPDG = PhotonMother->GetPdgCode();
		if (particleGrandMotherLabel > -1){
			PhotonGrandMother = (AliAODMCParticle*) AODMCTrackArray->At(PhotonMother->GetMother());
			particleGrandMotherPDG = PhotonGrandMother->GetPdgCode();
			particleGrandMotherNDaugthers = PhotonGrandMother->GetNDaughters();
		}	
	}
		
	// Check wether the first contribution was photon	
	if(abs(Photon->GetPdgCode()) == 22){
		isPhoton=kTRUE;
		// did it decay via the dalitz channel
		if (particleMotherLabel > -1 && particleMotherNDaugthers == 3) isDalitz = kTRUE;
		// Test wether particle stems from a shower or radiation
		if (abs(particleMotherPDG) == 11){									// check whether photon stems from electron
			isPhotonWithElecMother = kTRUE;
			if (particleGrandMotherLabel > -1){	 								// test whether first particle has a grandmother
				if (abs(particleGrandMotherPDG) == 22 ) isShower = kTRUE;	// check whether grandmother is a photon (meaning this is most likely a shower)
			}	
		}			
	}
	// Check wether the first contribution was electron
	if(abs(Photon->GetPdgCode()) == 11 ){
		isElectron=kTRUE;	
		if (particleMotherLabel > -1) {
			// was it a conversion
			if (abs(particleMotherPDG) == 22) isConversion = kTRUE;
			// did it decay via the dalitz channel
			if (particleGrandMotherLabel > -1 && particleGrandMotherNDaugthers == 3 ) isDalitz = kTRUE;
		}
		if (particleGrandMotherLabel > -1){										// check whether electron has a grandmother
			if (abs(particleGrandMotherPDG) == 11 ||  abs(particleGrandMotherPDG) == 22){	// test wether electron has photon or electron as grandmother (meaning will most likely be a shower)
				isShower = kTRUE; 
			}	
		}
	}
	
	// check wether there were other contributions to the cluster
	if (fNCaloPhotonMCLabels>1){
		// largest contribution was from photon and is not from shower or electron mother
		if (isPhoton && (!isShower || !isPhotonWithElecMother )){
// 			cout << "largest contribution is photon, with mother: " <<  particleMotherLabel <<"(" << particleMotherPDG<< ")" << "\t with " << particleMotherNDaugthers << " daugthers and grand mother: " << particleGrandMotherLabel  <<"(" << particleGrandMotherPDG << ")" << endl;
			
			AliAODMCParticle* dummyPart =NULL;
			AliAODMCParticle* dummyPartMother = NULL;
			AliAODMCParticle* dummyPartGrandMother = NULL;

			for (Int_t i = 1; i< fNCaloPhotonMCLabels; i++){
				if (i > 49) continue;													// abort if more than 20 entries to the cluster have been checked (more are not stored in these objects)
				dummyPart = (AliAODMCParticle*) AODMCTrackArray->At(GetCaloPhotonMCLabel(i));
				Int_t dummyPartMotherLabel = dummyPart->GetMother();
				Int_t dummyPartGrandMotherLabel = -1;
				Int_t dummyPartMotherPDG = -1; 
				Int_t dummyPartGrandMotherPDG = -1; 
					
				// check wether this particle has a mother & obtain the pdg code
				if (dummyPartMotherLabel > -1){
					dummyPartMother = (AliAODMCParticle*) AODMCTrackArray->At(dummyPart->GetMother());
					dummyPartGrandMotherLabel = dummyPartMother->GetMother();
					dummyPartMotherPDG = dummyPartMother->GetPdgCode(); 
					// check wether this particle has a grandmother & obtain its pdg code
					if (dummyPartGrandMotherLabel > -1){
						dummyPartGrandMother = (AliAODMCParticle*) AODMCTrackArray->At(dummyPartMother->GetMother());
						dummyPartGrandMotherPDG = dummyPartGrandMother->GetPdgCode();
					}
				}
				
				if (particleMotherLabel > -1){  												// test whether first particle has a mother
					if (dummyPartMotherLabel == particleMotherLabel) isMerged = kTRUE;			// test whether current and first particle have the same mother => i.e. other gamma from decay or dalitz electron
					if (dummyPartGrandMotherLabel > -1){	 									// test whether first particle has a grandmother
						// check wether particle is an electron from a conversion of a photon from the original mother
						if (abs(dummyPart->GetPdgCode()) == 11 && dummyPartGrandMotherLabel == particleMotherLabel ) isMergedPartConv = kTRUE; 	
						// check wether particle is an electron from a dalitz decay from the original mother
						if (abs(dummyPart->GetPdgCode()) == 11 && dummyPartMotherLabel == particleMotherLabel ) isDalitzMerged = kTRUE;
					}	
				}	
// 				cout << i << "\t" <<GetCaloPhotonMCLabel(i) << "\t" << dummyPart->GetPdgCode() << "\t mother is: " <<  dummyPartMotherLabel <<"(" << dummyPartMotherPDG<< ")"<< "\t grandmother is: " <<  dummyPartGrandMotherLabel<<"(" << dummyPartGrandMotherPDG<< ")"<<  endl;
				
			}	
		}
		// largest contribution was from electron & not a from a shower
		if (isElectron && !isShower){
// 			cout << "largest contribution is electron, with mother: " <<  particleMotherLabel <<"(" << particleMotherPDG<< ")" << "\t and grand mother: " << particleGrandMotherLabel  <<"(" << particleGrandMotherPDG << ")" << endl;
			AliAODMCParticle* dummyPart =NULL;
			AliAODMCParticle* dummyPartMother = NULL;
			AliAODMCParticle* dummyPartGrandMother = NULL;

			for (Int_t i = 1; i< fNCaloPhotonMCLabels; i++){
				if (i > 49) continue;
				
				dummyPart = (AliAODMCParticle*)AODMCTrackArray->At(GetCaloPhotonMCLabel(i));
				Int_t dummyPartMotherLabel = dummyPart->GetMother();
				Int_t dummyPartGrandMotherLabel = -1;
				Int_t dummyPartMotherPDG = -1; 
				Int_t dummyPartGrandMotherPDG = -1; 
				
				if (dummyPartMotherLabel > -1){
					dummyPartMother = (AliAODMCParticle*) AODMCTrackArray->At(dummyPart->GetMother());
					dummyPartGrandMotherLabel = dummyPartMother->GetMother();
					dummyPartMotherPDG = dummyPartMother->GetPdgCode(); 
					// check wether this particle has a grandmother & obtain its pdg code
					if (dummyPartGrandMotherLabel > -1){
						dummyPartGrandMother = (AliAODMCParticle*) AODMCTrackArray->At(dummyPartMother->GetMother());
						dummyPartGrandMotherPDG = dummyPartGrandMother->GetPdgCode();
					}
				}
				
				if (particleMotherLabel > -1){																			// test whether first particle has a mother
					if (isConversion && dummyPartMotherLabel == particleMotherLabel) isConversionFullyContained = kTRUE;		// test whether conversion is fully contained in cluster 

					if (particleGrandMotherLabel > -1){																	// test whether first particle has a grandmother															
						if (abs(dummyPart->GetPdgCode()) == 22){													// test whether this particle is a photon
							// check wether orginal electron and this photon stem from the same particle and electron stems from conversion
							if( dummyPartMotherLabel == particleGrandMotherLabel && (abs(dummyPartMotherPDG) != 11 ||  abs(dummyPartMotherPDG) != 22 )   ) isMergedPartConv = kTRUE;
							// check wether orginal electron and this photon stem from the same particle and electron originated in dalitz
							if( dummyPartMotherLabel == particleMotherLabel && (abs(dummyPartMotherPDG) != 11 ||  abs(dummyPartMotherPDG) != 22 )   ) isDalitzMerged = kTRUE;
						}
						if (abs(dummyPart->GetPdgCode()) == 11) {
							// check wether orginal electron and this electron stem from the same particle and electron stems from conversion
							if( dummyPartGrandMotherLabel == particleGrandMotherLabel &&  (abs(dummyPartGrandMotherPDG) != 11 ||  abs(dummyPartGrandMotherPDG) != 22 )   ) isMergedPartConv = kTRUE;
							// check wether orginal electron and this electron stem from the same particle and electron originated in dalitz decay
							if( dummyPartMotherLabel == particleMotherLabel && abs(particleMotherPDG) != 22 &&  (abs(dummyPartGrandMotherPDG) != 11 ||  abs(dummyPartGrandMotherPDG) != 22 )   ) isDalitzMerged = kTRUE;
						}
						
					}	
				}	
// 				cout << i << "\t" << GetCaloPhotonMCLabel(i) << "\t" << dummyPart->GetPdgCode() << "\t mother is: " <<  dummyPartMotherLabel <<"(" << dummyPartMotherPDG<< ")"<< "\t grandmother is: " <<  dummyPartGrandMotherLabel<<"(" << dummyPartGrandMotherPDG<< ")"<<  endl;
				
			}	
		}
	}		

	fCaloPhotonMCFlags = isPhoton *1 + isElectron *2 + isConversion*4+ isConversionFullyContained *8 + isMerged *16 + isMergedPartConv*32 + isDalitz *64 + isDalitzMerged *128 + isPhotonWithElecMother *256 + isShower * 512;		
// 	cout << "isPhoton: \t" << isPhoton << "\t isElectron: \t" << isElectron << "\t isConversion: \t" << isConversion <<  "\t isConversionFullyContained: \t" << isConversionFullyContained << "\t isMerged: \t" << isMerged << "\t isMergedPartConv: \t" << isMergedPartConv << "\t isPhotonWithElecMother: \t" << isPhotonWithElecMother<< "\t isDalitz: \t" << isDalitz << "\t isDalitzMerged: \t" << isDalitzMerged << "\t isShower: \t" << isShower << "\t"<< fCaloPhotonMCFlags<< endl;		
	
}