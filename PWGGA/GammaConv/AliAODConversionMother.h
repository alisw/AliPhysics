#ifndef ALIAODCONVERSIONMOTHER_H
#define ALIAODCONVERSIONMOTHER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

////////////////////////////////////////////////
//--------------------------------------------- 
// Class reconstructing the mother particle of conversion gammas
//---------------------------------------------
////////////////////////////////////////////////

//Author Daniel Lohner (Daniel.Lohner@cern.ch)

#include "TLorentzVector.h"
#include "AliAODConversionParticle.h"
#include "AliAODConversionPhoton.h"
#include "AliKFConversionMother.h"
#include "AliKFParticle.h"
#include "AliV0ParticleStrange.h"

class AliAODConversionMother : public AliAODConversionParticle{

	public: 

		//Default Constructor
		AliAODConversionMother();

		// Constructor for ESD to AOD Conversion
		AliAODConversionMother(const AliKFConversionMother *kf);

		//Constructor Decay Mother Particle
		AliAODConversionMother(const AliAODConversionPhoton *y1, const AliAODConversionPhoton *y2);
		// Constructor Mother particle from one photon and one meson
		AliAODConversionMother(const AliAODConversionMother *meson, const AliAODConversionPhoton *gamma);
		// Constructor Mother particle from two mesons
		AliAODConversionMother(const AliAODConversionMother *meson1, const AliAODConversionMother *meson2);
    // Contructor Mother particle from V0 and photon
    AliAODConversionMother(const AliV0ParticleStrange *v0, const AliAODConversionPhoton *gamma);

		
		//Destructor
		virtual ~AliAODConversionMother();

		// MC

		void SetMCLabel(Int_t i){fMCLabel=i;}
		Int_t GetMCLabel() const {return fMCLabel;}
        TParticle *GetMCParticle(AliMCEvent *mcEvent);
        Bool_t IsTrueMeson(AliMCEvent *mcEvent,Int_t pdgcode);

		///Set the Chi2 of reconstructed conversion gamma
		void SetChi2(Float_t chi2) {fChi2 = chi2;}

		//Get the Chi2 of particle
		Float_t Chi2() const {return fChi2;}

		///Set track or MC labels
		void SetLabel1(Int_t label){fLabel[0] = label;}
		void SetLabel2(Int_t label){fLabel[1] = label;}
		void SetLabel3(Int_t label){fLabel[2] = label;}
		void SetLabels(Int_t label1, Int_t label2, Int_t label3 = 0){fLabel[0] = label1; fLabel[1] = label2; fLabel[2] = label3;}

		Int_t GetLabel(Int_t i) const {return fLabel[i];}
		Int_t GetLabel1() const {return fLabel[0];}
		Int_t GetLabel2() const {return fLabel[1];}
		Int_t GetLabel3() const {return fLabel[2];}
			
		Double_t GetProductionRadius() const {return TMath::Sqrt(fProductionVtx[0]*fProductionVtx[0]+fProductionVtx[1]*fProductionVtx[1]);}
		Double_t GetProductionX() const {return fProductionVtx[0];}
		Double_t GetProductionY() const {return fProductionVtx[1];}
		Double_t GetProductionZ() const {return fProductionVtx[2];}

		void SetProductionX(Double_t x) {fProductionVtx[0]=x;}
		void SetProductionY(Double_t y) {fProductionVtx[1]=y;}
		void SetProductionZ(Double_t z) {fProductionVtx[2]=z;}
		void SetProductionPoint(Double_t* point){
			fProductionVtx[0] = point[0];
			fProductionVtx[1] = point[1];
			fProductionVtx[2] = point[2];
		}
			
		Float_t GetDCABetweenPhotons() const {return fdcaBetweenPhotons;}
		Float_t GetDCAZMotherPrimVtx() const {return fdcaZPrimVtx;}
		Float_t GetDCARMotherPrimVtx() const {return fdcaRPrimVtx;}
		UChar_t GetMesonQuality() const {return fQuality;}
			
		Double_t GetOpeningAngle() const { return fOpeningAngle;}

		Double_t GetAlpha() const { return fAlpha;}

		void SetWeight(Double_t weight) {fWeight=weight;}
		Double_t GetWeight() const {return fWeight;}

		Float_t CalculateDistanceBetweenPhotons(const AliAODConversionPhoton* y1, const AliAODConversionPhoton* y2 , Double_t prodPoint[3]);
		void CalculateDistanceOfClossetApproachToPrimVtx(const AliVVertex* primVertex);
		void DetermineMesonQuality(const AliAODConversionPhoton* y1, const AliAODConversionPhoton* y2);
		
		void SetTrueMesonValue(Int_t trueMeson) {fTrueMeson = trueMeson;}
		Int_t GetTrueMesonValue()const {return fTrueMeson;}
		
		
	private:
		Int_t fLabel[3]; 						// Labels of the decay photons
		Int_t fMCLabel; 						// MC Label
		Float_t fChi2; 							// Chi sq of reconstructed mother
		Double_t fOpeningAngle;
		Double_t fAlpha;
		Double_t fWeight; 						// Weight for BG Calculation
		Float_t fdcaBetweenPhotons; 			// dca between the two photons
		Double_t fProductionVtx[3]; 			// Production vertex
		Float_t fdcaZPrimVtx; 					// dca Z of meson to primary vertex
		Float_t fdcaRPrimVtx; 					// dca R of meson to primary vertex
		UChar_t fQuality; 						// Quality of the meson:
													// 0: garbage
													// 1: both photons quality 1
													// 2: 1 photon quality 1, 1 photon quality 2
													// 3: 1 photon quality 1, 1 photon quality 3
													// 4: both photons quality 2
													// 5: 1 photon quality 2, 1 photon quality 3
													// 6: both photons quality 3
		
		Int_t fTrueMeson;						// is true meson
													// 0 : no
													// 1 : pi0
													// 2 : eta
													// 3 : eta'
													// 4 : omega
		
	ClassDef(AliAODConversionMother,5)
};

#endif
