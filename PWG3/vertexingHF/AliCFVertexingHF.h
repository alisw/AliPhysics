#ifndef ALICFVERTEXINGHF_H
#define ALICFVERTEXINGHF_H


/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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

//-----------------------------------------------------------------------
// Class for HF corrections as a function of many variables and step 
// Author : C. Zampolli, CERN
// D. Caffarri, Univ & INFN Padova caffarri@pd.infn.it
// Base class for HF Unfolding - agrelli@uu.nl
//-----------------------------------------------------------------------

#include "AliCFContainer.h"
#include "AliAODRecoDecayHF.h"

class TH1I;
class TParticle ;
class TFile ;
class TClonesArray ;
class AliAODMCParticle;
class AliAODMCHeader;
class AliAODEvent;
class THnSparse;
class TClonesArray;
class AliESDtrackCuts;



class AliCFVertexingHF : public TObject {
	public:

        enum DecayChannel{kD0toKpi = 2, kDStartoKpipi = 21, kDplustoKpipi = 31, kLctopKpi = 32, kDstoKKpi = 33, kD0toKpipipi = 4};

	AliCFVertexingHF() ;
	AliCFVertexingHF(TClonesArray *mcArray, UShort_t originDselection);
	AliCFVertexingHF(const AliCFVertexingHF& c);
	AliCFVertexingHF& operator= (const AliCFVertexingHF& c);

	virtual ~AliCFVertexingHF();
		
	virtual Bool_t GetGeneratedValuesFromMCParticle(Double_t* /*vectorMC*/) {return kFALSE;} 
	virtual Bool_t GetRecoValuesFromCandidate(Double_t* /*vectorReco*/) const {return kFALSE;}
	virtual Bool_t CheckMCChannelDecay() const {return kFALSE;}
	virtual Bool_t SetRecoCandidateParam(AliAODRecoDecayHF* /*recoCand*/) {return kFALSE;}

        virtual void SetDecay3Prong(Int_t /*decay*/){};
	
	void   SetFillFromGenerated(Bool_t flag) {fFillFromGenerated = flag;}
	Bool_t GetFillFromGenerated() const {return fFillFromGenerated;}
		
	void  SetNVar(Int_t nVar) {fNVar = nVar;}  

	void  SetRecoPrimVertex (Double_t zPrimVertex) {fzPrimVertex = zPrimVertex;}
	void  SetMCPrimaryVertex (Double_t zMCVertex){fzMCVertex = zMCVertex;}
	void  SetMCLabel (Int_t mcLabel) {fmcLabel = mcLabel;}
	Int_t GetMCLabel () const {return  fmcLabel;}
		
	void   SetMCCandidateParam(Int_t label);

	Int_t  MCcquarkCounting(AliAODMCParticle* mcPart) const; 
	Bool_t CheckMCPartFamily(AliAODMCParticle */*mcPart*/, TClonesArray */*mcArray*/) const;
	//	Int_t  CheckOrigin(AliAODMCParticle* mcPart) const;
	Int_t  CheckOrigin() const;
  	Bool_t CheckMCDaughters() const;
	Bool_t FillMCContainer(Double_t *containerInputMC); 
	Bool_t FillRecoContainer(Double_t *containerInput); 
	Bool_t MCAcceptanceStep() const;
	Bool_t MCRefitStep(AliAODEvent *aodEvent, AliESDtrackCuts *trackCuts) const;
	Bool_t RecoStep();

	Double_t GetEtaProng(Int_t iProng) const;
	Double_t GetPtProng(Int_t iProng) const;

	Double_t GetPtCand() const {return fRecoCandidate->Pt();}
	Double_t GetYCand() const {return fRecoCandidate->Y();}

	Bool_t RecoAcceptStep(AliESDtrackCuts *trackCuts) const;
	
	Bool_t FillUnfoldingMatrix(Double_t fill[4]) const;
	
	void SetNProngs(Int_t nProngs){fProngs = nProngs;}
	void SetDselection(UShort_t originDselection); 
	UShort_t GetDselection() {return fOriginDselection;}; 
	Int_t CheckReflexion(Char_t isSign);
	Bool_t SetLabelArray();

	void SetCentralityValue(Float_t centValue) {fCentValue = centValue;}

	protected:
	
	TClonesArray      *fmcArray;               //mcArray candidate
	AliAODRecoDecayHF *fRecoCandidate;         // Reconstructed HF candidate 
	AliAODMCParticle  *fmcPartCandidate;
	Int_t fNDaughters;
	Int_t fNVar;                // get Number of variables for the container from the channel decay
	Double_t fzPrimVertex;       //Reco z primary vertex	
	Double_t fzMCVertex;         //MC z primary vertex
	
	Bool_t fFillFromGenerated;   //  flag to indicate whether data container should be filled  
	UShort_t fOriginDselection;      // flag to select D0 origins. 0 Only from charm 1 only from beauty 2 both from charm and beauty
	
	Bool_t fKeepDfromB;           //flag for the feed down from b quark decay. 			
	Bool_t fKeepDfromBOnly;       // flag to keep only the charm particles that comes from beauty decays
	Int_t fmcLabel;              // results of the MatchToMC()
	Int_t fProngs;               // n. of prongs	
	Int_t* fLabelArray;          //[fProngs]  array of labels

	Float_t fCentValue;         // centrality value

	ClassDef(AliCFVertexingHF, 3);
	
};

#endif
