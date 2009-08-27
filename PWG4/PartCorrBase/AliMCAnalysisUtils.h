#ifndef ALIMCANALYSISUTILS_H
#define ALIMCANALYSISUTILS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id:  $ */

//_________________________________________________________________________
// Class for analysis utils for MC data
// stored in stack or event header.
// Contains:
//  - method to check the origin of a given track/cluster
//  - method to obtain the generated jets
//
//*-- Author: Gustavo Conesa (INFN-LNF)

// --- ROOT system ---
#include <TObject.h> 
class TString ;
class TList ;

//--- AliRoot system ---
class AliCaloTrackReader ;
class AliStack ;

class AliMCAnalysisUtils : public TObject {
	
public: 
	
	AliMCAnalysisUtils() ; // ctor
	AliMCAnalysisUtils(const AliMCAnalysisUtils & g) ; // cpy ctor
	AliMCAnalysisUtils & operator = (const AliMCAnalysisUtils & g) ;//cpy assignment
	virtual ~AliMCAnalysisUtils() ;//virtual dtor

       //"Mostly" photon parent types on line 1,
       //                                                  
       //then "mostly" electron parent types on line 2, (e.g. electrons can
       //come from pi0 decay)                              
       //then charged particles on line 3,                                                                                    
       //followed by other and unknown on line 4                                                                              
       enum mcTypes {kMCPhoton, kMCPrompt, kMCFragmentation, kMCISR, kMCPi0Decay, kMCEtaDecay, kMCOtherDecay, kMCConversion,
		     kMCElectron, kMCEFromCFromB, kMCEFromC, kMCEFromB, kMCZDecay, kMCWDecay,
		     kMCMuon, kMCPion, kMCPi0, kMCKaon, kMCEta, kMCProton,
		     kMCOther, kMCUnknown};

	Int_t   CheckOrigin(const Int_t label, AliCaloTrackReader * reader, const Int_t input) ;
	Int_t   CheckOriginInStack(const Int_t label, AliStack * stack) ;
	Int_t   CheckOriginInAOD(const Int_t label, TClonesArray* mcparticles) ;
	
	TList * GetJets(AliCaloTrackReader * reader) ;
	
	void SetTagBit(Int_t &tag, const UInt_t set) {
		// Set bit of type set (mcTypes) in tag
		tag |= (1<<set) ; 
	} 
	
	Bool_t CheckTagBit(const Int_t tag, const UInt_t test) const {
		// Check if in tag the bit test (mcTypes) is set.
		if (tag & (1<<test) ) return  kTRUE ;    
		else return kFALSE ;
   }
  
	
	void Print(const Option_t * opt)const;
  	
	void SetDebug(Int_t deb) {fDebug=deb;}
	Int_t GetDebug() const {return fDebug;}	
	
	void SetMCGenerator(TString mcgen) {fMCGenerator=mcgen;}
	TString GetMCGenerator() const {return fMCGenerator;}	

private:
	Int_t   fCurrentEvent;        // Current Event
	Int_t	fDebug;               // Debug level
	TList * fJetsList;            // List of jets
	TString fMCGenerator;         // MC geneator used to generate data in simulation
	
	ClassDef(AliMCAnalysisUtils,3)
} ;


#endif //ALIMCANALYSISUTILS_H



