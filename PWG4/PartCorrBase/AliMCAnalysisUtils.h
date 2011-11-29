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
#include <TString.h> 
class TList ;
class TVector3;

//--- AliRoot system ---
class AliCaloTrackReader ;
class AliStack ;

class AliMCAnalysisUtils : public TObject {
	
 public: 
  AliMCAnalysisUtils() ; // ctor
  virtual ~AliMCAnalysisUtils() ;//virtual dtor
    
  //--------------------------------------
  //Enum with tag for origin of particles
  //--------------------------------------
  
  //"Mostly" photon parent types on line 1,
  //                                                  
  //then "mostly" electron parent types on line 2, (e.g. electrons can
  //come from pi0 decay)                              
  //then charged particles on line 3,                                                                                    
  //followed by other and unknown on line 4                                                                              
  enum mcTypes { kMCPhoton,     kMCPrompt,      kMCFragmentation, kMCISR,    
                 kMCPi0Decay,   kMCEtaDecay,    kMCOtherDecay,    kMCConversion,
                 kMCElectron,   kMCEFromCFromB, kMCEFromC,        kMCEFromB, kMCZDecay,   kMCWDecay,
                 kMCMuon,       kMCPion,        kMCPi0,           kMCKaon,   kMCEta,      kMCProton,   
                 kMCAntiProton, kMCNeutron,     kMCAntiNeutron,
                 kMCOther,      kMCUnknown,     kMCBadLabel                                                                                         } ;
  
  //--------------------------------------
  // Methods to check origin of clusters
  //--------------------------------------
  
  Int_t   CheckCommonAncestor(const Int_t index1, const Int_t index2, const AliCaloTrackReader* reader, 
			      Int_t & ancPDG, Int_t & ancStatus, TLorentzVector & momentum, TVector3 & v) ;
  Int_t   CheckOrigin(const Int_t label, const AliCaloTrackReader * reader, const Int_t input) ;
  //Check the label of the most significant particle but do checks on the rest of the contributing labels
  Int_t   CheckOrigin(const Int_t *label, const Int_t nlabels, const AliCaloTrackReader * reader, const Int_t input) ;
  
  Int_t   CheckOriginInStack(const Int_t *labels, const Int_t nlabels, AliStack * stack) ;
  Int_t   CheckOriginInAOD  (const Int_t *labels, const Int_t nlabels, const TClonesArray* mcparticles) ;
  
  void    CheckOverlapped2GammaDecay(const Int_t *labels, const Int_t nlabels, const Int_t mesonIndex, AliStack * stack, Int_t & tag);
  void    CheckOverlapped2GammaDecay(const Int_t *labels, const Int_t nlabels, const Int_t mesonIndex, const TClonesArray* mcparticles, Int_t & tag);
  
  //Check or set the bits produced in the above methods
  void    SetTagBit(Int_t &tag, const UInt_t set) const {
    // Set bit of type set (mcTypes) in tag
    tag |= (1<<set) ; 
  } 
  
  Bool_t  CheckTagBit(const Int_t tag, const UInt_t test) const {
    // Check if in tag the bit test (mcTypes) is set.
    if (tag & (1<<test) ) return  kTRUE ;    
    else return kFALSE ;
  }
  
  //--------------------------------------
  // Other methods
  //--------------------------------------
  
  // Method to recover MC jets stored in generator
  TList * GetJets(const AliCaloTrackReader * reader) ;
  
  void    SetDebug(Int_t deb)           { fDebug=deb           ; }
  Int_t   GetDebug()              const { return fDebug        ; }	
  
  void    SetMCGenerator(TString mcgen) { fMCGenerator = mcgen ; }
  TString GetMCGenerator()        const { return fMCGenerator  ; }	
  
  void    Print(const Option_t * opt) const;
  
 private:
  Int_t   fCurrentEvent;        // Current Event
  Int_t   fDebug;               // Debug level
  TList * fJetsList;            // List of jets
  TString fMCGenerator;         // MC geneator used to generate data in simulation
  
  AliMCAnalysisUtils & operator = (const AliMCAnalysisUtils & ) ; // cpy assignment
  AliMCAnalysisUtils(const AliMCAnalysisUtils & mcu) ;            // cpy ctor
  
  ClassDef(AliMCAnalysisUtils,3)

} ;

#endif //ALIMCANALYSISUTILS_H



