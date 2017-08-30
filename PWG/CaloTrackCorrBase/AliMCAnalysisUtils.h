#ifndef ALIMCANALYSISUTILS_H
#define ALIMCANALYSISUTILS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

//_________________________________________________________________________
/// \class AliMCAnalysisUtils
/// \ingroup CaloTrackCorrelationsBase
/// \brief Class with analysis utils for simulations.
///
/// Class containing utility methods for  analysis of simulations, ESD or AOD format.
/// Contains:
///  * method to check the origin of a given track/cluster, in several levels of ancestry
///  * methods to check if a cluster contains several particles, particularly decay photons pairs form same neutral meson.
///  * method to obtain the generated jets
///
/// More information can be found in this [twiki](https://twiki.cern.ch/twiki/bin/viewauth/ALICE/PhotonHadronCorrelations).
///
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, LPSC-IN2P3-CNRS
//_________________________________________________________________________

// --- ROOT system ---
#include <TObject.h>
#include <TString.h>
#include <TLorentzVector.h>
class TList ;
class TVector3;
class TClonesArray;

//--- AliRoot system ---
class AliMCEvent;
class AliGenEventHeader;

class AliMCAnalysisUtils : public TObject {
	
 public: 
  
  AliMCAnalysisUtils() ;          // ctor
  
  virtual ~AliMCAnalysisUtils() ; // virtual dtor
    
  //--------------------------------------
  // Enum with tag for origin of particles
  //--------------------------------------
  
  /// "Mostly" photon parent types on line 1,
  ///                                                  
  /// then "mostly" electron parent types on line 3, (e.g. electrons can
  /// come from pi0 decay)                              
  /// then charged particles on line 4,                                                                                    
  /// followed by other and unknown on line 5                                                                              
  enum mcTypes { kMCPhoton,     kMCPrompt,      kMCFragmentation, kMCISR,    
                 kMCPi0Decay,   kMCEtaDecay,    kMCOtherDecay,
                 kMCDecayPairLost, kMCDecayPairInCalo,
                 kMCDecayDalitz,kMCConversion,  kMCElectron,
                 kMCEFromCFromB,kMCEFromC,      kMCEFromB,        kMCZDecay, kMCWDecay,
                 kMCMuon,       kMCPion,        kMCPi0,           kMCKaon,   kMCEta,      kMCProton,   
                 kMCAntiProton, kMCNeutron,     kMCAntiNeutron,
                 kMCUnknown,    kMCBadLabel                                                         } ;
  
  //--------------------------------------
  // Methods to check origin of clusters
  //--------------------------------------
  
  Int_t   CheckCommonAncestor(Int_t index1, Int_t index2, const AliMCEvent* mcevent, 
			      Int_t & ancPDG, Int_t & ancStatus, TLorentzVector & momentum, TVector3 & prodVertex) ;
  
  Int_t   CheckOrigin(Int_t label, const AliMCEvent* mcevent) ;  
  Int_t   CheckOrigin(const Int_t *labels, Int_t nlabels, const AliMCEvent* mcevent, const TObjArray *arrayCluster = 0x0) ; 
  
  void    CheckOverlapped2GammaDecay(const Int_t *labels, Int_t nlabels, Int_t mesonIndex, const AliMCEvent* mcevent, Int_t & tag); 
  
  void    CheckLostDecayPair(const TObjArray *arrayCluster, Int_t iMom, Int_t iParent, const AliMCEvent* mcevent, Int_t & tag); 
  
  TLorentzVector GetMother     (Int_t label,const AliMCEvent* mcevent, Bool_t & ok);
  TLorentzVector GetMother     (Int_t label,const AliMCEvent* mcevent, Int_t & pdg, Int_t & status, Bool_t & ok);
  TLorentzVector GetMother     (Int_t label,const AliMCEvent* mcevent, Int_t & pdg, Int_t & status, Bool_t & ok, Int_t & momLabel);
  TLorentzVector GetGrandMother(Int_t label,const AliMCEvent* mcevent, Int_t & pdg, Int_t & status, Bool_t & ok, Int_t & grandMomLabel, Int_t & greatMomLabel);

  TLorentzVector GetMotherWithPDG(Int_t label, Int_t pdg,const AliMCEvent* mcevent, Bool_t & ok, Int_t & momLabel);
  
  void GetMCDecayAsymmetryAngleForPDG(Int_t label, Int_t pdg,const AliMCEvent* mcevent,
                                      Float_t & asy, Float_t & angle, Bool_t & ok);

  Int_t          GetNDaughters(Int_t label,const AliMCEvent* mcevent, Bool_t & ok);
  TLorentzVector GetDaughter  (Int_t daughter, Int_t label,const AliMCEvent* mcevent,
                               Int_t & pdg, Int_t & status, Bool_t & ok, Int_t & daugLabel, TVector3 & prodVertex);

  Int_t          GetNOverlaps(const Int_t * label, UInt_t nlabels,
                              Int_t mctag, Int_t mesonLabel,
                              AliMCEvent* mcevent,
                              Int_t *overpdg, Int_t *overlabel);
  
  //Check or set the bits produced in the above methods
  void    SetTagBit(Int_t &tag, UInt_t set) const {
    // Set bit of type set (mcTypes) in tag
    tag |= (1<<set) ; 
  } 
  
  Bool_t  CheckTagBit(Int_t tag, UInt_t test) const {
    // Check if in tag the bit test (mcTypes) is set.
    if (tag & (1<<test) ) return  kTRUE ;    
    else return kFALSE ;
  }
  
  //--------------------------------------
  // Other methods
  //--------------------------------------
    
  // Method to recover MC jets stored in generator
  TList * GetJets(AliMCEvent* mcevent, AliGenEventHeader * mcheader, Int_t eventNumber) ;
  
  void    SetDebug(Int_t deb)           { fDebug=deb           ; }
  Int_t   GetDebug()              const { return fDebug        ; }	
  
  enum generator {kPythia = 0, kHerwig = 1, kHijing = 2, kBoxLike = 3 } ;
  void    SetMCGenerator(Int_t   mcgen) ;
  void    SetMCGenerator(TString mcgen) ;
  Int_t   GetMCGenerator()        const { return fMCGenerator  ; }
  TString GetMCGeneratorString()  const { return fMCGeneratorString ; }
  
  void    Print(const Option_t * opt) const;
  void    PrintAncestry(AliMCEvent* mcevent, Int_t label, Int_t nGenerMax = 1000) const;
  void    PrintMCTag(Int_t tag) const;

 private:

  Int_t          fCurrentEvent;        ///<  Current Event number - GetJets()
  
  Int_t          fDebug;               ///<  Debug level
  
  TList        * fJetsList;            ///<  List of jets - GetJets()
  
  Int_t          fMCGenerator;         ///<  MC generator used to generate data in simulation
  
  TString        fMCGeneratorString;   ///<  MC generator used to generate data in simulation
  
  TLorentzVector fDaughMom;            //!<! particle momentum
  
  TLorentzVector fDaughMom2;           //!<! particle momentum
  
  TLorentzVector fMotherMom;           //!<! particle momentum
  
  TLorentzVector fGMotherMom;          //!<! particle momentum
  
  /// Copy constructor not implemented.
  AliMCAnalysisUtils & operator = (const AliMCAnalysisUtils & mcu) ; 
  
  /// Assignment operator not implemented.
  AliMCAnalysisUtils(              const AliMCAnalysisUtils & mcu) ; 
  
  /// \cond CLASSIMP
  ClassDef(AliMCAnalysisUtils,7) ;
  /// \endcond

} ;

#endif //ALIMCANALYSISUTILS_H



