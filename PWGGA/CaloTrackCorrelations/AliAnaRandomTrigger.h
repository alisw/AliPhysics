#ifndef ALIANARANDOMTRIGGER_H
#define ALIANARANDOMTRIGGER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

//_________________________________________________________________________
/// \class AliAnaRandomTrigger
/// \ingroup CaloTrackCorrelationsAnalysis 
/// \brief Gerenate a random trigger.
///
/// Generate a random trigger, input for other analysis
/// Set flat energy distribution over acceptance of EMCAL, PHOS or CTS
/// Be careful, correlate only with Min Bias events this random trigger particle
///
/// More information can be found in this [twiki](https://twiki.cern.ch/twiki/bin/viewauth/ALICE/PhotonHadronCorrelations)
/// and particularly in this [section](https://twiki.cern.ch/twiki/bin/viewauth/ALICE/PhotonHadronCorrelations#AliAnaRandomTrigger).
///
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, LPSC-IN2P3-CNRS
//_________________________________________________________________________

// Root system
class TH2F; 
#include <TRandom3.h>

// Analysis system
#include "AliAnaCaloTrackCorrBaseClass.h"
 
class AliAnaRandomTrigger : public AliAnaCaloTrackCorrBaseClass {
  
 public:
    
               AliAnaRandomTrigger() ;
    
  // Virtual destructor.
  virtual     ~AliAnaRandomTrigger() { ; }

  Bool_t       ExcludeDeadBadRegions(Float_t eta, Float_t phi);
  
  TObjString * GetAnalysisCuts();
  
  TList      * GetCreateOutputObjects();
    
  void         InitParameters();
    
  void         MakeAnalysisFillAOD()  ;
  
  void         MakeAnalysisFillHistograms() ; 
  
  void         Print(const Option_t * opt) const;
  
  void         SetEtaCut(Float_t min, Float_t max) { fEtaCut[0] = min ; fEtaCut[1] = max;}
  
  void         SetPhiCut(Float_t min, Float_t max) { fPhiCut[0] = min ; fPhiCut[1] = max;} // radians
  
  void         SetNumberOfRandomParticles(Int_t n) { fNRandom   = n   ; }
  
  void         SetTriggerDetector(TString det) ;
  void         SetTriggerDetector(Int_t   det) ;
  
 private:

  Int_t      fTriggerDetector ;       ///<  Detector : EMCAL, PHOS, CTS
  TString    fTriggerDetectorString ; ///<  Detector : EMCAL, PHOS, CTS
  Float_t    fEtaCut[2];              ///<  Eta acceptance
  Float_t    fPhiCut[2];              ///<  Phi acceptance, radians
  TRandom3   fRandom   ;              ///<  Random generator
  Int_t      fNRandom  ;              ///<  Number of random particles per event
  
  TLorentzVector fMomentum;           //!<! Avoid generating TLorentzVectors per event.
  
  //Constrol histograms 
  TH1F     * fhPt;                    //!<! pT distribution
  TH2F     * fhPhi;                   //!<! phi distribution vs pT, negative
  TH2F     * fhEta;                   //!<! eta distribution vs pT, negative
  TH2F     * fhEtaPhi;                //!<! eta vs phi distribution of positive charge
  
  /// Copy constructor not implemented.
  AliAnaRandomTrigger(              const AliAnaRandomTrigger & r) ;
    
  /// Assignment operator not implemented.
  AliAnaRandomTrigger & operator = (const AliAnaRandomTrigger & r) ;
  
  /// \cond CLASSIMP
  ClassDef(AliAnaRandomTrigger,5) ;
  /// \endcond

} ;

#endif //ALIANARANDOMTRIGGER_H



