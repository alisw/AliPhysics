#ifndef ALIANAPARTICLEPARTONCORRELATION_H
#define ALIANAPARTICLEPARTONCORRELATION_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

//_________________________________________________________________________
/// \class AliAnaParticlePartonCorrelation
/// \ingroup CaloTrackCorrelationsAnalysis 
/// \brief Particle-parton correlation
///
/// Class that contains the algorithm for the analysis of particle-parton correlation
/// Particle (for example direct gamma) must be found in a previous analysis
///
/// More information can be found in this [twiki](https://twiki.cern.ch/twiki/bin/viewauth/ALICE/PhotonHadronCorrelations)
/// and particularly in this [section](https://twiki.cern.ch/twiki/bin/viewauth/ALICE/PhotonHadronCorrelations#AliAnaParticlePartonCorrelation).
///
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, LPSC-IN2P3-CNRS
//_________________________________________________________________________

// --- ROOT ---
class TH2F ;

// --- ANALYSIS ---
#include "AliAnaCaloTrackCorrBaseClass.h"

class AliAnaParticlePartonCorrelation : public AliAnaCaloTrackCorrBaseClass {
       
public:
  
           AliAnaParticlePartonCorrelation() ;
    
  /// Virtual destructor.
  virtual ~AliAnaParticlePartonCorrelation() { ; }
  
  TList *  GetCreateOutputObjects();
  
  void     InitParameters();
    
  void     MakeAnalysisFillAOD()  ;
  
  void     MakeAnalysisFillHistograms() ;
  
  void     Print(const Option_t * opt) const;
  
private:
  
  TH2F * fhDeltaEtaNearParton; //!<! Difference of parton eta and prompt trigger particle eta
    
  TH2F * fhDeltaPhiNearParton; //!<! Difference of parton phi and prompt trigger particle phi
    
  TH2F * fhDeltaPtNearParton;  //!<! Difference of parton pT and prompt trigger particle pT
    
  TH2F * fhPtRatNearParton;    //!<! Ratio of parton pT and prompt trigger particle pT
  
  TH2F * fhDeltaEtaAwayParton; //!<! Difference of parton eta and prompt trigger particle eta
    
  TH2F * fhDeltaPhiAwayParton; //!<! Difference of parton phi and prompt trigger particle phi
    
  TH2F * fhDeltaPtAwayParton;  //!<! Difference of parton pT and prompt trigger particle pT
    
  TH2F * fhPtRatAwayParton;    //!<! Ratio of parton pT and prompt trigger particle pT
  
  /// Copy constructor not implemented.
  AliAnaParticlePartonCorrelation              (const AliAnaParticlePartonCorrelation & g) ;
    
  /// Assignment operator not implemented.
  AliAnaParticlePartonCorrelation & operator = (const AliAnaParticlePartonCorrelation & g) ;
  
  /// \cond CLASSIMP
  ClassDef(AliAnaParticlePartonCorrelation,1) ;
  /// \endcond

} ;

#endif //ALIANAPARTICLEPARTONCORRELATION_H



