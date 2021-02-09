/// \class AliAnalysisTaskJetShapeDeriv
/// \brief Background fluctuation studies: dMdpT spectrum for PYTHIA and single track embedding for derivative subtraction method
///
/// \ingroup PWGJEUSER
/// Most of the methods are in the class AliAnalysisTaskJetShapeBase
/// The Run() method is adapted to the area subtracted jets
/// 
///
/// \author Marta Werveij, updates Chiara Bianchin
/// \date Modified Nov 2015

#ifndef ALIANALYSISTASKJETSHAPEDERIV_H
#define ALIANALYSISTASKJETSHAPEDERIV_H

#include "AliAnalysisTaskJetShapeBase.h"

class AliAnalysisTaskJetShapeDeriv : public AliAnalysisTaskJetShapeBase {
 public:

  AliAnalysisTaskJetShapeDeriv();
  AliAnalysisTaskJetShapeDeriv(const char *name);
  virtual ~AliAnalysisTaskJetShapeDeriv();

  void                                UserCreateOutputObjects();

  //Setters
  void SetPartialExclusion(Bool_t b)                            { fPartialExclusion  = b   ; }
  
 protected:
  //Bool_t                              Run();
  Bool_t                              FillHistograms();

 private:
  Float_t         fM1st;                                           ///< 1st order subtracted jet mass
  Float_t         fM2nd;                                           ///< 2nd order subtracted jet mass
  Float_t         fDeriv1st;                                       ///< 1st derivative
  Float_t         fDeriv2nd;                                       ///< 2nd derivative
  Bool_t          fPartialExclusion;                               ///< randomly esclude areas according to Ncoll
  TH2F          **fh2MSubPtSubAll;                                 //!<! subtracted jet mass vs subtracted jet pT
  TH2F          **fh2PtTrueSubFacV1;                               //!<! true pT vs -(rho+rhom)*V1
  TH2F          **fh2PtRawSubFacV1;                                //!<! raw pT vs -(rho+rhom)*V1
  TH2F          **fh2PtCorrSubFacV1;                               //!<! subtracted pT vs -(rho+rhom)*V1
  TH2F          **fh2NConstSubFacV1;                               //!<! N constituents vs -(rho+rhom)*V1
  TH2F          **fh2PtTrueSubFacV2;                               //!<! true pT vs 0.5(rho+rhom)^2*V2
  TH2F          **fh2PtRawSubFacV2;                                //!<! raw pT vs 0.5(rho+rhom)^2*V2
  TH2F          **fh2PtCorrSubFacV2;                               //!<! subtracted pT vs 0.5(rho+rhom)^2*V2
  TH2F          **fh2NConstSubFacV2;                               //!<! N constituents vs 0.5(rho+rhom)^2*V2

  AliAnalysisTaskJetShapeDeriv(const AliAnalysisTaskJetShapeDeriv&);            // not implemented
  AliAnalysisTaskJetShapeDeriv &operator=(const AliAnalysisTaskJetShapeDeriv&); // not implemented

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskJetShapeDeriv, 13); //updated to the new histograms, -1 if don't commit
  /// \endcond
};
#endif



