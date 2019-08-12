/// \class AliAnalysisTaskJetShapeConst
/// \brief Background fluctuation studies: dMdpT spectrum for PYTHIA and single track embedding for Constituent subtraction method
///
/// \ingroup PWGJEUSER
/// Most of the methods are in the class AliAnalysisTaskJetShapeBase
/// The Run() method is adapted to the Constituent subtracted jets
/// 
///
/// \author Marta Werveij, updates Chiara Bianchin
/// \date Modified Nov 2015

#ifndef ALIANALYSISTASKJETSHAPECONST_H
#define ALIANALYSISTASKJETSHAPECONST_H


#include "AliAnalysisTaskJetShapeBase.h"

class TH2F;
class TH1F;

class AliAnalysisTaskJetShapeConst : public AliAnalysisTaskJetShapeBase {
 public:

  AliAnalysisTaskJetShapeConst();
  AliAnalysisTaskJetShapeConst(const char *name);
  virtual ~AliAnalysisTaskJetShapeConst();

  void                                UserCreateOutputObjects();
  //void                                Terminate(Option_t *option);
  void SetJetContainerSub(Int_t c)                              { fContainerSub      = c   ; }
  
 protected:

  Bool_t                              FillHistograms();

  
 private:
  
  TH1F 	        *fhptjetSMinusSingleTrack;                         //!<! pT distribution of jets subtracting the pT of the embedded track
  TH2F          *fhJet1vsJetTag;                                   //!<! N jet vs N jet tagged
  TH1F          *fhNconstit;                                       //!<! number of constituents of the matched jets
  TH1F          *fhAreaJet;                                        //!<! area of the matched jet

  
  AliAnalysisTaskJetShapeConst(const AliAnalysisTaskJetShapeConst&);            // not implemented
  AliAnalysisTaskJetShapeConst &operator=(const AliAnalysisTaskJetShapeConst&); // not implemented
  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskJetShapeConst, 14);
  /// \endcond
};
#endif



