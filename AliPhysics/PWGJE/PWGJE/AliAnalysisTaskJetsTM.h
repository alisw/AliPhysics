#ifndef AliAnalysisTaskJetsTM_cxx
#define AliAnalysisTaskJetsTM_cxx
//
// Thrust Major (TM) analysis of reconstructed jets.
// TM is the thrust in the plane perpendicular to the jet axis
// The present amalysis performs the following steps:
// (a) Construct to orthogonal unit vectors (e1, e2) in the plane perpendicular to the jet axis
// (b) Calculate the components of all particles with jT > 1 GeV with respect to e1, e2
// (c) Construct the sphericity matrix
// (d) Find the two orthogonal eigenvectors of the spericity matrix
// (e) Caluclate the components of all particles with jT < 1 GeV in the reference frame spanned by the eigenvectors
// (f) Calculate the azimuthal angle in this frame
//
//
// Author: andreas.morsch@cern.ch



class TH1F;
class TH2F;
class TList;
class TProfile;


#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskJetsTM : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskJetsTM(const char *name = "AliAnalysisTaskJetsTM");
  virtual ~AliAnalysisTaskJetsTM() {}
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
 private:

  TList       *fHists;          // List of histos
  TH1F        *fPtH;            // pT of reconstructed Jets
  TH1F        *fPtTH;           // pT of reconstructed tracks
  TH1F        *fPhiM;           // Phi Major distribtion
  TH2F        *fPhiMPt;         // Phi Major distribtion vs pt 
  TH2F        *fPhiMPtJ;        // Phi Major distribtion vs pt jet
  TH2F        *fPtSum;          // pT of reconstructed Jets
  
  AliAnalysisTaskJetsTM(const AliAnalysisTaskJetsTM&);            // not implemented
  AliAnalysisTaskJetsTM& operator=(const AliAnalysisTaskJetsTM&); // not implemented
  
  ClassDef(AliAnalysisTaskJetsTM, 1); // Jet Thrust Major Analysis
};

#endif
