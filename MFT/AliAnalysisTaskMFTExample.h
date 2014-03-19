#ifndef AliAnalysisTaskMFTExample_H
#define AliAnalysisTaskMFTExample_H

#include "AliAnalysisTaskSE.h"
#include "TH1D.h"
#include "TList.h"

//====================================================================================================================================================

class  AliAnalysisTaskMFTExample : public AliAnalysisTaskSE {

public:

  enum {kGenerated, kReconstructed};
 
  AliAnalysisTaskMFTExample();
  AliAnalysisTaskMFTExample(const char *name);

  virtual ~AliAnalysisTaskMFTExample() {
    delete fHistPtSingleMuons;
    delete fHistPtSingleMuonsFromJpsi;
    delete fHistPtDimuonsOS;
    delete fHistMassDimuonsOS;
    delete fHistPtDimuonsJpsi;
    delete fHistMassDimuonsJpsi;
    delete fHistResidualXVtxJpsi;
    delete fHistResidualYVtxJpsi;
    delete fHistResidualZVtxJpsi;
  }
  
  void SetVertexMode(Int_t vertexMode) { fVertexMode = vertexMode; }
  void SetVtxResolutionITS(Double_t sigmaX, Double_t sigmaY, Double_t sigmaZ) {
    fVtxResolutionITS[0] = sigmaX;
    fVtxResolutionITS[1] = sigmaY;
    fVtxResolutionITS[2] = sigmaZ;
  }

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);

private:

  Double_t fPrimaryVertex[3], fVtxResolutionITS[3];
  Int_t fVertexMode;

  TList *fHistogramList;

  TH1D *fHistPtSingleMuons, *fHistPtSingleMuonsFromJpsi;
  TH1D *fHistPtDimuonsOS, *fHistMassDimuonsOS;
  TH1D *fHistPtDimuonsJpsi, *fHistMassDimuonsJpsi;

  TH1D *fHistResidualXVtxJpsi, *fHistResidualYVtxJpsi, *fHistResidualZVtxJpsi;

  ClassDef(AliAnalysisTaskMFTExample, 1) // example of analysis

};

//====================================================================================================================================================

#endif
