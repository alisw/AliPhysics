#ifndef AliAnalysisTaskESDStrange_cxx
#define AliAnalysisTaskESDStrange_cxx

// macro to study V0s
// Author: H.Ricaud, Helene.Ricaud@IReS.in2p3.fr

class TH1F;
class TH2F;
class TList;
class AliESDEvent;

Double_t myRap(Double_t,Double_t);

#include "AliAnalysisTask.h"

class AliAnalysisTaskESDStrange : public AliAnalysisTask {
 public:
  AliAnalysisTaskESDStrange(const char *name = "AliAnalysisTaskESDStrange");
  virtual ~AliAnalysisTaskESDStrange() {}
  
  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
 private:
  AliESDEvent *fESD;    //ESD object
  TList       *fListHist;

  TH1F        *fHistTrackPerEvent;
  
  TH1F        *fHistPrimaryVertexX;
  TH1F        *fHistPrimaryVertexY;
  TH1F        *fHistPrimaryVertexZ;

  TH2F        *fHistDcaPosToPrimVertex;
  TH2F        *fHistDcaNegToPrimVertex;
  TH2F        *fHistDcaPosToPrimVertexZoom;
  TH2F        *fHistDcaNegToPrimVertexZoom;
  TH2F        *fHistRadiusV0;
  TH2F        *fHistDecayLengthV0;
  TH2F        *fHistDcaV0Daughters;
  TH2F        *fHistChi2;
  TH2F        *fHistCosPointAngle;
  TH2F        *fHistCosPointAngleZoom;
  TH2F        *fHistPtVsYK0s;
  TH2F        *fHistPtVsYK0sMI;
  TH2F        *fHistPtVsYLambda;
  TH2F        *fHistPtVsYLambdaMI;
  TH2F        *fHistPtVsYAntiLambda;
  TH2F        *fHistPtVsYAntiLambdaMI;

  TH1F        *fHistMassK0;
  TH1F        *fHistMassK0MI;
  TH1F        *fHistMassLambda;
  TH1F        *fHistMassLambdaMI;
  TH1F        *fHistMassAntiLambda;
  TH1F        *fHistMassAntiLambdaMI;

  TH2F        *fHistArmenterosPodolanski;
  TH2F        *fHistArmenterosPodolanskiMI;
   
  AliAnalysisTaskESDStrange(const AliAnalysisTaskESDStrange&); 
  AliAnalysisTaskESDStrange& operator=(const AliAnalysisTaskESDStrange&); 

  ClassDef(AliAnalysisTaskESDStrange, 1); 
};

#endif
