#ifndef AliAnalysisTaskAngularRatiosCorrelationEfficiency_H
#define AliAnalysisTaskAngularRatiosCorrelationEfficiency_H
class AliPIDResponse;

#include "AliAnalysisTaskSE.h"
#include "AliPIDResponse.h"
#include "AliEventCuts.h"
class AliAnalysisTaskAngularRatiosCorrelationEfficiency : public AliAnalysisTaskSE
{
public:
  AliAnalysisTaskAngularRatiosCorrelationEfficiency();
  AliAnalysisTaskAngularRatiosCorrelationEfficiency(const char *name);
  virtual ~AliAnalysisTaskAngularRatiosCorrelationEfficiency();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);
  //virtual void            SetFilterBit(UInt_t* FilterBit);
  virtual void SetFilterBit(UInt_t Bit)
  {
    filterBit = Bit;
  }
  virtual void SetParams(int f_iTask, int f_nPhiBins, int f_nVertexBins, int f_nPBins, int f_minCent, int f_maxCent,
                         Float_t f_minP, Float_t f_maxP, Float_t f_Vertexmin, Float_t f_Vertexmax,
                         Bool_t f_MCGen, Bool_t f_pbpb, int f_nSigma, int f_nCrossedRows)
  {
    iTask = f_iTask;
    nPhiBins = f_nPhiBins;
    nVertexBins = f_nVertexBins;
    nPBins = f_nPBins;
    minCent = f_minCent;
    maxCent = f_maxCent;
    minP = f_minP;
    maxP = f_maxP;
    Vertexmin = f_Vertexmin;
    Vertexmax = f_Vertexmax;
    MCGen = f_MCGen;
    pbpb = f_pbpb;
    nSigma = f_nSigma;
    nCrossedRows = f_nCrossedRows;
  }

private:
  AliAODEvent *fAOD;            //! input event
  AliPIDResponse *fPIDResponse; //! pid response object
  TList *fOutputList;           //! output list
  TH1F *fHistEventsCut;         //!
  TH1F *fHistTracksCut;         //!
  UInt_t filterBit;             //
  Bool_t MCGen, pbpb;
  AliEventCuts *fAliEventCuts; //!
  int iTask, nPhiBins, nVertexBins, nPBins, minCent, maxCent, nPhiWindows=16, nSigma, nCrossedRows;
  static const int nCentrClasses = 4, nEtaClasses = 16, nSorts = 8;
  Float_t minP, maxP, Vertexmin, Vertexmax;
  Float_t PtCut[3] = {0.2, 0.5, 0.5};
  Float_t nSigmaBoundary[3] = {0.5, 0.32, 0.7};

  //Efficiency Maps
  //Tracking:
  TH3D *NGenTracks[nCentrClasses][nEtaClasses][nSorts];            //!
  TH3D *EfficiencyTracking[nCentrClasses][nEtaClasses][nSorts];    //!
  TH3D *ContaminationTracking[nCentrClasses][nEtaClasses][nSorts]; //!
  TH3D *NReconstTracks[nCentrClasses][nEtaClasses][nSorts];        //!
  //PID:
  TH3D *NTrueTracks[nCentrClasses][nEtaClasses][nSorts];      //!
  TH3D *EfficiencyPID[nCentrClasses][nEtaClasses][nSorts];    //!
  TH3D *ContaminationPID[nCentrClasses][nEtaClasses][nSorts]; //!
  TH3D *NTracksInCut[nCentrClasses][nEtaClasses][nSorts];     //!

  TH2D *fDeDx;         //!
  TH2D *fDeDxSorts[4]; //!
  TH2D *fTOF;          //!
  TH2D *fTOFSorts[4];  //!

  TH1D *purityAll[nSorts];      //!
  TH1D *purity[nSorts][nSorts]; //!

  TH2D *fHistQASPDTrackletsvsV0MCent; //!

  AliAnalysisTaskAngularRatiosCorrelationEfficiency(const AliAnalysisTaskAngularRatiosCorrelationEfficiency &);            // not implemented
  AliAnalysisTaskAngularRatiosCorrelationEfficiency &operator=(const AliAnalysisTaskAngularRatiosCorrelationEfficiency &); // not implemented

  ClassDef(AliAnalysisTaskAngularRatiosCorrelationEfficiency, 1);
};

#endif
