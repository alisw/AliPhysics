#ifndef AliAnalysisTaskParticleYieldRatioCorrelationsEfficiency_H
#define AliAnalysisTaskParticleYieldRatioCorrelationsEfficiency_H

class TList;
class TH1F;
class TH1D;
class TH2D;
class TH3D;
class TF1;
class AliPIDResponse;

#include "AliAnalysisTaskSE.h"
#include "AliPIDResponse.h"
#include "AliEventCuts.h"
class AliAnalysisTaskParticleYieldRatioCorrelationsEfficiency : public AliAnalysisTaskSE
{
public:
  AliAnalysisTaskParticleYieldRatioCorrelationsEfficiency();
  AliAnalysisTaskParticleYieldRatioCorrelationsEfficiency(const char *name);
  virtual ~AliAnalysisTaskParticleYieldRatioCorrelationsEfficiency();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);

  //Analysis setters
  virtual void SetFilterBit(UInt_t Bit)
  {
    filterBit = Bit;
  }
  virtual void SetNumbersOfBins(int f_nPhiBins, int f_nVertexBins, int f_nPBins)
  {
    nPhiBins = f_nPhiBins;
    nVertexBins = f_nVertexBins;
    nPBins = f_nPBins;
  }
  virtual void SetCentLim(int f_minCent, int f_maxCent)
  {
    minCent = f_minCent;
    maxCent = f_maxCent;
  }
  virtual void SetPLim(Float_t f_minP, Float_t f_maxP)
  {
    minP = f_minP;
    maxP = f_maxP;
  }
  virtual void SetVertLim(Float_t f_Vertexmin, Float_t f_Vertexmax)
  {
    Vertexmin = f_Vertexmin;
    Vertexmax = f_Vertexmax;
  }
  virtual void SetParams(Bool_t f_IsMC, Bool_t f_pbpb, Float_t f_nSigma, int f_nCrossedRows, int f_movePhi)
  {
    IsMC = f_IsMC;
    pbpb = f_pbpb; // pp if false
    nSigma = f_nSigma;
    nCrossedRows = f_nCrossedRows;
    movePhi = f_movePhi;
  }
  virtual void SetCentralities(int f_nCentrClassesUsed, Float_t cent0, Float_t cent1, Float_t cent2, Float_t cent3, Float_t cent4, Float_t cent5)
    {
        nCentrClassesUsed = f_nCentrClassesUsed;
        Float_t f_CentrPercentiles[] = {cent0, cent1, cent2, cent3, cent4, cent5};
        for (int i = 0; i < nCentrClassesUsed+1; i++)
        {
        CentrPercentiles[i] = f_CentrPercentiles[i];
        }
    }
    virtual void SetCuts(Bool_t f_SPDvsV0MCut, Bool_t f_LargeTPCCut)
    {
        SPDvsV0MCut = f_SPDvsV0MCut;
        LargeTPCCut = f_LargeTPCCut;
    }
    
private:
  AliAODEvent *fAOD;            //! input event
  AliPIDResponse *fPIDResponse; //! pid response object
  TList *fOutputList;           //! output list
  TH1F *fHistEventsCut;         //!
  TH1F *fHistTracksCut;         //!
  UInt_t filterBit;             //
  Bool_t pbpb, IsMC;
  Bool_t SPDvsV0MCut, LargeTPCCut;
  AliEventCuts *fAliEventCuts; //!
  int nPhiBins, nVertexBins, nPBins, minCent, maxCent, nCrossedRows, movePhi;
  Float_t minP, maxP, Vertexmin, Vertexmax, nSigma;
  Float_t CentrPercentiles[10];
  int nCentrClassesUsed;
  static const int nCentrClasses = 10, nEtaClasses = 16, nSorts = 8;
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

  TH1D *purityAll[8]; //!
  TH1D *purity[8][8]; //!

  TH2D *fHistQASPDTrackletsvsV0MCent[3]; //!
  TH2D *fHistQAMultTPCvsESD[2]; //!
  TH2D *fHistQAMultTPCvsV0[2]; //!
  TH2D *fHistQAMultTrkvsMultTrkTOF[2]; //!

  AliAnalysisTaskParticleYieldRatioCorrelationsEfficiency(const AliAnalysisTaskParticleYieldRatioCorrelationsEfficiency &);            // not implemented
  AliAnalysisTaskParticleYieldRatioCorrelationsEfficiency &operator=(const AliAnalysisTaskParticleYieldRatioCorrelationsEfficiency &); // not implemented

  ClassDef(AliAnalysisTaskParticleYieldRatioCorrelationsEfficiency, 1);
};

#endif
