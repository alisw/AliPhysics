#ifndef AliAnalysisTaskParticleYieldRatioCorrelations_H
#define AliAnalysisTaskParticleYieldRatioCorrelations_H

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
class AliAnalysisTaskParticleYieldRatioCorrelations : public AliAnalysisTaskSE
{
public:
    AliAnalysisTaskParticleYieldRatioCorrelations();
    AliAnalysisTaskParticleYieldRatioCorrelations(const char *name);
    virtual ~AliAnalysisTaskParticleYieldRatioCorrelations();

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
    virtual void SetCorrections(Bool_t f_PIDCorr, Bool_t f_TrackingCorr)
    {
        f_PIDCorr = PIDCorr;
        f_TrackingCorr = TrackingCorr;
    }

private:
    AliAODEvent *fAOD;            //! input event
    AliPIDResponse *fPIDResponse; //! pid response object
    TList *fOutputList;           //! output list
    TList *fInputList;           //! output list
    TH1F *fHistEventsCut;         //! QA events
    TH1F *fHistTracksCut;         //! QA tracks
    UInt_t filterBit;             //
    Bool_t IsMC, pbpb;           // is MC or not; is Pb-Pb or pp
    Bool_t SPDvsV0MCut, LargeTPCCut,  PIDCorr, TrackingCorr;
    AliEventCuts *fAliEventCuts;  //!
    int nPhiBins, nVertexBins, nPBins, minCent, maxCent, nCrossedRows, movePhi;
    static const int nCentrClasses = 10, nEtaClasses = 16, nSorts = 8, nSubsamples = 20, nPhiWindows = 16;
    Float_t minP, maxP, Vertexmin, Vertexmax, nSigma;
    Float_t CentrPercentiles[10];
    int nCentrClassesUsed;
    static const int SortPairs = 6 * (6 + 1) / 2;

    //Declare hists of efficiency maps
    TH3D *EfficiencyPID[nCentrClasses][nEtaClasses][nSorts];         //!
    TH3D *EfficiencyTracking[nCentrClasses][nEtaClasses][nSorts];    //!
    TH3D *ContaminationPID[nCentrClasses][nEtaClasses][nSorts];      //!
    TH3D *ContaminationTracking[nCentrClasses][nEtaClasses][nSorts]; //!

    TH1D *NEvents;                 //!
    TH1D *NEventsSub[nSubsamples]; //!
    TH2D *fDeDx;                   //!
    TH2D *fDeDxSorts[4];           //!
    TH2D *fTOF;                    //!
    TH2D *fTOFSorts[4];            //!

    TH1D *fHistMomentsInEtaRec[nCentrClasses][nSorts]; //!
    TH1D *fHistMomentsInEtaGen[nCentrClasses][nSorts]; //!
    TH1D *fHistMomentsInEtaReg[nCentrClasses][nSorts]; //!

    TH1D *fHistMomentsInEtaRecSub[nSubsamples][nCentrClasses][nSorts]; //!
    TH1D *fHistMomentsInEtaGenSub[nSubsamples][nCentrClasses][nSorts]; //!
    TH1D *fHistMomentsInEtaRegSub[nSubsamples][nCentrClasses][nSorts]; //!

    TH1D *fHistCrossMomentsInEtaRec[nCentrClasses][SortPairs]; //!
    TH1D *fHistCrossMomentsInEtaGen[nCentrClasses][SortPairs]; //!
    TH1D *fHistCrossMomentsInEtaReg[nCentrClasses][SortPairs]; //!

    TH1D *fHistCrossMomentsInEtaRecSub[nSubsamples][nCentrClasses][SortPairs]; //!
    TH1D *fHistCrossMomentsInEtaGenSub[nSubsamples][nCentrClasses][SortPairs]; //!
    TH1D *fHistCrossMomentsInEtaRegSub[nSubsamples][nCentrClasses][SortPairs]; //!

    TH1D *fHistMomentsInPhiRec[nCentrClasses][nSorts]; //!
    TH1D *fHistMomentsInPhiGen[nCentrClasses][nSorts]; //!
    TH1D *fHistMomentsInPhiReg[nCentrClasses][nSorts]; //!

    TH1D *fHistMomentsInPhiRecSub[nSubsamples][nCentrClasses][nSorts]; //!
    TH1D *fHistMomentsInPhiGenSub[nSubsamples][nCentrClasses][nSorts]; //!
    TH1D *fHistMomentsInPhiRegSub[nSubsamples][nCentrClasses][nSorts]; //!

    TH1D *fHistCrossMomentsInPhiRec[nCentrClasses][SortPairs]; //!
    TH1D *fHistCrossMomentsInPhiGen[nCentrClasses][SortPairs]; //!
    TH1D *fHistCrossMomentsInPhiReg[nCentrClasses][SortPairs]; //!

    TH1D *fHistCrossMomentsInPhiRecSub[nSubsamples][nCentrClasses][SortPairs]; //!
    TH1D *fHistCrossMomentsInPhiGenSub[nSubsamples][nCentrClasses][SortPairs]; //!
    TH1D *fHistCrossMomentsInPhiRegSub[nSubsamples][nCentrClasses][SortPairs]; //!

    TH1D *fHistMomentsInAllAccRec[nCentrClasses]; //!
    TH1D *fHistMomentsInAllAccGen[nCentrClasses]; //!
    TH1D *fHistMomentsInAllAccReg[nCentrClasses]; //!

    TH1D *fHistMomentsInAllAccRecSub[nSubsamples][nCentrClasses]; //!
    TH1D *fHistMomentsInAllAccGenSub[nSubsamples][nCentrClasses]; //!
    TH1D *fHistMomentsInAllAccRegSub[nSubsamples][nCentrClasses]; //!

    //QA hists
    TH1D *fHistQAVx;                    //! Vx hist
    TH1D *fHistQAVy;                    //! Vy hist
    TH1D *fHistQAVz;                    //! Vz hist
    TH1F *fHistQAPt;                    //! Overal Pt spectrum
    TH1F *fHistQAEta;                   //! Overal Eta spectrum
    TH1F *fHistQAPhi;                   //! Overal Phi spectrum
    TH1D *fHistQAClustersTPC;           //! Number of TPC clusters distribution
    TH1D *fHistQACrossedRowsTPC;        //! Number of TPC crossed rows
    TH1D *fHistQAChi2perNDF;            //!
    TH1D *fHistQAClustersITS;           //! Number of ITS clusters distribution
                                        //
    TH2D *fHistQAEtaPhi;                //!
    TH2D *fHistQAMomPt;                //!
    TH2D *fHistQASPDTrackletsvsV0MCent[3]; //!
    TH2D *fHistQAMultTPCvsESD[2]; //!
    TH2D *fHistQAMultTPCvsV0[2]; //!
    TH2D *fHistQAMultTrkvsMultTrkTOF[2]; //!


    AliAnalysisTaskParticleYieldRatioCorrelations(const AliAnalysisTaskParticleYieldRatioCorrelations &);            // not implemented
    AliAnalysisTaskParticleYieldRatioCorrelations &operator=(const AliAnalysisTaskParticleYieldRatioCorrelations &); // not implemented

    ClassDef(AliAnalysisTaskParticleYieldRatioCorrelations, 1);
};

#endif
