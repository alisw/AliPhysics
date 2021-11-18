#ifndef AliAnalysisTaskAngularRatiosCorrelation_H
#define AliAnalysisTaskAngularRatiosCorrelation_H
class AliPIDResponse;

#include "AliAnalysisTaskSE.h"
#include "AliPIDResponse.h"
#include "AliEventCuts.h"
class AliAnalysisTaskAngularRatiosCorrelation : public AliAnalysisTaskSE
{
public:
    AliAnalysisTaskAngularRatiosCorrelation();
    AliAnalysisTaskAngularRatiosCorrelation(const char *name);
    virtual ~AliAnalysisTaskAngularRatiosCorrelation();

    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t *option);
    virtual void Terminate(Option_t *option);
    virtual void SetFilterBit(UInt_t Bit)
    {
        filterBit = Bit;
    }
    virtual void SetParams(int f_iTask, int f_nPhiBins, int f_nVertexBins, int f_nPBins, int f_minCent, int f_maxCent,
                           Float_t f_minP, Float_t f_maxP, Float_t f_Vertexmin, Float_t f_Vertexmax,
                           Bool_t f_MCGen, Bool_t f_pbpb, int f_nSigma, int f_nCrossedRows, int f_movePhi)
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
        movePhi = f_movePhi;
    }

private:
    AliAODEvent *fAOD;            //! input event
    AliPIDResponse *fPIDResponse; //! pid response object
    TList *fOutputList;           //! output list
    TH1F *fHistEventsCut;         //! QA events
    TH1F *fHistTracksCut;         //! QA tracks
    UInt_t filterBit;             //
    Bool_t MCGen, pbpb;           // is MC or not; is Pb-Pb or pp
    AliEventCuts *fAliEventCuts;  //!
    int iTask, nPhiBins, nVertexBins, nPBins, minCent, maxCent, nSigma, nCrossedRows, movePhi;
    static const int nCentrClasses = 4, nEtaClasses = 16, nSorts = 8, nSubsamples = 20, nPhiWindows = 16;
    Float_t minP, maxP, Vertexmin, Vertexmax;
    static const int SortPairs = 6 * (6 + 1) / 2;
    Float_t PtCut[3] = {0.2, 0.5, 0.5};
    Float_t nSigmaBoundary[3] = {0.5, 0.32, 0.7};

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
    TH2D *fHistQASPDTrackletsvsV0MCent; //!

    AliAnalysisTaskAngularRatiosCorrelation(const AliAnalysisTaskAngularRatiosCorrelation &);            // not implemented
    AliAnalysisTaskAngularRatiosCorrelation &operator=(const AliAnalysisTaskAngularRatiosCorrelation &); // not implemented

    ClassDef(AliAnalysisTaskAngularRatiosCorrelation, 1);
};

#endif
