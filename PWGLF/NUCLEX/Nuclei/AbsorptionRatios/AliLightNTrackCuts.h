#ifndef ALILIGHTNTRACKCUTS_H
#define ALILIGHTNTRACKCUTS_H

/*
 * AliLightNTrackCuts.h
 *
 *  Created on: Nov 14, 2017
 *      Author: gu74req
 */

#include "Rtypes.h"
#include "AliLightNTrack.h"
#include "AliLightNTrackMCHist.h"
#include "AliLightNTrackHist.h"
#include "TMath.h"
#include <iostream>

class AliLightNTrackCuts {
public:
    AliLightNTrackCuts();
    virtual ~AliLightNTrackCuts();
    static AliLightNTrackCuts *PrimProtonCuts(bool isMC, bool DCAPlots, bool CombSigma, bool ContribSplitting);
    static AliLightNTrackCuts *PrimDeuteronCuts(bool isMC, bool DCAPlots, bool CombSigma, bool ContribSplitting);
    static AliLightNTrackCuts *DecayProtonCuts(bool isMC,bool ContribSplitting);
    static AliLightNTrackCuts *DecayPionCuts(bool isMC,bool ContribSplitting);
    static AliLightNTrackCuts *Xiv0PionCuts(bool isMC,bool ContribSplitting);
    static AliLightNTrackCuts *Xiv0ProtonCuts(bool isMC,bool ContribSplitting);
    static AliLightNTrackCuts *XiBachPionCuts(bool isMC,bool ContribSplitting);
    //  static AliLightNTrackCuts *OmegaKaonCuts(bool isMC,bool ContribSplitting);
    
    //Setters for Plots
    void SetPlotDCADist(bool plot) {fDCAPlots=plot;};
    void SetPlotCombSigma(bool plot) {fCombSigma=plot;};
    void SetPlot3DPMassDcaxy(bool plot) {f3DPlotPmass2dca= kFALSE;};
    void SetPlotContrib(bool plot) {fContribSplitting=plot;};
    void SetIsMonteCarlo(bool isMC) {fMCData=isMC;};
    void SetFillQALater(bool FillLater){fFillQALater=FillLater;};
    bool GetIsMonteCarlo() {return fMCData;};
    //Setters for the Track Cuts
    void SetCheckFilterBit(bool check){fCheckFilterBit = check;};
    void SetFilterBit(UInt_t FilterBit){fFilterBit = FilterBit; fCheckFilterBit = kTRUE;};
    void SetPtRange(double pmin, double pmax){fpTmin = pmin; fpTmax = pmax; fcutPt = kTRUE;};
    void SetEtaRange(double etamin, double etamax){fetamin=etamin; fetamax=etamax; fcutEta = kTRUE;};
    double GetEtaMin() {return fetamin;}
    double GetEtaMax() {return fetamax;}
    void SetChi2perNDFCut(double Chi2){fChi2perNDF = Chi2; fCutChi2 = kTRUE;};
    void SetTPCRatioCut(double Ratio){fTPCRatio = Ratio;};
    void SetCutCharge(int charge){fcutCharge = kTRUE; fCharge = charge;};
    void SetCheckPileUpITS(bool check){fCheckPileUpITS=check;};
    void SetCheckPileUpTOF(bool check){fCheckPileUpTOF=check;};
    void SetCheckPileUp(bool check){fCheckPileUp=check;};
    void SetNClsTPC(int nCls){fnTPCCls = nCls; fcutnTPCCls = kTRUE;};
    void SetNClsITS(int nCls){fnITSCls = nCls; fcutnITSCls = kTRUE;};
    void SetDCAReCalculation(bool which){fDCAProp = which;};
    void SetDCAVtxXY(double dcaXY){fDCAToVertexXY = dcaXY; fCutDCAToVtxXY = kTRUE;};
    void SetCutDCAVtxXY(bool cutdcaXY){fCutDCAToVtxXY = cutdcaXY;};
    void SetDCAVtxZ(double dcaZ){fDCAToVertexZ = dcaZ; fCutDCAToVtxZ = kTRUE;};
    void SetCutDCAVtxZ(bool cutdcaZ){fCutDCAToVtxZ = cutdcaZ;};
    void SetCutSharedCls(bool cutit){fCutSharedCls = cutit;};
    void SetCheckTPCRefit(bool cutit){fCheckTPCRefit=cutit;};
    void SetCutTPCCrossedRows(bool cutit){fCutTPCCrossedRows = cutit;};
    void SetTPCCrossedRowsCut(double CrossedRows){fTPCCrossedRowsCut = CrossedRows;};
    void SetPID(AliPID::EParticleType pid, double pTPChresh, double sigValTPC, double sigValTOF)
    {fParticleID = pid; fPIDPTPCThreshold = pTPChresh; fNSigValueTPC = sigValTPC; fNSigValueTOF = sigValTOF; fCutPID = kTRUE;};
    void SetCutITSPID(double sigValITSmin, double sigValITSmax, bool cutit){fNSigValueITSmin = sigValITSmin,fNSigValueITSmax = sigValITSmax, fdoITSnSigmaCut = cutit;};
    void SetRapidityRange(double Ymin, double Ymax){fRapMin = Ymin; fRapMax = Ymax; fCutRapidity = kTRUE;};
    void SetRejLowPtPionsTOF(bool use){fRejectPions = use;};
    void SetCutSmallestSig(bool cutit){fCutHighPtSig = cutit;};
    void SetMassCut_ForDCA(double minMass, double maxMass){fMinMass = minMass; fMaxMass = maxMass;};
    //selection Methods
    bool isSelected(AliLightNTrack *Track);
    void BookQA(AliLightNTrack *Track);
    void BookMC(AliLightNTrack *Track);
    //  void FillSharedClusterQA(AliLightNTrack *Track);
    //Histogram things
    void Init();
    TList *GetQAHists() {return fHists->GetHistList();};
    TList *GetMCQAHists() {return fMCHists->GetHistList();};
    TString ClassName() {return "AliLightNTrackCuts";};
    void SetName(TString name){fHists->SetName(name.Data());};
    void FillStackGenerated(float p) {if(fMCHists)fMCHists->FillStackGen(p);};
    void FillStackGeneratedPrimary(float p) {if(fMCHists)fMCHists->FillStackGenPrimary(p);};
private:
    bool TrackingCuts(AliLightNTrack *Track);
    bool TPCPIDAODCuts(AliLightNTrack *Track);
    bool ITSPIDAODCuts(AliLightNTrack *Track);
    bool PIDAODCuts(AliLightNTrack *Track);
    bool SmallestNSig(AliLightNTrack *Track);
    bool DCACuts(AliLightNTrack *Track);
    bool MassCut_ForDCA(AliLightNTrack *Track);
    void BookTrackCuts();
    void FillMCContributions(AliLightNTrack *Track);
    AliLightNTrackMCHist *fMCHists; //!
    AliLightNTrackHist *fHists;     //!
    bool fMCData;                       //
    bool fDCAPlots;                     //
    bool fCombSigma;                    //
    bool f3DPlotPmass2dca;                    //
    bool fContribSplitting;             //
    bool fFillQALater;                  //
    bool fCheckFilterBit;               //
    bool fCheckPileUpITS;               //
    bool fCheckPileUpTOF;               //
    bool fCheckPileUp;                  //  Should only be used for Daughters of v0s
    UInt_t fFilterBit;                  //
    double fpTmin;                      //
    double fpTmax;                      //
    bool fcutPt;                        //
    double fetamin;                     //
    double fetamax;                     //
    bool fCutChi2;                      //
    double fChi2perNDF;                 //
    double fTPCRatio;                   //
    bool fcutEta;			      //
    double fRapMin;		      //
    double fRapMax;		      //
    bool fCutRapidity;		      //
    bool fcutCharge;                    //
    int fCharge;                        //
    int fnTPCCls;                       //
    bool fcutnTPCCls;                   //
    int fnITSCls;                       //
    bool fcutnITSCls;                   //
    bool fDCAProp;                      //  kTRUE means that the DCA gets recalculated by PropagateToDCA, kFALSE just uses the info stored in the AOD
    double fDCAToVertexXY;              //
    bool fCutDCAToVtxXY;                //
    bool fdoITSnSigmaCut;               //
    double fDCAToVertexZ;               //
    bool fCutDCAToVtxZ;                 //
    double fMinMass;                    //
    double fMaxMass;                    //
    bool fCutSharedCls;                 //
    bool fCheckTPCRefit;                //
    bool fCutTPCCrossedRows;            //
    double fTPCCrossedRowsCut;            //
    bool fCutPID;                       //
    bool fCutHighPtSig;                 // Reject tracks which have a lower Sigma for other particles (implemented for electrons, pion, kaons and protons)
    AliPID::EParticleType fParticleID;  //
    double fNSigValueTPC;               // defaults to 3
    double fNSigValueTOF;
    double fNSigValueITSmin;
    double fNSigValueITSmax;
    double fPIDPTPCThreshold;           // defaults to 0
    bool fRejectPions;                  // Supress Pions at low pT with the TOF, if information is available
    ClassDef(AliLightNTrackCuts,1);
};

#endif /* ALILIGHTNTRACKCUTS_H */
