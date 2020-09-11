#ifndef ALIANALYSISTASKSECHARMHADRONVN_H
#define ALIANALYSISTASKSECHARMHADRONVN_H

/* Copyright(c) 1998-2019, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//*************************************************************************************
// \class AliAnalysisTaskSECharmHadronvn
// \brief task for the analysis of D-meson vn
// \authors:
// F. Grosa, fabrizio.grosa@cern.ch
// F. Catalano, fabio.catalano@cern.ch
// A. Festanti, andrea.festanti@cern.ch
// G. Luparello, grazia.luparello@cern.ch
// F. Prino, prino@to.infn.it
// A. Rossi, andrea.rossi@cern.ch
// S. Trogolo, stefano.trogolo@cern.ch
///////////////////////////////////////////////////////////////////////////////////////

#include <THnSparse.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TString.h>

#include "AliAnalysisTaskSE.h"
#include "AliRDHFCuts.h"
#include "AliAODEvent.h"
#include "AliHFQnVectorHandler.h"
#include "AliAnalysisVertexingHF.h"
#include "AliHFMLResponse.h"

class AliAnalysisTaskSECharmHadronvn : public AliAnalysisTaskSE
{
    public:

    enum DecChannel{kDplustoKpipi,kD0toKpi,kDstartoKpipi,kDstoKKpi,kLctopK0S}; //more particles can be added
    enum EventPlaneMeth{kTPC,kTPCVZERO,kVZERO,kVZEROA,kVZEROC,kPosTPCVZERO,kNegTPCVZERO}; //Event plane to be calculated in the task
    enum FlowMethod{kEP,kSP,kEvShapeEP,kEvShapeSP,kEPVsMass,kEvShapeEPVsMass}; // Event Plane, Scalar Product or Event Shape Engeneering methods
    enum q2Method{kq2TPC,kq2PosTPC,kq2NegTPC,kq2VZERO,kq2VZEROA,kq2VZEROC}; // qn for Event Shape to be calculated in the task
    enum EventPlaneDet{kNone=-1,kFullTPC,kPosTPC,kNegTPC,kFullV0,kV0A,kV0C};

    AliAnalysisTaskSECharmHadronvn();
    AliAnalysisTaskSECharmHadronvn(const char *name, AliRDHFCuts *rdCuts, int decaychannel);

    virtual ~AliAnalysisTaskSECharmHadronvn();

    void SetTenderTaskName(TString name)                              {fTenderTaskName = name;}
    void SetAODMismatchProtection(int opt=1)                          {fAODProtection = opt;}

    void SetHarmonic(int harmonic)                                    {fHarmonic = harmonic;}
    void SetCalibrationType(int calibtype)                            {fCalibType = calibtype;}
    void SetNormMethod(int norm)                                      {fNormMethod = norm;}
    void SetOADBFileName(TString filename)                            {fOADBFileName = filename;}

    void SetFlowMethod(int meth)                                      {fFlowMethod = meth;}
    void SetEventPlaneDetector(int det)                               {fEvPlaneDet = det;}
    void SetSubEventDetectors(int detsubA, int detsubB)               {fSubEvDetA = detsubA; fSubEvDetB = detsubB;}
    void SetQnVectorDetConf(int detconf);
    void SetTPCEPOnly()                                               {SetQnVectorDetConf(kTPC);}
    void SetVZEROEP()                                                 {SetQnVectorDetConf(kVZERO);}
    void SetVZEROAEP()                                                {SetQnVectorDetConf(kVZEROA);}
    void SetVZEROCEP()                                                {SetQnVectorDetConf(kVZEROC);}
    void SetTPCEP()                                                   {SetQnVectorDetConf(kTPCVZERO);}
    void SetqnMethod(int qnmethod)                                    {fqnMeth = qnmethod;}
    void SetScalProdLimit(double limit=0.3)                           {limit < 1. ? fScalProdLimit = limit : fScalProdLimit = 1.;}
    void SetqnPercentileSelection(TString splinesfilepath)            {fPercentileqn=true; fqnSplineFileName=splinesfilepath;}

    void SetMinCentrality(float mincentr)                             {fMinCentr = mincentr;}
    void SetMaxCentrality(float maxcentr)                             {fMaxCentr = maxcentr;}

    void SetDecayChannel(int decaychannel)                            {fDecChannel = decaychannel;}
    void SetMassLimits(float range,int pdg);
    void SetMassLimits(float lowlimit, float uplimit);
    void SetNMassBins(int nbins)                                      {fNMassBins = nbins;}
    float GetUpperMassLimit() const                                   {return fUpmasslimit;}
    float GetLowerMassLimit() const                                   {return fLowmasslimit;}
    int GetNMassBins() const                                          {return fNMassBins;}

    void SetTPCHalvesEtaGap(double etagap = 0.2)                      {fEtaGapInTPCHalves=etagap;}
    void RemoveDauTracksFromqn(int removedau=1, bool remsoftpi=false) {fRemoveDauFromqn=removedau; fRemoveSoftPion=remsoftpi;}
    void SetRandomDownsamplFromqn(double fractokeep = 0.5)            {fEnableDownsamplqn=true; fFracToKeepDownSamplqn=fractokeep;}

    // methods for ML application
    void SetDoMLApplication(bool flag = kTRUE) {fApplyML = flag;}
    void SetMLConfigFile(TString path = ""){fConfigPath = path;}
    void SetMLBinsForSparse(int nbins = 300, double min = 0.85, double max = 1.) { fNMLBins = nbins; fMLOutputMin = min; fMLOutputMax = max;}

    // Implementation of interface methods
    virtual void UserCreateOutputObjects();
    virtual void LocalInit();
    virtual void UserExec(Option_t *option);

    private:

    double GetPhiInRange(double phi);
    double GetDeltaPsiSubInRange(double psi1, double psi2);
    void CalculateInvMasses(AliAODRecoDecayHF* d,float* &masses,int& nmasses);
    void GetMainQnVectorInfo(double &mainPsin, double &mainMultQn, double mainQn[2], double &SubAPsin, double &SubAMultQn, double SubAQn[2], double &SubBPsin, double &SubBMultQn, double SubBQn[2], AliHFQnVectorHandler* HFQnVectorHandler);
    void GetDaughterTracksToRemove(AliAODRecoDecayHF* d, int nDau, vector<AliAODTrack*> &trackstoremove);
    int IsCandidateSelected(AliAODRecoDecayHF *&d, int nDau, int absPdgMom, AliAnalysisVertexingHF *vHF, AliAODRecoDecayHF2Prong *dD0, double modelPred[2]);
    bool LoadSplinesForqnPercentile();

    static const int kVarForSparse = 10;

    AliAODEvent* fAOD;                      /// AOD event

    TList* fOutput;                         //!<! list send on output slot 2
    TH1F* fHistNEvents;                     //!<! histogram send on output slot 1
    TH1F* fHistCentrality[3];               //!<! hist. for cent distr (all,sel ev,out of cent)
    TH2F* fHistCandVsCent;                  //!<! hist. for number of selected candidates vs. centrality
    TH3F* fHistEvPlaneQncorr[6];            //!<! histograms for EP angle vs. centrality vs. qn
    TH2F* fHistqnVsCentrPercCalib[6];       //!<! histograms for qn percentile calibration with qn fine binning
    TH3F* fHistEPResolVsCentrVsqn[3];       //!<! histos needed to compute EP resolution vs. centrality vs. qn
    TH3F* fHistPercqnVsqnVsCentr;           //!<! histo of qn percentile vs. qn vs. centrality
    TH3F* fHistNtrklVsqnVsCentr;            //!<! histo of Ntracklets vs. qn vs. centrality
    THnSparseF* fHistMassPtPhiqnCentr;      //!<! THnSparse for the analysis of vn
    AliRDHFCuts* fRDCuts;                   /// cut values (saved in slot 3)

    TString fTenderTaskName;                /// name of tender task needed to get the calibrated Qn vectors

    float fMinCentr;                        /// minimum centrality
    float fMaxCentr;                        /// maximum centrality
    int fAODProtection;                     /// flag to activate protection against AOD-dAOD mismatch.
                                            /// -1: no protection,  0: check AOD/dAOD nEvents only,  1: check AOD/dAOD nEvents + TProcessID names

    int fHarmonic;                          /// flow harmonic
    int fCalibType;                         /// calibration type
    int fNormMethod;                        /// Qn-vector normalisation method
    TString fOADBFileName;                  /// AODB file name for calibrations (if Qn-framework not used)

    int fFlowMethod;                        /// method used to compute vn
    int fEvPlaneDet;                        /// detector for event plane
    int fSubEvDetA;                         /// detector for 1st subevent
    int fSubEvDetB;                         /// detector for 2nd subevent
    int fqnMeth;                            /// flag to select qn method
    double fScalProdLimit;                  /// max value for the scalar product histograms
    bool fPercentileqn;                     /// flag to replace qn with its percentile in the histograms
    TString fqnSplineFileName;              /// qn spline file name for qn percentile calibrations
    bool fLoadedSplines;                    /// flag to know if splines loaded with task or read from tender task
    TList* fqnSplinesList[6];               /// lists of splines used to compute the qn percentile

    int fDecChannel;                        /// decay channel identifier
    float fLowmasslimit;                    /// lower inv mass limit for histos
    float fUpmasslimit;                     /// upper inv mass limit for histos
    int fNMassBins;                         /// number of bins in the mass histograms

    double fEtaGapInTPCHalves;              /// eta gap between TPC subevents
    int fRemoveDauFromqn;                   /// flag to enable removal of D-meson daughter tracks from qn (1->remove single cand, 2->remove all cand of the analysed event)
    bool fRemoveSoftPion;                   /// flag to enable removal of soft pion too (only D*)
    bool fEnableDownsamplqn;                /// flag to enable random downsampling for qn
    double fFracToKeepDownSamplqn;          /// fraction of tracks to keep in qn with random downsampling

    /// variables for ML application
    bool fApplyML;                          /// flag to enable ML application
    TString fConfigPath;                    /// path to ML config file
    AliHFMLResponse* fMLResponse;           //!<! object to handle ML response
    int fNMLBins;                           /// number of bins for ML axis in thnsparse
    double fMLOutputMin;                    /// minimum value of ML ouptut in thnsparse axis
    double fMLOutputMax;                    /// maximum value of ML ouptut in thnsparse axis

    ClassDef(AliAnalysisTaskSECharmHadronvn,6); // AliAnalysisTaskSE for the HF vn analysis
};

#endif
