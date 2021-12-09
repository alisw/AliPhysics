#ifndef ALIANALYSISTASKJETQNVECTORS_H
#define ALIANALYSISTASKJETQNVECTORS_H

/* Copyright(c) 1998-2019, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//*************************************************************************************************
// \class AliAnalysisTaskJetQnVectors
// \brief task used to load the Qn calibrations and get the calibrated Qn vectors for HF analyses
// \authors:
// F. Grosa, fabrizio.grosa@cern.ch
// F. Catalano, fabio.catalano@cern.ch
// A. Dobrin, alexandru.dobrin@cern.ch
// A. Festanti, andrea.festanti@cern.ch
// G. Luparello, grazia.luparello@cern.ch
// F. Prino, prino@to.infn.it
// A. Rossi, andrea.rossi@cern.ch
// S. Trogolo, stefano.trogolo@cern.ch
///////////////////////////////////////////////////////////////////////////////////////////////////

#include <TH1F.h>
#include <TList.h>
#include <TString.h>
#include <TMath.h>

#include "AliAnalysisTaskSE.h"
#include "AliVEvent.h"
#include "AliAODEvent.h"
#include "AliJEQnVectorHandler.h"


class AliAnalysisTaskJetQnVectors : public AliAnalysisTaskSE {

public:

    enum Det{kFullTPC,kPosTPC,kNegTPC,kFullV0,kV0A,kV0C};

    AliAnalysisTaskJetQnVectors();
    AliAnalysisTaskJetQnVectors(const char *name, int harmonic, int calibType, TString oadbFileName);
    virtual ~AliAnalysisTaskJetQnVectors();

    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t *option);

    AliJEQnVectorHandler* GetQnVectorHandler() const                                                     {return fHFQnVecHandler;}
    TList* GetSplineForqnPercentileList(int det=kFullTPC) const;
    void SetUseAODBCalibrations(TString oadbFileName)                                                    {fOADBFileName = oadbFileName; fCalibType = AliJEQnVectorHandler::kQnCalib;}
    void SetUseQnFrameworkCalibrations()                                                                 {fCalibType = AliJEQnVectorHandler::kQnFrameworkCalib;}
    void SetNormalisationMethod(int normmethod)                                                          {fNormMethod = normmethod;}
    void SetTriggerInfo(TString trigClass, unsigned long long mask=0)                                    {fTriggerClass = trigClass; fTriggerMask = mask;}
    void LoadSplinesForqnPercentile(TString splinesfilepath);

    void EnableTPCPhiVsCentrDistrHistosVsRun()                                                           {fEnableTPCPhiVsCentrDistr=true;}
    void EnableQVecTPCVsCentrDistrHistosVsRun()                                                          {fEnableQvecTPCVsCentrDistr=true;}

    double Getq2V0()                                                                                     {return fq2V0;}

private:

    double GetEventPlaneAngle(double Qx, double Qy) const {return (TMath::Pi()+TMath::ATan2(-Qy,-Qx))/2;}

    TList *fOutputList;                              //!<! output list for histograms

    TH1F *fHistNEvents;                              //!<! histo with number of events
    TH1F *fHistCentrality;                           //!<! histo with centrality
    TH1F *fHistEventPlaneTPC[3];                     //!<! histos of TPC (Full, PosEta, NegEta) EP angle
    TH1F *fHistEventPlaneV0[3];                      //!<! histos of V0 (Full, V0A, V0C) EP angle
    TH2F *fHistqnVsCentrTPC[3];                    	 //!<! histos of q2TPC (Full, PosEta, NegEta) vs centrality (for spline calibration)
    TH2F *fHistqnVsCentrV0[3];                       //!<! histos of q2V0 (Full, V0A, V0C) vs centrality (for spline calibration)

    TH2F* fTPCPhiVsCentrDistr[2];                    //!<! histos of phi vs. centr of selected TPC tracks in eta>0 and eta<0
    TH2F* fQvecTPCVsCentrDistr[3];                   //!<! histos of TPC Q-vector vs. centr for tracks with eta>0 and eta<0

    bool fEnableTPCPhiVsCentrDistr;                  /// flag to enable histos of phi vs. centr
    bool fEnableQvecTPCVsCentrDistr;                 /// flag to enable histos of TPC Q-vector vs. centr

    AliJEQnVectorHandler* fHFQnVecHandler;           /// Qn-vector handler
    int fHarmonic;                                   /// Qn-vector harmonic
    int fCalibType;                                  /// type of calibrations used by handler
    int fNormMethod;                                 /// normalisation of Q vector

    TString fOADBFileName;                           /// OADB input file name

    TList* fSplineListqnPercTPC[3];                  /// Splines for qn percentile calibration for TPC
    TList* fSplineListqnPercV0[3];                   /// Splines for qn percentile calibration for V0

    AliAODEvent* fAOD;                               /// AOD event
    int fPrevEventRun;                               /// run number of event previously analysed

    TString fTriggerClass;                           /// trigger class
    unsigned long long fTriggerMask;                 /// trigger mask

    double fq2V0;                                    /// q2 vector from the V0     

    ClassDef(AliAnalysisTaskJetQnVectors, 4);
};

#endif
