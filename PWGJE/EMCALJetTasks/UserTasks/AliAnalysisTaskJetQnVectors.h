#ifndef ALIANALYSISTASKJETQNVECTORS_H
#define ALIANALYSISTASKJETQNVECTORS_H

/* Copyright(c) 1998-2019, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//*************************************************************************************************
// \class AliAnalysisTaskJetQnVectors
// \brief task used to load the Qn calibrations and get the calibrated Qn vectors for JE analyses
// \authors:
// C. Beattie, caitlin.beattie@cern.ch
// M. Sas,     mike.sas@cern.ch
///////////////////////////////////////////////////////////////////////////////////////////////////

#include <TH1F.h>
#include <TList.h>
#include <TString.h>
#include <TMath.h>

#include "AliAnalysisTaskSE.h"
#include "AliVEvent.h"
#include "AliAODEvent.h"
#include "AliJEQnVectorHandler.h"
#include "AliEventCuts.h"

class AliAnalysisTaskJetQnVectors : public AliAnalysisTaskSE {

public:

    enum Det{kFullTPC,kPosTPC,kNegTPC,kFullV0,kV0A,kV0C};

    AliAnalysisTaskJetQnVectors();
    AliAnalysisTaskJetQnVectors(const char *name, int harmonic, int calibType, TString oadbFileName1, TString oadbFileName2);
    virtual ~AliAnalysisTaskJetQnVectors();

    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t *option);

    void CreateQnVectorHandlers(); // Create the QnVector handlers, including loading the calibration files
    TDirectoryFile*  GetSplineForqnPercentileList(int det=kFullTPC) const;
    void SetUseQnFrameworkCalibrations()                                                                 {fCalibType = AliJEQnVectorHandler::kQnFrameworkCalib;}
    void SetNormalisationMethod(int normmethod)                                                          {fNormMethod = normmethod;}
    void SetTriggerInfo(TString trigClass, unsigned long long mask=0)                                    {fTriggerClass = trigClass; fTriggerMask = mask;}
    void SetRejectTPCPileup(bool val)                                                                    {fRejectTPCPileup = val;}
    void LoadSplinesForqnPercentile(TString splinesfilepath);

    void EnableTPCPhiVsCentrDistrHistosVsRun()                                                           {fEnableTPCPhiVsCentrDistr=true;}
    void EnableQVecTPCVsCentrDistrHistosVsRun()                                                          {fEnableQvecTPCVsCentrDistr=true;}

    double Getq2V0M()                                                                                    {return fq2V0M;}
    double Getq2V0A()                                                                                    {return fq2V0A;}
    double Getq2V0C()                                                                                    {return fq2V0C;}
    double Getq2TPC()                                                                                    {return fq2TPC;}
    double GetEPangleFullTPC()                                                                           {return fEPangleFullTPC;}
    double GetEPanglePosTPC()                                                                            {return fEPanglePosTPC;}
    double GetEPangleNegTPC()                                                                            {return fEPangleNegTPC;}
    double GetEPangleV0M()                                                                               {return fEPangleV0M;}
    double GetEPangleV0A()                                                                               {return fEPangleV0A;}
    double GetEPangleV0C()                                                                               {return fEPangleV0C;}

protected:
    AliEventCuts fEventCuts;                         //!<! Event selection  

private:

    double GetEventPlaneAngle(double Qx, double Qy) const {return (TMath::Pi()+TMath::ATan2(-Qy,-Qx))/2;}

    TList *fOutputList;                              //!<! output list for histograms

    TH1F *fHistNEvents;                              //!<! histo with number of events
    TH1F *fHistCentrality;                           //!<! histo with centrality
    TH3F *fHistResolution_epV0AV0C_qV0M;             //!<! histo with detector resolution
    TH3F *fHistResolution_epV0CTPC_qV0M;             //!<! histo with detector resolution
    TH3F *fHistResolution_epV0ATPC_qV0M;             //!<! histo with detector resolution
    TH3F *fHistResolution_epTPCpTPCn_qV0M;           //!<! histo with detector resolution
    TH3F *fHistResolution_epV0MTPCp_qV0M;            //!<! histo with detector resolution
    TH3F *fHistResolution_epV0MTPCn_qV0M;            //!<! histo with detector resolution
    TH3F *fHistResolution_epV0AV0C_qV0A;             //!<! histo with detector resolution
    TH3F *fHistResolution_epV0CTPC_qV0A;             //!<! histo with detector resolution
    TH3F *fHistResolution_epV0ATPC_qV0A;             //!<! histo with detector resolution
    TH3F *fHistResolution_epV0AV0C_qV0C;             //!<! histo with detector resolution
    TH3F *fHistResolution_epV0CTPC_qV0C;             //!<! histo with detector resolution
    TH3F *fHistResolution_epV0ATPC_qV0C;             //!<! histo with detector resolution
    TH1F *fHistEventPlaneTPC[3];                     //!<! histos of TPC (Full, PosEta, NegEta) EP angle
    TH1F *fHistEventPlaneV0[3];                      //!<! histos of V0 (Full, V0A, V0C) EP angle
    TH2F *fHistqnVsCentrTPC[3];                      //!<! histos of q2TPC (Full, PosEta, NegEta) vs centrality (for spline calibration)
    TH2F *fHistqnVsCentrV0[3];                       //!<! histos of q2V0 (Full, V0A, V0C) vs centrality (for spline calibration)

    TH2F* fTPCPhiVsCentrDistr[2];                    //!<! histos of phi vs. centr of selected TPC tracks in eta>0 and eta<0
    TH2F* fQvecTPCVsCentrDistr[3];                   //!<! histos of TPC Q-vector vs. centr for tracks with eta>0 and eta<0


    bool fEnableTPCPhiVsCentrDistr;                  /// flag to enable histos of phi vs. centr
    bool fEnableQvecTPCVsCentrDistr;                 /// flag to enable histos of TPC Q-vector vs. centr

    AliJEQnVectorHandler* fJEQnVecHandler1;          /// Qn-vector handler
    AliJEQnVectorHandler* fJEQnVecHandler2;          /// Qn-vector handler
    int fHarmonic;                                   /// Qn-vector harmonic
    int fCalibType;                                  /// type of calibrations used by handler
    int fNormMethod;                                 /// normalisation of Q vector

    TString fOADBFileName1;                           /// OADB input file name
    TString fOADBFileName2;                           /// OADB input file name

    TDirectoryFile* fSplineListqnPercTPC[3];         /// Splines for qn percentile calibration for TPC
    TDirectoryFile* fSplineListqnPercV0[3];          /// Splines for qn percentile calibration for V0

    AliAODEvent* fAOD;                               /// AOD event
    int fPrevEventRun;                               /// run number of event previously analysed

    TString fTriggerClass;                           /// trigger class
    unsigned long long fTriggerMask;                 /// trigger mask
    bool fRejectTPCPileup;                           /// TPC pileup rejection

    double fq2V0M;                                   /// q2 vector from the V0M   
    double fq2V0A;                                   /// q2 vector from the V0A    
    double fq2V0C;                                   /// q2 vector from the V0C     
    double fq2TPC;                                   /// q2 vector from the TPC
    double fEPangleFullTPC;                          /// EP Angle with calibrations from TPC Full
    double fEPanglePosTPC;                           /// EP Angle with calibrations from TPC eta>0
    double fEPangleNegTPC;                           /// EP Angle with calibrations from TPC eta<0
    double fEPangleV0M;                              /// EP Angle with calibrations from V0M
    double fEPangleV0C;                              /// EP Angle with calibrations from V0A
    double fEPangleV0A;                              /// EP Angle with calibrations from V0C

    ClassDef(AliAnalysisTaskJetQnVectors, 12);
};

#endif
