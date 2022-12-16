#ifndef ALIJEQNVECTORHANDLER_H
#define ALIJEQNVECTORHANDLER_H

/* Copyright(c) 1998-2019, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//*************************************************************************************
// \class AliJEQnVectorHandler
// \brief helper class to handle the Qn-vectors computed with different calibrations
// \adapted from HF task
// \authors:
// C. Beattie
// M. Sas
///////////////////////////////////////////////////////////////////////////////////////

#include <TMath.h>
#include <TFile.h>
#include <TBits.h>
#include <TH1D.h>
#include <TString.h>

#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliOADBContainer.h"
#include "AliAODVZERO.h"
#include "AliAnalysisTaskFlowVectorCorrections.h"
#include "AliQnCorrectionsManager.h"



class AliJEQnVectorHandler : public TObject
{
  public:
  
    enum QvecNorm {
        kQoverQlength,
        kQoverM,
        kQoverSqrtM,
        kNone  
    };

    enum CalibType {
        kQnCalib,
        kQnFrameworkCalib
    };

    AliJEQnVectorHandler();
    AliJEQnVectorHandler(int calibType, int normMeth, int harmonic, TString OADBfileName);

    virtual ~AliJEQnVectorHandler();

    bool SetAODEvent(AliAODEvent* event); 
    void ResetAODEvent(); 

    void SetHarmonic(int harmonic) {fHarmonic = harmonic;}
    void SetCalibrationType(int calibType) {fCalibType = calibType;}
    void SetNormalisationMethod(int normMeth) {fNormMethod = normMeth;}
    void SetTPCEtaLimits(double absetamin=0., double absetamx=0.8) {fEtaMinTPC=absetamin; fEtaMaxTPC=absetamx;}
    void SetTPCPtLimits(double ptmin=0.2, double ptmax=5) {fPtMinTPC=ptmin; fPtMaxTPC=ptmax;}
    void SetFractionOfTPCtracksToUse(double fracToKeep) {fFractionOfTracksForQnTPC = fracToKeep;}
    void SetCalibrationsOADBFileName(TString OADBfileName) {fOADBFileName = OADBfileName;}
    bool LoadOADBCalibrations();
    bool ComputeCalibratedQnVectorTPC(bool forceCalc = false);
    bool ComputeCalibratedQnVectorV0(bool forceCalc = false);

    int GetHarmonic() const {return fHarmonic;}
    int GetCalibrationType() const {return fCalibType;}
    int GetNormalisationMethod() const {return fNormMethod;}
    TString GetCalibrationsOADBFileName() const {return fOADBFileName;}

    void GetQnVecTPC(double QnVecFullTPC[2], double QnVecPosTPC[2], double QnVecNegTPC[2]);
    void GetQnVecV0(double QnVecFullV0[2], double QnVecV0A[2], double QnVecV0C[2]);
    void GetUnNormQnVecTPC(double QnVecFullTPC[2], double QnVecPosTPC[2], double QnVecNegTPC[2]);
    void GetUnNormQnVecV0(double QnVecFullV0[2], double QnVecV0A[2], double QnVecV0C[2]);

    void GetqnTPC(double &qnFullTPC, double &qnPosTPC, double &qnNegTPC);
    void GetqnV0(double &qnFullV0, double &qnV0A, double &qnV0C);

    void GetMultQnVecTPC(double &MultFullTPC, double &MultPosTPC, double &MultNegTPC);
    void GetMultQnVecV0(double &MultFullV0, double &MultV0A, double &MultV0C);

    void GetEventPlaneAngleTPC(double &EvPlaneFullTPC, double &EvPlanePosTPC, double &EvPlaneNegTPC);
    void GetEventPlaneAngleV0(double &EvPlaneFullV0, double &EvPlaneV0A, double &EvPlaneV0C);

    void RemoveTracksFromQnTPC(std::vector<AliAODTrack*> trToRem, double QnVecFullTPC[2], double QnVecPosTPC[2], double QnVecNegTPC[2], double &multFullTPC, double &multPosTPC, double &multNegTPC, bool getUnNormalised=false); 

    void EnablePhiDistrHistos();
    TH2F* GetPhiDistrHistosTPCPosEta() const {return fPhiVsCentrTPC[0];}
    TH2F* GetPhiDistrHistosTPCNegEta() const {return fPhiVsCentrTPC[1];}

  private:  

    //utils methods
    void ComputeQvecQnFrameworkTPC();
    void ComputeQvecQnFrameworkV0();
    void ComputeQvecTPC();
    void ComputeQvecV0();
    double ComputeEventPlaneAngle(double Qx, double Qy) const {return (TMath::Pi()+TMath::ATan2(-Qy,-Qx))/fHarmonic;}
    
    short GetVertexZbin() const;
    short GetCentBin() const;
    bool OpenInfoCalbration();

    bool IsTrackSelected(AliAODTrack* track);

    //data members
    unsigned int fHarmonic;                                    /// Harmonic of Q vector
    int fNormMethod;                                           /// Option for Q-vector normalisation
    int fCalibType;                                            /// Type of Q-vector calibrations to be applied
    
    double fEtaMinTPC;                                         /// Absolute minimum value of eta for TPC Q vector
    double fEtaMaxTPC;                                         /// Absolute maximum value of eta for TPC Q vector
    double fPtMinTPC;                                          /// Minimum value of pt for TPC Q vector
    double fPtMaxTPC;                                          /// Maximum value of pt for TPC Q vector

    double fQnVecFullTPC[2];                                   /// Qn vector of tracks in full TPC volume
    double fQnVecPosTPC[2];                                    /// Qn vector of tracks with positive pseudorapidity in TPC 
    double fQnVecNegTPC[2];                                    /// Qn vector of tracks with negative pseudorapidity in TPC 

    double fQnVecFullV0[2];                                    /// Qn vector of V0A+V0C sectors 
    double fQnVecV0A[2];                                       /// Qn vector of V0A sectors 
    double fQnVecV0C[2];                                       /// Qn vector of V0C sectors 

    double fMultFullTPC;                                       /// Multiplicity for Q vector with full TPC volume
    double fMultPosTPC;                                        /// Multiplicity for Q vector with positive pseudorapidity in TPC
    double fMultNegTPC;                                        /// Multiplicity for Q vector with negative pseudorapidity in TPC 

    double fMultFullV0;                                        /// Multiplicity for Q vector with V0A+V0C sectors
    double fMultV0A;                                           /// Multiplicity for Q vector with V0A sectors
    double fMultV0C;                                           /// Multiplicity for Q vector with V0C sectors

    double fQnVecNormFullTPC;                                  /// Norm of Q vector with full TPC volume before normalisation
    double fQnVecNormPosTPC;                                   /// Norm of Q vector with positive pseudorapidity in TPC before normalisation
    double fQnVecNormNegTPC;                                   /// Norm of Q vector with negative pseudorapidity in TPC before normalisation 

    double fQnVecNormFullV0;                                   /// Norm of Q vector with V0A+V0C sectors before normalisation
    double fQnVecNormV0A;                                      /// Norm of Q vector with V0A sectors before normalisation
    double fQnVecNormV0C;                                      /// Norm of Q vector with V0C sectors before normalisation

    TBits fUsedTrackPosIDs;                                    /// IDs of tracks (with positive ID) used for the computation of the Q vector (TPC)
    TBits fUsedTrackNegIDs;                                    /// IDs of tracks (with negative ID) used for the computation of the Q vector (TPC)

    //random rejection 
    double fFractionOfTracksForQnTPC;                          /// Fraction of tracks to keep when rejection enabled

    //data members needed for calibrations
    AliAODEvent* fAODEvent;                                    /// AOD event
    AliAODVZERO* fV0;                                          /// V0 info for the considered event
    int fRun;                                                  /// Run number
    double fZvtx;                                              /// Primary vertx Z 
    double fCentrality;                                        /// Event centrality

    TFile* fOADBFile;                                          /// OADB input file
    TString fOADBFileName;                                     /// OADB input file name
    bool fIsOADBFileOpen;                                      /// Flag to test whether the OADB file is open

    // OADB Objects for streaming
    AliOADBContainer * fMultV0BefCorPfpx;                      /// OADB object containing the hMultV0BefCorPfpx histograms

    TObjArray *fOADBzArray_contQx2am;           /// Array of OADB contQx2am object, index is z-vertex bin
    TObjArray *fOADBzArray_contQy2am;           /// Array of OADB contQy2am object, index is z-vertex bin
    TObjArray *fOADBzArray_contQx2as;           /// Array of OADB contQx2as object, index is z-vertex bin
    TObjArray *fOADBzArray_contQy2as;           /// Array of OADB contQy2as object, index is z-vertex bin

    TObjArray *fOADBzArray_contQx2cm;           /// Array of OADB contQx2cm object, index is z-vertex bin
    TObjArray *fOADBzArray_contQy2cm;           /// Array of OADB contQy2cm object, index is z-vertex bin
    TObjArray *fOADBzArray_contQx2cs;           /// Array of OADB contQx2cs object, index is z-vertex bin
    TObjArray *fOADBzArray_contQy2cs;           /// Array of OADB contQy2cs object, index is z-vertex bin

    TObjArray *fOADBcentArray_contTPCposEta;    /// Array of OADB contTPCposEta, index is cent bin
    TObjArray *fOADBcentArray_contTPCnegEta;    /// Array of OADB contTPCnegEta, index is cent bin

    int fCalibObjRun;                                          /// Run of loaded calibration objects

    TH1D* fHistMultV0;                                         /// Profile from V0 multiplicity

    TH1D* fQx2mV0A[14];                                        /// <Qxn> V0A
    TH1D* fQy2mV0A[14];                                        /// <Qyn> V0A
    TH1D* fQx2sV0A[14];                                        /// sigma Qxn V0A
    TH1D* fQy2sV0A[14];                                        /// sigma Qyn V0A
    
    TH1D* fQx2mV0C[14];                                        /// <Qxn> V0C
    TH1D* fQy2mV0C[14];                                        /// <Qyn> V0C
    TH1D* fQx2sV0C[14];                                        /// sigma Qxn V0C
    TH1D* fQy2sV0C[14];                                        /// sigma Qyn V0C

    bool fV0CalibZvtxDiff;                                     /// flag to properly manage Zvtx differential V0 calibrations

    TH1D* fWeightsTPCPosEta[9];                                /// Weights for TPC tracks with eta > 0
    TH1D* fWeightsTPCNegEta[9];                                /// Weights for TPC tracks with eta < 0
    bool fEnablePhiDistrHistos;                                /// Enable phi distribution histos
    TH2F* fPhiVsCentrTPC[2];                                   /// Phi vs. centr TH2 of selected TPC tracks in eta>0 and eta<0

    //data members needed for Qn-framework calibrations
    AliAnalysisTaskFlowVectorCorrections *fQnVectorTask;       /// Qn-framework task
    AliQnCorrectionsManager *fQnVectorMgr;                     /// Qn-framework manager

  /// \cond CLASSIMP
  ClassDef(AliJEQnVectorHandler, 6); ///
  /// \endcond
};

#endif
