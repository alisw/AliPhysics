
#ifndef AliAnalysisTaskRawJetWithEP_H
#define AliAnalysisTaskRawJetWithEP_H
/**
 * \file AliAnalysisTaskRawJetWithEP.h
 * \brief Declaration of class AliAnalysisTaskRawJetWithEP
 *
 * In this header file the class AliAnalysisTaskRawJetWithEP is declared.
 * This is a sample task that shows how to write a simple user analysis task
 * using the EMCal jet framework. It is also used to do automatic benchmark
 * tests of the software.
 *
 * \author Salvatore Aiola <salvatore.aiola@cern.ch>, Yale University
 * \date Apr 27, 2016
 */

/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TString.h>

#include "AliAnalysisTaskEmcalJet.h"
#include "THistManager.h"
#include "AliVEvent.h"
#include "AliAODEvent.h"
#include "AliEventCuts.h"
// #include "AliMultSelection.h"

#include "AliJEQnVectorHandler.h"
#include "AliAnalysisTaskJetQnVectors.h"
// #include "AliAnalysisTaskFlowVectorCorrections.h"
#include "AliQnCorrectionsManager.h"

class TClonesArray;
class TF1;

class AliJEQnVectorHandler;

class AliRhoParameter;
class AliLocalRhoParameter;


/**
 * \class AliAnalysisTaskRawJetWithEP
 * \brief Implementation of a sample jet analysis task.
 *
 * This class in an implementation of a sample task for EMCal jet analysis.
 * It derives from AliAnalysisTaskEmcalJet.
 * It performs a simple analysis, producing track, cluster and jet spectra.
 * It also performs a QA of the cluster-track matching.
 * Note: if jets are not used this class can be simplified by deriving
 * from AliAnalysisTaskEmcal and removing the functions DoJetLoop()
 * and AllocateJetHistograms().
 */
class AliAnalysisTaskRawJetWithEP : public AliAnalysisTaskEmcalJet {
  public:

    AliAnalysisTaskRawJetWithEP();
    AliAnalysisTaskRawJetWithEP(const char *name);
    virtual ~AliAnalysisTaskRawJetWithEP();

    void UserCreateOutputObjects();
    void Terminate(Option_t *option);

    // AliJEQnVectorHandler* GetQnVectorHandler() const {return fQ2VecHandler;}
    TList* GetSplineForqnPercentileList(int det=kFullTPC) const;
    // void SetUseAODBCalibrations(TString oadbFileName){fOADBFileName = oadbFileName; fCalibType = AliJEQnVectorHandler::kQnCalib;}
    // void SetUseQnFrameworkCalibrations(){fCalibType = AliJEQnVectorHandler::kQnFrameworkCalib;}

    static AliAnalysisTaskRawJetWithEP* AddTaskRawJetWithEP(
        const char *ntracks            = "usedefault",
        const char *nclusters          = "usedefault",
        const char* ncells             = "usedefault",
        const char *suffix             = "");

    enum runModeType{kLocal, kGrid};
    enum fitModulationType  { kNoFit, kV2, kV3, kCombined, kFourierSeries, \
                              kIntegratedFlow, kQC2, kQC4 }; // fit type
    
    enum DetFlev{kFullTPC,kPosTPC,kNegTPC,kFullV0,kV0A,kV0C};
    enum detectorType{ kTPC, kVZEROA, kVZEROC, kVZEROComb, kFixedEP};  // detector that was used for event plane
    enum qnVCalibType{kOrig, kJeHand};

    
    AliAnalysisTaskJetQnVectors* fV0Q2VectTask; ///< Reader for the Qn vector
    Double_t  fV0Q2Vector; ///< Calibrated q2 value 
    Double_t  fV0Ep2Angle; ///< Calibrated event-plane angle 

    //== s == Getter Prepare  ################################################
    TString GetLocalRhoName() const {return fLocalRhoName;}
    AliLocalRhoParameter* GetLocalRhoParameter() const {return fLocalRho;}
    //== e == Getter Prepare  ################################################
    
    //== s == Setter Prepare  ################################################
    void SetQnOadbFile(TString oadbFileName){fOADBFileName = oadbFileName;}
    void SetSplineFile(TString splineFileName){fSplinesFileName = splineFileName;}
    void SetCalibRefFile(TString calibRefFileName){fCalibRefFileName = calibRefFileName;}
    
    void SetModulationFitType(fitModulationType type) {fFitModulationType = type; }

    void SetEventQA(Int_t qaEventNum){fQaEventNum = qaEventNum;}
    void SetV0Combine(Bool_t bV0Combin){fV0Combin = bV0Combin;}
    void SetQnCalibType(qnVCalibType iQnVCalibType){fQnVCalibType = iQnVCalibType;}
    //== e == Setter Prepare  ################################################

  protected:
    AliEventCuts fEventCuts; //

    void       ExecOnce();
    Bool_t     Run();
    
    THistManager fHistManager;///< Histogram manager
    
    // ### protected #############################################################
  
  // ### private ###############################################################
  private:
    
    AliAODEvent* fAOD; /// AOD event
    TString fOADBFileName; /// OADB input file name
    TFile*  fOADBFile;
    TFile*  fCalibRefFile;
    TList*  fCalibRefObjList;
    TList*  fCalibV0Ref;
    Bool_t  fCalibQA;
    
    AliJEQnVectorHandler* fQ2VecHandler;/// Qn-vector handler
    AliJEQnVectorHandler* fQ3VecHandler;/// Qn-vector handler

    TList     *fOutputList; //!<! output list for histograms
    int       fCalibType;   /// type of calibrations used by handler
    int       fNormMethod;  /// normalisation of Q vector
    Bool_t    fV0Combin = kFALSE;
    
    TString   fSplinesFileName; /// Sprine input file name
    TString   fCalibRefFileName;

    void      LoadSpliForqnPerce();

    void       VzeroGainCalibQA();

    // AliMultSelection      *fMultSelection;      //! For Centrality 
    void       SetupPileUpRemovalFunctions();
    Bool_t     CheckEventIsPileUp2018(AliAODEvent *faod);

    void       AllocateEventPlaneHistograms();
    void       AllocateJetHistograms();
    void       AllocateTrackHistograms();
    
    void       MeasureTpcEPQA();
    
    void       SetModulationRhoFit();
    void       MeasureBkg();
    void       BkgFitEvaluation();
    
    void       DoEventPlane();
    void       DoJetLoop();
    void       DoTrackLoop();

    

    Bool_t     QnJEHandlarEPGet();
    Bool_t     QnGainCalibration();
    Bool_t     QnRecenteringCalibration();

    Double_t CalcEPAngle(double Qx, double Qy) const {return (TMath::Pi()+TMath::ATan2(-Qy,-Qx))/2;}
    Double_t CalcEPReso(Int_t n, Double_t &psiA, Double_t &psiB, Double_t &psiC);

    void  CalcRandomCone(Double_t &pt, Double_t &eta, Double_t &phi, \
      Double_t &leadingJetEta, Double_t &leadingJetPhi, Double_t &jetR
    ) const;

    TH1F*   GetResoFromOutputFile(detectorType det, Int_t h, TArrayD* cen);
    Double_t CalculateEventPlaneChi(Double_t res);

    static Double_t ChiSquarePDF(Int_t ndf, Double_t x) {
      Double_t n(ndf/2.), denom(TMath::Power(2, n)*TMath::Gamma(n));
      if (denom!=0)  return ((1./denom)*TMath::Power(x, n-1)*TMath::Exp(-x/2.)); 
      return -999; 
    }

    // note that the cdf of the chisquare distribution is the normalized lower incomplete gamma function
    static Double_t ChiSquareCDF(Int_t ndf, Double_t x) { return TMath::Gamma(ndf/2., x/2.); }

    static Double_t ChiSquare(TH1& histo, TF1* func) {
      // evaluate the chi2 using a poissonian error estimate on bins
      Double_t chi2(0.);
      for(Int_t i(0); i < histo.GetXaxis()->GetNbins(); i++) {
          if(histo.GetBinContent(i+1) <= 0.) continue;
          chi2 += TMath::Power((histo.GetBinContent(i+1) \
            - func->Eval(histo.GetXaxis()->GetBinCenter(1+i))), 2)\
              /histo.GetBinContent(i+1);
      }
      return chi2;
    }

    // static Double_t CalcEPChi(Double_t res)
    // {
    //   // return chi for given resolution to combine event plane estimates from two subevents
    //   // see Phys. Rev. C no. CS6346 (http://arxiv.org/abs/nucl-ex/9805001)
    //   Double_t chi(2.), delta(1.), con((TMath::Sqrt(TMath::Pi()))/(2.*TMath::Sqrt(2)));
    //   for (Int_t i(0); i < 15; i++) {
    //       chi = ((con*chi*TMath::Exp(-chi*chi/4.)*(TMath::BesselI0(chi*chi/4.)+TMath::BesselI1(chi*chi/4.))) < res) ? chi + delta : chi - delta;
    //       delta = delta / 2.;
    //   }
    //   return chi;
    // }
    
    

    qnVCalibType fQnVCalibType = kOrig;

    Int_t CheckRunNum = 0;
    Int_t fPreRunNum  = -1; /// run number of event previously analysed
    Int_t fRunOrder   = -1;
    
    /// Functions for Pile Up Event Removal:
    TF1   *fV0CutPU;      //!
    TF1   *fSPDCutPU;     //!
    TF1   *fMultCutPU;    //!
    TF1   *fCenCutLowPU;  //!
    TF1   *fCenCutHighPU; //!

    Int_t fQaEventNum = 0;
    Int_t chcekDrawBkg = 0;

    //qnVector 0:x, 1:y
    Double_t q2VecV0M[2]; //C+A
    Double_t q2VecV0C[2];
    Double_t q2VecV0A[2];
    Double_t q2V0[3];   //psi2 vzero 0:combin, 1:eta negative (C side), 2:eta positive (A side)
    Double_t q3VecV0M[2]; //C+A
    Double_t q3VecV0C[2];
    Double_t q3VecV0A[2];
    Double_t q3V0[3]; //psi2 vzero 0:combin, 1:eta negative (C side), 2:eta positive (A side)
    
    Double_t psi2V0[3]; //psi2 vzero 0:combin, 1:eta negative (C side), 2:eta positive (A side)
    Double_t psi3V0[3]; //psi3 vzero 0:combin, 1:eta negative (C side), 2:eta positive (A side)


    //qnVector 0:eta posi, 1:eta nega
    Double_t q2VecTpcM[2]; //posi + nega
    Double_t q2VecTpcP[2];
    Double_t q2VecTpcN[2];
    Double_t q2Tpc[3]; //psi2 TPC 0:eta negative (C side), 1:eta positive (A side)side)
    Double_t q3VecTpcM[2]; //posi + nega
    Double_t q3VecTpcP[2];
    Double_t q3VecTpcN[2];
    Double_t q3Tpc[3]; //psi2 TPC 0:eta negative (C side), 1:eta positive (A side)

    Double_t psi2Tpc[3]; //psi2 TPC 0:eta negative (C side), 1:eta positive (A side)
    Double_t psi3Tpc[3]; //psi3 TPC 0:eta negative (C side), 1:eta positive (A side)

    Double_t fV2ResoV0;
    Double_t fV3ResoV0;

    TH2F *fHCorrV0ChWeghts;     //!
    TH1D *fHCorrQ2xV0C; //!
    TH1D *fHCorrQ2yV0C; //!
    TH1D *fHCorrQ2xV0A; //!
    TH1D *fHCorrQ2yV0A; //!
    TH1D *fHCorrQ3xV0C; //!
    TH1D *fHCorrQ3yV0C; //!
    TH1D *fHCorrQ3xV0A; //!
    TH1D *fHCorrQ3yV0A; //! 


    // TFile*  fOADBFile;

    fitModulationType  fFitModulationType;     // fit modulation type
    
    TF1*               fFitModulation;         //-> modulation fit for rho
    TH1F*              hBkgTracks;         //-> modulation fit for rho
    

    AliAnalysisTaskRawJetWithEP(const AliAnalysisTaskRawJetWithEP&)           ; // not implemented
    AliAnalysisTaskRawJetWithEP &operator=(const AliAnalysisTaskRawJetWithEP&); // not implemented


    /// \cond CLASSIMP
    ClassDef(AliAnalysisTaskRawJetWithEP, 7);
    /// \endcond
  // ### private ###############################################################
};
#endif
