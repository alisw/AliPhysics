
#ifndef AliAnalysisTaskEmbeddingJetWithEP_H
#define AliAnalysisTaskEmbeddingJetWithEP_H
/**
 * \file AliAnalysisTaskEmbeddingJetWithEP.h
 * \brief Declaration of class AliAnalysisTaskEmbeddingJetWithEP
 *
 * In this header file the class AliAnalysisTaskEmbeddingJetWithEP is declared.
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
#include "AliYAMLConfiguration.h"
// #include "AliMultSelection.h"

#include "AliJEQnVectorHandler.h"
#include "AliAnalysisTaskJetQnVectors.h"
// #include "AliAnalysisTaskFlowVectorCorrections.h"
#include "AliQnCorrectionsManager.h"

#include "AliEmcalEmbeddingQA.h"

class TClonesArray;
class TF1;

class AliJEQnVectorHandler;

class AliRhoParameter;
class AliLocalRhoParameter;


/**
 * \class AliAnalysisTaskEmbeddingJetWithEP
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
class AliAnalysisTaskEmbeddingJetWithEP : public AliAnalysisTaskEmcalJet {
  public:

    AliAnalysisTaskEmbeddingJetWithEP();
    AliAnalysisTaskEmbeddingJetWithEP(const char *name);
    virtual ~AliAnalysisTaskEmbeddingJetWithEP();

    void UserCreateOutputObjects();
    void Terminate(Option_t *option);

    void SetRunList(bool removeDummyTask = false);

    // AliJEQnVectorHandler* GetQnVectorHandler() const {return fQ2VecHandler;}
    TList* GetSplineForqnPercentileList(int det=kFullTPC) const;

    static AliAnalysisTaskEmbeddingJetWithEP* AddTaskEmbeddingJetWithEP(
        const char *ntracks            = "usedefault",
        const char *nclusters          = "usedefault",
        const char* ncells             = "usedefault",
        const char *suffix             = "");

    enum runModeType{
      kLocal,           // run location (local)
      kGrid             // run location (grid)
    };

    enum fitModulationType  { 
      kNoFit,           // background fit type plane
      kV2,              // background fit type v2
      kV3,              // background fit type v3
      kCombined,        // background fit type v2 and v3
      kFourierSeries,   // background fit type FourierSeries
      kIntegratedFlow,  // background fit type integrated flow
      kQC2,             // background fit type qumrant2
      kQC4              // background fit type qumrant4
    };
    
    enum DetFlev{
      kFullTPC,         // use all TPC Qn Vector
      kPosTPC,          // use positive eta TPC Qn Vector
      kNegTPC,          // use negative eta TPC Qn Vector
      kFullV0,          // use A and C V0 Qn Vector
      kV0A,             // use C V0 Qn Vector
      kV0C              // use A V0 Qn Vector
    };

    enum detectorType{
      kTPC,             // use all TPC for Event Plane
      kVZEROA,          // use V0 A for Event Plane
      kVZEROC,          // use V0 C for Event Plane
      kVZEROComb,       // use V0 C and A for Event Plane
      kFixedEP          // use Event Plane
    };  // detector that was used for event plane
    
    
    
    PWG::Tools::AliYAMLConfiguration & GetYAMLConfiguration(){return fYAMLConfig;}

    //== s == Getter Prepare  ################################################
    TString GetLocalRhoName() const {return fLocalRhoName;}
    AliLocalRhoParameter* GetLocalRhoParameter() const {return fLocalRho;}
    //== e == Getter Prepare  ################################################
    
    //== s == Setter Prepare  ################################################
    void SetRunListFileName(std::string fileName) {fRunListFileName = fileName;}

    void SetCalibRefFile(TString calibRefFileName){fCalibRefFileName = calibRefFileName;}
    
    void SetModulationFitType(fitModulationType type) {fFitModulationType = type; }

    void SetPileupCut(Bool_t bPileupCut){fPileupCut =  bPileupCut;} 

    void SetPileupCutQA(Bool_t bPileupCutQA){fPileupCutQA =  bPileupCutQA;}
    void SetEventQA(Int_t qaEventNum){fQaEventNum = qaEventNum;}
    void SetV0Combine(Bool_t bV0Combin){fV0Combin = bV0Combin;}

    void SetBkgQA(Bool_t bBkgQA){fBkgQA = bBkgQA;}

    void SetMinPt(Double_t d){fMinPt = d;}
    void SetMaxPt(Double_t d){fMaxPt = d;}


    void SetDoJetMatchingGeom(Bool_t b)             { fDoJetMatchingGeom = b; }
    void SetDoDifferentialRM(Bool_t b)              { fDoDifferentialRM = b; }
    void SetDoJetMatchingMCFraction(Bool_t b)       { fDoJetMatchingMCFraction = b; }
    void SetRequireMatchedJetAcc(Bool_t b)          { fRequireMatchedJetAcc = b; }
    void SetJetMatchingR(Double_t r)                { fJetMatchingR = r; }
    void SetMinSharedPFraction(double d)            { fMinSharedPFraction = d; }
    void SetMCJetMinMatchingPt(Double_t min)        { fMCJetMinMatchingPt = min; }
    void SetDetJetMinMatchingPt(Double_t min)       { fDetJetMinMatchingPt = min; }
    //== e == Setter Prepare  ################################################

  protected:
    void       ExecOnce();
    Bool_t     Run();

    AliEventCuts fEventCuts;          ///< event selection utility
    PWG::Tools::AliYAMLConfiguration fYAMLConfig;  //!<! run List the YAML config
    
    THistManager fHistManager;        ///< Histogram manager
    std::vector<std::string> fUseRunList; //!<! run list vector
    
    // Embedding parameters
    AliEmcalEmbeddingQA fEmbeddingQA;  //!<! QA hists for embedding (will only be added if embedding)
    Bool_t  fDoJetMatchingGeom; ///< Do geometrical matching between det-level and truth-level jet container
    Bool_t  fDoJetMatchingMCFraction;  ///< Do MC-fraction based matching using PbPb det-level, pp det-level, and pp truth-level jet containers
    Bool_t  fDoDifferentialRM;         ///< This allows a differential RM depenent on the angularity of the particle lvl jet
    
    AliJetContainer* fMCJetContainer;    //!<!Pointer to jet container of truth-level jets
    AliJetContainer* fDetJetContainer;   //!<!Pointer to jet container of det-level jets
    AliJetContainer* fDetJetContainerPPIntermediate; //!<!Pointer to jet container of intermediate pp det-level jets, if MC-fraction matching
    Bool_t  fRequireMatchedJetAcc;  ///< Flag to require matched truth jet to be accepted (aside from geometrical acceptance)
    Double_t fJetMatchingR;         ///< Jet matching R threshold
    Double_t fMinSharedPFraction;    ///< Minimum shared momentum (pp det-level track pT in combined jet) / (pp det-level track pT)
    Double_t fMCJetMinMatchingPt;         ///< Min jet pT for MC jets being matched, for when container criteria is not applied
    Double_t  fDetJetMinMatchingPt;       ///< Min jet pT for Det jets being matched, for when container criteria is not applied
    Double_t  fPlotJetMatchCandThresh;    ///< Threshold for jet R to count candidates, affects plotting only

    // Plotting parameters
    Double_t  fMinPt;  ///< Histogram min pT limit
    Double_t  fMaxPt;  ///< Histogram max pT limit

    Bool_t  bUseJetCont2Acc = kFALSE; ///<
    // ### protected #############################################################
  
  // ### private ###############################################################
  private:
    
    AliAODEvent* fAOD;                //!<! AOD event
    TList*  fCalibV0Ref;              //!<!
    TFile*  fCalibRefFile;            //!<!
    TList*  fCalibRefObjList;         //!<!

    Bool_t  fPileupCut = kFALSE;      ///<
    Bool_t  fTPCQnMeasure = kFALSE;   ///<
    
    Bool_t  fPileupCutQA = kFALSE;    ///<
    Bool_t  fCalibQA = kFALSE;        ///<
    Bool_t  fGainCalibQA = kFALSE;    ///<
    Bool_t  fReCentCalibQA = kFALSE;  ///<
    Bool_t  fEPQA = kFALSE;           ///<
    Bool_t  fTrackQA = kFALSE;        ///<
    Bool_t  fBkgQA = kFALSE;          ///<
    
    
    
    AliJEQnVectorHandler* fQ2VecHandler;    ///< Qn-vector handler
    AliJEQnVectorHandler* fQ3VecHandler;    ///< Qn-vector handler

    TList         *fOutputList;             //!<! output list for histograms
    int           fCalibType;               ///< type of calibrations used by handler
    int           fNormMethod;              ///< normalisation of Q vector
    Bool_t        fV0Combin = kFALSE;       ///<
    int           fQaEventNum = -1;          ///<
    
    std::string   fRunListFileName;         ///< Run list file Name
    TString       fCalibRefFileName;        ///< Calibration input file name

    void          LoadSpliForqnPerce();
    void          VzeroGainCalibQA();

    // AliMultSelection      *fMultSelection;      //! For Centrality 
    void       SetupPileUpRemovalFunctions();
    Bool_t     CheckEventIsPileUp2018();

    void       AllocateJetHistograms();
    void       AllocateMatchedJetHistograms();

    void       MeasureTpcEPQA();
    
    void       SetModulationRhoFit();
    void       MeasureBkg();
    
    void       DoEventPlane();
    void       DoJetLoop();
    void       FillMatchedJetHistograms();

    void       ComputeJetMatches(AliJetContainer* jetCont1, AliJetContainer* jetCont2, \
      Bool_t bUseJetCont2Acc);
    void       SetJetClosestCandidate(AliEmcalJet* jet1, AliEmcalJet* jet2);
    const AliEmcalJet* GetMatchedPartLevelJet(const AliEmcalJet* detJet, Double_t detJetPt);
    Double_t   GetAngularity(const AliEmcalJet* jet);
    Double_t   GetRelativePhi(Double_t mphi, Double_t vphi);
    
    Bool_t     QnGainCalibration();
    Bool_t     QnRecenteringCalibration();

    Double_t CalcEPAngle(double Qx,double Qy) const {return (TMath::Pi()+TMath::ATan2(-Qy,-Qx))/2;}
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
    
    
    
    /// Functions for Pile Up Event Removal:
    TF1   *fV0CutPU;        //!<!
    TF1   *fSPDCutPU;       //!<!
    TF1   *fMultCutPU;      //!<!
    TF1   *fCenCutLowPU;    //!<!
    TF1   *fCenCutHighPU;   //!<!

    //qnVector 0:x, 1:y
    Double_t q2VecV0M[2];   //!<! Q2 V0 C+A vector(x,y)
    Double_t q2VecV0C[2];   //!<! Q2 V0 C+A vector(x,y)
    Double_t q2VecV0A[2];   //!<! Q2 V0 C+A vector(x,y)
    Double_t q2V0[3];       //!<! psi2 V0 0:combin, 1:eta negative (C side), 2:eta positive (A side)
    Double_t q3VecV0M[2];   //!<! Q3 V0 C+A vector(x,y)
    Double_t q3VecV0C[2];   //!<! Q3 V0 C vector(x,y)
    Double_t q3VecV0A[2];   //!<! Q3 V0 A vector(x,y)
    Double_t q3V0[3];       //!<! Q3 V0 0:combin, 1:eta negative (C side), 2:eta positive (A side)
    
    Double_t psi2V0[3];     //!<! psi2 V0 0:combin, 1:eta negative (C side), 2:eta positive (A side)
    Double_t psi3V0[3];     //!<! psi3 V0 0:combin, 1:eta negative (C side), 2:eta positive (A side)


    //qnVector 0:eta posi, 1:eta nega
    Double_t q2VecTpcM[2];  //!<!  Q2 TPC posi + nega
    Double_t q2VecTpcP[2];  //!<!  Q2 TPC nega
    Double_t q2VecTpcN[2];  //!<!  Q2 TPC nega
    Double_t q2Tpc[3];      //!<!  psi2 TPC 0:eta negative (C side), 1:eta positive (A side)side)
    Double_t q3VecTpcM[2];  //!<!  Q3 TPC posi + nega
    Double_t q3VecTpcP[2];  //!<!  Q3 TPC posi
    Double_t q3VecTpcN[2];  //!<!  Q3 TPC nega
    Double_t q3Tpc[3];      //!<!  psi2 TPC 0:eta negative (C side), 1:eta positive (A side)

    Double_t psi2Tpc[3];    //!<!  psi2 TPC 0:eta negative (C side), 1:eta positive (A side)
    Double_t psi3Tpc[3];    //!<!  psi3 TPC 0:eta negative (C side), 1:eta positive (A side)

    Double_t fV2ResoV0;     //!<!  V2 resolution value
    Double_t fV3ResoV0;     //!<!  V3 resolution value

    TH2F *fHCorrV0ChWeghts;   //!<!
    TH1D *fHCorrQ2xV0C;       //!<!
    TH1D *fHCorrQ2yV0C;       //!<!
    TH1D *fHCorrQ2xV0A;       //!<!
    TH1D *fHCorrQ2yV0A;       //!<!
    TH1D *fHCorrQ3xV0C;       //!<!
    TH1D *fHCorrQ3yV0C;       //!<!
    TH1D *fHCorrQ3xV0A;       //!<!
    TH1D *fHCorrQ3yV0A;       //!<! 

    fitModulationType  fFitModulationType;      ///< fit modulation type
    
    TF1*               fFitModulation;          //-> modulation fit for rho
    TH1F*              hBkgTracks;              //-> modulation fit for rho
    

    AliAnalysisTaskEmbeddingJetWithEP(const AliAnalysisTaskEmbeddingJetWithEP&)           ; // not implemented
    AliAnalysisTaskEmbeddingJetWithEP &operator=(const AliAnalysisTaskEmbeddingJetWithEP&); // not implemented


    /// \cond CLASSIMP
    ClassDef(AliAnalysisTaskEmbeddingJetWithEP, 32);
    /// \endcond
  // ### private ###############################################################
};
#endif
