
#ifndef AliAnalysisTaskEPCalibForJet_H
#define AliAnalysisTaskEPCalibForJet_H
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
#include "AliYAMLConfiguration.h"
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
class AliAnalysisTaskEPCalibForJet : public AliAnalysisTaskEmcalJet {
  public:

    AliAnalysisTaskEPCalibForJet();
    AliAnalysisTaskEPCalibForJet(const char *name);
    virtual ~AliAnalysisTaskEPCalibForJet();

    void UserCreateOutputObjects();
    void Terminate(Option_t *option);

    void SetRunList(bool removeDummyTask = false);

    // AliJEQnVectorHandler* GetQnVectorHandler() const {return fQ2VecHandler;}
    TList* GetSplineForqnPercentileList(int det=kFullTPC) const;

    static AliAnalysisTaskEPCalibForJet* AddTaskEPCalibForJet(
        const char *ntracks            = "usedefault",
        const char *nclusters          = "usedefault",
        const char* ncells             = "usedefault",
        const char *suffix             = "");

    enum runModeType{kLocal = 0,                //!< local run mode
                    kGrid = 1                   //!< grid run mode
    };
    enum fitModulationType  { kNoFit = 0,           //!< plato fit
                              kV2 = 1,              //!< only v2 function fit
                              kV3 = 2,              //!< only v3 function fit
                              kCombined = 3,        //!< v2 + v3 function fit
                              kFourierSeries = 4,   //!< Fourier series fit
                              kIntegratedFlow = 5,  //!< Integrated flow fit
                              kQC2 = 6,             //!< qc2 fit
                              kQC4 = 7              //!< qc4 fit
    }; // fit type
    
    enum DetFlev{ kFullTPC = 0,                     //!< use both sides TPC
                  kPosTPC = 1,                      //!< use only an eta positive side TPC
                  kNegTPC = 2,                      //!< use only an eta negative side TPC
                  kFullV0 = 3,                      //!< use both sides V0
                  kV0A = 4,                         //!< use only V0A
                  kV0C = 5                          //!< use only V0A
    };
    enum detectorType{kTPC = 0,                     //!< detector type use TPC
                      kVZEROA = 1,                  //!< detector type use only V0A
                      kVZEROC = 2,                  //!< detector type use only V0C
                      kVZEROComb = 3,               //!< detector type use V0A and V0C
                      kFixedEP = 4                  //!< fixed EP
    };
    enum qnVCalibType{kOrig = 0,                    //!<  way baesed on Flavor group
                      kJeHand = 1                   //!<  way based on JeHandlar
    };
    
    PWG::Tools::AliYAMLConfiguration & GetYAMLConfiguration(){return fYAMLConfig;}

    //== s == Getter Prepare  ################################################
    TString GetLocalRhoName() const {return fLocalRhoName;}
    AliLocalRhoParameter* GetLocalRhoParameter() const {return fLocalRho;}
    //== e == Getter Prepare  ################################################
    
    //== s == Setter Prepare  ################################################
    void SetRunListFileName(std::string fileName) {fRunListFileName = fileName;}

    void SetCalibRefFile(TString calibRefFileName){fCalibRefFileName = calibRefFileName;}
    
    void SetModulationFitType(fitModulationType type) {fFitModulationType = type; }

    void SetEventQA(Int_t qaEventNum){fQaEventNum = qaEventNum;}
    void SetV0Combine(Bool_t bV0Combin){fV0Combin = bV0Combin;}
    void SetQnCalibType(qnVCalibType iQnVCalibType){fQnVCalibType = iQnVCalibType;}
    
    void SetPileupCut(Bool_t bPileupCut){fPileupCut =  bPileupCut;}
    void SetTPCQnMeasure(Bool_t bTPCQnMeasure){fTPCQnMeasure =  bTPCQnMeasure;}
    
    void SetPileupCutQA(Bool_t bPileupCutQA){fPileupCutQA =  bPileupCutQA;}
    void SetGainCalibQA(Bool_t bGainCalibQA){fGainCalibQA =  bGainCalibQA;}
    void SetReCentCalibQA(Bool_t bReCentCalibQA){fReCentCalibQA =  bReCentCalibQA;}
    void SetEPQA(Bool_t bEPQA){fEPQA = bEPQA;}
    void SetTrackQA(Bool_t bTrackQA){fTrackQA = bTrackQA;}
    void SetBkgQA(Bool_t bBkgQA){fBkgQA = bBkgQA;}
    //== e == Setter Prepare  ################################################

  protected:
    void       ExecOnce();
    Bool_t     Run();
    
    // ### protected #############################################################
    THistManager fHistManager;                     ///< Histogram manager
    AliEventCuts fEventCuts;                       ///< event selection

    PWG::Tools::AliYAMLConfiguration fYAMLConfig;  //!<! run List the YAML config
    std::vector<std::string> fUseRunList;          //!<! run list vector

    AliAODEvent*  fAOD;               //!<! AOD event
    TList         *fOutputList;       //!<! output list for histograms
    TString       fCalibRefFileName;  ///<
    TFile*        fCalibRefFile;      //!<!
    TList*        fCalibRefObjList;   //!<!
    TList*        fCalibV0Ref;        //!<!

    AliJEQnVectorHandler* fQ2VecHandler;/// Qn-vector handler
    AliJEQnVectorHandler* fQ3VecHandler;/// Qn-vector handler

    Bool_t  fPileupCut = kFALSE;      ///<
    Bool_t  fTPCQnMeasure = kFALSE;   ///<
    
    Bool_t  fPileupCutQA = kFALSE;    ///<
    Bool_t  fCalibQA = kFALSE;        ///<
    Bool_t  fGainCalibQA = kFALSE;    ///<
    Bool_t  fReCentCalibQA = kFALSE;  ///<
    Bool_t  fEPQA = kFALSE;           ///<
    Bool_t  fTrackQA = kFALSE;        ///<
    Bool_t  fBkgQA = kFALSE;          ///<
    Bool_t  fV0Combin = kFALSE;       ///<
    
    
    // fRunEventList
    
    int           fCalibType;   /// type of calibrations used by handler
    qnVCalibType  fQnVCalibType = kOrig; ///< fCalibration Type
    
    int           fNormMethod;  /// normalisation of Q vector
    
    std::string   fRunListFileName;   ///<
    TString       fSplinesFileName;   ///< Sprine input file name
    

    void       VzeroGainCalibQA();

    // AliMultSelection      *fMultSelection;      //! For Centrality 
    void       SetupPileUpRemovalFunctions();
    Bool_t     CheckEventIsPileUp2018();

    void       AllocatePileupCutHistograms();
    void       AllocateGainCalibHistograms();
    void       AllocateReCentCalibHistograms();
    void       AllocateEventPlaneHistograms();
    void       AllocateTrackHistograms();
    void       AllocateBkgHistograms();
    
    void       MeasureTpcEPQA();
    void       MeasureQnTPC();
    
    void       SetModulationRhoFit();
    void       MeasureBkg();
    void       BkgFitEvaluation();
    
    void       DoMeasureChGainDiff();
    void       DoEventPlane();
    void       DoTrackLoop();

    Bool_t     QnV0GainCalibration();
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

    fitModulationType  fFitModulationType;     ///< fit modulation type
    
    TF1*               fFitModulation;          //-> modulation fit for rho
    TH1F*              hBkgTracks;              //-> modulation fit for rho
    
    Int_t CheckRunNum = -1;         //!<! check run number for local 
    Int_t fQaEventNum = -1;         //!<! qa check run number 

  private:
    AliAnalysisTaskEPCalibForJet(const AliAnalysisTaskEPCalibForJet&)           ; // not implemented
    AliAnalysisTaskEPCalibForJet &operator=(const AliAnalysisTaskEPCalibForJet&); // not implemented


    /// \cond CLASSIMP
    ClassDef(AliAnalysisTaskEPCalibForJet, 22);
    /// \endcond
  // ### private ###############################################################
};
#endif
