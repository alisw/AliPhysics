
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

class TClonesArray;
class TF1;

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
        TString EPCailbType = "JeHand",
        TString EPCalibJEHandRefFileName = "alien:///alice/cern.ch/user/t/tkumaoka/calibV0TPCRun2Vtx10P118qPass3.root",
        TString EPCalibOrigRefFileName = "alien:///alice/cern.ch/user/t/tkumaoka/CalibV0GainCorrectionLHC18q_Oct2021.root",
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
    enum CalibType {
        kQnCalib,
        kQnFrameworkCalib
    };

    PWG::Tools::AliYAMLConfiguration & GetYAMLConfiguration(){return fYAMLConfig;}

    //== s == Getter Prepare  ################################################
    TString GetLocalRhoName() const {return fLocalRhoName;}
    AliLocalRhoParameter* GetLocalRhoParameter() const {return fLocalRho;}
    //== e == Getter Prepare  ################################################
    
    //== s == Setter Prepare  ################################################
    void SetRunListFileName(std::string fileName) {fRunListFileName = fileName;}
    
    void SetModulationFitType(fitModulationType type) {fFitModulationType = type; }

    void SetEventQA(Int_t qaEventNum){fQaEventNum = qaEventNum;}
    void SetV0Combine(Bool_t bV0Combin){fV0Combin = bV0Combin;}
    void SetQnCalibType(TString iQnVCalibType){fQnVCalibType = iQnVCalibType;}
    
    void SetPileupCut(Bool_t bPileupCut){fPileupCut =  bPileupCut;}
    void SetRejectTPCPileup(Bool_t bPileupCut){fRejectTPCPileup = bPileupCut;}
    void SetTPCQnMeasure(Bool_t bTPCQnMeasure){fTPCQnMeasure =  bTPCQnMeasure;}
    
    void SetPileupCutQA(Bool_t bPileupCutQA){fPileupCutQA =  bPileupCutQA;}
    void SetGainCalibQA(Bool_t bGainCalibQA){fGainCalibQA =  bGainCalibQA;}
    void SetReCentCalibQA(Bool_t bReCentCalibQA){fReCentCalibQA =  bReCentCalibQA;}
    void SetDoReCentCalib(Bool_t bDoReCentCalib){fDoReCentCalib =  bDoReCentCalib;}
    void SetEPQA(Bool_t bEPQA){fEPQA = bEPQA;}
    void SetEPResoQA(Bool_t bEPResoQA){fEPResoQA = bEPResoQA;}
    void SetTrackQA(Bool_t bTrackQA){fTrackQA = bTrackQA;}
    void SetBkgQA(Bool_t bBkgQA){fBkgQA = bBkgQA;}



    void SetCalibOrigRefObjList(TList *hWgtsV0ZDC){this->fCalibRefObjList = (TList *) hWgtsV0ZDC->Clone();}

    void SetLRefMultV0BefCorPfpx(AliOADBContainer *hList){this->fMultV0BefCorPfpx =   (AliOADBContainer *) hList->Clone();}
    void SetLRefQx2am(TObjArray *hList){this->fOADBzArray_contQx2am = (TObjArray *)hList->Clone();}
    void SetLRefQy2am(TObjArray *hList){this->fOADBzArray_contQy2am = (TObjArray *)hList->Clone();}
    void SetLRefQx2as(TObjArray *hList){this->fOADBzArray_contQx2as = (TObjArray *)hList->Clone();}
    void SetLRefQy2as(TObjArray *hList){this->fOADBzArray_contQy2as = (TObjArray *)hList->Clone();}
    void SetLRefQx3am(TObjArray *hList){this->fOADBzArray_contQx3am = (TObjArray *)hList->Clone();}
    void SetLRefQy3am(TObjArray *hList){this->fOADBzArray_contQy3am = (TObjArray *)hList->Clone();}
    void SetLRefQx3as(TObjArray *hList){this->fOADBzArray_contQx3as = (TObjArray *)hList->Clone();}
    void SetLRefQy3as(TObjArray *hList){this->fOADBzArray_contQy3as = (TObjArray *)hList->Clone();}
    void SetLRefQx2cm(TObjArray *hList){this->fOADBzArray_contQx2cm = (TObjArray *)hList->Clone();}
    void SetLRefQy2cm(TObjArray *hList){this->fOADBzArray_contQy2cm = (TObjArray *)hList->Clone();}
    void SetLRefQx2cs(TObjArray *hList){this->fOADBzArray_contQx2cs = (TObjArray *)hList->Clone();}
    void SetLRefQy2cs(TObjArray *hList){this->fOADBzArray_contQy2cs = (TObjArray *)hList->Clone();}
    void SetLRefQx3cm(TObjArray *hList){this->fOADBzArray_contQx3cm = (TObjArray *)hList->Clone();}
    void SetLRefQy3cm(TObjArray *hList){this->fOADBzArray_contQy3cm = (TObjArray *)hList->Clone();} 
    void SetLRefQx3cs(TObjArray *hList){this->fOADBzArray_contQx3cs = (TObjArray *)hList->Clone();} 
    void SetLRefQy3cs(TObjArray *hList){this->fOADBzArray_contQy3cs = (TObjArray *)hList->Clone();}
    void SetLRefTPCposEta(TObjArray *hList){this->fOADBcentArray_contTPCposEta = (TObjArray *) hList->Clone();}
    void SetLRefTPCnegEta(TObjArray *hList){this->fOADBcentArray_contTPCnegEta = (TObjArray *) hList->Clone();}

    bool ExtractRecentPara(TFile *RefFile, TObjArray *lRefQx2am, TObjArray *lRefQy2am, TObjArray *lRefQx2as, TObjArray *lRefQy2as, TObjArray *lRefQx3am, TObjArray *lRefQy3am, TObjArray *lRefQx3as, TObjArray *lRefQy3as, TObjArray *lRefQx2cm, TObjArray *lRefQy2cm, TObjArray *lRefQx2cs, TObjArray *lRefQy2cs, TObjArray *lRefQx3cm,TObjArray *lRefQy3cm, TObjArray *lRefQx3cs, TObjArray *lRefQy3cs, TObjArray *lRefTPCposEta, TObjArray *lRefTPCnegEta);

    bool SetAODEvent(AliAODEvent* event); 
    void ResetAODEvent(); 
    
    void SetCalibrationType(int calibType) {fCalibType = calibType;}
    void SetTPCEtaLimits(double absetamin=0., double absetamx=0.8) {fEtaMinTPC=absetamin; fEtaMaxTPC=absetamx;}
    void SetTPCPtLimits(double ptmin=0.2, double ptmax=5) {fPtMinTPC=ptmin; fPtMaxTPC=ptmax;}
    void SetFractionOfTPCtracksToUse(double fracToKeep) {fFractionOfTracksForQnTPC = fracToKeep;}
    // void SetCalibrationsOADBFileName(TString OADBfileName) {fOADBFileName = OADBfileName;}


    int GetCalibrationType() const {return fCalibType;}
    int GetNormalisationMethod() const {return fNormMethod;}
    // TString GetCalibrationsOADBFileName() const {return fOADBFileName;}
    
    void GetQnVec(double QnVecFull[2], double QnVecA[2], double QnVecC[2], Double_t QnNorm[3], Double_t Multi[3]);
    void Getqn(Double_t Qn[3], Double_t QnNorm[3], Double_t Multi[3]);


    void EnablePhiDistrHistos();
    TH2F* GetPhiDistrHistosTPCPosEta() const {return fPhiVsCentrTPC[0];}
    TH2F* GetPhiDistrHistosTPCNegEta() const {return fPhiVsCentrTPC[1];}
    //== e == Setter Prepare  ################################################




  protected:
    void       ExecOnce();
    Bool_t     Run();
    
    // ### protected #############################################################
    THistManager fHistManager;                     ///< Histogram manager
    AliEventCuts fEventCuts;                       ///< event selection

    PWG::Tools::AliYAMLConfiguration fYAMLConfig;  //!<! run List the YAML config
    std::vector<std::string> fUseRunList;          //!<! run list vector

    AliAODEvent*  fAOD;                 //< AOD event
    TList         *fOutputList;         //!<! output list for histograms
    TList*        fCalibV0Ref;          ///<

    Bool_t  fPileupCut       = kFALSE;  ///<
    Bool_t  fRejectTPCPileup = kFALSE;  ///<
    Bool_t  fTPCQnMeasure    = kFALSE;  ///<
    
    Bool_t  fPileupCutQA     = kFALSE;  ///<
    Bool_t  fCalibQA         = kFALSE;  ///<
    Bool_t  fGainCalibQA     = kFALSE;  ///<
    Bool_t  fDoReCentCalib   = kFALSE;  ///<
    Bool_t  fReCentCalibQA   = kFALSE;  ///<
    Bool_t  fEPQA            = kFALSE;  ///<
    Bool_t  fEPResoQA        = kFALSE;  ///<
    Bool_t  fTrackQA         = kFALSE;  ///<
    Bool_t  fBkgQA           = kFALSE;  ///<
    Bool_t  fV0Combin        = kFALSE;  ///<
    
    
    // fRunEventList
    
    int           fCalibType;   /// type of calibrations used by handler
    TString       fQnVCalibType = "kOrig"; ///< fCalibration Type
    
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
    void       AllocateEPResoHistograms();
    void       AllocateTrackHistograms();
    void       AllocateBkgHistograms();
    
    void       ReSetValuable();

    void       MeasureTpcEPQA();
    void       MeasureQnTPC();
    
    void       SetModulationRhoFit();
    void       MeasureBkg();
    void       BkgFitEvaluation();
    
    void       DoMeasureChGainDiff();
    Bool_t     DoEventPlane();
    void       DoTrackLoop();

    Bool_t     QnJEHandlarEPGet();
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
    
    
    //qnVector 0:x, 1:y
    Double_t befGainCalibQ2VecV0M[2];   //!<! Q2 V0 C+A vector(x,y)
    Double_t befGainCalibQ2VecV0C[2];   //!<! Q2 V0 C+A vector(x,y)
    Double_t befGainCalibQ2VecV0A[2];   //!<! Q2 V0 C+A vector(x,y)
    Double_t befGainCalibQ2V0[3];       //!<! psi2 V0 0:combin, 1:eta negative (C side), 2:eta positive (A side)
    Double_t befGainCalibQ3VecV0M[2];   //!<! Q3 V0 C+A vector(x,y)
    Double_t befGainCalibQ3VecV0C[2];   //!<! Q3 V0 C vector(x,y)
    Double_t befGainCalibQ3VecV0A[2];   //!<! Q3 V0 A vector(x,y)
    Double_t befGainCalibQ3V0[3];       //!<! Q3 V0 0:combin, 1:eta negative (C side), 2:eta positive (A side)

    Double_t V0Mult2[3];    /// For q2 V0 0:combin, 1:eta negative (C side), 2:eta positive (A side)
    Double_t V0Mult3[3];    /// For q3 V0 0:combin, 1:eta negative (C side), 2:eta positive (A side)
    
    Double_t q2VecV0M[2];   //!<! Q2 V0 C+A vector(x,y)
    Double_t q2VecV0C[2];   //!<! Q2 V0 C+A vector(x,y)
    Double_t q2VecV0A[2];   //!<! Q2 V0 C+A vector(x,y)
    Double_t q2V0[3];       //!<! psi2 V0 0:combin, 1:eta negative (C side), 2:eta positive (A side)
    Double_t q2NormV0[3];   ///< Q2 V0 0:combin, 1:eta negative (C side), 2:eta positive (A side)
    Double_t q3VecV0M[2];   //!<! Q3 V0 C+A vector(x,y)
    Double_t q3VecV0C[2];   //!<! Q3 V0 C vector(x,y)
    Double_t q3VecV0A[2];   //!<! Q3 V0 A vector(x,y)
    Double_t q3V0[3];       //!<! Q3 V0 0:combin, 1:eta negative (C side), 2:eta positive (A side)
    Double_t q3NormV0[3];   ///< Q3 V0 0:combin, 1:eta negative (C side), 2:eta positive (A side)

    Double_t psi2V0[3];     //!<! psi2 V0 0:combin, 1:eta negative (C side), 2:eta positive (A side)
    Double_t psi3V0[3];     //!<! psi3 V0 0:combin, 1:eta negative (C side), 2:eta positive (A side)


    //qnVector 0:eta posi, 1:eta nega
    Double_t TpcMult2[3];  ///< For q2 V0 0:combin, 1:eta negative (C side), 2:eta positive (A side)
    Double_t TpcMult3[3];  ///< For q3 V0 0:combin, 1:eta negative (C side), 2:eta positive (A side)
    Double_t q2VecTpcM[2];  //!<!  Q2 TPC posi + nega
    Double_t q2VecTpcP[2];  //!<!  Q2 TPC nega
    Double_t q2VecTpcN[2];  //!<!  Q2 TPC nega
    Double_t q2Tpc[3];      //!<!  psi2 TPC 0:eta negative (C side), 1:eta positive (A side)side)
    Double_t q2NormTpc[3];  ///<  psi2 TPC 0:eta negative (C side), 1:eta positive (A side)side)
    Double_t q3VecTpcM[2];  //!<!  Q3 TPC posi + nega
    Double_t q3VecTpcP[2];  //!<!  Q3 TPC posi
    Double_t q3VecTpcN[2];  //!<!  Q3 TPC nega
    Double_t q3Tpc[3];      //!<!  psi2 TPC 0:eta negative (C side), 1:eta positive (A side)
    Double_t q3NormTpc[3];  ///<  psi2 TPC 0:eta negative (C side), 1:eta positive (A side)side)

    Double_t psi2Tpc[3];    //!<!  psi2 TPC 0:eta negative (C side), 1:eta positive (A side)
    Double_t psi3Tpc[3];    //!<!  psi3 TPC 0:eta negative (C side), 1:eta positive (A side)

    Double_t fV2ResoV0;     //!<!  V2 resolution value
    Double_t fV3ResoV0;     //!<!  V3 resolution value

    /// Functions for Pile Up Event Removal:
    TF1   *fV0CutPU;        //!<!
    TF1   *fSPDCutPU;       //!<!
    TF1   *fMultCutPU;      //!<!
    TF1   *fCenCutLowPU;    //!<!
    TF1   *fCenCutHighPU;   //!<!

    TList *fCalibRefObjList;   ///<
    TH2F *fHCorrV0ChWeghts;   //!<!
    TH1D *fHCorrQ2xV0M;       //!<!
    TH1D *fHCorrQ2yV0M;       //!<!
    TH1D *fHCorrQ2xV0C;       //!<!
    TH1D *fHCorrQ2yV0C;       //!<!
    TH1D *fHCorrQ2xV0A;       //!<!
    TH1D *fHCorrQ2yV0A;       //!<!

    TH1D *fHCorrQ3xV0M;       //!<!
    TH1D *fHCorrQ3yV0M;       //!<!
    TH1D *fHCorrQ3xV0C;       //!<!
    TH1D *fHCorrQ3yV0C;       //!<!
    TH1D *fHCorrQ3xV0A;       //!<!
    TH1D *fHCorrQ3yV0A;       //!<! 

    fitModulationType  fFitModulationType;     ///< fit modulation type
    
    TF1*               fFitModulation;          //-> modulation fit for rho
    TH1F*              hBkgTracks;              //-> modulation fit for rho
    

    /// ========================================================================================
    //utils methods
    void ComputeQvecV0(Double_t QnVecV0M[2],Double_t QnVecV0C[2],Double_t QnVecV0A[2], Double_t QnNorm[3], Double_t Multi[3], unsigned int harmonic);
    void ComputeQvecTpc(Double_t QnVecTpcM[2],Double_t QnVecTpcN[2],Double_t QnVecTpcP[2], Double_t QnNorm[3], Double_t Multi[3], unsigned int harmonic);
    Double_t ComputeEventPlaneAngle(Double_t QnVec[2], Double_t harmonic) const {return (TMath::Pi()+TMath::ATan2(-QnVec[1],-QnVec[0]))/harmonic;};
    
    short GetVertexZbin() const;
    Int_t GetCentBin();
    bool OpenInfoCalbration();

    bool IsTrackSelected(AliAODTrack* track);

    //data members
    double fEtaMinTPC;         ///< Absolute minimum value of eta for TPC Q vector
    double fEtaMaxTPC;         ///< Absolute maximum value of eta for TPC Q vector
    double fPtMinTPC;          ///< Minimum value of pt for TPC Q vector
    double fPtMaxTPC;          ///< Maximum value of pt for TPC Q vector

    TBits fUsedTrackPosIDs;    ///< IDs of tracks (with positive ID) used for the computation of the Q vector (TPC)
    TBits fUsedTrackNegIDs;    ///< IDs of tracks (with negative ID) used for the computation of the Q vector (TPC)

    //random rejection 
    double fFractionOfTracksForQnTPC;   ///< Fraction of tracks to keep when rejection enabled

    //data members needed for calibrations
  
    AliAODVZERO* fV0;                        ///< V0 info for the considered event
    int fRun;                                ///< Run number
    double fZvtx;                            ///< Primary vertx Z 
    double fCentrality;                      ///< Event centrality

    bool fIsOADBFileOpen;                    ///< Flag to test whether the OADB file is open

    // OADB Objects for streaming
    AliOADBContainer * fMultV0BefCorPfpx;    ///< OADB object containing the hMultV0BefCorPfpx histograms
    TObjArray *fOADBzArray_contQx2am; ///< Array of OADB contQx2am object, index is z-vertex bin
    TObjArray *fOADBzArray_contQy2am; ///< Array of OADB contQy2am object, index is z-vertex bin
    TObjArray *fOADBzArray_contQx2as; ///< Array of OADB contQx2as object, index is z-vertex bin
    TObjArray *fOADBzArray_contQy2as; ///< Array of OADB contQy2as object, index is z-vertex bin

    TObjArray *fOADBzArray_contQx2cm; ///< Array of OADB contQx2cm object, index is z-vertex bin
    TObjArray *fOADBzArray_contQy2cm; ///< Array of OADB contQy2cm object, index is z-vertex bin
    TObjArray *fOADBzArray_contQx2cs; ///< Array of OADB contQx2cs object, index is z-vertex bin
    TObjArray *fOADBzArray_contQy2cs; ///< Array of OADB contQy2cs object, index is z-vertex bin

    TObjArray *fOADBzArray_contQx3am; ///< Array of OADB contQx3am object, index is z-vertex bin
    TObjArray *fOADBzArray_contQy3am; ///< Array of OADB contQy3am object, index is z-vertex bin
    TObjArray *fOADBzArray_contQx3as; ///< Array of OADB contQx3as object, index is z-vertex bin
    TObjArray *fOADBzArray_contQy3as; ///< Array of OADB contQy3as object, index is z-vertex bin

    TObjArray *fOADBzArray_contQx3cm; ///< Array of OADB contQx3cm object, index is z-vertex bin
    TObjArray *fOADBzArray_contQy3cm; ///< Array of OADB contQy3cm object, index is z-vertex bin
    TObjArray *fOADBzArray_contQx3cs; ///< Array of OADB contQx3cs object, index is z-vertex bin
    TObjArray *fOADBzArray_contQy3cs; ///< Array of OADB contQy3cs object, index is z-vertex bin

    TObjArray *fOADBcentArray_contTPCposEta; ///< Array of OADB contTPCposEta, index is cent bin
    TObjArray *fOADBcentArray_contTPCnegEta; ///< Array of OADB contTPCnegEta, index is cent bin

    int fCalibObjRun;            ///< Run of loaded calibration objects

    TH1D* fHistMultV0;           ///< Profile from V0 multiplicity

    TH1D* fQx2mV0A[14];          ///< <Qxn> V0A
    TH1D* fQy2mV0A[14];          ///< <Qyn> V0A
    TH1D* fQx2sV0A[14];          ///< sigma Qxn V0A
    TH1D* fQy2sV0A[14];          ///< sigma Qyn V0A
    
    TH1D* fQx2mV0C[14];          ///< <Qxn> V0C
    TH1D* fQy2mV0C[14];          ///< <Qyn> V0C
    TH1D* fQx2sV0C[14];          ///< sigma Qxn V0C
    TH1D* fQy2sV0C[14];          ///< sigma Qyn V0C

    TH1D* fQx3mV0A[14];          ///< <Qxn> V0A
    TH1D* fQy3mV0A[14];          ///< <Qyn> V0A
    TH1D* fQx3sV0A[14];          ///< sigma Qxn V0A
    TH1D* fQy3sV0A[14];          ///< sigma Qyn V0A
    
    TH1D* fQx3mV0C[14];          ///< <Qxn> V0C
    TH1D* fQy3mV0C[14];          ///< <Qyn> V0C
    TH1D* fQx3sV0C[14];          ///< sigma Qxn V0C
    TH1D* fQy3sV0C[14];          ///< sigma Qyn V0C

    bool fV0CalibZvtxDiff;       //< flag to properly manage Zvtx differential V0 calibrations

    TH1D* fWeightsTPCPosEta[11];  ///< Weights for TPC tracks with eta > 0
    TH1D* fWeightsTPCNegEta[11];  ///< Weights for TPC tracks with eta < 0
    bool fEnablePhiDistrHistos;  ///< Enable phi distribution histos
    TH2F* fPhiVsCentrTPC[2];     ///< Phi vs. centr TH2 of selected TPC tracks in eta>0 and eta<0
    



    Int_t CheckRunNum = -1;         //!<! check run number for local 
    Int_t fQaEventNum = -1;         //!<! qa check run number 

  private:
    AliAnalysisTaskEPCalibForJet(const AliAnalysisTaskEPCalibForJet&)           ; // not implemented
    AliAnalysisTaskEPCalibForJet &operator=(const AliAnalysisTaskEPCalibForJet&); // not implemented


    /// \cond CLASSIMP
    ClassDef(AliAnalysisTaskEPCalibForJet, 80);
    /// \endcond
  // ### private ###############################################################
};
#endif

