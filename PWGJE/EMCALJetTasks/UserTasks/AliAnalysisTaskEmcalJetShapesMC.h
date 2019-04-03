#ifndef ALIANALYSISTASKEMCALJETSHAPESMC_H
#define ALIANALYSISTASKEMCALJETSHAPESMC_H


/// \class AliAnalysisTaskEmcalJetShapesMC
/// \brief Task to store and correlate the MC shapes
///
/// Task to store in a tree a certain number of jet shapes at generated level.
/// This can be used both for MC studies and as input for TMVA analysis.
/// No subtraction is implemented for the moment even the fglag are in place.
///
/// \author Alice Davide Caffarri <davide.caffarri@cern.ch>, CERN
/// \author Alice Leticia Cunqueiro Mendez <leticia.cunqueiro.mendez@cern.ch>, Westfaelische Wilhelms-Universitaet Muenster (DE)
/// \date Apr 21, 2016

class TH1;
class TH2;
class TH3;
class TH3F;
class TTree;
class THnSparse;
class TClonesArray;
class TArrayI;
class TRandom3;
class TF1;
class AliAnalysisManager;
class AliJetContainer;
class AliEmcalJetFinder;
class AliFJWrapper;



#include "AliAnalysisTaskEmcalJet.h"
#include "AliFJWrapper.h"
#include "AliClusterContainer.h"


class AliAnalysisTaskEmcalJetShapesMC : public AliAnalysisTaskEmcalJet {
 public:

  enum JetShapeType {
      kGenShapes = 0
  };
  enum JetShapeSub {
    kNoSub = 0,
    kConstSub = 1,
    kDerivSub = 2
  };
  enum JetSelectionType {
    kInclusive = 0,
    kRecoil = 1
  };

  enum DerivSubtrOrder {
    kSecondOrder = 0,
    kFirstOrder = 1
  };

  AliAnalysisTaskEmcalJetShapesMC();
  AliAnalysisTaskEmcalJetShapesMC(const char *name);
  virtual ~AliAnalysisTaskEmcalJetShapesMC();

  void                                UserCreateOutputObjects();
  void                                Terminate(Option_t *option);

  //Setters
  void SetJetContainer(Int_t c)                             { fContainer     = c   ; }
  void SetMinFractionShared(Double_t f)                     { fMinFractionShared = f   ; }
  void SetJetShapeType(JetShapeType t)                      { fJetShapeType       = t   ; }
  void SetJetShapeSub(JetShapeSub t)                        { fJetShapeSub     = t   ; }
  void SetJetSelection(JetSelectionType t)                  { fJetSelection    = t   ; }
  void SetSwitchSDkT(Int_t c)                               { fSwitchSDKtNSub   = c   ; }
  void SetSwitchSDMin(Int_t c)                              { fSwitchSDMinNSub   = c   ; }
  void SetSwitchkT(Int_t c)                                 { fSwitchKtNSub   = c   ; }
  void SetSwitchMin(Int_t c)                                { fSwitchMinNSub   = c   ; }
  void SetSwitchAkT(Int_t c)                                { fSwitchAktNSub   = c   ; }
  void SetAdditionalTracks(Int_t c)                         { fAdditionalTracks =c   ;}
  void SetHardCutoff(Float_t c)                             { fHardCutoff =c   ;}
  void SetOptionalPartonInfo(Int_t c)                       { fOptionalPartonInfo = c; }
  void SetJetPtThreshold(Float_t f)                         { fPtThreshold     = f   ; }
  void SetRMatching(Float_t f)                              { fRMatching = f ;}
  void SetJetRadius(Float_t f)                              { fJetRadius = f ;}
  void SetSubjetRadius(Float_t f)                           { fSubjetRadius = f ;}
  void SetSelectShapes(Int_t c)                             { fSelectedShapes = c;}
  void SetPtTriggerSelections(Float_t minpT, Float_t maxpT) { fminpTTrig = minpT; fmaxpTTrig = maxpT; }
  void SetAngularWindowRecoilJet (Float_t t)                {fangWindowRecoil = t; }
  Float_t GetMinPtTriggerSelection()                        {return fminpTTrig;}
  Float_t GetMaxPtTriggerSelection()                        {return fmaxpTTrig;}
  void SetMediumParameters(Float_t t, Float_t c)            {fqhat=t; fxlength=c;}
  void SetCentralitySelectionOn(Bool_t t)                   { fCentSelectOn = t;}
  void SetOneConstSelectionOn(Bool_t t)                     { fOneConstSelectOn =t;}
  void SetMinCentrality(Float_t t)                          { fCentMin = t ; }
  void SetMaxCentrality(Float_t t)                          { fCentMax = t ; }
  void SetSemigoodCorrect(Int_t yesno)                      {fSemigoodCorrect=yesno;}
  void SetHolePos(Float_t poshole)                          { fHolePos = poshole; }
  void SetHoleWidth(Float_t holewidth)                      { fHoleWidth = holewidth; }
  void SetDerivativeSubtractionOrder(Int_t c)               {fDerivSubtrOrder = c;}
  void DoJetELoss(Bool_t b, Float_t f)                      {fScaleELoss = b; xfraction = f;}
  void AddMedScat(Bool_t b, Float_t f, Int_t n)             {fAddMedScat = b; fAddMedScatPtFrac = f; fAddMedScatN = n;}

  static AliAnalysisTaskEmcalJetShapesMC* AddTaskJetShapesMC(
	 const char * njetsBase,
         const Double_t jetradius,
         const Double_t subjetradius,
         const char *ntracksPartLevel,
         const char *type,
         const char *CentEst,
         Int_t       pSel,
         TString     trigClass      = "",
         TString     kEmcalTriggers = "",
         TString     tag            = "",
         const char *rhoName = "",
         AliAnalysisTaskEmcalJetShapesMC::JetShapeType jetShapeType=AliAnalysisTaskEmcalJetShapesMC::kGenShapes,
         AliAnalysisTaskEmcalJetShapesMC::JetShapeSub jetShapeSub=AliAnalysisTaskEmcalJetShapesMC::kNoSub,
         AliAnalysisTaskEmcalJetShapesMC::JetSelectionType jetSelection=AliAnalysisTaskEmcalJetShapesMC::kInclusive,
         Float_t minpTHTrigger =0.,
	 Float_t maxpTHTrigger =0.,
         AliAnalysisTaskEmcalJetShapesMC::DerivSubtrOrder derivSubtrOrder = AliAnalysisTaskEmcalJetShapesMC::kSecondOrder );

 protected:

  Bool_t                              RetrieveEventObjects();
  Bool_t                              Run();
  Bool_t                              FillHistograms();

  Float_t                            GetJetMass(AliEmcalJet *jet,Int_t jetContNb=0);
  Float_t                            Angularity(AliEmcalJet *jet, Int_t jetContNb=0);
  Float_t                            GetJetAngularity(AliEmcalJet *jet, Int_t jetContNb=0);
  Float_t                            PTD(AliEmcalJet *jet, Int_t jetContNb=0);
  Float_t                            GetJetpTD(AliEmcalJet *jet, Int_t jetContNb=0);
  Float_t                            Circularity(AliEmcalJet *jet, Int_t jetContNb=0);
  Float_t                            GetJetCircularity(AliEmcalJet *jet, Int_t jetContNb=0);
  Float_t                            LeSub(AliEmcalJet *jet, Int_t jetContNb=0);
  Float_t                            GetJetLeSub(AliEmcalJet *jet, Int_t jetContNb=0);
  Float_t                            GetJetNumberOfConstituents(AliEmcalJet *jet,Int_t jetContNb=0);
  Float_t                            GetSigma2(AliEmcalJet *jet, Int_t jetContNb=0);
  Float_t                            Sigma2(AliEmcalJet *jet, Int_t jetContNb=0);
  void                               NTValues(AliEmcalJet *jet, Int_t jetContNb, Float_t* nTFractions);
  void                                SoftDrop(AliEmcalJet *fJet,AliJetContainer *fJetCont, Double_t zcut, Double_t beta, Int_t ReclusterAlgo);
  void                                RecursiveParents(AliEmcalJet *fJet,AliJetContainer *fJetCont, Int_t ReclusterAlgo,Float_t PartonFlavor);
  AliEmcalJetFinder*                 Recluster(AliEmcalJet *Jet, Int_t JetContNb, Double_t JetRadius, Double_t SubJetRadius, Double_t SubJetMinPt, Int_t Algorithm, const char* Name);

  //Double_t                           NSubJettiness(AliEmcalJet *Jet, Int_t JetContNb,  AliEmcalJetFinder *Reclusterer, Int_t N, Int_t A, Int_t B);
  Double_t                           SubJetOrdering(AliEmcalJet *Jet, AliEmcalJetFinder *Reclusterer, Int_t N, Int_t Type, Bool_t Index);

  Double_t                           GetSubjetFraction(AliEmcalJet *Jet, Int_t JetContNb, Double_t JetRadius,  AliEmcalJetFinder *Reclusterer);

  Float_t                            CoreFrac(AliEmcalJet *jet, Int_t jetContNb=0);
  Float_t                            GetJetCoreFrac(AliEmcalJet *jet, Int_t jetContNb=0);


  Double_t                           FjNSubJettiness(AliEmcalJet *Jet, Int_t JetContNb, Int_t N, Int_t Algorithm, Double_t Beta, Int_t Option, Double_t Beta_SD=0.0, Double_t ZCut=0.1, Int_t SoftDropOn=0);


  Int_t                              SelectTrigger(Float_t minpT, Float_t maxpT);
  Double_t                           RelativePhi(Double_t mphi, Double_t vphi);

  Int_t                               fContainer;              // jets to be analyzed 0 for Base, 1 for subtracted.
  Float_t                             fMinFractionShared;          // only fill histos for jets if shared fraction larger than X
  JetShapeType                        fJetShapeType;               // jet type to be used
  JetShapeSub                         fJetShapeSub;                // jet subtraction to be used
  JetSelectionType                    fJetSelection;               // Jet selection: inclusive/recoil jet
  Float_t                             fShapesVar[13];              //     jet shapes used for the tagging
  Float_t                             fPtThreshold;
  Float_t                             fRMatching;
  Float_t                             fJetRadius;
  Float_t                             fSubjetRadius;
  Int_t                               fSelectedShapes;                //chose set of shapes
  Int_t                               fSwitchKtNSub;
  Int_t                               fSwitchMinNSub;
  Int_t                               fSwitchAktNSub;
  Int_t                               fSwitchSDKtNSub;
  Int_t                               fSwitchSDMinNSub;
  Int_t                               fAdditionalTracks;             //number of extra tracks to stress the grooming
  Float_t                             fHardCutoff;                  //hard cutoff in the iterative procedure
  Int_t                               fOptionalPartonInfo;          //parton info only valid for pythia6
  Float_t                             fminpTTrig;                   //min - max pT for trigger particle in case of recoil jet
  Float_t                             fmaxpTTrig;
  Float_t                             fangWindowRecoil;             //angular window for btb recoil analysis
  Int_t                               fSemigoodCorrect;             //if==1 we run over semigood runs
  Float_t                             fHolePos;                          //position in radians of the bad TPC sector
  Float_t                             fHoleWidth;                       //width of the hole in radians
 TRandom3                             *fRandom;                     //! Random number generator
 Float_t                              fqhat;                        //qhat
 Float_t                              fxlength;                     //medium length
 Bool_t                              fCentSelectOn;                // switch on/off centrality selection
  Float_t                             fCentMin;                     // min centrality value
  Float_t                             fCentMax;                     // max centrality value
  Bool_t                              fOneConstSelectOn;                // switch on/off one constituent selection
  Int_t                               fDerivSubtrOrder;
  Bool_t                              fScaleELoss;
  Double_t                            xfraction;
  Bool_t                              fAddMedScat;
  Double_t                            fAddMedScatPtFrac;
  Double_t                            fAddMedScatN;

  TH2F                                *fPhiJetCorr6;//
  TH2F                                *fPhiJetCorr7;//
  TH2F                                *fEtaJetCorr6;//
  TH2F                                *fEtaJetCorr7;//
  TH2F                                *fPtJetCorr;//
  TH1F                                *fPtJet;//
  TH2F                                *fhpTjetpT; //
  TH1F                                *fhPt;//
  TH1F                                *fhPhi;//
  THnSparse                           *fHLundIterative;//
  THnSparse                           *fHLundIterative_ktaxis;//
  THnSparse                           *fHLundIterativeInject;//
  TH2F                                *fNbOfConstvspT;//

  TTree           *fTreeObservableTagging;//!<! Tree with tagging variables subtracted MC or true MC or raw
  TF1             *fTf1Omega;//! to generate omega according to BDMPS tail
  TF1             *fTf1Kt; //!  to generate kT according to BDMS tail

  TF1             *fTf1SoftOmega;//! to generate omega for soft background
  TF1             *fTf1SoftKt; //!  to generate kT for soft background
 private:
  AliAnalysisTaskEmcalJetShapesMC(const AliAnalysisTaskEmcalJetShapesMC&);            // not implemented
  AliAnalysisTaskEmcalJetShapesMC &operator=(const AliAnalysisTaskEmcalJetShapesMC&); // not implemented

  ClassDef(AliAnalysisTaskEmcalJetShapesMC, 2);
};
#endif
