#ifndef AliAnalysisTaskCharm_H
#define AliAnalysisTaskCharm_H

//########################################################################
//#                                                                      #
//#   Class to produce charm cocktail (analysis task)                    #
//#                                                                      #
//#  Authors:                                                            #
//#   Raphaelle Bailhache, Uni Frankfurt / Raphaelle.Bailhache@cern.ch   #
//#   Carsten Klein, Uni Frankfurt / Carsten.Klein@cern.ch               #
//#   Jerome Jung, Uni Frankfurt / s8286523@uni-frankfurt.de             #
//#   Sebastian Scheid, Uni Frankfurt / s.scheid@cern.ch                 #
//#                                                                      #
//########################################################################

#include "AliAnalysisTaskSE.h"
// local files
#include "AliCocktailSmearing.h"
class TH1F;
class TH2F;
class TH3F;
class TList;
class TGraph;
class AliMCEvent;
class AliInputEventHandler;

class AliAnalysisTaskCharm : public AliAnalysisTaskSE, public AliCocktailSmearing {
  
public:
  AliAnalysisTaskCharm(); ///< default constructor probably needed for AnalysisManager or such...
  AliAnalysisTaskCharm(const Char_t* name); ///< named constructor which also creates input and output objects.
  virtual ~AliAnalysisTaskCharm();
  
  virtual void UserCreateOutputObjects();
  virtual void Init();
  virtual void LocalInit(){Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);
  
  void         SetProcessType(Int_t processType)            { fProcessType = processType; }
  void         SetPtCutHigh(Double_t max)                   { fPtCutHigh = max; }
  void         SetPtCutLow(Double_t min)                    { fPtCutLow = min; }
  void         SetEtaCutMinMax(Double_t min, Double_t max)  { fEtamin = min; fEtamax = max;}
  void         ScaleByRAA(Bool_t b)                         { fScaleByRAA = b; }
  void         ScaleByCNM(Bool_t b,TGraph *cnmgraph)        { fScaleByCNM = b; fgraphCNM = cnmgraph;}
  void         TakeptOfDCNM(Bool_t b)                       { fTakeptOfDCNM = b; }
  void         SetNbEvent(Int_t b)                          { fNbEvent = b;  }
  void         SetApplywm(Bool_t b)                         { fApplywm = b;  }
  void         SetApplyEventw(Bool_t b)                     { fEventWeight = b;  }
  void         Selectoneccbar(Bool_t b)                     { fSelectoneccbar = b;  }
  void         Selectcleanhistory(Bool_t b)                 { fSelectcleanhistory = b;  }
  
  
  
private:
  void           CreateHistos();
  Double_t       pt_cut200(Double_t pT);
  Double_t       pt_cut300(Double_t pT);
  Double_t       pt_cut400(Double_t pT);
  Double_t       pt_cutHigh(Double_t pT);
  Double_t       pt_cutLow(Double_t pT);
  Double_t       scale_RAA(Double_t pT);
  Double_t       scale_CNM(Double_t pT);
  Bool_t         Inphenixacc(Double_t phi, Double_t pt, Int_t pdg);
  
  
protected:
  AliMCEvent*              fMcEvent;    //! MC event
  AliInputEventHandler*    fMcHandler;  //! MC EventHandler
  Int_t                    fNbEvents;   // number of events
  Int_t                    fProcessType; // Select the process type 
  Double_t                 fPtCutHigh;
  Double_t                 fPtCutLow;
  Double_t                 fEtamin;
  Double_t                 fEtamax;
  Bool_t                   fScaleByRAA;
  Bool_t                   fScaleByCNM;
  TGraph                  *fgraphCNM;
  Bool_t                   fTakeptOfDCNM;
  Bool_t                   fApplywm;
  Bool_t                   fEventWeight;
  Bool_t                   fSelectoneccbar;
  Bool_t                   fSelectcleanhistory;
  Int_t                    fNbEvent;
  Int_t                    fNbEventCounter;
  // Histogram with 4 bins: for # of events 
  TH1F *hNEvents;                                      //!
  TH1F *hNEventsW;                                     //!
  // Histograms for mesons including charm quark
  TH2F *hQuarkMethod1;                                 //!
  TH2F *hQuarkMethod2;                                 //!
  TH2F *hCharm;                                        //!
  TH2F *hDpm;                                          //!
  TH2F *hD0;                                           //!  
  TH2F *hDs;                                           //!
  TH2F *hDstar;                                        //!
  TH2F *hLambdac;                                      //!
  // for comparisons to LHCb , 2.0 < y < 4.5 
  TH1F *hDeltaPhi_D0;                                  //! Delta phi distribution
  TH1F *hDeltaPhi_D0_LHCb;                             //! Delta phi distribution in the LHCb acceptance (3 > pT > 12 GeV/c and 2 < y < 4)
  TH1F *hDeltaPhi_D0_LHCb_e;                           //! Delta phi distribution in the LHCb acceptance (3 > pT > 12 GeV/c and 2 < y < 4)
  ///////////////////////////////////////////////////////
  // Histograms for rapidity distributions
  // rap and eta do not make any difference for electrons, nevertheless, wanted to have both cases in hRapElectron and hEtaElectron 
  // for comparisons to FONLL ,  pT > 0.2 and 0.5 GeV/c :
  TH1F *hRapElectronPt200;                             //!
  TH1F *hRapElectronPt500;                             //!
  ///////////////////////////////////////////////////////
  //
  TH1F *hPtCharmQuarkMethod1;                          //!
  TH1F *hPtCharmQuarkMethod1CNM;                       //!
  TH1F *hRapCharmQuarkMethod1;                         //!
  TH1F *hRapCharmQuarkMethod2;                         //!
  TH1F *hRapCharmQuarkMethod1Pair;                     //!
  TH1F *hRapCharmQuarkMethod2Pair;                     //!
  TH2F *hRapCharmQuark2DMethod1;                       //!
  TH2F *hRapCharmQuark2DMethod2;                       //!
  TH1F *hPhiCharmQuark1DMethod1;                       //!
  TH1F *hPhiCharmQuark1DMethod2;                       //!
  TH2F *hPtCharmQuark2DMethod1;                        //!
  TH2F *hPtCharmQuark2DMethod2;                        //!
  ///////////////////////////////////////////////////////
  TH2F *hRapCquarkDmeson;                              //!
  TH2F *hPtCquarkDmeson;                               //!
  TH2F *hRapDmesonElectron;                            //!
  TH2F *hPtDmesonElectron;                             //!
  // hRapDmesonElectron is empty, so, wanted to fill it in another way  
  TH2F *hRapDMom2e;                                    //!
  // generated pt vs eta of electrons from charm
  TH2F *hPtEtaElectron;                                //!
  TH2F *hPtEtaElectronE;                               //!
  TH2F *hPtEtaElectronP;                               //!
  // Correlation in eta of the positron and electron
  TH2F *hEtaPositronElectron;                          //!
  TH2F *hEtaPositronElectronDelta;                     //! 
  TH1F *hPhiPositronElectronDelta;                     //!
  TH1F *hPhiPositronElectronDeltaonehighDmeson;        //!
  TH1F *hPhiPositronElectronDeltabothhighDmeson;       //!
  // Histograms for Pt spectra : c-->e , cBar->e 
  TH1F *hPte_eta08;                                    //!
  TH1F *hPteP_eta08;                                   //!
  TH1F *hPteM_eta08;                                   //!
  TH1F *hPte_y08;                                      //!
  TH1F *hPteP_y08;                                     //!
  TH1F *hPteM_y08;                                     //!
  // Histograms for invariant mass spectra (ULS,LS): c-->e , cBar->e
  TH1F *hMee_ULS_simulated;                            //!
  TH1F *hMee_LS_simulated;                             //!
  TH2F *hMeePtee_ULS_eta05;                            //!
  TH2F *hMeePtee_LS_eta05;                             //!
  TH2F *hMeePtee_ULS_eta035;                           //!
  TH2F *hMeePtee_LS_eta035;                            //!
  TH3F *hPhiee_ULS_eta08_pt200;                        //!
  TH3F *hPhiee_LS_eta08_pt200;                         //!
  TH3F *hPhiee_ULS_eta08_pt400;                        //!
  TH3F *hPhiee_LS_eta08_pt400;                         //!
  TH2F *hMeePtee_ULS_eta035_phenixacc;                 //!
  TH2F *hMeePtee_LS_eta035_phenixacc;                  //!
  TH3F *hPhiee_ULS_eta08_highoneDmeson_pt200;          //!
  TH3F *hPhiee_LS_eta08_highoneDmeson_pt200;           //!
  TH3F *hPhiee_ULS_eta08_highbothDmeson_pt200;         //!
  TH3F *hPhiee_LS_eta08_highbothDmeson_pt200;          //!
  TH3F *hPhiee_ULS_eta08_highoneDmeson_pt400;          //!
  TH3F *hPhiee_LS_eta08_highoneDmeson_pt400;           //!
  TH3F *hPhiee_ULS_eta08_highbothDmeson_pt400;         //!
  TH3F *hPhiee_LS_eta08_highbothDmeson_pt400;          //!
  TH1F *hMee_ULS_eta08;                                //!
  TH1F *hMee_LS_eta08;                                 //!
  TH1F *hMee_ULS_eta08_pt200;                          //!
  TH1F *hMee_LS_eta08_pt200;                           //!
  // for 1D histograms, one pt-cut (200 MeV) is enough to do checks, others can be projected from 2D...
  //
  // Histograms for invariant mass vs pair pt (ULS,LS): c-->e , cBar->e
  TH2F *hMeePtee_ULS_eta1;                             //!
  TH2F *hMeePtee_LS_eta1;                              //!
  TH2F *hMeePtee_ULS_eta08;                            //!
  TH2F *hMeePtee_LS_eta08;                             //!
  TH2F *hMeePtee_ULS_eta08_pt200;                      //!
  TH2F *hMeePtee_LS_eta08_pt200;                       //!
  TH2F *hMeePtee_ULS_eta_pt;                           //!
  TH2F *hMeePtee_LS_eta_pt;                            //!
  TH2F *hMeePtee_ULS_eta08_pt300;                      //!
  TH2F *hMeePtee_LS_eta08_pt300;                       //!
  TH2F *hMeePtee_ULS_eta08_pt400;                      //!
  TH2F *hMeePtee_LS_eta08_pt400;                       //!
  TH2F *hMeePtee_ULS_eta08_pt200_opAngle50;            //!
  TH2F *hMeePtee_LS_eta08_pt200_opAngle50;             //!
  TH2F *hMeePtee_ULS_eta08_pt300_opAngle50;            //!
  TH2F *hMeePtee_LS_eta08_pt300_opAngle50;             //!
  TH2F *hMeePtee_ULS_eta08_pt400_opAngle50;            //!
  TH2F *hMeePtee_LS_eta08_pt400_opAngle50;             //!
  // opening angle
  TH2F *hMeeOpAngle_ULS_eta08_pt200;                   //!
  TH2F *hMeeOpAngle_LS_eta08_pt200;                    //!
  // Mother D meson pt
  TH3F *hMotherPt_ULS_eta08_pt200;                     //!
  TH3F *hMotherPt_LS_eta08_pt200;                      //!
  TH3F *hMotherPt_ULS_eta08_pt400;                     //!
  TH3F *hMotherPt_LS_eta08_pt400;                      //!
  // monitoring cnm weights
  TH2F *hweightcnmD;                                   //!
  TH2F *hweightcnmHFE;                                 //! 
  TH2F *hweightcnmee;                                  //!
  
  TList       *fOutputList;                            //! Output list

  AliAnalysisTaskCharm(const AliAnalysisTaskCharm &c); // not implemented
  AliAnalysisTaskCharm& operator= (const AliAnalysisTaskCharm &c); // not implemented
  
  ClassDef(AliAnalysisTaskCharm,3)
};

#endif

