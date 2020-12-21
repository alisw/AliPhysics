#ifndef AliAnalysisTaskBeauty_H
#define AliAnalysisTaskBeauty_H


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
class TList;
class AliMCEvent;
class AliInputEventHandler;

class AliAnalysisTaskBeauty : public AliAnalysisTaskSE, public AliCocktailSmearing {
  
public:
  AliAnalysisTaskBeauty(); ///< default constructor probably needed for AnalysisManager or such...
  AliAnalysisTaskBeauty(const Char_t* name); ///< named constructor which also creates input and output objects.
  virtual ~AliAnalysisTaskBeauty();
  
  virtual void UserCreateOutputObjects();
  virtual void Init();
  virtual void LocalInit(){Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);
  
  void         SetProcessType(Int_t processType)  { fProcessType = processType; }
  void         SetPtCutHigh(Double_t max)         { fPtCutHigh = max; }
  void         ScaleByRAA(Bool_t b)               { fScaleByRAA = b; }
  void         SetApplyEventw(Bool_t b)           { fEventWeight = b;  }
  
  
private:
  void           CreateHistos();
  Double_t       pt_cut200(double pT);
  Double_t       pt_cut300(double pT);
  Double_t       pt_cut400(double pT);
  Double_t       pt_cutHigh(Double_t pT);
  Double_t       scale_RAA(Double_t pT);
  Bool_t         Inphenixacc(Double_t phi, Double_t pt, Int_t pdg);
  
  
protected:
  AliMCEvent*              fMcEvent;    //! MC event
  AliInputEventHandler*    fMcHandler;  //! MC EventHandler
  Int_t                    fNbEvents;   // number of events
  Int_t                    fProcessType; // Select the process type
  Bool_t                   fEventWeight;
  Double_t                 fPtCutHigh;
  Bool_t                   fScaleByRAA;
  // Histogram with 4 bins: for # of events 
  TH1F *hNEvents;                                      //!
  TH1F *hNEventsW;                                     //!
  // Histograms for mesons including charm quark
  TH2F *hQuarkMethod1;                                 //!
  TH2F *hQuarkMethod2;                                 //!
  TH2F *hBeauty;                                       //!
  TH2F *hBpm;                                          //!
  TH2F *hB0;                                           //!
  TH2F *hBs;                                           //!
  TH2F *hLambdab;                                      //!
  ////  
  // for comparisons to FONLL ,  pT > 0.2 and 0.5 GeV/c :
  TH1F *hRapElectron_be_Pt200;                         //!
  TH1F *hRapElectron_be_Pt500;                         //!
  TH1F *hRapElectron_bce_Pt200;                        //!
  TH1F *hRapElectron_bce_Pt500;                        //!
  ///////////////////////////////////////////////////////
  TH1F *hRapBeautyQuarkMethod1;                        //!
  TH1F *hRapBeautyQuarkMethod2;                        //!
  TH2F *hRapBeautyQuark2DMethod1;                      //!
  TH2F *hRapBeautyQuark2DMethod2;                      //!
  ///////////////////////////////////////////////////////
  TH2F *hRapBMom2e;                                    //!
  TH2F *hRapBMom2c2e;                                  //!
  // generated pt vs eta of electrons from beauty
  TH2F *hPtEtaElectron_be;                             //!
  TH2F *hPtEtaElectron_bce;                            //!
  // Histograms for Pt spectra : b-->e , bBar->e
  TH1F *hPte_eta08_be;                                 //!
  TH1F *hPteP_eta08_be;                                //!
  TH1F *hPteM_eta08_be;                                //!
  TH1F *hPte_y08_be;                                   //!
  TH1F *hPteP_y08_be;                                  //!
  TH1F *hPteM_y08_be;                                  //!
  // Histograms for Pt spectra : b->c->e , bBar->cBar->e
  TH1F *hPte_eta08_bce;                                //!
  TH1F *hPteP_eta08_bce;                               //!
  TH1F *hPteM_eta08_bce;                               //!
  TH1F *hPte_y08_bce;                                  //!
  TH1F *hPteP_y08_bce;                                 //!
  TH1F *hPteM_y08_bce;                                 //!
  // Histograms (ULS,LS)for all combinations of  b->e and b->c->e
  TH1F *hMee_ULS_simulated;                            //!
  TH1F *hMee_LS_simulated;                             //!
  TH1F *hMee_ULS_eta05;                                //!
  TH1F *hMee_LS_eta05;                                 //!
  TH1F *hMee_ULS_eta08;                                //!
  TH1F *hMee_LS_eta08;                                 //!
  TH1F *hMee_ULS_eta035;                               //!
  TH1F *hMee_LS_eta035;                                //!
  TH1F *hMee_ULS_eta08_pt200;                          //!
  TH1F *hMee_LS_eta08_pt200;                           //!
  TH1F *hMee_ULS_eta08_pt400;                          //!
  TH1F *hMee_LS_eta08_pt400;                           //!
  TH1F *hMee_ULS_eta035_phenixacc;                     //!
  TH1F *hMee_LS_eta035_phenixacc;                      //!
  TH2F *hMeePtee_ULS_eta08;                            //!
  TH2F *hMeePtee_LS_eta08;                             //!
  TH2F *hMeePtee_ULS_eta08_pt200;                      //!
  TH2F *hMeePtee_LS_eta08_pt200;                       //!
  TH2F *hMeePtee_ULS_eta08_pt400;                      //!
  TH2F *hMeePtee_LS_eta08_pt400;                       //!
  TH2F *hMeePtee_ULS_eta08_pt200_opAngle50;            //!
  TH2F *hMeePtee_LS_eta08_pt200_opAngle50;             //!
  TH2F *hMeePtee_ULS_eta08_pt300_opAngle50;            //!
  TH2F *hMeePtee_LS_eta08_pt300_opAngle50;             //!
  TH2F *hMeePtee_ULS_eta08_pt400_opAngle50;            //!
  TH2F *hMeePtee_LS_eta08_pt400_opAngle50;             //!
  // Histograms (ULS,LS),  b-->e , bBar->e
  TH1F *hMee_ULS_simulated_be;                         //!
  TH1F *hMee_LS_simulated_be;                          //!
  TH1F *hMee_ULS_eta05_be;                             //!
  TH1F *hMee_LS_eta05_be;                              //!
  TH1F *hMee_ULS_eta08_be;                             //!
  TH1F *hMee_LS_eta08_be;                              //!
  TH1F *hMee_ULS_eta035_be;                            //!
  TH1F *hMee_LS_eta035_be;                             //!
  TH1F *hMee_ULS_eta08_pt200_be;                       //!
  TH1F *hMee_LS_eta08_pt200_be;                        //!
  TH1F *hMee_ULS_eta08_pt400_be;                       //!
  TH1F *hMee_LS_eta08_pt400_be;                        //!
  TH1F *hMee_ULS_eta035_phenixacc_be;                  //!
  TH1F *hMee_LS_eta035_phenixacc_be;                   //!
  TH2F *hMeePtee_ULS_eta08_be;                         //!
  TH2F *hMeePtee_LS_eta08_be;                          //!
  TH2F *hMeePtee_ULS_eta08_pt200_be;                   //!
  TH2F *hMeePtee_LS_eta08_pt200_be;                    //!
  TH2F *hMeePtee_ULS_eta08_pt400_be;                   //!
  TH2F *hMeePtee_LS_eta08_pt400_be;                    //!
  //Histograms (ULS,LS),  b->c->e , bBar->c->e
  TH1F *hMee_ULS_simulated_bce;                        //!
  TH1F *hMee_LS_simulated_bce;                         //!
  TH1F *hMee_ULS_eta05_bce;                            //!
  TH1F *hMee_LS_eta05_bce;                             //!
  TH1F *hMee_ULS_eta08_bce;                            //!
  TH1F *hMee_LS_eta08_bce;                             //!
  TH1F *hMee_ULS_eta035_bce;                           //!
  TH1F *hMee_LS_eta035_bce;                            //!
  TH1F *hMee_ULS_eta08_pt200_bce;                      //!
  TH1F *hMee_LS_eta08_pt200_bce;                       //!
  TH1F *hMee_ULS_eta08_pt400_bce;                      //!
  TH1F *hMee_LS_eta08_pt400_bce;                       //!
  TH1F *hMee_ULS_eta035_phenixacc_bce;                 //!
  TH1F *hMee_LS_eta035_phenixacc_bce;                  //!
  TH2F *hMeePtee_ULS_eta08_bce;                        //!
  TH2F *hMeePtee_LS_eta08_bce;                         //!
  TH2F *hMeePtee_ULS_eta08_pt200_bce;                  //!
  TH2F *hMeePtee_LS_eta08_pt200_bce;                   //!
  TH2F *hMeePtee_ULS_eta08_pt400_bce;                  //!
  TH2F *hMeePtee_LS_eta08_pt400_bce;                   //!
  // opening angle
  TH2F *hMeeOpAngle_ULS_eta08_pt200;                   //!
  TH2F *hMeeOpAngle_LS_eta08_pt200;                    //!
  
  TList       *fOutputList;                            //! Output list

  AliAnalysisTaskBeauty(const AliAnalysisTaskBeauty &c); // not implemented
  AliAnalysisTaskBeauty& operator= (const AliAnalysisTaskBeauty &c); // not implemented

  
  ClassDef(AliAnalysisTaskBeauty,1)
};

#endif

