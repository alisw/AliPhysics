#ifndef AliEbyEPidQATaskOnFlyKine_cxx
#define AliEbyEPidQATaskOnFlyKine_cxx


//=========================================================================//
//             AliEbyE OnFLy QA Tasks for Charge and PID                   //
//                         For Testing Only                                //
//                   Satyajit Jena | sjena@cern.ch                         //
//=========================================================================//


class TH1D;
class TH2F;
class TH3F;
class TString;
class TList;


#include "AliAnalysisTaskSE.h"


class AliEbyEPidQATaskOnFlyKine: public AliAnalysisTaskSE {
 public:
  AliEbyEPidQATaskOnFlyKine( const char *name = "HigherMomentAnalysis");
  virtual ~AliEbyEPidQATaskOnFlyKine();

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

    
  void SetKinematicCut(Double_t eta, Double_t pt, Double_t vx) {
    fEtaCut = eta;
    fPtCut  = pt;
    fVtxZ   = vx;
  }

  static const Char_t* fgkPidName[4];
  static const Char_t* fgkPidLatex[4][2];
    
 private:

  TList *fThnList;         //!
  Int_t fCentrality; //
  Double_t fEtaCut;       //!
  Double_t fPtCut;        //!
  Double_t fVtxZ;
  TH2F *fHistImpNpart; //
  TH2F *fHistImpMult; //
  TH2F *fHistNpartMult;  //
  TH1F *fHistStat; //
  
  TH2F *fHistPt[4][2]; //!
  TH2F *fHistEtaY[4]; //!
  TH2F *fHistPhi[4]; //!
  TH2F *fHistPhiPt[4]; //!
  TH2F *fHistMult[4][2]; //!
  TH2F *fHistMultTot[4][2]; //!

  AliEbyEPidQATaskOnFlyKine(const AliEbyEPidQATaskOnFlyKine&);
  AliEbyEPidQATaskOnFlyKine& operator = (const AliEbyEPidQATaskOnFlyKine&);//Not implimented..
  ClassDef(AliEbyEPidQATaskOnFlyKine, 1);

};

#endif

 
