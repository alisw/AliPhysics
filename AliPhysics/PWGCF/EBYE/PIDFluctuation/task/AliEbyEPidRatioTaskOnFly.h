#ifndef AliEbyEPidRatioTaskOnFly_cxx
#define AliEbyEPidRatioTaskOnFly_cxx

//=========================================================================//
//                                                                         //
//           Analysis Task for Particle Ratio Fluctuaions                  //
//              Author: Satyajit Jena || Deepika Jena                      //
//                      sjena@cern.ch || drathee@cern.ch                   //
//                                                                         //
//=========================================================================//


class TH1D;
class TH2F;
class TH3F;
class TString;
class TList;


#include "AliAnalysisTaskSE.h"


class AliEbyEPidRatioTaskOnFly: public AliAnalysisTaskSE {
 public:
  AliEbyEPidRatioTaskOnFly( const char *name = "HigherMomentAnalysis");
  virtual ~AliEbyEPidRatioTaskOnFly();

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

    
  void SetKinematicsCutsAOD(Double_t ptl, Double_t pth, Double_t eta){
    fPtLowerLimit   = ptl;
    fPtHigherLimit  = pth;
    fEtaLowerLimit  = -eta;
    fEtaHigherLimit = eta;
  }

  static const Char_t* fgkPidName[4];
  static const Char_t* fgkPidLatex[4][2];
  static const Char_t* fgkPidTitles[4][2];

  void FillHistSetCent();

 private:

  TList *fThnList;//!
  
 
  Double_t fPtLowerLimit; //
  Double_t fPtHigherLimit;//
  Double_t fEtaLowerLimit;//
  Double_t fEtaHigherLimit;//
  Int_t fCentrality; //

  Int_t      fOrder;                 //  Max order of higher order distributions
  Double_t **fRedFactp;              //!  Array of particle/anti-particle reduced factorial
  Int_t      fNp[4][2];                    //  Array of particle/anti-particle counts
  
  
  AliEbyEPidRatioTaskOnFly(const AliEbyEPidRatioTaskOnFly&);
  AliEbyEPidRatioTaskOnFly& operator = (const AliEbyEPidRatioTaskOnFly&);//Not implimented..
  ClassDef(AliEbyEPidRatioTaskOnFly, 1);

};

#endif

 
