#ifndef AliEbyEPidTaskFastGen_cxx
#define AliEbyEPidTaskFastGen_cxx

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
#include "THnSparse.h"



class AliEbyEPidTaskFastGen: public AliAnalysisTaskSE {
 public:
  AliEbyEPidTaskFastGen( const char *name = "HigherMomentAnalysis");
  virtual ~AliEbyEPidTaskFastGen();

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

  void SetVertexDiamond(Double_t vx, Double_t vy, Double_t vz) {fVxMax = vx;fVyMax = vy; fVzMax = vz;}
  
  
  void SetKinematicsCutsAOD(Double_t ptl, Double_t pth, Double_t eta){
    fPtLowerLimit   = ptl;
    fPtHigherLimit  = pth;
    fEtaLowerLimit  = -eta;
    fEtaHigherLimit = eta;
  }

   
 enum ESparseData_t {
    kAll_part=0,
    kCent_imp=1,
    kN_part=2,
    kN_ch=3,
    kNch_plus,
    kNch_minus,
    kN_pi,
    kNpi_plus,
    kNpi_minus,
    kN_k,
    kNka_plus,
    kNka_minus,
    kN_pr,
    kNpr_plus,
    kNpr_minus,
    kNSparseData
  };
  
 private:

  TList *fThnList;
  
  Double_t fVxMax;               //vxmax
  Double_t fVyMax;//vymax
  Double_t fVzMax;//vzmax
  Double_t fPtLowerLimit;
  Double_t fPtHigherLimit;
  Double_t fEtaLowerLimit;
  Double_t fEtaHigherLimit;
    
  THnSparseI *fHistoCorrelationMC; 

  AliEbyEPidTaskFastGen(const AliEbyEPidTaskFastGen&);
  AliEbyEPidTaskFastGen& operator = (const AliEbyEPidTaskFastGen&);//Not implimented..
  ClassDef(AliEbyEPidTaskFastGen, 1);

};

#endif

 
