#ifndef ALICOMPARISONRES_H
#define ALICOMPARISONRES_H

//------------------------------------------------------------------------------
// Class to keep information from comparison of 
// reconstructed and MC particle tracks (TPC resolution).   
// 
// Author: J.Otwinowski 04/02/2008 
//------------------------------------------------------------------------------

class TFile;
class AliMCInfo;
class AliESDRecInfo;
class AliESDEvent; 
class AliESD;
class AliESDfriend;
class AliMCInfoCuts;
class AliRecInfoCuts;
class TH1I;
class TH3F;
class TH3;
class TProfile;
class TProfile2D;
class AliESDVertex;

#include "TNamed.h"

class AliComparisonRes : public TNamed {
public :
  AliComparisonRes(); 
  virtual ~AliComparisonRes();
  void      InitHisto();
  void      InitCuts();
  void      Exec(AliMCInfo* infoMC, AliESDRecInfo *infoRC);
  void      ProcessConstrained(AliMCInfo* infoMC, AliESDRecInfo *infoRC);
  void      Process(AliMCInfo* infoMC, AliESDRecInfo *infoRC);

  // Selection cuts
  void SetAliRecInfoCuts(AliRecInfoCuts* cuts=0) {fCutsRC = cuts;}   
  void SetAliMCInfoCuts(AliMCInfoCuts* cuts=0) {fCutsMC = cuts;}  
   
  AliRecInfoCuts*  GetAliRecInfoCuts() const {return fCutsRC;}  
  AliMCInfoCuts*   GetAliMCInfoCuts()  const {return fCutsMC;}  

  // Merge output objects (needed by PROOF) 
  virtual Long64_t Merge(TCollection* list);

  // Analyse output histograms
  void Analyse();
  static TH1F*       MakeResol(TH2F * his, Int_t integ, Bool_t type); 


private:
  //
  // Control histograms
  //
  TH2F* fPtResolLPT;        //-> pt resolution - low pt
  TH2F* fPtResolHPT;        //-> pt resolution - high pt 
  TH2F* fPtPullLPT;         //-> pt resolution - low pt
  TH2F* fPtPullHPT;         //-> pt resolution - high pt 

  //
  // Resolution constrained param
  //
  TH2F   *fCPhiResolTan;   //-> angular resolution -  constrained
  TH2F   *fCTanResolTan;   //-> angular resolution -  constrained
  TH2F   *fCPtResolTan;    //-> pt resolution      -  constrained
  TH2F   *fCPhiPullTan;    //-> angular resolution -  constrained
  TH2F   *fCTanPullTan;    //-> angular resolution -  constrained
  TH2F   *fCPtPullTan;     //-> pt resolution      -  constrained

  //
  // Histograms for track resolution parameterisation
  //

  TH2F* f1Pt2Resol1PtTPC;      //-> (1/mcpt-1/pt)/(1+1/pt)^2 vs 1/pt (TPC)
  TH2F* f1Pt2Resol1PtTPCITS;   //-> (1/mcpt-1/pt)/(1+1/pt)^2 vs 1/pt (TPC+ITS)
  TH2F* fYResol1PtTPC;         //-> (mcy-y)/(0.2+1/pt) vs 1/pt (TPC) 
  TH2F* fYResol1PtTPCITS;      //-> (mcy-y)/(0.2+1/pt) vs 1/pt (TPC + ITS) 
  TH2F* fZResol1PtTPC;         //-> (mcz-z)/(0.2+1/pt) vs 1/pt (TPC)
  TH2F* fZResol1PtTPCITS;      //-> (mcz-z)/(0.2+1/pt) vs 1/pt (TPC+ITS)
  TH2F* fPhiResol1PtTPC;       //-> (mcphi-phi)/(0.1+1/pt) vs 1/pt (TPC)
  TH2F* fPhiResol1PtTPCITS;    //-> (mcphi-phi)/(0.1+1/pt) vs 1/pt (TPC+ITS)
  TH2F* fThetaResol1PtTPC;     //-> (mctheta-theta)/(0.1+1/pt) vs 1/pt (TPC)
  TH2F* fThetaResol1PtTPCITS;  //-> (mctheta-theta)/(0.1+1/pt) vs 1/pt (TPC+ITS)
  
  // constrained
  TH2F* fC1Pt2Resol1PtTPC;      //-> (1/mcpt-1/pt)/(1+1/pt)^2 vs 1/pt (TPC)
  TH2F* fC1Pt2Resol1PtTPCITS;   //-> (1/mcpt-1/pt)/(1+1/pt)^2 vs 1/pt (TPC+ITS)
  TH2F* fCYResol1PtTPC;         //-> (mcy-y)/(0.2+1/pt) vs 1/pt (TPC) 
  TH2F* fCYResol1PtTPCITS;      //-> (mcy-y)/(0.2+1/pt) vs 1/pt (TPC + ITS) 
  TH2F* fCZResol1PtTPC;         //-> (mcz-z)/(0.2+1/pt) vs 1/pt (TPC)
  TH2F* fCZResol1PtTPCITS;      //-> (mcz-z)/(0.2+1/pt) vs 1/pt (TPC+ITS)
  TH2F* fCPhiResol1PtTPC;       //-> (mcphi-phi)/(0.1+1/pt) vs 1/pt (TPC)
  TH2F* fCPhiResol1PtTPCITS;    //-> (mcphi-phi)/(0.1+1/pt) vs 1/pt (TPC+ITS)
  TH2F* fCThetaResol1PtTPC;     //-> (mctheta-theta)/(0.1+1/pt) vs 1/pt (TPC)
  TH2F* fCThetaResol1PtTPCITS;  //-> (mctheta-theta)/(0.1+1/pt) vs 1/pt (TPC+ITS)

  AliESDVertex *fVertex;  //! 

  // Global cuts objects
  AliRecInfoCuts*  fCutsRC;      // selection cuts for reconstructed tracks
  AliMCInfoCuts*  fCutsMC;       // selection cuts for MC tracks

  AliComparisonRes(const AliComparisonRes&); // not implemented
  AliComparisonRes& operator=(const AliComparisonRes&); // not implemented

  ClassDef(AliComparisonRes,1);
};

#endif
