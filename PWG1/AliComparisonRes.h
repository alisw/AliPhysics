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
class TString;
class AliESDVertex;

#include "TNamed.h"
#include "AliComparisonObject.h"

class AliComparisonRes : public AliComparisonObject {
public :
  AliComparisonRes(); 
  virtual ~AliComparisonRes();

  // Init data members
  virtual void     Init();

  // Execute analysis
  virtual void      Exec(AliMCInfo* infoMC, AliESDRecInfo *infoRC);

  // Merge output objects (needed by PROOF) 
  virtual Long64_t Merge(TCollection* list);

  // Analyse output histograms
  virtual void Analyse();

  // Get analysis folder
  virtual TFolder* GetAnalysisFolder() {return fAnalysisFolder;}

  // Process events
  void      ProcessConstrained(AliMCInfo* infoMC, AliESDRecInfo *infoRC);
  void      Process(AliMCInfo* infoMC, AliESDRecInfo *infoRC);

  // Create folder for analysed histograms
  TFolder *CreateFolder(TString folder = "folderRes",TString title = "Analysed Resolution histograms");

  // Export objects to folder
  TFolder *ExportToFolder(TObjArray * array=0);

  // Selection cuts
  void SetAliRecInfoCuts(AliRecInfoCuts* cuts=0) {fCutsRC = cuts;}   
  void SetAliMCInfoCuts(AliMCInfoCuts* cuts=0) {fCutsMC = cuts;}  
   
  AliRecInfoCuts*  GetAliRecInfoCuts() const {return fCutsRC;}  
  AliMCInfoCuts*   GetAliMCInfoCuts()  const {return fCutsMC;}  

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

  TH2F* f1Pt2ResolS1PtTPC;      //-> (1/mcpt-1/pt)/(1+1/pt)^2 vs sqrt(1/pt) (TPC)
  TH2F* f1Pt2ResolS1PtTPCITS;   //-> (1/mcpt-1/pt)/(1+1/pt)^2 vs sqrt(1/pt) (TPC+ITS)
  TH2F* fYResolS1PtTPC;         //-> (mcy-y)/(0.2+1/pt) vs sqrt(1/pt) (TPC) 
  TH2F* fYResolS1PtTPCITS;      //-> (mcy-y)/(0.2+1/pt) vs sqrt(1/pt) (TPC + ITS) 
  TH2F* fZResolS1PtTPC;         //-> (mcz-z)/(0.2+1/pt) vs sqrt(1/pt) (TPC)
  TH2F* fZResolS1PtTPCITS;      //-> (mcz-z)/(0.2+1/pt) vs sqrt(1/pt) (TPC+ITS)
  TH2F* fPhiResolS1PtTPC;       //-> (mcphi-phi)/(0.1+1/pt) vs sqrt(1/pt) (TPC)
  TH2F* fPhiResolS1PtTPCITS;    //-> (mcphi-phi)/(0.1+1/pt) vs sqrt(1/pt) (TPC+ITS)
  TH2F* fThetaResolS1PtTPC;     //-> (mctheta-theta)/(0.1+1/pt) vs sqrt(1/pt) (TPC)
  TH2F* fThetaResolS1PtTPCITS;  //-> (mctheta-theta)/(0.1+1/pt) vs sqrt(1/pt) (TPC+ITS)
  
  // constrained
  TH2F* fC1Pt2ResolS1PtTPC;      //-> (1/mcpt-1/pt)/(1+1/pt)^2 vs sqrt(1/pt) (TPC)
  TH2F* fC1Pt2ResolS1PtTPCITS;   //-> (1/mcpt-1/pt)/(1+1/pt)^2 vs sqrt(1/pt) (TPC+ITS)
  TH2F* fCYResolS1PtTPC;         //-> (mcy-y)/(0.2+1/pt) vs sqrt(1/pt) (TPC) 
  TH2F* fCYResolS1PtTPCITS;      //-> (mcy-y)/(0.2+1/pt) vs sqrt(1/pt) (TPC + ITS) 
  TH2F* fCZResolS1PtTPC;         //-> (mcz-z)/(0.2+1/pt) vs sqrt(1/pt) (TPC)
  TH2F* fCZResolS1PtTPCITS;      //-> (mcz-z)/(0.2+1/pt) vs sqrt(1/pt) (TPC+ITS)
  TH2F* fCPhiResolS1PtTPC;       //-> (mcphi-phi)/(0.1+1/pt) vs sqrt(1/pt) (TPC)
  TH2F* fCPhiResolS1PtTPCITS;    //-> (mcphi-phi)/(0.1+1/pt) vs sqrt(1/pt) (TPC+ITS)
  TH2F* fCThetaResolS1PtTPC;     //-> (mctheta-theta)/(0.1+1/pt) vs sqrt(1/pt) (TPC)
  TH2F* fCThetaResolS1PtTPCITS;  //-> (mctheta-theta)/(0.1+1/pt) vs sqrt(1/pt) (TPC+ITS)

  AliESDVertex *fVertex;  //! 

  // Global cuts objects
  AliRecInfoCuts*  fCutsRC;      // selection cuts for reconstructed tracks
  AliMCInfoCuts*  fCutsMC;       // selection cuts for MC tracks

  // analysis folder 
  TFolder *fAnalysisFolder; // folder for analysed histograms

  AliComparisonRes(const AliComparisonRes&); // not implemented
  AliComparisonRes& operator=(const AliComparisonRes&); // not implemented

  ClassDef(AliComparisonRes,1);
};

#endif
