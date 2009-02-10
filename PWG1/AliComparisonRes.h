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

  // getters
  TH3F *GetMCVertex()  { return fMCVertex; }
  TH3F *GetRecVertex() { return fRecVertex; } 

  TH3F *GetPhiTanPtTPC() { return fPhiTanPtTPC; } 
  TH3F *GetPhiTanPtTPCITS() { return fPhiTanPtTPCITS; } 

  TH2F *GetPtResolTPC()     { return fPtResolTPC; } 
  TH2F *GetPtPullTPC()      { return fPtPullTPC; } 
  TH2F *GetPhiResolTanTPC() { return fPhiResolTanTPC; } 
  TH2F *GetTanResolTanTPC() { return fTanResolTanTPC; } 
  TH2F *GetPhiPullTanTPC()  { return fPhiPullTanTPC; } 
  TH2F *GetTanPullTanTPC()  { return fTanPullTanTPC; } 

  TH2F *GetPtResolTPCITS()     { return fPtResolTPCITS; } 
  TH2F *GetPtPullTPCITS()      { return fPtPullTPCITS; } 
  TH2F *GetPhiResolTanTPCITS() { return fPhiResolTanTPCITS; } 
  TH2F *GetTanResolTanTPCITS() { return fTanResolTanTPCITS; } 
  TH2F *GetPhiPullTanTPCITS()  { return fPhiPullTanTPCITS; } 
  TH2F *GetTanPullTanTPCITS()  { return fTanPullTanTPCITS; } 

  //
  // Resolution constrained param
  //
  TH2F *GetCPhiResolTan() { return fCPhiResolTan; } 
  TH2F *GetCTanResolTan() { return fCTanResolTan; } 
  TH2F *GetCPtResolTan()  { return fCPtResolTan; } 
  TH2F *GetCPhiPullTan()  { return fCPhiPullTan; } 
  TH2F *GetCTanPullTan()  { return fCTanPullTan; } 
  TH2F *GetCPtPullTan()   { return fCPtPullTan; } 

  //
  // Histograms for track resolution parameterisation
  //
  TH2F *Get1Pt2ResolS1PtTPC()    { return f1Pt2ResolS1PtTPC; } 
  TH2F *Get1Pt2ResolS1PtTPCITS() { return f1Pt2ResolS1PtTPCITS; } 
  TH2F *GetYResolS1PtTPC()       { return fYResolS1PtTPC; } 
  TH2F *GetYResolS1PtTPCITS()    { return fYResolS1PtTPCITS; } 
  TH2F *GetZResolS1PtTPC()       { return fZResolS1PtTPC; } 
  TH2F *GetZResolS1PtTPCITS()    { return fZResolS1PtTPCITS; } 
  TH2F *GetPhiResolS1PtTPC()     { return fPhiResolS1PtTPC; } 
  TH2F *GetPhiResolS1PtTPCITS()  { return fPhiResolS1PtTPCITS; } 
  TH2F *GetThetaResolS1PtTPC()   { return fThetaResolS1PtTPC; } 
  TH2F *GetThetaResolS1PtTPCITS(){ return fThetaResolS1PtTPCITS; } 

  // constrained
  //
  TH2F *GetC1Pt2ResolS1PtTPC()    { return fC1Pt2ResolS1PtTPC; } 
  TH2F *GetC1Pt2ResolS1PtTPCITS() { return fC1Pt2ResolS1PtTPCITS; } 
  TH2F *GetCYResolS1PtTPC()       { return fCYResolS1PtTPC; } 
  TH2F *GetCYResolS1PtTPCITS()    { return fCYResolS1PtTPCITS; } 
  TH2F *GetCZResolS1PtTPC()       { return fCZResolS1PtTPC; } 
  TH2F *GetCZResolS1PtTPCITS()    { return fCZResolS1PtTPCITS; } 
  TH2F *GetCPhiResolS1PtTPC()     { return fCPhiResolS1PtTPC; } 
  TH2F *GetCPhiResolS1PtTPCITS()  { return fCPhiResolS1PtTPCITS; } 
  TH2F *GetCThetaResolS1PtTPC()   { return fCThetaResolS1PtTPC; } 
  TH2F *GetCThetaResolS1PtTPCITS(){ return fCThetaResolS1PtTPCITS; } 

private:
  //
  // Control histograms
  //

  TH3F *fMCVertex;  //-> MC primary vertex 
  TH3F *fRecVertex; //-> Reconstructed primary vertex

  TH3F *fPhiTanPtTPC; //-> phi vs tantheta vs pt
  TH3F *fPhiTanPtTPCITS; //-> phi vs tantheta vs pt

  // TPC only
  TH2F *fPtResolTPC;        //-> pt resolution
  TH2F *fPtPullTPC;         //-> pt pull
  TH2F *fPhiResolTanTPC;       //-> angular resolution 
  TH2F *fTanResolTanTPC;       //-> angular resolution
  TH2F *fPhiPullTanTPC;        //-> angular resolution
  TH2F *fTanPullTanTPC;        //-> angular resolution

  // TPC+ITS
  TH2F *fPtResolTPCITS;        //-> pt resolution
  TH2F *fPtPullTPCITS;         //-> pt pull
  TH2F *fPhiResolTanTPCITS;       //-> angular resolution 
  TH2F *fTanResolTanTPCITS;       //-> angular resolution
  TH2F *fPhiPullTanTPCITS;        //-> angular resolution
  TH2F *fTanPullTanTPCITS;        //-> angular resolution

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
