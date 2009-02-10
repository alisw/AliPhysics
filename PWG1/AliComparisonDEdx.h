#ifndef ALICOMPARISONDEdx_H
#define ALICOMPARISONDEdx_H

//------------------------------------------------------------------------------
// Class to keep information from comparison of 
// reconstructed and MC particle tracks (TPC dE/dx).   
// 
// Author: J.Otwinowski 04/02/2008 
//------------------------------------------------------------------------------

class TFile;
class AliMCInfo;
class AliESDRecInfo;
class AliESDEvent; 
class AliESD;
class AliESDfriend;
class AliRecInfoCuts;
class AliMCInfoCuts;
class TH1I;
class TH3F;
class TH3;
class TProfile;
class TProfile2D;
class TGraph2D;
class TGraph; 
class TString;

#include "TNamed.h"
#include "AliComparisonObject.h"

//class AliComparisonDEdx : public TNamed {
class AliComparisonDEdx : public AliComparisonObject {
public :
  AliComparisonDEdx(); 
  ~AliComparisonDEdx();

  // Init data members
  virtual void Init();

  // Execute analysis
  virtual void Exec(AliMCInfo* infoMC, AliESDRecInfo *infoRC);

  // Merge output objects (needed by PROOF) 
  virtual Long64_t Merge(TCollection* list);

  // Analyse output histograms
  virtual void Analyse();

  // Get analysis folder
  virtual TFolder* GetAnalysisFolder() {return fAnalysisFolder;}

  // Create folder for analysed histograms
  TFolder *CreateFolder(TString folder = "folderDEdx",TString title = "Analysed DEdx histograms");

  // Export objects to folder
  TFolder *ExportToFolder(TObjArray * array=0);

  // Process events
  void      Process(AliMCInfo* infoMC, AliESDRecInfo *infoRC);

  // Selection cuts
  void SetAliRecInfoCuts(AliRecInfoCuts* cuts=0) {fCutsRC = cuts;}
  void SetAliMCInfoCuts(AliMCInfoCuts* cuts=0)   {fCutsMC = cuts;} 

  void SetMCPtMin(const Float_t cuts=0) {fMCPtMin = cuts;} 
  void SetMCAbsTanThetaMax(const Float_t cuts=1e99) {fMCAbsTanThetaMax = cuts;} 
  void SetMCPdgCode(const Int_t cuts=0) {fMCPdgCode = cuts;} 

  AliRecInfoCuts*  GetAliRecInfoCuts() const {return fCutsRC;}      
  AliMCInfoCuts*   GetAliMCInfoCuts()  const {return fCutsMC;}     
  Float_t GetMCPtMin() const {return fMCPtMin;}
  Float_t GetMCAbsTanThetaMax() const {return fMCAbsTanThetaMax;}
  Int_t GetMCPdgCode() const {return fMCPdgCode;} 

  static TH1F*     MakeResol(TH2F * his, Int_t integ, Bool_t type); 

  //
  // TPC dE/dx 
  TH2F* GetTPCSignalNormTan() {return fTPCSignalNormTan;}
  TH2F* GetTPCSignalNormSPhi() {return fTPCSignalNormSPhi;}
  TH2F* GetTPCSignalNormTPhi() {return fTPCSignalNormTPhi;}
  //
  TH3F* GetTPCSignalNormTanSPhi() {return fTPCSignalNormTanSPhi;}
  TH3F* GetTPCSignalNormTanTPhi() {return fTPCSignalNormTanTPhi;}
  TH3F* GetTPCSignalNormTanSPt() {return fTPCSignalNormTanSPt;}
  

private:

  // TPC dE/dx 
  TH2F* fTPCSignalNormTan;    //-> TPC signal normalized to the calculated MC signal 
  TH2F* fTPCSignalNormSPhi;   //-> TPC signal normalized to the calculated MC signal
  TH2F* fTPCSignalNormTPhi;   //-> TPC signal normalized to the calculated MC signal
  //
  TH3F* fTPCSignalNormTanSPhi;   //-> TPC signal normalized to the calculated MC signal
  TH3F* fTPCSignalNormTanTPhi;   //-> TPC signal normalized to the calculated MC signal
  TH3F* fTPCSignalNormTanSPt;    //-> TPC signal normalized to the calculated MC signal
  
  // Selection cuts
  AliRecInfoCuts*  fCutsRC; // selection cuts for reconstructed tracks
  AliMCInfoCuts*   fCutsMC; // selection cuts for MC tracks

  Float_t fMCPtMin;               // min. MC pt cut
  Float_t fMCAbsTanThetaMax;      // max. MC abs[tan(theta)] cut
  Int_t fMCPdgCode;               // selected particle with Pdg code

  // analysis folder 
  TFolder *fAnalysisFolder; // folder for analysed histograms

  AliComparisonDEdx(const AliComparisonDEdx&); // not implemented
  AliComparisonDEdx& operator=(const AliComparisonDEdx&); // not implemented

  ClassDef(AliComparisonDEdx,1);
};

#endif
