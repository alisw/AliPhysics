#ifndef ALICOMPARISONEFF_H
#define ALICOMPARISONEFF_H

//------------------------------------------------------------------------------
// Class to keep information from comparison of 
// reconstructed and MC particle tracks (TPC efficiency).   
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
class TGeoManager; 
class TString;
class TStatToolkit; 
class AliMagF;
class AliESDVertex;

#include "TNamed.h"
#include "AliComparisonObject.h"

class AliComparisonEff : public AliComparisonObject {
public :
  AliComparisonEff(); 
  ~AliComparisonEff();

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
  TFolder *CreateFolder(TString folder = "folderEff",TString title = "Analysed Efficiency histograms");

  // Export objects to folder
  TFolder *ExportToFolder(TObjArray * array=0);

  // Process events
  void Process(AliMCInfo* infoMC, AliESDRecInfo *infoRC);

  // Selection cuts
  void SetAliRecInfoCuts(AliRecInfoCuts* cuts=0) {fCutsRC = cuts;}
  void SetAliMCInfoCuts(AliMCInfoCuts* cuts=0) {fCutsMC = cuts;} 
  
  // Getters
  AliRecInfoCuts*  GetAliRecInfoCuts() const {return fCutsRC;} 
  AliMCInfoCuts*   GetAliMCInfoCuts()  const {return fCutsMC;}

   

private:

  // Control histograms
  TH1F *fMCPt;
  TH1F *fMCRecPt;
  TH1F *fMCRecPrimPt;
  TH1F *fMCRecSecPt;
  
  TProfile* fEffTPCPt;      //->TPC efficiency as function of Pt (tan+-1)
  TProfile* fEffTPCPtMC;    //->MC -TPC efficiency as function of Pt (tan+-1)
  TProfile* fEffTPCPtF;     //->efficiency for findable tracks

  TProfile* fEffTPCPt_P;    //->TPC efficiency as function of Pt (tan+-1) - Protons
  TProfile* fEffTPCPtMC_P;  //->MC -TPC efficiency as function of Pt (tan+-1) - Protons
  TProfile* fEffTPCPtF_P;   //->efficiency for findable tracks - Protons

  TProfile* fEffTPCPt_Pi;   //->TPC efficiency as function of Pt (tan+-1) - Pions
  TProfile* fEffTPCPtMC_Pi; //->MC -TPC efficiency as function of Pt (tan+-1) - Pions
  TProfile* fEffTPCPtF_Pi;  //->efficiency for findable tracks - Pions

  TProfile* fEffTPCPt_K;    //->TPC efficiency as function of Pt (tan+-1) - Kaons
  TProfile* fEffTPCPtMC_K;  //->MC -TPC efficiency as function of Pt (tan+-1) - Kaons
  TProfile* fEffTPCPtF_K;   //->efficiency for findable tracks - Kaons

  //
  TProfile* fEffTPCTan;     //->TPC efficiency as function of Tan (pt>0.15
  TProfile* fEffTPCTanMC;   //->MC -TPC efficiency as function of Tan (pt>0.15)
  TProfile* fEffTPCTanF;    //->efficiency for findable tracks Tan (pt>0.15)
  //
  TProfile2D* fEffTPCPtTan;    //->TPC efficiency as function of Pt and tan
  TProfile2D* fEffTPCPtTanMC;  //->MC -TPC efficiency as function of Pt and tan
  TProfile2D* fEffTPCPtTanF;   //->TPC efficiency as function of Pt and tan

  // idx - 0 (isPrim), idx - 1 (isPrim && infoRC->GetStatus(1)==3)
  // idx - 2 (infoRC->GetStatus(1)==3),  idx - 3 (infoRC->GetStatus(1)==3 && !isPrim )
  //
  
  TH2F* fTPCPtDCASigmaIdeal[4]; //->TPC efficiency vs Pt vs DCA/Sigma (tan+-1)
  TH2F* fTPCPtDCASigmaFull[4];  //->TPC efficiency vs Pt vs DCA/Sigma (tan+-1, full systematics)
  TH2F* fTPCPtDCASigmaDay0[4];  //->TPC efficiency vs Pt vs DCA/Sigma (tan+-1, goofie systematics)

  TH2F* fTPCPtDCAXY[4];     //->TPC efficiency as Pt vs DCA_XY (tan+-1)
  TH2F* fTPCPtDCAZ[4];      //->TPC efficiency as Pt vs DCA_Z (tan+-1)

  // Pid = 0 - electrons,  1 - muons, 2 - kaons, 3 - pions, 4 - protons   
  TH3F* fTPCPtDCASigmaIdealPid[4]; //->TPC efficiency vs Pt vs DCA/Sigma (tan+-1)
  TH3F* fTPCPtDCASigmaFullPid[4];  //->TPC efficiency vs Pt vs DCA/Sigma (tan+-1, full systematics)
  TH3F* fTPCPtDCASigmaDay0Pid[4];  //->TPC efficiency vs Pt vs DCA/Sigma (tan+-1, goofie systematics)
  TH3F* fTPCPtDCAXYPid[4];     //->TPC efficiency vs Pt vs DCA_XY (tan+-1)
  TH3F* fTPCPtDCAZPid[4];      //->TPC efficiency vs Pt vs DCA_Z (tan+-1)

  // Global cuts objects
  AliRecInfoCuts* fCutsRC;     // selection cuts for reconstructed tracks
  AliMCInfoCuts*  fCutsMC;     // selection cuts for MC tracks

  // Magnet (needed for DCA calculations) 
  AliESDVertex* fVertex;  //! 
  
  // analysis folder 
  TFolder *fAnalysisFolder; // folder for analysed histograms

  AliComparisonEff(const AliComparisonEff&); // not implemented
  AliComparisonEff& operator=(const AliComparisonEff&); // not implemented

  ClassDef(AliComparisonEff,1);
};

#endif
