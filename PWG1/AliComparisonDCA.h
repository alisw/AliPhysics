#ifndef ALICOMPARISONDCA_H
#define ALICOMPARISONDCA_H

//------------------------------------------------------------------------------
// Class to keep information from comparison of 
// reconstructed and MC particle tracks (DCA - Distance of Closest Approach 
// to the vertex).   
// 
// Author: J.Otwinowski 04/02/2008 
//------------------------------------------------------------------------------

class AliMCInfo;
class AliESDRecInfo;
class AliESDEvent; 
class AliRecInfoCuts;
class AliMCInfoCuts;
class AliESDVertex;
class TH3F;
class TH3;
class TString;
class TNamed;

#include "THnSparse.h"
#include "AliComparisonObject.h"

class AliComparisonDCA : public AliComparisonObject {
public :
  AliComparisonDCA(); 
  AliComparisonDCA(Char_t* name, Char_t* title, Int_t analysisMode, Bool_t hptGenerator);
  ~AliComparisonDCA();

  // Init data members
  virtual void Init();

  // Execute analysis
  virtual void Exec(AliMCInfo* const infoMC, AliESDRecInfo *const infoRC);

  // Merge output objects (needed by PROOF) 
  virtual Long64_t Merge(TCollection* const list);

  // Analyse output histograms
  virtual void Analyse();

  // Get analysis folder
  virtual TFolder* GetAnalysisFolder() const {return fAnalysisFolder;}

  // Create folder for analysed histograms
  TFolder *CreateFolder(TString folder = "folderDCA",TString title = "Analysed DCA histograms");

  // Export objects to folder
  TFolder *ExportToFolder(TObjArray * array=0);

  // Process events
  void  ProcessTPC(AliMCInfo* const infoMC, AliESDRecInfo* const infoRC);
  void  ProcessTPCITS(AliMCInfo* const infoMC, AliESDRecInfo* const infoRC);
  void  ProcessConstrained(AliMCInfo* const infoMC, AliESDRecInfo* const infoRC); // not implemented

  // Selection cuts
  void SetAliRecInfoCuts(AliRecInfoCuts* const cuts=0) {fCutsRC = cuts;}
  void SetAliMCInfoCuts(AliMCInfoCuts* const cuts=0) {fCutsMC = cuts;}  

  AliRecInfoCuts*  GetAliRecInfoCuts() const {return fCutsRC;}
  AliMCInfoCuts*   GetAliMCInfoCuts()  const {return fCutsMC;}

  // getters
  /*
  TH3F  *GetD0TanSPtTPCITS() const {return fD0TanSPtTPCITS;}
  TH3F  *GetD1TanSPtTPCITS() const {return fD1TanSPtTPCITS;}
  TH3F  *GetD0TanSPt() const {return fD0TanSPt;}
  TH3F  *GetD1TanSPt() const {return fD1TanSPt;}
  TH3F  *GetD0TanSPtTPC() const {return fD0TanSPtTPC;}
  TH3F  *GetD1TanSPtTPC() const {return fD1TanSPtTPC;}
  */

  // DCA
  THnSparse* GetDCAHisto() const {return fDCAHisto;}

private:

  // DCA histograms
  THnSparseF *fDCAHisto; //-> dca_r:dca_z:eta:pt 
 
  /*
  TH3F  *fD0TanSPtTPCITS; //-> distance to vertex y (TPC+ITS clusters) 
  TH3F  *fD1TanSPtTPCITS; //-> distance to vertex z (TPC+ITS clusters) 
  TH3F  *fD0TanSPt;     //-> distance to vertex y  
  TH3F  *fD1TanSPt;     //-> distance to vertex z 
  TH3F  *fD0TanSPtTPC;  //-> distance to vertex y (only TPC track parameters) 
  TH3F  *fD1TanSPtTPC;  //-> distance to vertex z (only TPC track parameters)
  */

  // Global cuts objects
  AliRecInfoCuts*  fCutsRC; // selection cuts for reconstructed tracks
  AliMCInfoCuts*  fCutsMC;  // selection cuts for MC tracks

  // analysis folder 
  TFolder *fAnalysisFolder; // folder for analysed histograms

  AliComparisonDCA(const AliComparisonDCA&); // not implemented
  AliComparisonDCA& operator=(const AliComparisonDCA&); // not implemented

  ClassDef(AliComparisonDCA,1);
};

#endif
