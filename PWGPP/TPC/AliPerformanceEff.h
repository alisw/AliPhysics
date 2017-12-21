#ifndef ALIPERFORMANCEEFF_H
#define ALIPERFORMANCEEFF_H

//------------------------------------------------------------------------------
// Class to keep information from comparison of 
// reconstructed and MC particle tracks (TPC efficiency).   
// 
// Author: J.Otwinowski 04/02/2008 
// Changes by J.Salzwedel 10/30/2014
//------------------------------------------------------------------------------

class TFile;
class TParticle;
class TString;
class TNamed;
class THnSparse;
class AliMCInfo;
class AliESDRecInfo;
class AliVEvent; 
class AliVfriendEvent; 
class AliMCEvent; 
class AliMCParticle; 
class AliVtrack;
class TRootIOCtor;

#include "AliPerformanceObject.h"

class AliPerformanceEff : public AliPerformanceObject {
public :
  AliPerformanceEff(TRootIOCtor*);
  AliPerformanceEff(const Char_t* name="AliPerformanceEff",const Char_t*title="AliPerformanceEff",Int_t analysisMode=0, Bool_t hptGenerator=kFALSE);
  virtual ~AliPerformanceEff();

  // Init data members
  virtual void Init();

  // Execute analysis 
  virtual void  Exec(AliMCEvent* const mcEvent, AliVEvent *const vEvent, AliVfriendEvent *const vFriendEvent, const Bool_t bUseMC, const Bool_t bUseVfriend);

  // Merge output objects (needed by PROOF) 
  virtual Long64_t Merge(TCollection* list);

  // Analyse output histograms 
  virtual void Analyse();

  // Get analysis folder
  virtual TFolder* GetAnalysisFolder() const {return fAnalysisFolder;}

  // Create folder for analysed histograms
  TFolder *CreateFolder(TString folder = "folderEff",TString title = "Analysed Efficiency histograms");

  // Export objects to folder
  TFolder *ExportToFolder(TObjArray *array=0);

  // Process events
  void ProcessTPC(AliMCEvent* const mcEvent, AliVEvent *const vEvent);
  void ProcessTPCITS(AliMCEvent* const mcEvent, AliVEvent *const vEvent);
  void ProcessConstrained(AliMCEvent* const mcEvent, AliVEvent *const vEvent);
  void ProcessTPCSec(AliMCEvent* const mcEvent, AliVEvent *const vEvent);

  Bool_t IsRecTPC(AliVTrack *track);
  Bool_t IsRecTPCITS(AliVTrack *track);
  Bool_t IsRecConstrained(AliVTrack *track);

  Bool_t IsFindable(const AliMCEvent *mcEvent, Int_t label);
  Bool_t HasTPCReference(const AliMCEvent *mcEvent, Int_t label);
  Int_t TransformToPID(TParticle *mcPart);

  THnSparseF* GetEffHisto() const {return fEffHisto;}
  THnSparseF* GetEffSecHisto() const {return fEffSecHisto;}
  
  static void SetfReadNClsTree(bool v) {fReadNClsTree = v;}

private:

  static const Int_t fgkMaxClones = 3, fgkMaxFakes = 3;
  
  // Helper Method
  TH1D* AddHistoEff(Int_t axis, const Char_t *name, const Char_t* vsTitle, const Int_t type, const Int_t secondary = 0);
  TH1D* WeightedProjection(THnSparseF* src, Int_t axis, Int_t nWeights, Int_t* weightCoords);

  // Control histograms
  THnSparseF *fEffHisto; //-> mceta:mcphi:mcpt:pid:isPrim:recStatus:findable:charge
  THnSparseF *fEffSecHisto; //-> mceta:mcphi:mcpt:pid:isPrim:recStatus:findable:mcR:mother_phi:mother_eta:charge
  
  TH1D* fTrackPtNCls; //
  TH2D* fTrackNClsFound; //

  // analysis folder 
  TFolder *fAnalysisFolder; // folder for analysed histograms
  
  int fNEvent; //!
  static bool fReadNClsTree; //
  TFile* fNClsTreeFile; //!
  TTree* fNClsTree; //!
  std::vector<short>* fNClsVec; //!

  AliPerformanceEff(const AliPerformanceEff&); // not implemented
  AliPerformanceEff& operator=(const AliPerformanceEff&); // not implemented

  ClassDef(AliPerformanceEff,4);
};

#endif
