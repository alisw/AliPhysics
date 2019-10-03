#ifndef AliAnalysisTaskMuonAODfromGeneral_H
#define AliAnalysisTaskMuonAODfromGeneral_H

/* $Id$ */ 

/* 19 Nov 2007
   Class declaration for the specific muon AOD generation
   Extracts only muon tracks from a general AOD and builds dimuons
   Livio Bianchi, Universita' di Torino
*/
#include "TTree.h"
#include "TH1.h"
#include "TChain.h"
#include "AliAODEvent.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "AliAODEventInfo.h"

class AliAnalysisTaskMuonAODfromGeneral : public AliAnalysisTask {
 public:
  AliAnalysisTaskMuonAODfromGeneral() : AliAnalysisTask(), fInfos(0), fDimuons(0), fChain(0), fOrgAOD(0), fNewAOD(0), ft(0), fBeamEnergy(0) {}
  AliAnalysisTaskMuonAODfromGeneral(const char *name, Double_t BeamEnergy);
  virtual ~AliAnalysisTaskMuonAODfromGeneral() {}
  
  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();				
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);			
  
  void SetBeamEnergy(Double_t BeamEnergy){fBeamEnergy=BeamEnergy;}
  Double_t GetBeamEnergy(){return fBeamEnergy;}
  
 private:
  AliAnalysisTaskMuonAODfromGeneral(const AliAnalysisTaskMuonAODfromGeneral&); // Not implemented
  AliAnalysisTaskMuonAODfromGeneral& operator=(const AliAnalysisTaskMuonAODfromGeneral&); // Not implemented

  AliAODEventInfo	*fInfos;
/*  TClonesArray 	*fInfos;*/
  TClonesArray 	*fDimuons;
  TChain 	*fChain;
  AliAODEvent 	*fOrgAOD;
  AliAODEvent   *fNewAOD;
  TTree 	*ft; // Output Tree
  Double_t 	fBeamEnergy;
  ClassDef(AliAnalysisTaskMuonAODfromGeneral, 1); // example of analysis
};
#endif

