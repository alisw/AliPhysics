#ifndef ALIANALYSISTASKSOH_cxx
#define ALIANALYSISTASKSOH_cxx

// $Id$

class TList;
class TH1F;
class TH2F;
class TArrayI;
class AliESDEvent;
class AliMCEvent;
class AliESDtrack;
class AliESDCaloCluster;
class AliESDtrackCuts;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskSOH : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskSOH();
  AliAnalysisTaskSOH(const char *name);
  virtual ~AliAnalysisTaskSOH();
  
  void                UserCreateOutputObjects();
  void                UserExec(Option_t *option);
  void                Terminate(Option_t *);
  
  void                SetEsdTrackCuts(AliESDtrackCuts *cuts)         { fEsdTrackCuts     = cuts ; }
  void                SetHybridTrackCuts1(AliESDtrackCuts *cuts)     { fHybridTrackCuts1 = cuts ; }
  void                SetHybridTrackCuts2(AliESDtrackCuts *cuts)     { fHybridTrackCuts2 = cuts ; }

 private:

  AliESDtrack        *GetAcceptTrack(AliESDtrack *esdtrack) ;
  Bool_t              IsGoodCluster(AliESDCaloCluster *cluster)             ;
  Bool_t              IsGoodMcParticle(AliVParticle* vParticle, Int_t ipart);
  void                ProcessTrack()                       ;
  void                ProcessMc()                          ;

  AliESDEvent        *fESD;                      //!esd event
  AliMCEvent         *fMC;                       //!mv event
  AliESDtrackCuts    *fEsdTrackCuts;             // esd track cuts
  AliESDtrackCuts    *fHybridTrackCuts1;         // hybrid track cuts
  AliESDtrackCuts    *fHybridTrackCuts2;         // hybrid track cuts
  TArrayI            *fTrackIndices;             //!selected track index
  TList              *fOutputList;               //!output list
  TH1F               *fHEventStat;               //!statistics histo
  TH1F               *fHTrkEffParGenPt;          //!mc truth pt spectrum
  TH1F               *fHTrkEffDetGenPt;          //!mc detector level pt spectrum
  TH1F               *fHTrkEffDetRecPt;          //!reconstructed detector level pt spectrum
  TH2F               *fHEOverPVsPt;              //!(cluster energy over reconstructed track p) vs. track pt
  TH2F               *fHEMCalResponsePion;       //!same as above for pions 
  TH2F               *fHEMCalResponseElec;       //!same as above for electrons
  TH2F               *fHEMCalResponseProton;     //!same as above for protons
  TH2F               *fHEMCalRecPhiEtaClus;      //!EMCal cluster phi vs. eta
  TH2F               *fHEMCalRecPhiEtaTrk;       //!EMCal cluster label matched track phi vs. eta
  TH2F               *fHClsPhiEta;               //!EMCal cluster phi vs. eta
  TH2F               *fHEMCalRecdPhidEta;        //!(EMCal cluster phi - track phi) vs. (EMCal cluster eta - track eta)
  TH2F               *fHEMCalRecdPhidEtaP;       //!same as above for positive charge tracks
  TH2F               *fHEMCalRecdPhidEtaM;       //!same as above for negative charge tracks
  TH2F               *fHEMCalRecdPhidEta_Truth;  //!same as above with mc truth matching
  TH2F               *fHEMCalRecdPhidEtaP_Truth; //!same as above with positive truth charge matching
  TH2F               *fHEMCalRecdPhidEtaM_Truth; //!same as above with negative truth charge matching
  TH2F               *fHEMCalRecdPhidEtaposEta;  //!same as above for positive eta
  TH2F               *fHEMCalRecdPhidEtanegEta;  //!same as above for negative eta

  AliAnalysisTaskSOH(const AliAnalysisTaskSOH&); // not implemented
  AliAnalysisTaskSOH& operator=(const AliAnalysisTaskSOH&); // not implemented
  
  ClassDef(AliAnalysisTaskSOH, 1); // Analysis task Saehanseul Oh
};
#endif
