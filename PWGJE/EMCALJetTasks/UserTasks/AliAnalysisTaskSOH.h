#ifndef ALIANALYSISTASKSOH_H
#define ALIANALYSISTASKSOH_H

// $Id$

class TList;
class TH1F;
class TH2F;
class TH3F;
class THnSparse;
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
  void                ProcessCluster()                     ;
  void                ProcessMc()                          ;
  void                ProcessScaleFactor()                 ;

  AliESDEvent        *fESD;                      //!esd event
  AliMCEvent         *fMC;                       //!mv event
  AliESDtrackCuts    *fEsdTrackCuts;             // esd track cuts
  AliESDtrackCuts    *fHybridTrackCuts1;         // hybrid track cuts
  AliESDtrackCuts    *fHybridTrackCuts2;         // hybrid track cuts
  TArrayI            *fTrackIndices;             //!selected track index
  TArrayI            *fClusterIndices;           //!cluster with two matched MC track index
  TObjArray          *fClusterArray;             //!selected cluster array
  
  TList              *fOutputList;               //!output list
  TH1F               *fHEventStat;               //!statistics histo
  TH3F               *fHTrkEffParGenPtEtaPhi;    //!charged pion mc truth pt spectrum
  TH3F               *fHTrkEffDetGenPtEtaPhi;    //!charged pion mc detector level pt spectrum
  TH3F               *fHTrkEffDetRecPtEtaPhi;    //!charged pion reconstructed detector level pt spectrum
  TH3F               *fHTrkEffDetRecFakePtEtaPhi;//!fake and secondary tracks pt spectrum
  TH1F               *fHTrkEffParGenPt_rmInj;    //!mc truth without injected signal pt spectrum
  TH1F               *fHTrkEffDetGenPt_rmInj;    //!generated detector level particle without injected signal pt spectrum 
  TH1F               *fHTrkEffParGenPt_PiRmInj;    //!mc truth(pion) without injected signal pt spectrum
  TH1F               *fHTrkEffDetGenPt_PiRmInj;    //!generated detector level particle(pion) without injected signal pt spectrum 
  TH1F               *fHTrkEffParGenPt_all;      //!mc truth pt spectrum
  TH1F               *fHTrkEffDetGenPt_all;      //!generated detector level pt spectrum
  TH1F               *fHTrkEffParGenPt_Dch;      //!charged D meson truth pt spectrum
  TH1F               *fHTrkEffDetGenPt_Dch;      //!charged D meson generated detector level pt spectrum
  TH1F               *fHTrkEffParGenPt_Dn;       //!neutral D meson truth pt spectrum
  TH1F               *fHTrkEffDetGenPt_Dn;       //!neutral D meson generated detector level pt spectrum
  TH1F               *fHTrkEffParGenPt_Ds;       //!strange D meson truth pt spectrum
  TH1F               *fHTrkEffDetGenPt_Ds;       //!strange D meson generated detector level pt spectrum
  TH1F               *fHTrkEffParGenPt_cL;       //!charmed lambda truth pt spectrum
  TH1F               *fHTrkEffDetGenPt_cL;       //!charmed lambda generated detector level pt spectrum
  TH1F               *fHTrkEffParGenPt_JPsi;     //!J/Psi truth pt spectrum
  TH1F               *fHTrkEffDetGenPt_JPsi;     //!J/Psi reconstructed detector level pt spectrum
  TH1F               *fHScaleFactor;             //!scale factor spectrum
  TH1F               *fHScaleFactor100HC;        //!scale factor with 100% HC spectrum
  TH2F               *fHEOverPVsPt;              //!(cluster energy over reconstructed track p) vs. track pt
  TH2F               *fHEMCalResponsePion;       //!same as above for pions 
  TH2F               *fHEMCalResponseElec;       //!same as above for electrons
  TH2F               *fHEMCalResponseProton;     //!same as above for protons
  TH2F               *fHEMCalRecdPhidEta;        //!(EMCal cluster phi - track phi) vs. (EMCal cluster eta - track eta)
  TH2F               *fHEMCalRecdPhidEtaP;       //!same as above for positive charge tracks
  TH2F               *fHEMCalRecdPhidEtaM;       //!same as above for negative charge tracks
  TH2F               *fHEMCalRecdPhidEta_Truth;  //!same as above with mc truth matching
  TH2F               *fHEMCalRecdPhidEtaP_Truth; //!same as above with positive truth charge matching
  TH2F               *fHEMCalRecdPhidEtaM_Truth; //!same as above with negative truth charge matching
  TH2F               *fHEMCalRecdPhidEtaposEta;  //!same as above for positive eta
  TH2F               *fHEMCalRecdPhidEtanegEta;  //!same as above for negative eta
  TH2F               *fHPhotonEdiff100HC;        //!(truth E - calculated E in 100% HC)/truth E vs. truth E with photon
  TH2F               *fHPhotonEdiff70HC;         //!(truth E - calculated E in 70% HC)/truth E vs. truth E with photon
  TH2F               *fHPhotonEdiff30HC;         //!(truth E - calculated E in 30% HC)/truth E vs. truth E with photon
  TH2F               *fHPhotonEdiff0HC;          //!(truth E - cluster E)/truth E vs. truth E with photon
  TH2F               *fHPhotonEVsClsE;           //!cluster E vs. truth photon E
  TH2F               *fHistEsub1Pch;             //!(subtracted E in 100% HC) vs. total track P, clusters with 1 matching track
  TH2F               *fHistEsub2Pch;             //!(subtracted E in 100% HC) vs. total track P, clusters with 2 matching tracks
  TH2F               *fHistEsub1PchRat;          //!(subtracted E in 100% HC)/total track P vs. total track P, clusters with 1 matching track
  TH2F               *fHistEsub2PchRat;          //!(subtracted E in 100% HC)/total track P vs. total track P, clusters with 2 matching tracks
  THnSparse          *fHClsEoverMcE_All;         //!cluster E/MC particle E, cluster with only one matching particle
  THnSparse          *fHClsEoverMcE_Photon;      //!above for photon
  THnSparse          *fHClsEoverMcE_Elec;        //!above for electron
  THnSparse          *fHClsEoverMcE_Pion;        //!above for pion

  AliAnalysisTaskSOH(const AliAnalysisTaskSOH&); // not implemented
  AliAnalysisTaskSOH& operator=(const AliAnalysisTaskSOH&); // not implemented
  
  ClassDef(AliAnalysisTaskSOH, 9); // Analysis task Saehanseul Oh
};
#endif
