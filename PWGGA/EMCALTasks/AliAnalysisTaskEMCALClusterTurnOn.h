#ifndef ALIANALYSISTASKEMCALCLUSTERTURNON_H
#define ALIANALYSISTASKEMCALCLUSTERTURNON_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

  ///////////////////////////////////////////////////////////////////////////
  ///\class AliAnalysisTaskEMCALClusterTurnOn
  ///
  ///
  ///
  /// 
  ///
  /// 
  /// 
  /// 
  ///////////////////////////////////////////////////////////////////////////

  //ROOT System

class TH1D;
class TH2D;
class TH3D;
class TF1;
class THnSparse;
class TList;
class TObjArray;
class AliEMCALGeometry;
class AliESDCaloCells;
class AliESDEvent;
class AliESDtrack;
class TClonesArray;
class TList;
class TString;
class AliVCluster;
class AliVParticle;
class AliESDtrackCuts;
class AliAODEvent;
class AliAODCaloCells;
class AliVCluster;
class AliMCEvent;
class AliStack;
class TParticle;
class AliClusterContainer;
class AliParticleContainer;
class AliTrackContainer;
class AliEmcalParticle;
  //AliRoot Syste
class AliEMCALTrack;
  //class AliMagF;
class AliEMCALRecoUtils;
  //class AliAnalysisFilter;
class AliAODTrack;
class AliAODCaloCluster;
class AliESDCaloCluster;
class AliVCaloCells;
class AliAODMCParticle;
class AliGenPythiaEventHeader;
  //class AliEventPoolManager;

#include "AliAnalysisTaskEmcal.h"
#include <vector>
using std::vector;

class AliAnalysisTaskEMCALClusterTurnOn: public AliAnalysisTaskEmcal {
public:
  AliAnalysisTaskEMCALClusterTurnOn();
  AliAnalysisTaskEMCALClusterTurnOn(const char *name, Bool_t histo=kFALSE);
  virtual ~AliAnalysisTaskEMCALClusterTurnOn();
  
  void                     UserCreateOutputObjects();
  
  void                     SetIsoConeRadius(Float_t r)                                     { fIsoConeRadius = r ;}
  void                     SetCTMdeltaEta (Float_t r)                                      { fdetacut = r ;}
  void                     SetCTMdeltaPhi (Float_t r)                                      { fdphicut = r ;}
  void                     SetQA (Bool_t QA)                                               { fQA = QA;}
  void                     SetThn (Bool_t Th)                                              { fThn = Th;}
  void                     SetOnlyRecalc (Bool_t rec)                                              { fRecalc = rec;}
  void                     SetNLMCut (Bool_t isNLMCut, Int_t NLMCut, Int_t NLMmin)                       { fIsNLMCut = isNLMCut;
    fNLMCut = NLMCut; fNLMmin = NLMmin;}
  void                     SetTMClusterRejection (Bool_t tm)                               { fTMClusterRejected = tm;}
  void                     SetTMClusterRejectioninCone (Bool_t tm)                         { fTMClusterInConeRejected = tm;}
  void			   SetPtBinning(vector<Double_t> binedges)			   { fBinsPt = binedges; }
  void			   SetPtClBinning(vector<Double_t> binedges)			   { fBinsPtCl = binedges; }
  void			   SetRejectionBinning(vector<Double_t> binedges)			   { fBinsRejection = binedges; }
  void			   SetEnergyBinning(vector<Double_t> binedges)			   { fBinsEnergy = binedges; }
  void			   SetEtaBinning(vector<Double_t> binedges)			   { fBinsEta = binedges; }
  void			   SetPhiBinning(vector<Double_t> binedges)			   { fBinsPhi = binedges; }
  void			   SetEtaClBinning(vector<Double_t> binedges)			   { fBinsClEta = binedges; }
  void			   SetPhiClBinning(vector<Double_t> binedges)			   { fBinsClPhi = binedges; }
  void                     SetM02cut (Bool_t M02)                                          { fM02cut = M02;}
  void                     SetFastOrMasking(Bool_t mask)                                   { fMaskFastOrCells = mask; }
  void                     SetFastOrPath(TString path)                                     { fFastOrPath = path;}

protected:
  
  
  Bool_t                   ClustTrackMatching(AliVCluster *emccluster);

  Int_t                    GetNLM(AliVCluster *coi, AliVCaloCells* cells);
  Int_t                    GetNLM(AliVCluster* coi, AliVCaloCells* cells, Int_t *absIdList, Float_t *maxEList);
  Bool_t                   AreNeighbours(Int_t abscell1, Int_t abscell2) const;
  Float_t                  RecalEnClust(AliVCluster* cluster, AliVCaloCells * cells);
  void                     RecalAmpCell(Float_t  & amp, Int_t absId) const ;
  
  Bool_t                   CheckBoundaries(TLorentzVector vecCOI);
  
  Double_t*                 GenerateFixedBinArray(Int_t n, Double_t min, Double_t max) const;
  void                     ExecOnce();
  Bool_t                   Run();
  void                     FillTHnSparse(AliEMCALGeometry* geom, AliVCluster *coi, TLorentzVector vecCOI, Double_t RejectedAt);
  
  AliAODEvent *fAODEvent;  //!<!
  AliVEvent *fVevent;      //!<! AliVEvent
  
  AliEMCALRecoUtils          *fEMCALRecoUtils; //!<!  EMCAL utils for cluster rereconstruction.
  
  
    //    TList       *fOutputList; //!<! Output list
    //    TGeoHMatrix *fGeomMatrix[12];//!<! Geometry misalignment matrices for EMCal
  Float_t     fIsoConeRadius;                  // Radius for the Isolation Cont
  Double_t    fdetacut;                        // cut on deta between track and cluster
  Double_t    fdphicut;                        // cut on dphi between track and cluster
  Double_t    fM02mincut;                      // lambda0^2 minimum cut
  Double_t    fM02maxcut;                      // lambda0^2 maximum cut
  Bool_t      fQA;                             // Flag for few further QA plots wrt the ones already done in the EMCALTask
  Bool_t      fThn;                            // Flag for ThnSparse
  Bool_t      fRecalc;                         // only analyze events with L1 recalc patch above 2012 threshold
  Int_t       fNDimensions;                    //!<!number of Dimensions for the THnSPARSE Reconstruction
  Int_t       fClusDimensions;                 //!<!number of Dimensions for the THnSparse with cluster info
  Bool_t      fIsNLMCut;                       // NLM cut available
  Int_t       fNLMCut;                         // number of NLM cut
  Int_t       fNLMmin;                         // minimum number of NLM
  Bool_t      fTMClusterRejected;              // able/disable TM cluster rejection
  Bool_t      fTMClusterInConeRejected;        // able/disable TM cluster rejection in isolation cone
  TString     fFastOrPath;
  
    // Initialization of variables for THnSparse

  std::vector<Double_t>    fBinsPt;
  std::vector<Double_t>    fBinsPtCl;
  std::vector<Double_t>    fBinsRejection;
  std::vector<Double_t>    fBinsEnergy;
  std::vector<Double_t>    fBinsEta;
  std::vector<Double_t>    fBinsPhi;
  std::vector<Double_t>    fBinsClEta;
  std::vector<Double_t>    fBinsClPhi;

  

    //IMPLEMENT ALL THE HISTOGRAMS AND ALL THE OUTPUT OBJECTS WE WANT!!!
  TH2D        *fEtaPhiClus;                     ///< EMCAL Cluster Distribution EtaPhi ---QA
  TH1D        *fPT;                             //!<! Pt distribution
  TH1D        *fE;                              //!<! E distribution
  TH2D        *hEt_M02;                         //!<! E vs M02 distribution
  TH2D        *fNLM;                            //!<! NLM distribution
  TH1D        *fVz;                             //!<! Veretex Z distribution
  TH1D        *fEvents;                         //!<! Number of Events
  TH1D        *fPtaftTime;                      //!<! E distribution for clusters after cluster time cut
  TH1D        *fCelldist;                       //!<! Cell distribution after time cut
  TH1D        *fPtaftCell;                      //!<! Pt distribution for clusters after NCells cut
  TH1D        *fNLMdist;                        //!<! NLM distribution after nCells cut
  TH1D        *fPtaftNLM;                       //!<! Pt distribution for clusters after NLM cut
  TH1D        *fPtaftTM;                        //!<! E distribution for neutral clusters
  TH1D        *fDTBC;                           //!<! DTBC distribution after TM
  TH1D        *fPtaftDTBC;                      //!<! E distribution for NC after DistanceToBadChannel cut
  TH1D        *fPtaftFC;                        //!<! E distribution for clusters after fiducial cut
  TH2D        *fTriggerbit;                     //!<! fired Trigger Bits 
//  TH1D        *hL0Amplitude;                    //!<! ADC Amplitudes of EMC L0
  TH2D        *hmaxADC;                         //!<! max L1 ADC/Event
  TH2D        *hADCpos0;                        //!<! ADC of patches at (0,0)
  TH1D        *hSmearedE;                       //!<! smeared patch E distribution
  TH1D        *hPatchE;                         //!<! patch energy
//  TH2D        *hmaxL0ADC;                       //!<! max L0 ADC/Event
  TH2D        *hL1PatchPosition;                //!<! position of L1 trigger patch
  TH2D        *hFastOrPatchE;                   //!<! FastOr# E distribution
//  TH1D        *fL0triggered;                    //!<! max cluster energy of L0 triggered events
  TH2D        *fEventsover10;                   //!<! # of Events over threshold
  TH1D        *fL1triggered;                    //!<! max cluster energy of L1 triggered events
  TH1D        *fClusTime;                       //!<! Time distribution for clusters
  TH2D        *fPt_trig;                        //!<! Et - fired trigger class
  Bool_t      fM02cut;                          
  Bool_t      fMaskFastOrCells;                 
  
  THnSparse   *fOutputTHnS;                    //!<! pT,Rejection,cell info
  TH2D        *hFastOrIndexLeadingCluster;     //!<! leading Cluster E per FastOr
  THnSparse   *fOutTHnS_Clust;                 //!<! pT,Rejection,cluster info
  
  std::vector<Int_t>       MaskedFastOrs;      
  
private:
  AliAnalysisTaskEMCALClusterTurnOn(const AliAnalysisTaskEMCALClusterTurnOn&);            // not implemented
  AliAnalysisTaskEMCALClusterTurnOn&operator=(const AliAnalysisTaskEMCALClusterTurnOn&); // not implemented
  
    /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskEMCALClusterTurnOn, 10);    //EMCAL Neutrals base analysis task
                                                       /// \endcond
};
#endif
