#ifndef ALIANALYSISTASKFULLPPJET_H
#define ALIANALYSISTASKFULLPPJET_H

// **************************************
// Task used for the analysis of full pp jets
// In additional, functions needed for systematic
// uncertainties are also included
// *******************************************

#include <vector>

class TFormula;
class TH2F;
class TH1F;
class TF1;
class THnSparse;
class TRandom3;
class TObjArray;
class TClonesArray;
class TObject;
class TString;
class TProfile2D;
class AliAODEvent;
class AliESDEvent;
class AliMCEvent;
class AliStack;
class AliESDtrack;
class AliESDtrackCuts;
class AliEMCALGeometry;
class AliEMCALRecoUtils;
class AliESDCaloCluster;
class AliFJWrapper;
class AliAODJet;

#include "AliAnalysisTaskSE.h"
#include "fastjet/JetDefinition.hh"
#include "fastjet/PseudoJet.hh"

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskFullppJet : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskFullppJet();
  AliAnalysisTaskFullppJet(const char *name);
  virtual ~AliAnalysisTaskFullppJet();

  Bool_t Notify();
  void UserCreateOutputObjects();
  void UserExec(Option_t *option);
  void Terminate(Option_t *);
  void PrintConfig();
  
  enum {kTPCOnlyVtx = 1<<3,
        kTrigger = 1<<4,
	kHighZ = 1<<5,
	kSuspicious =   1<<6};

  void SetMCAna(Bool_t mc)                        { fIsMC=mc;                       }
  void SetMCStandalone(Bool_t mc)                 { fMCStandalone=mc;               }
  void SetChargedMC(Bool_t mc)                    { fChargedMC=mc;                  }
  void SetRunPeriod(char *p)                      { fPeriod=p;                      }
  void SetOfflineTrigger(Bool_t t)                { fOfflineTrigger=t;              }
  void SetXsec(Float_t Xsec)                      { fXsecScale=Xsec;                }
  void SetConstrainChargedInEMCal(Bool_t con)     { fConstrainChInEMCal=con;        }
  void SetRejectNK(Bool_t reject)                 { fRejectNK=reject;               }
  void SetRejectWD(Bool_t reject)                 { fRejectWD=reject;               }
  void SetVerbosity(Int_t v)                      { fVerbosity = v;                 }
  void SetTrackCutsType(Int_t type)               { fTrackCutsType = type;          }
  void SetEsdTrackCuts(AliESDtrackCuts *cuts)     { fEsdTrackCuts=cuts;             }
  void SetHybridTrackCuts1(AliESDtrackCuts *cuts) { fHybridTrackCuts1=cuts;         }
  void SetHybridTrackCuts2(AliESDtrackCuts *cuts) { fHybridTrackCuts2=cuts;         }
  void SetKinCutType(Int_t type)                  { fKinCutType = type;             }
  void SetZvtx(Double_t zvtx)                     { fZVtxMax = zvtx;                }
  void SetEtaMax(Double_t eta)                    { fTrkEtaMax = eta;               }
  void SetdEdxRange(Double_t min, Double_t max)   { fdEdxMin=min; fdEdxMax=max;     }
  void SetEoverPRange(Double_t min, Double_t max) { fEoverPMin=min; fEoverPMax=max; }
  void SetMatchType(Int_t type)                   { fMatchType=type;                }
  void SetRejectExoticCluster(Bool_t reject)      { fRejectExoticCluster=reject;    }
  void SetRemoveProblematicSM4(Bool_t remove)     { fRemoveBadChannel=remove;       }
  void SetUseGoodSM(Bool_t good)                  { fUseGoodSM=good;                }
  void SetStudySubEInHC(Bool_t study)             { fStudySubEInHC=study;           }
  void SetStudyMcOverSubE(Bool_t study)           { fStudyMcOverSubE=study;         }
  void SetRejectElectron(Bool_t reject)           { fElectronRejection=reject;      }
  void SetCorrectHadron(Bool_t correct)           { fHadronicCorrection=correct;    }
  void SetHCFraction(Float_t fraction)            { fFractionHC=fraction;           } 
  void SetHCLowerPtCutMIP(Double_t pt)            { fHCLowerPtCutMIP=pt;            }
  void SetNonStdBranch(char* name)                { fNonStdBranch=name;             }
  void SetNonStdFile(char* name)                  { fNonStdBranch=name;             }
  void SetAlgorithm(char *algo)                   { fAlgorithm=algo;                }
  void SetRadius(char *r)                         { fRadius=r;                      }
  void SetRecombinationScheme(Int_t scheme)       { fRecombinationScheme=scheme;    }
  void SetSpotGoodJet(Bool_t spot)                { fSpotGoodJet=spot;              }
  void SetFindChargedOnlyJet(Bool_t ch)           { fFindChargedOnlyJet=ch;         }
  void SetFindNeutralOnlyJet(Bool_t ne)           { fFindNeutralOnlyJet=ne;         }
  void SetCheckTrkEffCorr(Bool_t check)           { fCheckTrkEffCorr=check;         }
  void SetTrkEffCorrCutZ(Double_t zcut)           { fTrkEffCorrCutZ=zcut;           }    
  void SetSmearMC(Double_t smear)                 { fSmearMC=smear;                 }
  void SetRunUE(Bool_t run)                       { fRunUE=run;                     }

  //--------------------------------
  // Kinematic cut
  //--------------------------------
  void SetPtRange(Double_t minMB, Double_t maxMB, Double_t minHT, Double_t maxHT) 
  { fTrkPtMin[0]=minMB; fTrkPtMax[0]=maxMB; fTrkPtMin[1]=minHT; fTrkPtMax[1]=maxHT; }

  void SetEtRange(Double_t minMB, Double_t maxMB, Double_t minHT, Double_t maxHT)
  { fClsEtMin[0]=minMB; fClsEtMax[0]=maxMB; fClsEtMin[1]=minHT; fClsEtMax[1]=maxHT;  }

  //--------------------------------
  // Jet quality cut
  //--------------------------------
  void SetJetNEFCut(Double_t min, Double_t max)
  { fJetNEFMin=min; fJetNEFMax=max; }

  //---------------------------------
  // Systematic study
  //---------------------------------
  void SetSaveQAHistos(Bool_t save)                 { fSaveQAHistos=save;           }
  void SetSysJetTrigEff(Bool_t sys)                 { fSysJetTrigEff=sys;           }
  void SetVaryJetTrigEff(Double_t vary)             { fVaryJetTrigEff=vary;         }
  void SetSysTrkPtRes(Bool_t sys)                   { fSysTrkPtRes=sys;             }
  void SetVaryTrkPtRes(Double_t vary)               { fVaryTrkPtRes=vary;           }
  void SetSysTrkEff(Bool_t sys)                     { fSysTrkEff=sys;               }
  void SetVaryTrkEff(Double_t vary)                 { fVaryTrkEff=vary;             }
  void SetSysTrkClsMth(Bool_t sys)                  { fSysTrkClsMth=sys;            }
  void SetSysTrkClsCut(Double_t deta, Double_t dphi){ fCutdEta=deta; fCutdPhi=dphi; }
  void SetSysNonLinearity(Bool_t sys)               { fSysNonLinearity=sys;         }
  void SetSysClusterEScale(Bool_t sys)              { fSysClusterEScale=sys;        }
  void SetVaryClusterEScale(Double_t vary)          { fVaryClusterEScale=vary;      }
  void SetSysClusterERes(Bool_t sys)                { fSysClusterERes=sys;          }
  void SetVaryClusterERes(Double_t vary)            { fVaryClusterERes=vary;        }
    

protected:
  AliESDtrack  *GetAcceptTrack(AliESDtrack *esdtrack);
  Int_t        RunOfflineTrigger();
  Double_t     GetOfflineTriggerProbability(AliESDCaloCluster *cluster);
  Int_t        GetClusterSuperModule(AliESDCaloCluster *cluster);
  void         ProcessMC(const Int_t r=0);
  void         GetMCInfo();
  Bool_t       IsGoodMcPartilce(const AliVParticle* vParticle, const Int_t ipart);
  Int_t        FindSpatialMatchedJet(fastjet::PseudoJet jet, AliFJWrapper *jetFinder, const Double_t radius);
  Int_t        FindEnergyMatchedJet(AliFJWrapper *jetFinder1, const Int_t index1, AliFJWrapper *jetFinder2, const Double_t fraction=0.5);
  Bool_t       IsPrimaryVertexOk(const Double_t trueVtxZ) const;
  Bool_t       IsTPCOnlyVtx() const;
  Bool_t       IsLEDEvent() const;
  void         CheckExoticEvent();
  void         CheckEventTriggerBit();
  void         BookHistos();
  void         GetESDTrax();
  Bool_t       IsElectron(AliESDtrack *track, Double_t clsE) const;
  void         GetESDEMCalClusters();
  Bool_t       IsGoodCluster(AliESDCaloCluster *cluster);
  Double_t     SubtractClusterEnergy(AliESDCaloCluster *cluster, Double_t &eRes, Double_t &MIPE, Double_t &MCsubE);
  Double_t     GetClusterEnergyResolution(AliESDCaloCluster *cluster);
  void         FindDetJets(const Int_t s=0, const Int_t a=0, const Int_t r=0);
  void         FillAODJets(TClonesArray *fJetArray, AliFJWrapper *jetFinder, const Bool_t isTruth = 0);
  void         AnalyzeJets(AliFJWrapper *jetFinder, const Int_t type, const Int_t r);
  void         RunAnalyzeUE(AliFJWrapper *jetFinder);
  Bool_t       IsGoodJet(fastjet::PseudoJet jet, Double_t radius);
  Bool_t       IsGoodJet(AliAODJet *jet, Double_t radius);
  Double_t     GetLeadingZ(const Int_t jetIndex, AliFJWrapper *jetFinder);
  Double_t     GetZ(const Double_t trkPx, const Double_t trkPy, const Double_t trkPz, const Double_t jetPx, const Double_t jetPy, const Double_t jetPz) const;
  Double_t     GetMeasuredJetPtResolution(const Int_t jetIndex, AliFJWrapper *jetFinder);
  Double_t     GetTrkEff(Double_t inPt);
  Double_t     GetJetMissPtDueToTrackingEfficiency(const Int_t jetIndex, AliFJWrapper *jetFinder, const Int_t radiusIndex);
  Double_t     GetExoticEnergyFraction(AliESDCaloCluster *cluster);
  Double_t     GetSmearedTrackPt(AliESDtrack *track);
  Double_t     GetAdditionalSmearing(AliESDtrack *track);

  enum { kNBins = 20 };

 private:
  Int_t             fVerbosity;                         //  Control output
  Int_t             fEDSFileCounter;                    //  Keep track of the ESD file inside a chain
  Int_t             fNTracksPerChunk;                   //  Number of tracks per ESD file; used for debugging
  TString           fPeriod;                            //  Data period
  AliESDEvent       *fESD;                              //! ESD object
  AliAODEvent       *fAOD;                              //! AOD object 
  AliMCEvent        *fMC;                               //! MC object
  AliStack          *fStack;                            //! MC stack
  TObjArray         *fTrackArray;                       //! Array of input tracks
  TObjArray         *fClusterArray;                     //! Array of input clusters
  TArrayI           *fMcPartArray;                      //! Array of MC particles
  Bool_t            fIsMC;                              //  Flag if analyzing MC data
  Bool_t            fMCStandalone;                      //  Flag if only analyze Particle-Level MC
  Bool_t            fChargedMC;                         //  Flag if finding only charged MC jets
  Float_t           fXsecScale;                         //  Corss section
  Double_t          fCentrality;                        //! V0M for current event
  Double_t          fZVtxMax;                           //  Max vertex z cut
  Int_t             fTriggerType;                       //  0-MB, 1-EMC
  Bool_t            fIsTPCOnlyVtx;                      //  Flag of events with ONLY TPC vertex
  Bool_t            fIsExoticEvent3GeV;                 //  Flag of events with exotic cluster above 3 GeV
  Bool_t            fIsExoticEvent5GeV;                 //  Flag of events with exotic cluster above 5 GeV
  Bool_t            fIsEventTriggerBit;                 //  Flag of triggered events with valid trigger bit
  Bool_t            fOfflineTrigger;                    //  Run offline trigger
  TH2F              *fTriggerMask;                      //! Offline trigger mask
  TH1D              *fTriggerCurve[10];                 //! Trigger turn-on curves of EMCal clusters
  TF1               *fTriggerEfficiency[10];            //! Fit of trigger turn-on curves for EMCal clusters above 4-5 GeV
  AliEMCALGeometry  *fGeom;                             //! EMCal goemetry utility
  AliEMCALRecoUtils *fRecoUtil;                         //! Reco utility
  AliESDtrackCuts   *fEsdTrackCuts;                     //! Track cuts for good tracks
  AliESDtrackCuts   *fHybridTrackCuts1;                 //! Track cuts for tracks without SPD hit
  AliESDtrackCuts   *fHybridTrackCuts2;                 //! Track cuts for tracks witout SPD hit or ITS refit
  Int_t             fTrackCutsType;                     //  0-Global track, 1-TPCOnly track
  Int_t             fKinCutType;                        //  0-cut on track before jet finding, 1-cut on jet with high-pt tracks
  Double_t          fTrkEtaMax;                         //  Max |eta| cut
  Double_t          fTrkPtMin[2];                       //  Min pt cut
  Double_t          fTrkPtMax[2];                       //  Max pt cut
  Double_t          fdEdxMin;                           //  Min dE/dx cut
  Double_t          fdEdxMax;                           //  Max dE/dx cut
  Double_t          fEoverPMin;                         //  Min E/P cut
  Double_t          fEoverPMax;                         //  Max E/P cut
  Double_t          fClsEtMin[2];                       //  Min et cut
  Double_t          fClsEtMax[2];                       //  Max et cut
  Int_t             fMatchType;                         //  0-extrapolation, 1-extrapolation + MC label
  Bool_t            fRejectExoticCluster;               //  Flag to reject exotic cluster
  Bool_t            fRemoveBadChannel;                  //  Flag to remove problematic region in SM4
  Bool_t            fUseGoodSM;                         //  Flag to not use trigger bit in SM2,3,4,5
  Bool_t            fStudySubEInHC;                     //  If true, the hadronic correction will be ingored. For physics, it should be set to kFALSE
  Bool_t            fStudyMcOverSubE;                   //  Study the over-subtraction of hadronic correction in simualtion. 
  Bool_t            fElectronRejection;                 //  Switches on electron correction to avoid double counting of electrons 
  Bool_t            fHadronicCorrection;                //  switches on hadronic correction to avoid double counting of hadrons
  Float_t           fFractionHC;                        //  fraction of hadronic correction
  Double_t          fHCLowerPtCutMIP;                   //  Lower track pt cut for MIP correction    
  TF1               *fClusterEResolution;               //! Parameterization of cluster energy resolution from test beam results
  TF1               *fMCNonLin;                         //! Non-linearity of simualtion
  Double_t          fJetNEFMin;                         //  Min jet NEF cut
  Double_t          fJetNEFMax;                         //  Max jet NEF cut
  Bool_t            fSpotGoodJet;                       //  Good jet catching
  Bool_t            fFindChargedOnlyJet;                //  Find jets with TPC tracks
  Bool_t            fFindNeutralOnlyJet;                //  Find jets with EMCal clusters
  Bool_t            fCheckTrkEffCorr;                   //  Check the procedure of tracking efficiency correction
  Double_t          fTrkEffCorrCutZ;                    //  Cut on the tracks that are added back Z < 0.3
  TF1               *fTrkEffFunc[3];                    //! Fit function of tracking efficiency
  TH1F              *fhCorrTrkEffPtBin[2][2];           //! Number of tracks per jet pt bin, used to correct the tracking efficiency explicitly
  TH1F              *fhCorrTrkEffSample[2][2][kNBins];  //! Tracking efficiency estimated from simulation
  TRandom3          *fRandomGen;                        //! Random number generator
  Bool_t            fRunUE;                             //  Run analysis of underlying event
  
  Bool_t            fSysJetTrigEff;                     //  Flag of systematic uncertainty of jet trigger efficiency
  Double_t          fVaryJetTrigEff;                    //  Variation of cluster E-scale for systematic uncertainty of jet trigger efficiency
  Bool_t            fSysTrkPtRes;                       //  Flag of systematic uncertainty of tracking momentum resolution
  Double_t          fVaryTrkPtRes;                      //  Variation of tracking momentum resolution
  Bool_t            fSysTrkEff;                         //  Flag of systematic uncertainty of tracking efficiency
  Double_t          fVaryTrkEff;                        //  Variation of tracking efficiency
  Bool_t            fSysTrkClsMth;                      //  Flag of systematic uncertainty of track-cluster matching
  Double_t          fCutdEta;                           //  Variation of dEta cut
  Double_t          fCutdPhi;                           //  Variation of dPhi cut
  Bool_t            fSysNonLinearity;                   //  Flag of systematic uncertainty of EMCal non-linearity
  Bool_t            fSysClusterEScale;                  //  Flag of systematic uncertainty of EMCal energy scale
  Double_t          fVaryClusterEScale;                 //  Variation of EMCal energy scale
  Bool_t            fSysClusterERes;                    //  Flag of systematic uncertainty of EMCal energy resolution
  Double_t          fVaryClusterERes;                   //  Variation of EMCal energy resolution
  
  TString           fNonStdBranch;                      //! Non-std branch name for AOD jets
  TString           fNonStdFile;                        //! Name of optional file that the non-std branch is written to
  TString           fAlgorithm;                         //! name of algorithm
  TString           fRadius;                            //! Jet cone radius
  Int_t             fRecombinationScheme;               //  Recombination scheme of jet finding
  AliFJWrapper      *fDetJetFinder[3][2][2];            //! Jet finder
  TClonesArray      *fJetTCA[3][2][2];	                //! TCA of jets: in - akt - 0.4
  Bool_t            fConstrainChInEMCal;                //  Constain charged particle to be in EMCal acceptance
  Bool_t            fRejectNK;                          //  Reject neutron and K_L
  Bool_t            fRejectWD;                          //  Reject primaries, mainly k^0_S,  that decay weakly
  Bool_t            fSmearMC;                           //  Flag of smearing tracking resolution in MC to match data. Obselete.
  TF1               *fTrkPtResData;                     //! Parameterazation of momentum resolution estimated from data
  AliFJWrapper      *fTrueJetFinder[2];                 //! Jet finder for particle jets
  TClonesArray      *fMcTruthAntikt[2];                 //! TCA of MC truth anti-kt jets

  TList             *fOutputList;                       //! Output list
  Bool_t            fSaveQAHistos;                      //  Flag of saving QA histograms
  TH1F              *fhJetEventStat;                    //! Event counts used for jet analysis
  TH1F              *fhEventCheck;                      //! Event statistics for vertex types
  TH2F              *fhTrkPhiEtaCheck;                  //! Phi vs Eta of tracks in events with only TPC vertex
  TH1F              *fhChunkQA;                         //! Check if the chunk is corrupted
  TH1F              *fVertexGenZ[2];                    //! Generated event vertex z
  TH1F              *fEventZ[2];                        //! reconstructed event vertex z
  TH1F              *fhNTrials[2];                      //! # of trials
  TH1F              *fhNMatchedTrack[2];                //! # of matched tracks per cluster
  TH2F              *fhSubEVsTrkPt[2][4];               //! Subtracted energy due to hadronic correction
  TH2F              *fhNeutralPtInJet[3][2];            //! pt of neutral constituents in jet
  TH2F              *fhChargedPtInJet[3][2];            //! pt of charged constituents in jet
  TH2F              *fhLeadNePtInJet[3][2];             //! pt of leading neutral constituents in jet
  TH2F              *fhLeadChPtInJet[3][2];             //! pt of leading charged constituents in jet
  TH2F              *fhChLeadZVsJetPt[2][2];            //! Leading charged constituent Z vs jet pt
  TH3F              *fhJetPtVsZ[3][2];                  //! Jet pt vs constituent Z vs constituent type
  TH3F              *fRelTrkCon[2][2];                  //! Jet pt vs track pt contribution vs track class
  TH2F              *fhJetPtWithTrkThres[2][2];         //! pt of jets containing tracks above certian threshold
  TH2F              *fhJetPtWithClsThres[2][2];         //! pt of jets containing clusters above certian threshold
  TH2F              *fhJetPtVsLowPtCons[2][2];          //! Contribution of low pt particles to jet energy
  THnSparse         *fJetEnergyFraction[3][2];          //! Jet energy fraction
  THnSparse         *fJetNPartFraction[3][2];           //! Jet NPart fraction
  TH1F              *fJetCount[3][2];                   //! pT distribution of pions detected 
  TH2F              *fhSubClsEVsJetPt[2][2][5];         //! f*subtracted cluster energy vs jet pt
  TH2F              *fhHCTrkPtClean[2][2][5];           //! Cleanly subtracted charged pt
  TH2F              *fhHCTrkPtAmbig[2][2][5];           //! Ambiguously subtracted charged pt
  TH2F              *fHCOverSubE[2][5];                 //! Error made by hadronic correction assessed by using particle jet
  TH2F              *fHCOverSubEFrac[2][5];             //! Error made by hadronic correction assessed by using particle jet
  TH3F              *fhFcrossVsZleading[2][2];          //! Jet pt vs Fcross vs Zleading
  TH1F              *fhJetPtInExoticEvent[2][2];        //! Jet pt in exotic events
  TH2F              *fhJetPtVsUE[2];                    //! Underlying event contribution
  TH1F              *fhClsE[2];                         //! Cluster energy distribution

  AliAnalysisTaskFullppJet(const AliAnalysisTaskFullppJet&);            // not implemented
  AliAnalysisTaskFullppJet &operator=(const AliAnalysisTaskFullppJet&); // not implemented

  ClassDef(AliAnalysisTaskFullppJet, 1);
};

#endif
