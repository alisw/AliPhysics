#ifndef ALIANALYSISTASKHJETDPHI_H
#define ALIANALYSISTASKHJETDPHI_H

#include <vector>

class TH2F;
class TH1F;
class TF1;
class THnSparse;
class TClonesArray;
class TObject;
class TString;
class AliAODEvent;
class AliESDEvent;
class AliAODExtension;
class AliMCEvent;
class AliRhoParameter;
class TRandom3;
class AliEmcalJet;
class AliVTrack;
class AliNamedArrayI;
class AliAODTrack;
class AliESDtrackCuts;
class AliAODJetEventBackground;
class AliNamedString;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskHJetDphi : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskHJetDphi();
  AliAnalysisTaskHJetDphi(const char *name);
  virtual ~AliAnalysisTaskHJetDphi();

  void UserCreateOutputObjects();
  Bool_t UserNotify();
  void UserExec(Option_t *option);
  void Terminate(Option_t *);
  void PrintConfig();

  void SetVerbosity(Int_t i)                                   { fVerbosity = i;                 }
  void SetIsEmbedding(Bool_t b)                                { fIsEmbedding = b;               }
  void SetAnaType(Int_t i)                                     { fAnaType = i;                   }  
  void SetRunPeriod(char *p)                                   { fPeriod = p;                    }
  void SetCollisionSystem(char *s)                             { fCollisionSystem = s;           }
  void SetIsMC(Bool_t mc)                                      { fIsMC = mc;                     }
  void SetAnalyzeMCTruth(Bool_t mc)                            { fAnalyzeMCTruth = mc;           }
  void SetOfflineTrgMask(AliVEvent::EOfflineTriggerTypes mask) { fOfflineTrgMask = mask;         }
  void SetMaxVtxZ(Double_t z)                                  { fMaxVtxZ = z;                   }
  void SetFilterMask(UInt_t filter)                            { fFilterMask = filter;           }
  void SetRequireITSRefit(Bool_t r)                            { fRequireITSRefit = r;           }
  void SetRequireSharedClsCut(Bool_t r)                        { fRequireSharedClsCut = r;       }
  void SetNonStdFile(char* s)                                  { fNonStdFile = s;                }
  void SetMcParticleArrName(char *s)                           { fMcParticleArrName = s;         }
  void SetEmbTrkArrName(char *s)                               { fEmbTrkArrName = s;             }
  void SetTrackArrName(char *s)                                { fTrackArrName=s;                }
  void SetSwitchOnAvoidTpcHole(Bool_t cut)                     { fSwitchOnAvoidTpcHole=cut;      }
  void SetCutTPCBoundary(Bool_t cut)                           { fCutTPCBoundary=cut;            }
  void SetDistToTPCBoundary(Double_t dist)                     { fDistToTPCBoundary=dist;        }
  void SetTrkPtRange(Double_t min, Double_t max)               { fMinTrkPt=min; fMaxTrkPt=max;   }
  void SetTrkPhiRange(Double_t min, Double_t max)              { fMinTrkPhi=min; fMaxTrkPhi=max; }
  void SetTrkEtaRange(Double_t min, Double_t max)              { fMinTrkEta=min; fMaxTrkEta=max; }
  void SetRadius(Double_t rad)                                 { fRadius=rad;                    }
  void SetJetArrName(char *s)                                  { fJetArrName=s;                  }
  void SetPLJetArrName(char *s)                                { fPLJetArrName = s;              }
  void SetDLJetArrName(char *s)                                { fDLJetArrName = s;              }
  void SetRhoName(char *s)                                     { fRhoName=s;                     }
  void SetRunTrkQA(Bool_t run)                                 { fRunTrkQA=run;                  }
  void SetRunJetQA(Bool_t run)                                 { fRunJetQA=run;                  }
  void SetRunSingleInclHJet(Bool_t run)                        { fRunSingleInclHJet=run;         }
  void SetTTtype(Bool_t type)                                  { fTTtype=type;                   }
  void SetTTRange(Double_t min, Double_t max)                  { fTTMinPt=min; fTTMaxPt=max;     }
  void SetJetPtMin(Double_t min)                               { fJetPtMin = min;                }
  void SetRunPLHJet(Bool_t run)                                { fRunPLHJet = run;               }
  void SetRunDLHJet(Bool_t run)                                { fRunDLHJet = run;               }
  void SetRunLeadTrkQA(Bool_t run)                             { fRunLeadTrkQA=run;              }
  void SetStudyKtEffects(Bool_t study)                         { fStudyKtEffects=study;          }
  void SetKtValue(Double_t kt)                                 { fKtValue=kt;                    }
  void SetAODfilterBits(Int_t b0 = 0, Int_t b1 = 0)            { fAODfilterBits[0] = b0  ; fAODfilterBits[1] = b1  ; }
  void SetRunBkgFlow(Bool_t run)                               { fRunBkgFlow=run;                }


protected:
  Bool_t                            RetrieveArraies();
  void                              RunTrackQA();
  void                              RunJetQA(const TClonesArray *jetArray, const Double_t rho, THnSparse *hJetPt, THnSparse *hJetArea, THnSparse *hJetQA);
  Int_t                             FindSingleIncTrigger(const TClonesArray *trackArray, Double_t &trigPt, Double_t &trigPhi, Double_t &trigEta, const Int_t arrayType);
  void                              RunSingleInclHJetCorr(Double_t trigPt, Double_t trigPhi, Double_t trigEta, const TClonesArray *jetArray, Double_t rho, THnSparse *hTT, THnSparse *hn);
  void                              RunLeadTrkQA();
  void                              StudyKtEffects();
  Bool_t                            AcceptTrack(AliVParticle *track);
  Bool_t                            IsGoodAODtrack(AliVParticle *track);
  Bool_t                            IsGoodJet(Double_t jetEta);
  Double_t                          GetLeadingPt(const Int_t jetIndex);
  Double_t                          GetZ(const Double_t trkPx, const Double_t trkPy, const Double_t trkPz, const Double_t jetPx, const Double_t jetPy, const Double_t jetPz);
  Double_t                          CalculateDPhi(const Double_t phi1, const Double_t phi2);
  Double_t                          CalculatePhi(const Double_t py, const Double_t px);
  Int_t                             LocateToTPCHole(const Double_t phi);
  Int_t                             GetParticleType(Int_t pdg_input);
  Double_t                          GetAODTrackPtRes(AliAODTrack *track);
  Int_t                             GetPtHardBin(Double_t ptHard);

 private:
  Int_t                             fVerbosity;                  //  Control output
  Bool_t                            fIsEmbedding;                //  Flag of embedding trains
  Int_t                             fAnaType;                    //  0-EMCal jet framework; 1-charged jet framework
  TString                           fPeriod;                     //  Run period
  TString                           fCollisionSystem;            //  Collision system
  Bool_t                            fIsMC;                       //  Flag if MC data
  Bool_t                            fAnalyzeMCTruth;             //  Analyze MC truth
  AliMCEvent                        *fMC;                        //! MC events
  AliVEvent                         *fEvent;                     //! Input event 
  AliESDEvent                       *fESD;                       //! ESD event
  AliAODEvent                       *fAODIn;                     //! Input AOD event
  AliAODEvent                       *fAODOut;                    //! Output AOD event
  AliAODExtension                   *fAODExtension;              //! where we take the jets from can be input or output AOD
  AliVEvent::EOfflineTriggerTypes   fOfflineTrgMask;             //  Mask of offline triggers to accept
  Int_t                             fTriggerType;                //! Trigger type of the event
  Double_t                          fCentrality;                 //! V0M for current event
  Double_t                          fMaxVtxZ;                    //  Maximum z of vertices
  AliESDtrackCuts                   *fEsdTrkCut;                 //! Track cuts for global tracks in ESD
  AliESDtrackCuts                   *fEsdHybCut;                 //! Track cuts for complementary tracks in ESD
  UInt_t                            fFilterMask;                 //  Filter mask to select AOD tracks
  Bool_t                            fRequireITSRefit;            //  Flag to require ITS refit for AOD tracks
  Bool_t                            fRequireSharedClsCut;        //  Flag to require cut on fraction of shared TPC clusters
  Bool_t                            fIsInit;                     //! Flag if all the arraies are successfully retrieved
  TString                           fNonStdFile;                 //  Name of delta aod file to catch the extension
  TString                           fMcParticleArrName;          //  Name of the input mc particles
  TClonesArray                      *fMcParticleArray;           //! Array of input mc particles
  AliNamedArrayI                    *fMcParticlelMap;            //! Array of mc map for EMCal train
  TString                           fEmbTrkArrName;              //  Name of PbPb tracks + PYTHIA tracks
  TClonesArray                      *fEmbTrkArray;               //! Array of PbPb tracks + PYTHIA tracks
  TString                           fTrackArrName;               //  Name of the input track array
  TClonesArray                      *fTrackArray;                //! Array of input tracks
  Int_t                             fTriggerTrkIndex;            //! Index of the trigger track in the event 
  Double_t                          fTriggerTrkPt;               //! Trigger track pt
  Bool_t                            fSwitchOnAvoidTpcHole;       //  Switch on to avoid TPC hole for the recoil jets
  Int_t                             fAvoidTpcHole;               //  Flag if TPC hole is present
  Bool_t                            fCutTPCBoundary;             //  Flag of reqiring trigger tracks stay away from TPC boundary
  Double_t                          fDistToTPCBoundary;          //  Distance to TPC boundary
  Double_t                          fMinTrkPt;                   //  Minimum pt for tracks
  Double_t                          fMaxTrkPt;                   //  Maximum pt for tracks
  Double_t                          fMinTrkEta;                  //  Minimum eta for tracks
  Double_t                          fMaxTrkEta;                  //  Maximum eta for tracks
  Double_t                          fMinTrkPhi;                  //  Minimum phi for tracks
  Double_t                          fMaxTrkPhi;                  //  Maximum phi for tracks
  Double_t                          fRadius;                     //  Jet radius
  TString                           fJetArrName;                 //  Name of the found jet array
  TString                           fPLJetArrName;               //  Name of the embedded PYTHIA jet array on particle level
  TString                           fDLJetArrName;               //  Name of the embedded PYTHIA jet array on detector level
  TClonesArray                      *fJetArray;                  //! Array of the found jets
  TClonesArray                      *fPLJetArray;                //! Array of the embedded PYTHIA jet array on particle level
  TClonesArray                      *fDLJetArray;                //! Array of the embedded PYTHIA jet array on detector level
  TString                           fRhoName;                    //  Name of the rho parameter
  AliRhoParameter                   *fRho;                       //! Rho parameter
  Double_t                          fRhoValue;                   //! Value of the rho parameter
  AliAODJetEventBackground          *fEvtBkg;                    //! Event background for LEGO train
  AliNamedString                    *fPtHardBinName;             //! Pt hard bin param
  Int_t                             fPtHardBin;                  //! Pt hard bin  
  Bool_t                            fRunTrkQA;                   //  Flag to run track QA
  Bool_t                            fRunJetQA;                   //  Flag to run jet QA
  Bool_t                            fRunSingleInclHJet;          //  Flag to run h+jet
  Int_t                             fTTtype;                     //  0-single inclusive; 1-leading (not implemented yet)
  Double_t                          fTTMinPt;                    //  Minimum pt for TT
  Double_t                          fTTMaxPt;                    //  Maximum pt for TT
  Double_t                          fJetPtMin;                   //  Minimum pt for jets
  Bool_t                            fRunPLHJet;                  //  Run h+jet for detector-level jets
  Bool_t                            fRunDLHJet;                  //  Run h+jet for particle-level jets
  Bool_t                            fRunLeadTrkQA;               //  Run QA for trigger hadron
  Bool_t                            fStudyKtEffects;             //  Study kt effects
  Double_t                          fKtValue;                    //  Value of input kt 
  TRandom3                          *fRandom;                    //! Random number generator
  Int_t                             fAODfilterBits[2];           //  AOD track filter bit map
  Bool_t                            fRunBkgFlow;                 //  Vary rho for recoil jets

  // List of histograms
  TList                             *fOutputList;                //! Output list
  TH1F                              *fhEventStat;                //!
  TH1F                              *fhNTrials;                  //!
  TH1F                              *fhPtHardBins;               //!

  // Event properties
  TH1F                              *fhVtxZ[4];                  //!
  TH1F                              *fhCentrality[4];            //!
  TH2F                              *fhRhoVsCent[4];             //!

  // track QA
  TH2F                              *fhTrkPt[4];                 //!
  THnSparse                         *fhTrkQA[4];                 //!
  THnSparse                         *fhTrkPtRes[4];              //!
  THnSparse                         *fhTrkPhiRes[4];             //!

  // jet QA
  THnSparse                         *fhJetPt[4][3];              //!
  THnSparse                         *fhJetArea[4][3];            //!
  THnSparse                         *fhJetQA[4][3];              //!

  // h+jet analysis
  TH1F                              *fhNumberOfTT[4];            //!
  THnSparse                         *fhTTPt[4][3];               //!
  THnSparse                         *fHJetPhiCorr[4][3];         //!
  THnSparse                         *fHJetPhiCorrUp[4];          //!
  THnSparse                         *fHJetPhiCorrDown[4];        //!

  // additional studies
  THnSparse                         *fhLeadTrkQA[4];             //!
  THnSparse                         *fhKtEffects[4];             //!

  AliAnalysisTaskHJetDphi(const AliAnalysisTaskHJetDphi&);            // not implemented
  AliAnalysisTaskHJetDphi &operator=(const AliAnalysisTaskHJetDphi&); // not implemented

  ClassDef(AliAnalysisTaskHJetDphi, 3);
};

#endif
