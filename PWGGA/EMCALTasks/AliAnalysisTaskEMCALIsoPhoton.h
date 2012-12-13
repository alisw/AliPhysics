#ifndef AliAnalysisTaskEMCALIsoPhoton_h
#define AliAnalysisTaskEMCALIsoPhoton_h

// $Id$

class TH1F;
class TH2F;
class TH3F;
class THnSparse;
class TList;
class TObjArray;
class AliEMCALGeometry;
class AliESDCaloCells;
class AliESDEvent;
class AliESDtrack;
class AliESDtrackCuts;
class AliVCluster;
class AliMCEvent;
class AliStack;
class TParticle;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskEMCALIsoPhoton : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskEMCALIsoPhoton();
  AliAnalysisTaskEMCALIsoPhoton(const char *name);
  virtual ~AliAnalysisTaskEMCALIsoPhoton() {}

  void                   UserCreateOutputObjects();
  void                   UserExec(Option_t *option);
  void                   Terminate(Option_t *);

  void                   GetCeIso(TVector3 vec, Float_t &iso, Float_t &phiband, Float_t &core);
  Double_t               GetCrossEnergy(const AliVCluster *cluster, Short_t &idmax);
  Double_t               GetMaxCellEnergy(const AliVCluster *cluster, Short_t &id) const; 
  void                   GetTrIso(TVector3 vec, Float_t &iso, Float_t &phiband, Float_t &core);
  void                   FillClusHists();
  void                   FillMcHists();
  Float_t                GetClusSource(const AliVCluster *cluster);
  void                   FollowGamma();
  void                   GetDaughtersInfo(int firstd, int lastd, int selfid, const char *indputindent);
  Float_t                GetMcPtSumInCone(Float_t etaclus, Float_t phiclus, Float_t R);
  void                   SetExotCut(Double_t c)                 { fExoticCut          = c;       }
  void                   SetGeoName(const char *n)              { fGeoName            = n;       }
  void                   SetIsoConeR(Double_t r)                { fIsoConeR           = r;       }
  void                   SetPeriod(const char *n)               { fPeriod             = n;       }
  void                   SetTriggerBit(const char *tb)          { fTrigBit            = tb;      }
  void                   SetPrimTrackCuts(AliESDtrackCuts *c)   { fPrTrCuts           = c;       }
  void                   SetTrainMode(Bool_t t)                 { fIsTrain            = t;       }
  void                   SetMcMode(Bool_t mc)                   { fIsMc               = mc;      }
  void                   SetDebugOn(Bool_t d)                   { fDebug              = d;       }
  void                   SetPathStringSelect(char *p)           { fPathStrOpt         = p;       }
  void                   SetEtCut(Double_t ec)                  { fECut               = ec;      }
  
 protected:
  TRefArray             *fCaloClusters;          //!pointer to EMCal clusters
  TObjArray             *fSelPrimTracks;         //!pointer to ESD primary tracks
  TClonesArray          *fTracks;                //!track input array
  AliESDCaloCells       *fEMCalCells;            //!pointer to EMCal cells
  AliESDtrackCuts       *fPrTrCuts;              //pointer to hold the prim track cuts
  AliEMCALGeometry      *fGeom;                  // geometry utils
  TString                fGeoName;               // geometry name (def = EMCAL_FIRSTYEARV1)
  TString                fPeriod;                // string to the LHC period
  TString                fTrigBit;               // string to the trigger bit name
  Bool_t                 fIsTrain;               // variable to set train mode
  Bool_t                 fIsMc;                  // variable to set mc mode
  Bool_t                 fDebug;                 // variable to set on/off debugging printouts
  TString                fPathStrOpt;            // variable to set the name of files to be analyzed (MC only)
  Double_t               fExoticCut;             // variable to set the cut on exotic clusters
  Double_t               fIsoConeR;              // variable to set the isolation cone radius
  Int_t                  fNDimensions;           // variable to set the number of dimensions of n-sparse
  Double_t               fECut;                  // variable to set the minimum E of a cluster
  Int_t                  fTrackMult;             // global variable with the event multiplicity        
  TString                fMcIdFamily;            // string that holds the ids of all particles originated from the prompt photon
  Int_t                  fNClusForDirPho;        // number of clusters from prompt photon per event
  Float_t                fDirPhoPt;              // prompt photon pt (assumes only one per event)
  Float_t                fHigherPtCone;          // higher pt inside the cone around the candidate

  
 private:
  AliESDEvent *fESD;      //! ESD object
  AliMCEvent  *fMCEvent;  //! MC event object
  AliStack    *fStack;    //!MC particles stack object

  TList       *fOutputList; //! Output list
  //histograms for events with 1+ track pt>1
  TH1F        *fEvtSel;                    //!evt selection counter: 0=all trg, 1=pv cut 
  TH1F        *fNClusEt10;                 //!number of clusters w/ Et>10 in the event
  TH1F        *fRecoPV;                    //!histogram to record if an event has a prim. vert.
  TH1F        *fPVtxZ;                     //!primary vertex Z before cut
  TH1F        *fTrMultDist;                //!track multiplicity distribution
  TH3F        *fMCDirPhotonPtEtaPhi;       //!direct produced photon pt
  TH1F        *fDecayPhotonPtMC;           //!decay photon pt
  TH2F        *fCellAbsIdVsAmpl;           //!cell abs id vs cell amplitude (energy)
  TH2F        *fNClusHighClusE;            //!total number of clusters vs. highest clus energy in the event
  TH2F        *fHigherPtConeM02;           //!M02 vs. the higher pt of a track inside the cone
  TH2F        *fClusEtMcPt;                //!cluster et x mc-pt
  TH2F        *fClusMcDetaDphi;            //!delta-eta x delta-phi(reco-mc)
  TH2F        *fNClusPerPho;               //!delta-eta x delta-phi(reco-mc)
  TH2F        *fMcPtInConeBG;              //!sum of mc-pt of "primary" particles inside de cone, as a function of NET-ISO in BG template
  TH2F        *fMcPtInConeSBG;             //!sum of mc-pt of "primary" particles inside de cone, as a function of NET-ISO in SBG range
  TH2F        *fMcPtInConeBGnoUE;          //!sum of mc-pt of "primary" particles inside de cone, as a function of NET-ISO in BG template no UE sub
  TH2F        *fMcPtInConeSBGnoUE;         //!sum of mc-pt of "primary" particles inside de cone, as a function of NET-ISO in SBG range no UE sub
  TH3F        *fMCDirPhotonPtEtaPhiNoClus; //!pt x eta x phi for prompt photons that didn't produce clusters
  THnSparse   *fHnOutput;                  //!Output matrix with 7 dimensions

  AliAnalysisTaskEMCALIsoPhoton(const AliAnalysisTaskEMCALIsoPhoton&); // not implemented
  AliAnalysisTaskEMCALIsoPhoton& operator=(const AliAnalysisTaskEMCALIsoPhoton&); // not implemented
  
  ClassDef(AliAnalysisTaskEMCALIsoPhoton, 1); // Class to analyse isolated photons
};
#endif
