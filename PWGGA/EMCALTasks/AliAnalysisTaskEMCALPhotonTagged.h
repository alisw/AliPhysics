#ifndef AliAnalysisTaskEMCALPhotonTagged_h
#define AliAnalysisTaskEMCALPhotonTagged_h

// $Id$

class TH1F;
class TH2F;
class TH3F;
class THnSparse;
class TList;
class TObjArray;
class AliEMCALGeometry;
class AliOADBContainer;
class AliESDCaloCells;
class AliESDEvent;
class AliESDtrack;
class AliESDtrackCuts;
class AliAODEvent;
class AliAODCaloCells;
class AliVCluster;
class AliVCaloCells;
class AliMCEvent;
class AliStack;
class TParticle;
class AliAODMCParticle;
class TGeoHMatrix;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskEMCALPhotonTagged : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskEMCALPhotonTagged();
  AliAnalysisTaskEMCALPhotonTagged(const char *name);
  virtual ~AliAnalysisTaskEMCALPhotonTagged() {}

  void                   UserCreateOutputObjects();
  void                   UserExec(Option_t *option);
  void                   Terminate(Option_t *);

  Double_t               GetCrossEnergy(const AliVCluster *cluster, Short_t &idmax);
  Double_t               GetMaxCellEnergy(const AliVCluster *cluster, Short_t &id) const; 
  Double_t               GetTrackMatchedPt(Int_t matchIndex);
  void                   FillClusHists();
  void                   FillMcHists();
  void                   FillQA();
  Float_t                GetClusSource(const AliVCluster *cluster); //returns the pt of the prompt-photon; if doesn't descend from pro-pho returns -0.1
  void                   FollowGamma();
  void                   GetDaughtersInfo(int firstd, int lastd, int selfid, const char *indputindent);
  bool                   IsMcPDG(Int_t mclabel, Int_t PDG);
  bool                   IsExotic(AliVCluster *c);
  void                   CheckTriggerPatch();
  Double_t               SmearM02(Double_t m02);
  Int_t                  GetNumberOfLocalMaxima(AliVCluster* cluster, AliVCaloCells* cells);
  Int_t                  GetNumberOfLocalMaxima(AliVCluster* cluster, AliVCaloCells* cells,
						Int_t *absIdList,     Float_t *maxEList);
  Bool_t                 AreNeighbours(Short_t absId1, Short_t absId2);
  Bool_t                 IsPi0M02(Double_t M02, Double_t Et);
  AliVCaloCells          *GetVCaloCells();
  //setters
  void                   SetExotCut(Double_t c)                 { fExoticCut          = c;       }
  void                   SetGeoName(const char *n)              { fGeoName            = n;       }
  void                   SetPeriod(const char *n)               { fPeriod             = n;       }
  void                   SetTriggerBit(const char *tb)          { fTrigBit            = tb;      }
  void                   SetPrimTrackCuts(AliESDtrackCuts *c)   { fPrTrCuts           = c;       }
  void                   SetComplTrackCuts(AliESDtrackCuts *c)  { fCompTrCuts         = c;       }
  void                   SetTrainMode(Bool_t t)                 { fIsTrain            = t;       }
  void                   SetMcMode(Bool_t mc)                   { fIsMc               = mc;      }
  void                   SetDebugOn(Bool_t d)                   { fDebug              = d;       }
  void                   SetPathStringSelect(char *p)           { fPathStrOpt         = p;       }
  void                   SetEtCut(Double_t ec)                  { fECut               = ec;      }
  void                   SetImportGeometryFromFile(Bool_t  im, 
                                           TString pa = "")     { fImportGeometryFromFile = im ; 
                                                                  fImportGeometryFilePath = pa ; }    
  void                  SetTrackFilterBit(ULong_t bit)          { fFilterBit = bit;  }
  void                  SetHybridOn()                           { fSelHybrid = kTRUE; }
  void                  SetFillQA()                             { fFillQA = kTRUE; }
  void                  SelectCPVFromTrack(Bool_t b)            { fCpvFromTrack = b; }
  void                  SetEtPtHistoBinning(Int_t n, 
					    Double_t lowx, 
					    Double_t highx)     { fNBinsPt = n; fPtBinLowEdge = lowx; fPtBinHighEdge = highx; }
  void                  SetPileUpRejSPD()                       { fPileUpRejSPD       = kTRUE;   }
  void                  SetDistanceToBadCh(Double_t d)          { fDistToBadChan      = d;       }
  void                  SetSmearM02On(Double_t s)               { fSigmaSmear         = s;       }
  void                  SetNLMCut(Int_t n)                      { fNLMCut             = n;       }
 protected:
  TObjArray             *fESDClusters;           //!pointer to EMCal clusters
  TObjArray             *fAODClusters;           //!pointer to EMCal clusters
  TObjArray             *fSelPrimTracks;         //!pointer to ESD primary tracks
  TClonesArray          *fTracks;                //!track input array
  TClonesArray          *fAODMCParticles;        //!MC particles array for AOD analysis
  AliESDCaloCells       *fESDCells;              //!pointer to EMCal cells, esd
  AliAODCaloCells       *fAODCells;              //!pointer to EMCal cells, aod  
  AliVCaloCells         *fVCells;                //!pointer to EMCal virtual cells
  AliESDtrackCuts       *fPrTrCuts;              //pointer to hold the prim track cuts
  AliESDtrackCuts       *fCompTrCuts;            //pointer to hold complementary track cuts (a la Gustavo)
  AliEMCALGeometry      *fGeom;                  // geometry utils
  TString                fGeoName;               // geometry name (def = EMCAL_FIRSTYEARV1)
  AliOADBContainer      *fOADBContainer;         //!OADB container used to load misalignment matrices
  TVector3               fVecPv;                 // vector to hold the event's primary vertex
  TString                fPeriod;                // string to the LHC period
  TString                fTrigBit;               // string to the trigger bit name
  Bool_t                 fIsTrain;               // variable to set train mode
  Bool_t                 fIsMc;                  // variable to set mc mode
  Bool_t                 fDebug;                 // variable to set on/off debugging printouts
  TString                fPathStrOpt;            // variable to set the name of files to be analyzed (MC only)
  Double_t               fExoticCut;             // variable to set the cut on exotic clusters
  Int_t                  fNDimensions;           // variable to set the number of dimensions of n-sparse
  Double_t               fECut;                  // variable to set the minimum E of a cluster
  Int_t                  fTrackMult;             // global variable with the event multiplicity        
  TString                fMcIdFamily;            // string that holds the ids of all particles originated from the prompt photon
  Int_t                  fNClusForDirPho;        // number of clusters from prompt photon per event
  Float_t                fDirPhoPt;              // prompt photon pt (assumes only one per event)
  Bool_t                 fImportGeometryFromFile;  // Import geometry settings in geometry.root file
  TString                fImportGeometryFilePath;  // path fo geometry.root file
  Double_t               fMaxPtTrack;            //track with highest pt in event
  Double_t               fMaxEClus;              //cluster with highest energy in event
  Int_t                  fNCells50;              // variable to keep the number of cells with E>50 MeV
  ULong_t                fFilterBit;             // Track selection bit, for AODs 
  Bool_t                 fSelHybrid;             // bool to select hybrid tracks
  Bool_t                 fFillQA;                // bool to fill the QA plots
  TString                fClusIdFromTracks;      // string to hold the list of cluster ids given by tracks
  Bool_t                 fCpvFromTrack;          // set the track-matching method to track->GetEMCALcluster()
  Int_t                  fNBinsPt;               // set the number of bins in axis of histograms filled with pt (or Et)
  Double_t               fPtBinLowEdge;          // low edge of the first pt (Et) bin
  Double_t               fPtBinHighEdge;         // high edge of the first pt (Et) bin
  Int_t                  fNCuts;                 // number of cuts (QA purposes)

  Bool_t                 fPileUpRejSPD;          // flag to set pile-up rejection via SPD (multiple vertices)
  Double_t               fDistToBadChan;         // distance to bad channel
  Double_t               fSigmaSmear;            //std dev of the gaussian smearing for the MC M02 values
  Int_t                  fNLMCut;                //maximum of local maxima in a cluster

 private:
  AliESDEvent *fESD;      //! ESD object
  AliAODEvent *fAOD;      //! AOD object
  AliVEvent   *fVEvent;   //! AliVEvent
  AliMCEvent  *fMCEvent;  //! MC event object
  AliStack    *fStack;    //!MC particles stack object
  TGeoHMatrix *fGeomMatrix[12];//! Geometry misalignment matrices for EMCal
  TList       *fOutputList; //! Output list
  //histograms for events with 1+ track pt>1
  TH1F        *fEvtSel;                    //!evt selection counter: 0=all trg, 1=pv cut 
  TH1F        *fNClusEt10;                 //!number of clusters w/ Et>10 in the event
  TH1F        *fClusArrayNames;            //!cluster array name for each analysed event (CaloClusters, EmcCaloClusters or Other)
  TH1F        *fRecoPV;                    //!histogram to record if an event has a prim. vert.
  TH1F        *fPVtxZ;                     //!primary vertex Z before cut
  TH1F        *fTrMultDist;                //!track multiplicity distribution
  TH3F        *fMCDirPhotonPtEtaPhi;       //!direct produced photon pt, eta, phi
  TH1F        *fDecayPhotonPtMC;           //!decay photon pt
  TH2F        *fCellAbsIdVsAmpl;           //!cell abs id vs cell amplitude (energy)
  TH2F        *fNClusHighClusE;            //!total number of clusters vs. highest clus energy in the event
  TH2F        *fClusEtMcPt;                //!cluster et x mc-pt
  TH2F        *fClusMcDetaDphi;            //!delta-eta x delta-phi(reco-mc)
  TH2F        *fNClusPerPho;               //!delta-eta x delta-phi(reco-mc)
  THnSparse   *fHnOutput;                  //!Output matrix with 7 dimensions

  //QA histos
  TList       *fQAList;           //!output list holding QA histos
  TH1F        *fNTracks;          //!number of tracks from Array->GetEntries()
  TH1F        *fEmcNCells;        //!number of emcal cells in the event
  TH1F        *fEmcNClus;         //!# of emcal clusters
  TH1F        *fEmcNClusCut;      //!# of clusters in an event with at least 1 clus with E > fECut ("triggered event")
  TH1F        *fNTracksECut;      //!number of tracks from Array->GetEntries() in "triggered event"
  TH1F        *fEmcNCellsCut;     //!number of emcal cells in a in "triggered event"
  TH1F        *fEmcClusETM1;      //!emcal track matched cluster energy (TracDx,z method)
  TH1F        *fEmcClusETM2;      //!emcal track matched cluster energy (track->GetEMCALcluster() method)
  TH1F        *fEmcClusNotExo;    //!cluster energy (exotics removed)
  TH2F        *fEmcClusEClusCuts; //!cluster E spectrum per cluster cut (none, exotic, exo+cpv1, exo+cpv1+time, exo+cpv1+time+m02)
  TH2F        *fEmcClusEPhi;      //!cluster E spectrum vs. phi
  TH2F        *fEmcClusEPhiCut;   //!cluster E spectrum vs. phi in "triggered event"
  TH2F        *fEmcClusEEta;      //!cluster E spectrum vs. eta
  TH2F        *fEmcClusEEtaCut;   //!cluster E spectrum vs. eta in "triggered event"
  TH2F        *fTrackPtPhi;       //!selected tracks pt vs. phi
  TH2F        *fTrackPtPhiCut;    //!selected tracks pt vs. phi in "triggered event"
  TH2F        *fTrackPtEta;       //!selected tracks pt vs. eta
  TH2F        *fTrackPtEtaCut;    //!selected tracks pt vs. eta in "triggered event"
  TH2F        *fMaxCellEPhi;      //!max cell energy vs. cell phi
  TH2F        *fDetaDphiFromTM;   //!dphi vs deta of track->GetEMCALcluster() clusters
  TH2F        *fEoverPvsE;        //!E/p for tracks with 80<TPCsignal<100 vs cluster E (check material)
  TH2F        *fTrackDEtaDPhiPho; //!dEta-dPhi of TM for inclusive true pi0 photons (MC only)
  TH2F        *fTrackDEtaDPhiPi;  //!dEta-dPhi of TM for inclusive true pions (MC only)
  TH2F        *fTrackDEtaDPhiPPho;//!dEta-dPhi of TM for inclusive true prompt photons (MC only)
  //trigger histos
  TH1F        *fETrigg;           //!energy returned by trigger patch info
  //shower shape studies



  AliAnalysisTaskEMCALPhotonTagged(const AliAnalysisTaskEMCALPhotonTagged&); // not implemented
  AliAnalysisTaskEMCALPhotonTagged& operator=(const AliAnalysisTaskEMCALPhotonTagged&); // not implemented
  
  ClassDef(AliAnalysisTaskEMCALPhotonTagged, 1); // Class to analyse isolated photons
};
#endif
