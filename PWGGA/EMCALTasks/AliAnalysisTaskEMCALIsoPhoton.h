#ifndef AliAnalysisTaskEMCALIsoPhoton_h
#define AliAnalysisTaskEMCALIsoPhoton_h

// $Id$

class TH1F;
class TH2F;
class TObjArray;
class AliEMCALGeometry;
class AliESDCaloCells;
class AliESDEvent;
class AliESDtrack;
class AliESDtrackCuts;
class AliVCluster;

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
  void                   SetExotCut(Double_t c)                 { fExoticCut          = c;       }
  void                   SetGeoName(const char *n)              { fGeoName            = n;       }
  void                   SetIsoConeR(Double_t r)                { fIsoConeR           = r;       }
  void                   SetPeriod(const char *n)               { fPeriod             = n;       }
  void                   SetTriggerBit(const char *tb)          { fTrigBit            = tb;      }
  void                   SetPrimTrackCuts(AliESDtrackCuts *c)   { fPrTrCuts           = c;       }
  void                   SetTrainMode(Bool_t t)                 { fIsTrain            = t;       }
  
 protected:
  TRefArray             *fCaloClusters;          //!pointer to EMCal clusters
  TObjArray             *fSelPrimTracks;         //!pointer to ESD primary tracks
  AliESDCaloCells       *fEMCalCells;            //!pointer to EMCal cells
  AliESDtrackCuts       *fPrTrCuts;              //!pointer to hold the prim track cuts
  AliEMCALGeometry      *fGeom;                  // geometry utils
  TString                fGeoName;               // geometry name (def = EMCAL_FIRSTYEARV1)
  TString                fPeriod;                // string to the LHC period
  TString                fTrigBit;               // string to the trigger bit name
  Bool_t                 fIsTrain;               // variable to set train mode
  Double_t               fExoticCut;             // variable to set the cut on exotic clusters
  Double_t               fIsoConeR;              // variable to set the isolation cone radius
  
 private:
  AliESDEvent *fESD;      //! ESD object
  TList       *fOutputList; //! Output list
  //histograms for events with 1+ track pt>1
  TH1F        *fEvtSel;                  //!evt selection counter: 0=all trg, 1=pv cut 
  TH1F        *fPVtxZ;                   //!primary vertex Z before cut
  TH2F        *fCellAbsIdVsAmpl;         //!cell abs id vs cell amplitude (energy)
  TH2F        *fNClusHighClusE;          //!total number of clusters vs. highest clus energy in the event
  TH2F        *fM02Et;                   //!M02 vs Et for all clusters
  TH2F        *fM02EtTM;                 //!M02 vs Et for clusters with track-match (dEta=0.01 && dPhi=0.025)
  TH2F        *fM02EtCeIso1;             //!M02 vs Et for clusters with isolation neutral Et<1GeV
  TH2F        *fM02EtCeIso2;             //!M02 vs Et for clusters with isolation neutral Et<2GeV
  TH2F        *fM02EtCeIso5;             //!M02 vs Et for clusters with isolation neutral Et<5GeV
  TH2F        *fM02EtTrIso1;             //!M02 vs Et for clusters with isolation charged Et<1GeV
  TH2F        *fM02EtTrIso2;             //!M02 vs Et for clusters with isolation charged Et<2GeV
  TH2F        *fM02EtTrIso5;             //!M02 vs Et for clusters with isolation charged Et<5GeV
  TH2F        *fM02EtAllIso1;            //!M02 vs Et for clusters with isolation total Et<1GeV
  TH2F        *fM02EtAllIso2;            //!M02 vs Et for clusters with isolation total Et<2GeV
  TH2F        *fM02EtAllIso5;            //!M02 vs Et for clusters with isolation total Et<5GeV
  TH2F        *fCeIsoVsEtPho;            //!Neutral isolation Et vs. cluster Et, 0.10<M02<0.30
  TH2F        *fTrIsoVsEtPho;            //!Charged isolation Et vs. cluster Et, 0.10<M02<0.30
  TH2F        *fAllIsoVsEtPho;           //!Total isolation Et vs. cluster Et, 0.10<M02<0.30
  TH2F        *fCeIsoVsEtPi0;            //!Neutral isolation Et vs. cluster Et, pi0 selection (BG)
  TH2F        *fTrIsoVsEtPi0;            //!Charged isolation Et vs. cluster Et, pi0 selection (BG)
  TH2F        *fAllIsoVsEtPi0;           //!Total isolation Et vs. cluster Et, pi0 selection (BG)
  //track matched stuff
  TH2F        *fM02EtCeIso1TM;           //!Track-matched M02 vs Et for clusters with isolation neutral Et<1GeV
  TH2F        *fM02EtCeIso2TM;           //!Track-matched M02 vs Et for clusters with isolation neutral Et<2GeV
  TH2F        *fM02EtCeIso5TM;           //!Track-matched M02 vs Et for clusters with isolation neutral Et<5GeV
  TH2F        *fM02EtTrIso1TM;           //!Track-matched M02 vs Et for clusters with isolation charged Et<1GeV
  TH2F        *fM02EtTrIso2TM;           //!Track-matched M02 vs Et for clusters with isolation charged Et<2GeV
  TH2F        *fM02EtTrIso5TM;           //!Track-matched M02 vs Et for clusters with isolation charged Et<5GeV
  TH2F        *fM02EtAllIso1TM;          //!Track-matched M02 vs Et for clusters with isolation total Et<1GeV
  TH2F        *fM02EtAllIso2TM;          //!Track-matched M02 vs Et for clusters with isolation total Et<2GeV
  TH2F        *fM02EtAllIso5TM;          //!Track-matched M02 vs Et for clusters with isolation total Et<5GeV
  TH2F        *fCeIsoVsEtPhoTM;          //!Track-matched Neutral isolation Et vs. cluster Et, 0.10<M02<0.30
  TH2F        *fTrIsoVsEtPhoTM;          //!Track-matched Charged isolation Et vs. cluster Et, 0.10<M02<0.30
  TH2F        *fAllIsoVsEtPhoTM;         //!Track-matched Total isolation Et vs. cluster Et, 0.10<M02<0.30
  TH2F        *fCeIsoVsEtPi0TM;          //!Track-matched Neutral isolation Et vs. cluster Et, pi0 selection (BG)
  TH2F        *fTrIsoVsEtPi0TM;          //!Track-matched Charged isolation Et vs. cluster Et, pi0 selection (BG)
  TH2F        *fAllIsoVsEtPi0TM;         //!Track-matched Total isolation Et vs. cluster Et, pi0 selection (BG)
   
  AliAnalysisTaskEMCALIsoPhoton(const AliAnalysisTaskEMCALIsoPhoton&); // not implemented
  AliAnalysisTaskEMCALIsoPhoton& operator=(const AliAnalysisTaskEMCALIsoPhoton&); // not implemented
  
  ClassDef(AliAnalysisTaskEMCALIsoPhoton, 1); // Class to analyse isolated photons
};
#endif
