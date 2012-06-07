#ifndef AliAnalysisTaskEMCALIsoPhoton_cxx
#define AliAnalysisTaskEMCALIsoPhoton_cxx

class TH1;
class TH2;
class TObjArray;
class AliESDEvent;
class AliESDtrack;
class AliESDtrackCuts;
class AliVCluster;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskEMCALIsoPhoton : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskEMCALIsoPhoton() : 
  AliAnalysisTaskSE(), 
  
    fCaloClusters(0),
    fSelPrimTracks(0),
    fEMCalCells(0),
    fPrTrCuts(0),

    fGeom(0x0),
  
    fESD(0),
  
    fOutputList(0),

    fEvtSel(0),

    fPVtxZ(0),                   //!primary vertex Z before cut
    fCellAbsIdVsAmpl(0),         //!cell abs id vs cell amplitude (energy)
    fNClusHighClusE(0),          //!total number of clusters vs. highest clus energy in the event
    fM02Et(0),                   //!M02 vs Et for all clusters
    fM02EtTM(0),                 //!M02 vs Et for clusters with track-match (dEta=0.01 && dPhi=0.025)
    fM02EtCeIso1(0),             //!M02 vs Et for clusters with isolation neutral Et<1GeV
    fM02EtCeIso2(0),             //!M02 vs Et for clusters with isolation neutral Et<2GeV
    fM02EtCeIso5(0),             //!M02 vs Et for clusters with isolation neutral Et<5GeV
    fM02EtTrIso1(0),             //!M02 vs Et for clusters with isolation charged Et<1GeV
    fM02EtTrIso2(0),             //!M02 vs Et for clusters with isolation charged Et<2GeV
    fM02EtTrIso5(0),             //!M02 vs Et for clusters with isolation charged Et<5GeV
    fM02EtAllIso1(0),            //!M02 vs Et for clusters with isolation total Et<1GeV
    fM02EtAllIso2(0),            //!M02 vs Et for clusters with isolation total Et<2GeV
    fM02EtAllIso5(0),            //!M02 vs Et for clusters with isolation total Et<5GeV
    fCeIsoVsEtPho(0),            //!Neutral isolation Et vs. cluster Et, 0.05<M02<0.30
    fTrIsoVsEtPho(0),            //!Charged isolation Et vs. cluster Et, 0.05<M02<0.30
    fAllIsoVsEtPho(0),           //!Total isolation Et vs. cluster Et, 0.05<M02<0.30
    //track matched stuff
    fM02EtCeIso1TM(0),           //!Track-matched M02 vs Et for clusters with isolation neutral Et<1GeV
    fM02EtCeIso2TM(0),           //!Track-matched M02 vs Et for clusters with isolation neutral Et<2GeV
    fM02EtCeIso5TM(0),           //!Track-matched M02 vs Et for clusters with isolation neutral Et<5GeV
    fM02EtTrIso1TM(0),           //!Track-matched M02 vs Et for clusters with isolation charged Et<1GeV
    fM02EtTrIso2TM(0),           //!Track-matched M02 vs Et for clusters with isolation charged Et<2GeV
    fM02EtTrIso5TM(0),           //!Track-matched M02 vs Et for clusters with isolation charged Et<5GeV
    fM02EtAllIso1TM(0),          //!Track-matched M02 vs Et for clusters with isolation total Et<1GeV
    fM02EtAllIso2TM(0),          //!Track-matched M02 vs Et for clusters with isolation total Et<2GeV
    fM02EtAllIso5TM(0),          //!Track-matched M02 vs Et for clusters with isolation total Et<5GeV
    fCeIsoVsEtPhoTM(0),          //!Track-matched Neutral isolation Et vs. cluster Et, 0.05<M02<0.30
    fTrIsoVsEtPhoTM(0),          //!Track-matched Charged isolation Et vs. cluster Et, 0.05<M02<0.30
    fAllIsoVsEtPhoTM(0)          //!Track-matched Total isolation Et vs. cluster Et, 0.05<M02<0.30


    
  
  
  {}
  AliAnalysisTaskEMCALIsoPhoton(const char *name);
  virtual ~AliAnalysisTaskEMCALIsoPhoton() {}

  void   UserCreateOutputObjects();
  void   UserExec(Option_t *option);
  void   Terminate(Option_t *);

  void         SetGeoName(const char *n)              { fGeoName            = n;       }
  void         SetPeriod(const char *n)               { fPeriod             = n;       }
  void         SetTrainMode(Bool_t t)                 { fIsTrain            = t;       }
  void         SetExotCut(Double_t c)                 { fExoticCut          = c;       }
  void         SetIsoConeR(Double_t r)                { fIsoConeR           = r;       }
  void         SetPrimTrackCuts(AliESDtrackCuts *c)   { fPrTrCuts           = c;       }
  void         FillClusHists();
  void         GetCeIso(TVector3 vec, Float_t &iso, Float_t &phiband, Float_t &core);
  void         GetTrIso(TVector3 vec, Float_t &iso, Float_t &phiband, Float_t &core);
  Double_t     GetCrossEnergy(const AliVCluster *cluster, Short_t &idmax);
  Double_t     GetMaxCellEnergy(const AliVCluster *cluster, Short_t &id) const; 




  
 protected:
  TRefArray             *fCaloClusters;          //!pointer to EMCal clusters
  TObjArray             *fSelPrimTracks;         //!pointer to ESD primary tracks
  AliESDCaloCells       *fEMCalCells;            //!pointer to EMCal cells
  AliESDtrackCuts       *fPrTrCuts;              //!pointer to hold the prim track cuts
  AliEMCALGeometry      *fGeom;                   // geometry utils
  TString               fGeoName;                // geometry name (def = EMCAL_FIRSTYEARV1)
  TString               fPeriod;                 // string to the LHC period
  Bool_t                fIsTrain;                //variable to set train mode
  Double_t              fExoticCut;              //variable to set the cut on exotic clusters
  Double_t              fIsoConeR;               //variable to set the isolation cone radius
  
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
  TH2F        *fCeIsoVsEtPho;            //!Neutral isolation Et vs. cluster Et, 0.05<M02<0.30
  TH2F        *fTrIsoVsEtPho;            //!Charged isolation Et vs. cluster Et, 0.05<M02<0.30
  TH2F        *fAllIsoVsEtPho;           //!Total isolation Et vs. cluster Et, 0.05<M02<0.30
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
  TH2F        *fCeIsoVsEtPhoTM;          //!Track-matched Neutral isolation Et vs. cluster Et, 0.05<M02<0.30
  TH2F        *fTrIsoVsEtPhoTM;          //!Track-matched Charged isolation Et vs. cluster Et, 0.05<M02<0.30
  TH2F        *fAllIsoVsEtPhoTM;         //!Track-matched Total isolation Et vs. cluster Et, 0.05<M02<0.30



   
  AliAnalysisTaskEMCALIsoPhoton(const AliAnalysisTaskEMCALIsoPhoton&); // not implemented
  AliAnalysisTaskEMCALIsoPhoton& operator=(const AliAnalysisTaskEMCALIsoPhoton&); // not implemented
  
  ClassDef(AliAnalysisTaskEMCALIsoPhoton, 1); // example of analysis
};

#endif
