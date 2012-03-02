#ifndef AliAnalysisTaskTrgContam_cxx
#define AliAnalysisTaskTrgContam_cxx

class TH1;
class TH2;
class TObjArray;
class AliESDEvent;
class AliESDtrack;
class AliESDtrackCuts;
class AliVCluster;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskTrgContam : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskTrgContam() : 
    AliAnalysisTaskSE(), 
    fCaloClusters(0),
    fEMCalCells(0),
    fGeom(0x0),
    fGeoName(),
    fPeriod(),
    fIsTrain(0),
    fTrigThresh(0),
    fExoticCut(0),
    fESD(0),
    fOutputList(0),
    fEvtSel(0),
    fClusEt(0),
    fClusEtTM(0),
    fClusEtLead(0),
    fClusEtSubLead(0),
    fClusEtLeadTM(0),
    fClusEtSubLeadTM(0),
    fClusEtExotic(0), 
    fClusEtExoticTM(0),
    fClusEtSingleExotic(0),
    fM02Et(0),
    fM02EtTM(0),
    fM02EtExot(0),
    fM02EtExotTM(0) {}
  AliAnalysisTaskTrgContam(const char *name);
  virtual ~AliAnalysisTaskTrgContam() {}

  void   UserCreateOutputObjects();
  void   UserExec(Option_t *option);
  void   Terminate(Option_t *);

  void         SetGeoName(const char *n)              { fGeoName            = n;       }
  void         SetPeriod(const char *n)               { fPeriod             = n;       }
  void         SetTrainMode(Bool_t t)                 { fIsTrain            = t;       }
  void         SetTrigThresh(Double_t t)              { fTrigThresh         = t;       }
  void         SetExotCut(Double_t c)                 { fExoticCut          = c;       }
  void         FillClusHists();
  Double_t     GetCrossEnergy(const AliVCluster *cluster, Short_t &idmax);
  Double_t     GetMaxCellEnergy(const AliVCluster *cluster, Short_t &id) const; 
  
 protected:
  TRefArray             *fCaloClusters;          //!pointer to EMCal clusters
  AliESDCaloCells       *fEMCalCells;            //!pointer to EMCal cells
  AliEMCALGeometry      *fGeom;                   // geometry utils
  TString               fGeoName;                // geometry name (def = EMCAL_FIRSTYEARV1)
  TString               fPeriod;                 // string to the LHC period
  Bool_t                fIsTrain;                //variable to set train mode
  Double_t              fTrigThresh;             //variable to set the trigger threshold
  Double_t              fExoticCut;              //variable to set the cut on exotic clusters
  
 private:
  AliESDEvent *fESD;      //! ESD object
  TList       *fOutputList; //! Output list
  //histograms for events with 1+ track pt>1
  TH1F        *fEvtSel;                  //!evt selection counter: 0=all trg, 1=pv cut 
  TH1F        *fClusEt;                  //!cluster Et spectrum
  TH1F        *fClusEtTM;                //!cluster(matched to a track) Et spectrum
  TH1F        *fClusEtLead;              //!leading trigger cluster Et
  TH1F        *fClusEtSubLead;           //!sub-leading trigger cluster Et
  TH1F        *fClusEtLeadTM;            //!leading trigger cluster (TM) Et
  TH1F        *fClusEtSubLeadTM;         //!subleading trigger cluster (TM) Et
  TH1F        *fClusEtExotic;            //!exotic trigger clusters Et
  TH1F        *fClusEtExoticTM;          //!exotic trigger clusters (TM) Et
  TH1F        *fClusEtSingleExotic;      //!exotic trigger only clusters Et 
  TH2F        *fM02Et;                   //!M02xEt for trigger clusters
  TH2F        *fM02EtTM;                 //!M02xEt for trigger clusters with track matched
  TH2F        *fM02EtExot;               //!M02xEt for trigger clusters of exotic
  TH2F        *fM02EtExotTM;             //!M02xEt for trigger TM clusters of exotic
   
  AliAnalysisTaskTrgContam(const AliAnalysisTaskTrgContam&); // not implemented
  AliAnalysisTaskTrgContam& operator=(const AliAnalysisTaskTrgContam&); // not implemented
  
  ClassDef(AliAnalysisTaskTrgContam, 1); // example of analysis
};
#endif
