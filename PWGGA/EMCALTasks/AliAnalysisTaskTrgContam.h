#ifndef AliAnalysisTaskTrgContam_h
#define AliAnalysisTaskTrgContam_h

// $Id$

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
  AliAnalysisTaskTrgContam(); 
  AliAnalysisTaskTrgContam(const char *name);
  virtual ~AliAnalysisTaskTrgContam() {}

  void         UserCreateOutputObjects();
  void         UserExec(Option_t *option);
  void         Terminate(Option_t *);

  Double_t     GetCrossEnergy(const AliVCluster *cluster, Short_t &idmax);
  Double_t     GetMaxCellEnergy(const AliVCluster *cluster, Short_t &id) const; 
  void         FillClusHists();
  void         SetExotCut(Double_t c)                 { fExoticCut          = c;       }
  void         SetGeoName(const char *n)              { fGeoName            = n;       }
  void         SetPeriod(const char *n)               { fPeriod             = n;       }
  void         SetTrainMode(Bool_t t)                 { fIsTrain            = t;       }
  void         SetTrigThresh(Double_t t)              { fTrigThresh         = t;       }
  
 protected:
  TRefArray            *fCaloClusters;           //!pointer to EMCal clusters
  AliESDCaloCells      *fEMCalCells;             //!pointer to EMCal cells
  AliEMCALGeometry     *fGeom;                   // geometry utils
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
  TH1F        *fCellEnergy;              //!cell energy spectrum (all)
  TH2F        *fM02Et;                   //!M02xEt for trigger clusters
  TH2F        *fM02EtTM;                 //!M02xEt for trigger clusters with track matched
  TH2F        *fM02EtExot;               //!M02xEt for trigger clusters of exotic
  TH2F        *fM02EtExotTM;             //!M02xEt for trigger TM clusters of exotic
   
  AliAnalysisTaskTrgContam(const AliAnalysisTaskTrgContam&); // not implemented
  AliAnalysisTaskTrgContam& operator=(const AliAnalysisTaskTrgContam&); // not implemented
  
  ClassDef(AliAnalysisTaskTrgContam, 1); // Trigger contamination class
};
#endif
