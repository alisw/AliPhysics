#ifndef AliAnalysisTaskShoShaTests_h
#define AliAnalysisTaskShoShaTests_h

// $Id$

class TH1;
class TH2;
class TObjArray;
class AliESDEvent;
class AliESDtrack;
class AliESDtrackCuts;
class AliVCluster;
class AliEMCALGeometry;
class AliOADBContainer;
class TGeoHMatrix;


#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskShoShaTests : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskShoShaTests(); 
  AliAnalysisTaskShoShaTests(const char *name);
  virtual ~AliAnalysisTaskShoShaTests() {}

  void                  UserCreateOutputObjects();
  void                  UserExec(Option_t *option);
  void                  Terminate(Option_t *);

  Double_t              GetCrossEnergy(const AliVCluster *cluster, Short_t &idmax);
  Double_t              GetMaxCellEnergy(const AliVCluster *cluster, Short_t &id) const; 
  Double_t              NeutClusPairInvMass(const AliVCluster *cl1, Int_t ic);
  Int_t                 GetAncestorPdg(const Int_t label);
  void                  FillClusHists();
  void                  SetExotCut(Double_t c)                 { fExoticCut          = c;       }
  void                  SetGeoName(const char *n)              { fGeoName            = n;       }
  void                  SetPeriod(const char *n)               { fPeriod             = n;       }
  void                  SetTrainMode(Bool_t t)                 { fIsTrain            = t;       }
  void                  SetTrigThresh(Double_t t)              { fTrigThresh         = t;       }
  
 protected:
  TRefArray            *fCaloClusters;            //!pointer to EMCal clusters
  AliESDCaloCells      *fEMCalCells;              //!pointer to EMCal cells
  AliEMCALGeometry     *fGeom;                    // geometry utils
  TString               fGeoName;                 // geometry name (def = EMCAL_FIRSTYEARV1)
  AliOADBContainer      *fOADBContainer;         //!OADB container used to load misalignment matrices
  TString               fPeriod;                  // string to the LHC period
  Bool_t                fIsTrain;                 // variable to set train mode
  Bool_t                fIsMC;                    // flag for MC events
  Double_t              fTrigThresh;              // variable to set the trigger threshold
  Double_t              fExoticCut;               // variable to set the cut on exotic clusters
  Double_t              fEClusCut;                // variable to set the minimum cluster E
  Double_t              fLowPi0MCut;              // variable to set the minimum pi0 peak cut
  Double_t              fHighPi0MCut;             // variable to set the maximum pi0 peak cut                                                                                                                

  
 private:
  AliESDEvent          *fESD;                     //!esd event
  AliMCEvent           *fMCEvent;                 //! MC event object 
  AliStack             *fStack;                   //!MC particles stack object   
  TGeoHMatrix          *fGeomMatrix[12];          //! Geometry misalignment matrices for EMCal
  TList                *fOutputList;              //!output list
  Double_t             *fPvPos;                   //!Prim Vertex position array of coord
  TH1F                 *fEvtSel;                  //!evt selection counter: 0=all trg, 1=pv cut 
  TH1F                 *fPVZ;                     //!primary vertex z distribution
  TH1F                 *fClusEt;                  //!cluster Et spectrum
  TH1F                 *fClusEtTM;                //!cluster(matched to a track) Et spectrum
  TH1F                 *fClusEtLead;              //!leading trigger cluster Et
  TH1F                 *fClusEtSubLead;           //!sub-leading trigger cluster Et
  TH1F                 *fClusEtLeadTM;            //!leading trigger cluster (TM) Et
  TH1F                 *fClusEtSubLeadTM;         //!subleading trigger cluster (TM) Et
  TH1F                 *fClusEtExotic;            //!exotic trigger clusters Et
  TH1F                 *fClusEtExoticTM;          //!exotic trigger clusters (TM) Et
  TH1F                 *fClusEtSingleExotic;      //!exotic trigger only clusters Et 
  TH1F                 *fCellEnergy;              //!cell energy spectrum (all)
  TH1F                 *fInvMassEMCNN;            //!inv mass of EMC neutral cluster pairs
  TH2F                 *fM02Et;                   //!M02xEt for trigger clusters
  TH2F                 *fM02EtPi0MassClCl;        //!M02xEt for trigger clusters
  TH2F                 *fM02EtPi0MassClClTruPi0;  //!M02xEt for trigger clusters from pi0 MC truth
  TH2F                 *fM02EtPi0MassClClTruPiC;  //!M02xEt for trigger clusters from pi+- MC truth
  TH2F                 *fM02EtPi0MassClClTruEta;  //!M02xEt for trigger clusters from eta MC truth
  TH2F                 *fM02EtPi0MassClClTruK_0;  //!M02xEt for trigger clusters from k0 MC truth
  TH2F                 *fM02EtPi0MassClClTruK_C;  //!M02xEt for trigger clusters from k+- MC truth
  TH2F                 *fM02EtPi0MassClClTruPro;  //!M02xEt for trigger clusters from proton MC truth
  TH2F                 *fM02EtPi0MassClClTruNeu;  //!M02xEt for trigger clusters from neutron MC truth
  TH2F                 *fM02EtTM;                 //!M02xEt for trigger clusters with track matched
  TH2F                 *fM02EtExot;               //!M02xEt for trigger clusters of exotic
  TH2F                 *fM02EtExotTM;             //!M02xEt for trigger TM clusters of exotic
   
  AliAnalysisTaskShoShaTests(const AliAnalysisTaskShoShaTests&); // not implemented
  AliAnalysisTaskShoShaTests& operator=(const AliAnalysisTaskShoShaTests&); // not implemented
  
  ClassDef(AliAnalysisTaskShoShaTests, 1); // Trigger contamination class
};
#endif
