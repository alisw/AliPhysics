#ifndef AliAnalysisTaskEMCALPi0V2ShSh_h
#define AliAnalysisTaskEMCALPi0V2ShSh_h

// $Id: AliAnalysisTaskEMCALPi0V2ShSh.h$

class TH1F;
class TH1D;
class TH2F;
class THnSparse;
class TList;
class TObjArray;
class AliOADBContainer;
class AliEMCALGeometry;
class AliESDEvent;
class AliESDtrack;
class AliESDCaloCells;
class AliAODEvent;
class AliAODCaloCells;
class AliVCluster;
class AliCentrality;

#include "AliAnalysisTaskSE.h"


class AliAnalysisTaskEMCALPi0V2ShSh : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskEMCALPi0V2ShSh();
  AliAnalysisTaskEMCALPi0V2ShSh(const char *name);
  virtual ~AliAnalysisTaskEMCALPi0V2ShSh() {}

  void                   UserCreateOutputObjects();
  void                   UserExec(Option_t *option);
  void                   VZEROEventPlane();
  void                   FillClusterHists();
  void                   FillTrackHists();
  void                   Terminate(Option_t *);

 protected:
  AliEventplane         *fEventPlane;
  Double_t               fCentralityV0M;
  TObjArray             *fESDClusters;           //!pointer to EMCal clusters
  TObjArray             *fAODClusters;           //!pointer to EMCal clusters
  AliESDCaloCells       *fESDCells;              //!pointer to EMCal cells, esd
  AliAODCaloCells       *fAODCells;              //!pointer to EMCal cells, aod  
  AliEMCALGeometry      *fGeom;                  // geometry utils
  TString                fGeoName;               // geometry name (def = EMCAL_FIRSTYEARV1)
  AliOADBContainer      *fOADBContainer;         //!OADB container used to load misalignment matrices
  
  
 private:
  AliESDEvent *fESD;                 //! ESD object
  AliAODEvent *fAOD;                 //! AOD object
  TList       *fOutputList;          //! General Output list
  TGeoHMatrix *fGeomMatrix[12];      //! Geometry misalignment matrices for EMCal

  Double_t    fEPTPC;
  Double_t    fEPTPCResolution;
  Double_t    fEPV0;
  Double_t    fEPV0A;
  Double_t    fEPV0C;
  Double_t    fEPV0Ar;
  Double_t    fEPV0Cr;
  Double_t    fEPV0r;
  Double_t    fEPV0A4r;
  Double_t    fEPV0A5r;
  Double_t    fEPV0A6r;
  Double_t    fEPV0A7r;
  Double_t    fEPV0C0r;
  Double_t    fEPV0C1r;
  Double_t    fEPV0C2r;
  Double_t    fEPV0C3r;

  //histograms
  TH1F        *fHistAllcentV0;
  TH1F        *fHistAllcentV0r;
  TH1F        *fHistAllcentV0A;
  TH1F        *fHistAllcentV0C;
  TH1F        *fHistAllcentTPC;

  TH2F        *fHistEPTPC;
  TH2F        *fHistEPTPCResolution;

  TH2F        *fHistEPV0;
  TH2F        *fHistEPV0A;
  TH2F        *fHistEPV0C;
  TH2F        *fHistEPV0Ar;
  TH2F        *fHistEPV0Cr;
  TH2F        *fHistEPV0r;
  TH2F        *fHistEPV0A4r;
  TH2F        *fHistEPV0A7r;
  TH2F        *fHistEPV0C0r;
  TH2F        *fHistEPV0C3r;

  TH2F        *fHistdifV0A_V0C0r;
  TH2F        *fHistdifV0A_V0C3r;
  TH2F        *fHistdifV0C0r_V0C3r;
  TH2F        *fHistdifV0C_V0A4r;
  TH2F        *fHistdifV0C_V0A7r;
  TH2F        *fHistdifV0A4r_V0A7r;
  TH2F        *fHistdifV0Ar_V0Cr;	

  TH1F        *fHistClusterEta;
  TH1F        *fHistClusterPhi;
  TH1F        *fHistClusterE;
  TH1F        *fHistClusterEt;
  TH1F        *fHistClusterN;
  TH1F        *fHistClusterM02;
  TH2F        *fHistClusterEN;
  TH2F        *fHistClusterEM02;
  TH2F        *fHistClusterPhiEta;
  TH2F	      *fHistClusterEtN;
  TH2F        *fHistClusterEtM02;
  TH1D        *fHistClusterdphiV0;

  TH1F        *fHistTrackPt;
  TH1F        *fHistTrackEta;
  TH1F        *fHistTrackPhi;
  TH2F        *fHistTrackPhiEta;

  THnSparse   *fClusterPbV0;
  THnSparse   *fClusterPbV0A;
  THnSparse   *fClusterPbV0C;
  THnSparse   *fClusterPbTPC;
 

  AliAnalysisTaskEMCALPi0V2ShSh(const AliAnalysisTaskEMCALPi0V2ShSh&); // not implemented
  AliAnalysisTaskEMCALPi0V2ShSh& operator=(const AliAnalysisTaskEMCALPi0V2ShSh&); // not implemented
  
  ClassDef(AliAnalysisTaskEMCALPi0V2ShSh, 1);
};
#endif
