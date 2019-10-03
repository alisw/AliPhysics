#ifndef AliAnalysisTaskNeutralMesons_cxx
#define AliAnalysisTaskNeutralMesons_cxx

#include "AliAnalysisTaskSE.h"

class TList;
class TF1;
class TH1D;
class TH2D;
class TFile;
class TObjArray;
class TArrayD;
class TString;
class TLorentzVector;
class TClonesArray;
class TParticle;
class AliMCEvent;
class AliStack;
class AliESDEvent;
class AliVEvent;
class AliVCluster;
class AliEMCALGeometry;
class AliEMCALRecoUtils;
class AliAODEvent;
class AliAODCaloCluster;
class AliOADBContainer;
class AliGenEventHeader;


class AliAnalysisTaskNeutralMesons : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskNeutralMesons();
  AliAnalysisTaskNeutralMesons(const char *name);
  virtual ~AliAnalysisTaskNeutralMesons();

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
  void           SetMCtype(Int_t mct)  {fMCtype=mct;}
  void           SetMCprt(Int_t mcprt)  {fMCpart=mcprt;}
  void           SetGeoName(TString geoname) {fGeoName=geoname;}
  void           SetLowEnergyCut(Double_t enrg) {fLowEnergCut=enrg;}
  void           SetNcellsPerClusterCut(Bool_t ncellFlag) {fNcellsInClusterFlag=ncellFlag;}
  void           SetNewClBranchName(TString branch) {fBranchName=branch;}
     
private:
  Bool_t       EsdVertexOk() const;
  Bool_t       IsGoodCluster(AliVCluster *cluster);
  void         MakeClusterCorrections(AliVCluster *cluster);
  void         MakeInvMass(const TLorentzVector& p1, const TLorentzVector& p2);
  void         MakeInvMassMixed(const TLorentzVector& pa, const TLorentzVector& pb);
  virtual Bool_t ComparePtHardAndJetPt() ;
  virtual Bool_t ComparePtHardAndClusterPt() ;
  virtual AliGenEventHeader* GetGenEventHeader() const ; 
  
  AliVEvent           *fESD;                  //! event object
  AliMCEvent          *fMC;                   // ! mc obj
  AliEMCALGeometry    *fGeom;                 //! EMCal geometry utility
  AliEMCALRecoUtils   *fRecoUtil;             //! Reco utility
  TList               *fOutputList;           //! Output list
  Bool_t              fRejectExoticCluster;   //
  Bool_t              fRemoveBadChannel;      //
  TH1D                *fhEclus;                //! cluster energy
  TH1D                *fhEclusSM0;               //! cluster energy sm0
  TH1D                *fhEclusSM1;               //! cluster energy sm1
  TH1D                *fhEclusSM2;               //! cluster energy sm2
  TH1D                *fhEclusSM3;               //! cluster energy sm3
  TH1D                *fhEclusSM4;               //! cluster energy sm4
  TH1D                *fhEclusSM5;               //! cluster energy sm5
  TH1D                *fhEclusSM6;               //! cluster energy sm6
  TH1D                *fhEclusSM7;               //! cluster energy sm7
  TH1D                *fhEclusSM8;               //! cluster energy sm8
  TH1D                *fhEclusSM9;               //! cluster energy sm9
  TH2D                *fhEVsTime;              //! cluster energy vs time
  TH2D                *fhEVsNcells;           //! cluster energy vs n cells
  TH1D                *fhZvertex;             //! event vertex distribution
  TFile               *f;                     //! file with trigger maps
  TF1                 *f1;                    // !first time cut
  TF1                 *f2;                    // !second time cut
  TF1                 *f3;                    // !non Linearity
  TH2D                *fPi0;                  //!pi0 mass
  TH2D                *fEta;                  //!eta mass
  TH1D                *fhZvertexAll;             //! event vertex distribution all
  TClonesArray        *fclusterList; // new clusters
  AliAODEvent         *fAOD;    // AOD object
  TParticle           *fparticle;    //! mc particle
  TH2D                *fhpi0Decay; //! sim
  TH2D                *fhpi0Inter;  //! sim
  TH2D                *fhRecDecay; //! sim
  TH2D                *fhRecInter; //! ssim
  TH1D                *fPrim; //! primaries
  TH1D                *fPrimNoTRD; //! primaries 
  TH1D                *fK0; //! K p
  TH2D                *fPi0mixed;                  //!pi0 mixed mass
  TH2D                *fEtamixed;                  //!eta mixed mass
  TObjArray           *fPool; //!pool of good clusters  
  TArrayD             *fClusterZvtx;     //!  for mixing
  TArrayD             *fNEvClusters;     //!  for mixing
  TH2D                *fhHitMap; //!hit map
  TH1D                *fNcl;           //! in list
  Int_t                fNeventCls;
  Int_t                fTrpi0;
  TH1D                *fK0Pt;           //! in list
  TFile               *fFileEff;           //! in list
  TH1D                *fhEff;           //! in list
  Double_t             fWeight;
  TH2D                *fhLowm;                  //! mass
  TH2D                *fhHighm;                  //! mass
  TH2D                *fPi0Conv;                  //!  mass
  TH2D                *fhLowmConv;                  //!  mass
  TH2D                *fhHighmConv;                  //!  mass
  TH2D                *fhPhotonEtaPhi ;             //! photon maps
  Int_t                fMCtype;   //1:Pythia, 2:jet=jet
  Int_t                fMCpart; //111 pi0, 221 eta
  TH1D                *fPrimOldGeo; //! primaries pointing to the first 4 SMs
  TString              fGeoName;  // geometry name
  Double_t             fLowEnergCut;  // low energy cluster cut
  Bool_t               fNcellsInClusterFlag;  // number of cells per cluster
  TString              fBranchName; // branch with the new clusters
    
  AliAnalysisTaskNeutralMesons(const AliAnalysisTaskNeutralMesons&); // not implemented
  AliAnalysisTaskNeutralMesons& operator=(const AliAnalysisTaskNeutralMesons&); // not implemented
  
  ClassDef(AliAnalysisTaskNeutralMesons, 1);
};

#endif
