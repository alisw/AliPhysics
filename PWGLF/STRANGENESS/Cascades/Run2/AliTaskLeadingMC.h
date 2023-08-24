#ifndef AliTaskLeadingMC_h
#define AliTaskLeadingMC_h

class AliESDEvent;
class TTree;
class TDatabasePDG;
class AliMCEvent;
class TH3F;
class TH2F;

#include "AliAnalysisTask.h"
#include "AliAnalysisFilter.h"

class AliTaskLeadingMC : public AliAnalysisTask {
 public:
  AliTaskLeadingMC(const char *name = "AliTaskLeadingMC"); 
  AliTaskLeadingMC(const AliTaskLeadingMC& source):AliAnalysisTask(source) {}
  virtual ~AliTaskLeadingMC() {}
  AliTaskLeadingMC& operator=(const AliTaskLeadingMC& source) {}
  
  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);

  void SetEtaThreshold(float val) { fEtaThreshold = val; }
  float GetEtaThreshold() const { return fEtaThreshold; }
  void SetEnergyThreshold(float val) { fEnergyThreshold = val; }
  float GetEnergyThreshold() const { return fEnergyThreshold; }
  void SetEtaBarrel(float val) { fEtaBarrel = val; }
  float GetEtaBarrel() const { return fEtaBarrel; }
  
  void fillZDCreco();
  void loopMC(AliMCEvent *mcEvent);
  void loopTrackRef(TTree *treeTR, AliMCEvent *mcEvent);
  void loopTrack(AliMCEvent *mcEvent);
  void  SetTrackCuts(AliAnalysisFilter* fTrackFilter);

  void SetZDCPGeo(float xmin=9.,float xmax=27.,float ymin=-7.,float ymax=7.,float zmin=10000.,float zmax=13000.);
  void SetZDCNGeo(float xmin=-4.,float xmax=4.,float ymin=-4.,float ymax=4.,float zmin=11000.,float zmax=12500.);
  
  void AskTrackRef(bool value=true) { fAskTrackRef = value; }


private: 
  // Notation
  // 1=(C)lockwise - 2=(A)nticlockwise

  static const int fgkDim = 50; // max array dimension
  Double_t ComputeSpherocity();
  AliAnalysisFilter   *fTrackFilter;       //!<! track filter for spherocity estimator 
  Float_t             fSpherocity;        ///< stores value of spherocity
  Int_t fNTracksSpherocity; //! number of tracks to calculate spherocity


  bool fAskTrackRef = false;

  AliESDEvent *fESD = nullptr; //!  ESD event
  TTree *fTree = nullptr;      //!
  TDatabasePDG *fDB = nullptr; //!

  TH3F * fH_SPD_VZERO; //!
  TH3F * fH_SPD_ZDC; //!
  TH2F * fH_SPD_VZERO_ev; //!
  TH2F * fH_SPD_ZDC_ev; //!

  float fV0Perc = 0.;
  float	fZdcPerc = 0.;
  float fZdcPercFired = 0.;
  float fMultRef5 = 0.;
  float fMultRef8 = 0.;
  float fMultSPDcl = 0.;
  float fMultSPDtr = 0.;
  int   fInelGT0 = 0;
  int   fSPDtracklets = 0;
  int   fSPDtrackletsA = 0;
  int   fSPDtrackletsC = 0;
  int   fTOFclusters = 0;
  int   fTOFclustersTrg = 0;

  int fIsTrackRef = 0;

// ZDC GEO PARAMETER (to select hit from TrackRef)
  float fZDCN_Zmin = 0.;
  float fZDCN_Zmax = 0.;
  float fZDCN_Ymin = 0.;
  float fZDCN_Ymax = 0.;
  float fZDCN_Xmin = 0.;
  float fZDCN_Xmax = 0.;

  float fZDCP_Zmin = 0.;
  float fZDCP_Zmax = 0.;
  float fZDCP_Ymin = 0.;
  float fZDCP_Ymax = 0.;
  float fZDCP_Xmin = 0.;
  float fZDCP_Xmax = 0.;

  // particles with eta > threshold (loop of MC stack)
  float fEtaThreshold = 8.;
  float fEnergyThreshold = 0.;
  float fEtaBarrel = 0.5;
  
  Int_t fP_cand_leadA = 0;         //! N charged candidate leading in A side
  Int_t fP_cand_leadC = 0;         //! N charged candidate leading in C side
  Int_t fN_cand_leadA = 0;         //! N neutral candidate leading in A side
  Int_t fN_cand_leadC = 0;         //! N neutral candidate leading in C side
  Int_t fPdg_cand_leadP2[fgkDim];  //! pdg of candidates
  Int_t fPdg_cand_leadP1[fgkDim];  //!
  Int_t fPdg_cand_leadN2[fgkDim];  //!
  Int_t fPdg_cand_leadN1[fgkDim];  //!
  Int_t fLabel_cand_P1[fgkDim];    //! label in the stack
  Int_t fLabel_cand_P2[fgkDim];    //!
  Int_t fLabel_cand_N1[fgkDim];    //!
  Int_t fLabel_cand_N2[fgkDim];    //!
  Int_t fPdgM_cand_leadP2[fgkDim]; // mother pdg of candidates
  Int_t fPdgM_cand_leadP1[fgkDim]; //!
  Int_t fPdgM_cand_leadN2[fgkDim]; //!
  Int_t fPdgM_cand_leadN1[fgkDim]; //!
  Int_t fLabelM_cand_P1[fgkDim];   //! label of mother in the stack
  Int_t fLabelM_cand_P2[fgkDim];   //!
  Int_t fLabelM_cand_N1[fgkDim];   //!
  Int_t fLabelM_cand_N2[fgkDim];   //!
  Float_t fE_p_cand_leadA[fgkDim]; //! candidate energy
  Float_t fE_p_cand_leadC[fgkDim];//!
  Float_t fE_n_cand_leadA[fgkDim]; //!
  Float_t fE_n_cand_leadC[fgkDim]; //!
  Float_t fMomCP1x[fgkDim];        //! candidate momentum
  Float_t fMomCP2x[fgkDim];        //!
  Float_t fMomCN1x[fgkDim];        //!
  Float_t fMomCN2x[fgkDim];        //!
  Float_t fMomCP1y[fgkDim];        //!
  Float_t fMomCP2y[fgkDim];        //!
  Float_t fMomCN1y[fgkDim];        //!
  Float_t fMomCN2y[fgkDim];        //!
  Float_t fMomCP1z[fgkDim];        //!
  Float_t fMomCP2z[fgkDim];        //!
  Float_t fMomCN1z[fgkDim];        //!
  Float_t fMomCN2z[fgkDim];        //!

  // ZDC recon infos
  Float_t fAdcZDCN1[5]; //! 0 common, 1-4 single towers
  Float_t fAdcZDCN2[5]; //!
  Float_t fAdcZDCP1[5]; //!
  Float_t fAdcZDCP2[5]; //!

  // MC particles with hits in the ZDC
  Int_t fP1hits = 0;            //! n hits in ZDCP C side
  Int_t fP2hits = 0;            //! n hits in ZDCP A side
  Int_t fN1hits = 0;            //! n hits in ZDCN C side
  Int_t fN2hits = 0;            //! n hits in ZDCN A side
  Float_t fE_p_leadA[fgkDim];   //! energy
  Float_t fE_p_leadC[fgkDim];   //!
  Float_t fE_n_leadA[fgkDim];   //!
  Float_t fE_n_leadC[fgkDim];   //!
  Int_t fPdg_leadingN1[fgkDim]; //! pdg code
  Int_t fPdg_leadingN2[fgkDim]; //!
  Int_t fPdg_leadingP1[fgkDim]; //!
  Int_t fPdg_leadingP2[fgkDim]; //!
  Int_t fLabel_leadP1[fgkDim];  //! label in the stack
  Int_t fLabel_leadP2[fgkDim];  //!
  Int_t fLabel_leadN1[fgkDim];  //!
  Int_t fLabel_leadN2[fgkDim];  //!
  Float_t fXP1[fgkDim];         //! coordinates of hits
  Float_t fXP2[fgkDim];         //!
  Float_t fXN1[fgkDim];         //!
  Float_t fXN2[fgkDim];         //!
  Float_t fYP1[fgkDim];         //!
  Float_t fYP2[fgkDim];         //!
  Float_t fYN1[fgkDim];         //!
  Float_t fYN2[fgkDim];         //!
  Float_t fZP1[fgkDim];         //!
  Float_t fZP2[fgkDim];         //!
  Float_t fZN1[fgkDim];         //!
  Float_t fZN2[fgkDim];         //!
  Float_t fMomP1x[fgkDim];      //! hit particle momentum
  Float_t fMomP2x[fgkDim];      //!
  Float_t fMomN1x[fgkDim];      //!
  Float_t fMomN2x[fgkDim];      //!
  Float_t fMomP1y[fgkDim];      //!
  Float_t fMomP2y[fgkDim];      //!
  Float_t fMomN1y[fgkDim];      //!
  Float_t fMomN2y[fgkDim];      //!
  Float_t fMomP1z[fgkDim];      //!
  Float_t fMomP2z[fgkDim];      //!
  Float_t fMomN1z[fgkDim];      //!
  Float_t fMomN2z[fgkDim];      //!
  Int_t fLabelMP1[fgkDim];      //! label mother
  Int_t fLabelMP2[fgkDim];      //!
  Int_t fLabelMN1[fgkDim];      //!
  Int_t fLabelMN2[fgkDim];      //!


  // other observables
  Int_t fNch = 0;            //! total charged multiplicity
  Int_t fNchEta = 0;         //! multiplicity in central region
  Int_t fNchEtaA = 0;        //! multiplicity in central region
  Int_t fNchEtaC = 0;        //! multiplicity in central region
  Float_t fEnergyEta = 0;    //! energy in central region
  Int_t fNMPI = 0;           //! number of multiparton interaction
  Int_t fNLambdaEta = 0;     //!
  Int_t fNXiEta = 0;         //!
  float fPtXiEta[100];       //!
  Int_t fNAntiXiEta = 0;     //!
  float fPtAntiXiEta[100];   //!
  Int_t fNOmegaEta = 0;      //!
  Int_t fNPiEta = 0;         //!
  Int_t fNPi0Eta = 0;        //!
  Int_t fNKchEta = 0;        //!
  Int_t fNK0Eta = 0;         //!
  float fSumPtLambdaEta = 0; //!
  float fSumPtXiEta = 0;     //!
  float fSumPtOmegaEta = 0;  //!
  float fSumPtPiEta = 0;     //!
  float fMaxChargePt = 0;    //!
  float fEffEnergy = 0;      //!
  Int_t fNXiEtaFrag = 0;     //!
  Int_t fNXiEtaUp = 0;       //!
  Int_t fNXiEtaDown = 0;     //!
  float fPtXiEtaFrag[100];   //!
  float fPtXiEtaUp[100];     //!
  float fPtXiEtaDown[100];   //!

  // temporary variables
  Float_t fPt = 0;    //!
  Float_t fEta = 0;   //!
  Float_t fE = 0;     //!
  Float_t fTheta = 0; //!
  Float_t fPdg = 0;   //!
  Float_t fSC = 0;    //!
  Float_t fPx = 0;    //!
  Float_t fPy = 0;    //!
  Float_t fPz = 0;    //!
  Float_t fTRpz = 0;  //!
  Float_t fX = 0;     //!
  Float_t fY = 0;     //!
  Float_t fZ = 0;     //! impact point (cm) in the ZDC
  Int_t fLabel = 0;   //!
	
  ClassDef(AliTaskLeadingMC, 2); // example of analysis
};

#endif
