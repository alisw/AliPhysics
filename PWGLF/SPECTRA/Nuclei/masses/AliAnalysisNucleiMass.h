#ifndef ALIANALYSISNUCLEIMASS_H
#define ALIANALYSISNUCLEIMASS_H

// ROOT includes
#include <TList.h>

// AliRoot includes
#include <AliAnalysisTaskSE.h>
#include <AliPIDResponse.h>

class AliAODEvent;
class AliESDEvent;
class AliVEvent;
class TH2F;
class TH2D;
class TH1F;
class TF1;
class TF2;
class TH2D;
class TGraph;
class AliESDtrackCuts;
class TProfile;
class TFile;
class TObject;

class AliAnalysisNucleiMass : public AliAnalysisTaskSE {
 public:
  AliAnalysisNucleiMass();
  AliAnalysisNucleiMass(const char *name);
  
  virtual ~AliAnalysisNucleiMass();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);


  //Cuts on the events
  void SetCentrality(Double_t CentMin=0., Double_t CentMax=100.) {Centrality[0]=CentMin; Centrality[1]=CentMax;};
  //Cuts on the tracks
  void SetFilterBit(Int_t TestFilterBit=16) {FilterBit=TestFilterBit;}
   //geometrical cuts
  void SetEtaLimit(Double_t etaMin=-0.8, Double_t etaMax=0.8) {EtaLimit[0]=etaMin;EtaLimit[1]=etaMax;}
  void SetDCACut(Double_t DCAxyCUT=0.1, Double_t DCAzCUT=1000.0) {DCAxyCut=DCAxyCUT; DCAzCut=DCAzCUT;}
   //other cuts 
  void SetNsigmaTPCCut(Double_t nSigmaTpcCut=2) {NsigmaTpcCut=nSigmaTpcCut;}
  void SetNminTPCcluster(Int_t nMinTPCcluster=0) {NminTpcCluster=nMinTPCcluster;}
  void SetTrdCut(Int_t kTRDcut=0) {iTrdCut=kTRDcut;}
  
  //Settings
  void SetisSignalCheck(Int_t IsSignalCheck=2) {kSignalCheck=IsSignalCheck;}
  void SetMtofMethod(Int_t iMtofMethod=1) {iMtof=iMtofMethod;}
  void SetPvtxNucleiCorrection(Int_t kMomVtxCorr=1) {kPvtxCorr=kMomVtxCorr;}
  //void SetiTriggerSel(Int_t ItriggerSel=-99) {iTriggerSel=ItriggerSel;}

 private:
  AliAnalysisNucleiMass(const AliAnalysisNucleiMass &old); 
  AliAnalysisNucleiMass& operator=(const AliAnalysisNucleiMass &source);
    
  static const Int_t nbin=46;                      // Number of pt bins in Tof Mass distributions
  static const Int_t nBconf=2;                     // Number of Magnetic field configuration (B++ and B--)
  static const Int_t nPart=9;                      // Number of particle type: e,mu,pi,K...
  static const Int_t nSpec=18;                     // Number of particle species: particles: e+,e-,mu+,mu-,...
    
  //Variables settings with public methods:
  Double_t Centrality[2];                          // Centrality bin (min and max)
  Int_t FilterBit;                                 // Filter Bit to be used
  Double_t EtaLimit[2];                            // Eta windows in analysis
  Double_t DCAxyCut;                               // Cut on DCA-xy
  Double_t DCAzCut;                                // Cut on DCA-z
  Double_t NsigmaTpcCut;                           // number of sigma Tpc Cut
  Int_t NminTpcCluster;                            // Number of minimum TPC clusters
  Int_t iTrdCut;                                   // iTrdCut==0-> No TRD cut; iTrdCut==1-> Yes TRD cut: yes TRD; iTrdCut==2->Yes TRD cut: no TRD; 
  Int_t kSignalCheck;                              // kSignalCheck==1->Fill all plots ; kSignalCheck==0->Fill only TH1 ; kSignalCheck==2-> Fill TH1 and some TH2 usefull in analysis
  Int_t iMtof;                                     // iMtof==1->m~pVtx ; iMtof==2->m~pExp ; iMtof==4->m~<p> (same correction for particle and antiparticle) ;
  Int_t kPvtxCorr;                                 // kPvtxCorr==1->Momentum at the primary vertex for (anti)nuclei is rescaled ; kPvtxCorr==0->no correction
  
  //other:
  Int_t iBconf;                                   //! If Magnetic Field configuration is down or up
  Bool_t kTOF;                                    //! kTOFout and kTIME required
  
  static const Int_t iTriggerSel=-99;                   // -99->no trigger required ; 0-> if kMB ; 16-> if kCentral ; 17-> if kSemiCentral ; -2 -> No MB, No Central and No SemiCentral  

  AliAODEvent* fAOD;                              //! AOD object
  AliESDEvent* fESD;                              //! ESD object
  AliVEvent* fEvent;                              //! general object
  AliPIDResponse *fPIDResponse;                   //! pointer to PID response
  TList *fList[nBconf];                           //! lists for slot
  
  TH1I *htriggerbits[nBconf][2];                  //! Trigger bits distribution
  TH1F *htemp[nBconf];                            //! Temp. plot: avoid a problem with the merge of the output when a TList is empty (of the opposite magnetic field configuration)
  TH1F *hCentrality[nBconf][2];                   //! Centrality of the selected and analyzed events
  TH1F *hZvertex[nBconf][2];                      //! z-vertex distribution before and after the cuts on the event
  
  TH1F *hEta[nBconf];                             //! Eta distribution of the tracks
  TH1F *hPhi[nBconf];                             //! Phi particle distribution
  TH2F *fEtaPhi[nBconf];                          //! Phi vs Eta particle distribution
  TH1F *hNTpcCluster[nBconf];                     //! # of the TPC clusters after the track cuts
  TH1F *hNTrdSlices[nBconf];                      //! Number of the TRD slices after the track cuts

  //TPC info:
  TH2F *fdEdxVSp[nBconf][2];                      //! dedx vs pTpc
  TProfile *hDeDxExp[nBconf][9];                  //! TPC spline used
  TH2F *fNsigmaTpc[nBconf][18];                    //! NsigmaTPC vs. pTpc
  TH2F *fNsigmaTpc_kTOF[nBconf][18];              //! NsigmaTPC vs. pt when kTOF is required 
  
  //TOF info:
  TH2F *fBetaTofVSp[nBconf][2];                   //! beta vs pVtx
  TProfile *hBetaExp[nBconf][9];                  //! TOF expected beta
  TH2F *fNsigmaTof[nBconf][9];                    //! NsigmaTOF vs. pT
  
  //TPC and TOF conbined
  TH2F *fM2vsP_NoTpcCut[nBconf][1][2];            //! M2 vs. P
  TH2F *fM2vsP[nBconf][1][18];                    //! M2 vs. P with NsigmaTpcCut for each particle species
  TH2F *fM2vsZ[nBconf][2];                        //! M2 vs. Z

  //DCA distributions
  TH1D *hDCAxy[nBconf][18][nbin];                 //! DCAxy distribution with NsigmaTpcCut for each particle species, in p bins
  TH1D *hDCAz[nBconf][18][nbin];                  //! DCAz distribution with NsigmaTpcCut for each particle species, in p bins
  //TH2F *h2DCA[nBconf][18][nbin];                //! DCAxy vs DCAz with NsigmaTpcCut for each particle species, in p bins
  TH2F *h2DCAap[nBconf][18];                      //! DCAxy vs DCAz with NsigmaTpcCut for each particle species
    
  //TOF mass distributions
  TH1D *hM2CutDCAxy[nBconf][18][nbin];            //! Tof m2 distribution with NsigmaTpcCut
  
  //...
  TH2F *fPmeanVsBetaGamma[nBconf][18];            //! <p>/p vs beta*gamma for pi,K,p
  TProfile *prPmeanVsBetaGamma[nBconf][18];       //! <p>/p vs beta*gamma for pi,K,p (profile)
  
  //Parameterizations:
  TF2 *fPvtxTrueVsReco[4];                        //! TF1 pVtx_True vs pVtx_Reco calculated with MC for d,t,He3,He4
  TProfile *prPvtxTrueVsReco[nBconf][4];          //! TProfile pVtx_True vs pVtx_Reco calculated with MC for d,t,He3,He4

  TF1 *fPmeanVsBGcorr[14];                        //! <p>/p as a function of beta*gamma for pi,K,p,d,t,He3,He4
  TProfile *prPmeanVsBGcorr[nBconf][14];          //! <p>/p vs beta*gamma for pi,K,p,d,t,He3,He4 as calculated from the parameterizations
 
  //------------------------------Methods----------------------------------------
  void MomVertexCorrection(Double_t p, Double_t *pC, Double_t eta, Int_t FlagPid);
  
  void FillDCAdist(Double_t DCAxy, Double_t DCAz, Double_t charge, Int_t FlagPid, Int_t stdFlagPid[9], Double_t *pC);
  
  void GetMassFromPvertex(Double_t beta, Double_t p, Double_t &M2);
  void GetZTpc(Double_t dedx, Double_t pTPC, Double_t M2, Double_t &Z2);

  void GetMassFromPvertexCorrected(Double_t beta, Double_t *pC, Double_t *Mass2);
  
  void GetMassFromExpTimes(Double_t beta, Double_t *IntTimes, Double_t *Mass2);
  void GetPmeanVsBetaGamma(Double_t *IntTimes, Double_t *pVtx, Int_t FlagPid, Int_t FlagPidTof, Double_t charge);

  void GetMassFromMeanMom(Double_t beta, Double_t *IntTimes, Double_t *pVtx, Double_t eta, Double_t charge, Double_t *Mass2, Int_t FlagPid, Int_t FlagPidTof);
  
  void SetPvtxCorrections();
  void SetPmeanCorrections();

  ClassDef(AliAnalysisNucleiMass, 4);
};

#endif
