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
  void SetAbsEtaLimit(Double_t etaMin=0., Double_t etaMax=0.8) {EtaLimit[0]=etaMin;EtaLimit[1]=etaMax;}
  void SetDCACut(Double_t DCAxyCUT=0.1, Double_t DCAzCUT=1000.0) {DCAxyCut=DCAxyCUT; DCAzCut=DCAzCUT;}
   //other cuts 
  void SetNsigmaTPCCut(Double_t nSigmaTpcCut=2) {NsigmaTpcCut=nSigmaTpcCut;}
  void SetNminTPCcluster(Int_t nMinTPCcluster=0) {NminTpcCluster=nMinTPCcluster;}
  void SetTrdCut(Int_t kTRDcut=0) {iTrdCut=kTRDcut;}
  //Settings
  void SetisSignalCheck(Int_t IsSignalCheck=2) {kSignalCheck=IsSignalCheck;}
  void SetMtofMethod(Int_t iMtofMethod=1) {iMtof=iMtofMethod;}

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
  Int_t iMtof;                                     // iMtof==1->m~pVtx ; iMtof==2->m~pExp ; iMtof==4->m~pExp(MCcorrected) for (d,He2); iMtof==4->m~pExp(MCcorrected) (p,d,He3)
  
  //other:
  Int_t iBconf;                                   //! If Magnetic Field configuration is down or up
  Bool_t kTOF;                                    //! kTOFout and kTIME required
  
  AliAODEvent* fAOD;                              //! AOD object
  AliESDEvent* fESD;                              //! ESD object
  AliVEvent* fEvent;                              //! general object
  AliPIDResponse *fPIDResponse;                   //! pointer to PID response
  TList *fList[nBconf];                           //! lists for slot
  
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
  TH2F *fNsigmaTpc[nBconf][9];                    //! NsigmaTPC vs. pTpc
  TH2F *fNsigmaTpc_kTOF[nBconf][18];              //! NsigmaTPC vs. pt when kTOF is required and in DCAxyCut 
  
  //TOF info:
  TH2F *fBetaTofVSp[nBconf][2];                   //! beta vs pVtx
  TProfile *hBetaExp[nBconf][9];                  //! TOF expected beta
  TH2F *fNsigmaTof[nBconf][9];                    //! NsigmaTOF vs. pT
  TH2F *fNsigmaTof_DcaCut[nBconf][18];                    //! NsigmaTOF vs. pT

  //TPC and TOF conbined
  TH2F *fM2vsPt_NoTpcCut[nBconf][2][2];           //! M2 vs. Pt w/o the DCAxyCut
  TH2F *fM2vsPt[nBconf][2][18];                   //! M2 vs. Pt with NsigmaTpcCut for each particle species, w/o the DCAxyCut
  TH2F *fM2vsZ[nBconf][10];                       //! M2 vs. Z in various pT bins

  //DCA distributions
  TH1D *hDCAxy[nBconf][18][nbin];                 //! DCAxy distribution with NsigmaTpcCut for each particle species, in pT bins
  TH1D *hDCAz[nBconf][18][nbin];                  //! DCAz distribution with NsigmaTpcCut for each particle species, in pT bins
  
  //TOF mass distributions
  TH1D *hM2CutDCAxy[nBconf][18][nbin];            //! Tof m2 distribution in DCAxyCut and with NsigmaTpcCut
  TH1D *hM2CutGroundDCAxy[nBconf][18][nbin];      //! Tof m2 distribution in the background of DCAxyCut (secondary nuclei selection) and with NsigmaTpcCut
 
  //Parametrizations
  TF1 *fPmeanVsPexp[3];                           //! Parameterization of (<p>-pExp)/pExp vs pExp for p,d,He3

  //------------------------------Methods----------------------------------------
  void GetMassFromPvertex(Double_t beta, Double_t p, Double_t &M2);
  void GetZTpc(Double_t dedx, Double_t pTPC, Double_t M2, Double_t &Z2);
  void GetMassFromExpTimes(Double_t beta, Double_t *IntTimes, Double_t *Mass2, Int_t iCorr=4);
 
  ClassDef(AliAnalysisNucleiMass, 1);
};

#endif
