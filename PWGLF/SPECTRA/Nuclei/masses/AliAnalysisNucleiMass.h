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

  void SetCentrality(Float_t *fCt) {fCentrality[0]=fCt[0];fCentrality[1]=fCt[1];};
  void SetFilterBit(Int_t TestFilterBit) {FilterBit=TestFilterBit;}
  void SetNTPCcluster(Int_t nTPCcluster) {NminTPCcluster=nTPCcluster;}
  void SetDCAzCut(Float_t fDCAzCut) {DCAzCUT =fDCAzCut;}
  void SetDCAxyCut(Float_t fDCAxyCut) {DCAxyCUT=fDCAxyCut;}
  void SetkTPCcut(Bool_t isTPCcut) {kTPCcut=isTPCcut;}
  void SetNsigmaTPCCut(Float_t NsigmaTpcCut) {NsigmaTPCCut=NsigmaTpcCut;}
  void SetisSignalCheck(Bool_t IsSignalCheck) {isSignalCheck=IsSignalCheck;}
  void SetMomBin(Int_t iMomBin) {MomType=iMomBin;}
  void SetAbsEtaLimit(Double_t *etaLimit) {EtaLimit[0]=etaLimit[0];EtaLimit[1]=etaLimit[1];}
  void SetTRDanalysis(Bool_t kTrdAnalysis=kFALSE, Int_t iTrd=1) {kTRDana=kTrdAnalysis;iTRD=iTrd;}

 private:
  AliAnalysisNucleiMass(const AliAnalysisNucleiMass &old); 
  AliAnalysisNucleiMass& operator=(const AliAnalysisNucleiMass &source);
    
  TFile *fmism;                     //! For load the mism time distr
  TH1F *hmism;                      //! The mism time distr
    
  TFile *fchDist;                   //! Load the tof chan dist from IP
  TH1D *hChDist;                    //! The tof chan dist from IP
  
  static const Int_t nbin = 46;     // number of pt bins

  Double_t EtaLimit[2];                 // Eta windows in analysis

  Int_t MomType;                    // type of momentum bins in analysis (7 are all ON): (Flag: 001(1)->pT 010(2)->p 100(3)->pTPC)
 
  Bool_t kTRDana;                    //TRD analysis: 0->No 1->Yes
  
  Int_t iTRD;                        //TRD: 2->No TRD, 4->Yes TRD, 1->indifferent
 
  Bool_t fMC;                       // if MC

  Float_t fCentrality[2];           // centrality bin (min and max)
  
  Int_t FilterBit;                  // filter be to be used

  Int_t NminTPCcluster;             // min TPC cluster number

  Float_t DCAzCUT;                  // cut on DCA-z
  Float_t DCAxyCUT;                 // cut on DCA-xy

  Bool_t kTPCcut;                   // to apply a TPC 2 sigma cut

  Bool_t kTPC;                      //! is > NminTPCcluster 
  Bool_t kTOF;                      //! kTOFout and kTIME required

  Int_t iBconf;                      //! if Magnetic Configuration is down or up 

  Bool_t isSignalCheck;               // if write with an appropriate binning the plots of the various signals (QA,...) 

  Float_t NsigmaTPCCut;              // number of sigma Tpc Cut

  AliAODEvent* fAOD;                //! AOD object
  
  AliESDEvent* fESD;                //! ESD object
  
  AliVEvent* fEvent;                //! general object
 
  AliPIDResponse *fPIDResponse;     //! pointer to PID response

  TList *fList1[2];                    //! lists for slot

  TH1F *hNeventSelected[2];            //! selected Event counter  

  TH1F *hNevent[2];                    //! analyzed Event counter
  
  TH1F *hZvertex[2];                   //! z-vertex distribution

  TH1F *hEtaDistribution[2][2];          //! Eta distribution of the tracks

  TH2F *fEtaSpecies[2][18];           //! Eta distribution of the each particle identified by the TPC
  
  TH1F *hPhi[2][6];                   //! Phi particle distribution

  TH2F *fEtaPhi[2][6];                 //! Phi vs Eta particle distribution

  TH2F *fPhiSpecies[2][18];          //! Phi vs Eta particle distribution identified by the TPC

  TH2F *fdEdxVSp[2][3];                //! dedx vs p plots

  TH2F *fBetaTofVSp[2];                //! beta vs p plots

  TH1F *hTOFSignalPion[2];             //! pion  TOF signal

  TH2F *fM2vsP_NoTpcCut[2][3];         //! M2 vs. P
  
  TH2F *fNsigmaTPC[2][9];              //! NsigmaTPC vs. pT
  
  TH2F *fNsigmaTOF[2][9];              //! NsigmaTOF vs. pT

  TH2F *fNsigmaTPCvsP_kTOFtrue[2][18]; //! NsigmaTPC vs. p with kTOFout && kTIME for provide TPC different cuts effect
 
  TProfile *hDeDxExp[2][9];            //! TPC spline used

  TProfile *hBetaExp[2][9];            //! TOF expected beta
  
  TH2F *fM2vsZ[2][15];                 //! M2 vs. Z in different pT range

  TH2F *fM2vsZwithTPC[2][15];          //! M2 vs. Z in different pT range with 2sigmaTPC cut

  TH2F *fM2vsP[2][18];                 //! M2 vs. P with 2 sigma TPC cut for each particle species

  TH1D *hDCAxy[2][18][nbin];           //! DCA distribution in 2 sigma TPC cut for each particle species, in pT bins

  TH1D *hM2CutDCAxy[2][18][nbin];      //! M^{2} IN DCA cut (in 2 sigma TPC cut), in pT bins

  TH1D *hDCAz[2][18][nbin];            //! DCAz distribution in 2 sigma TPC cut for each particle species, in pT bins

  TH1D *hM2CutGroundDCAxy[2][18][nbin];//! M^{2} OUT DCA cut (in 2 sigma TPC cut), in pT bins

  TH2F *fM2vsP_NoTpcCut_DCAxyCut[2][3];//! M^{2} vs. P with a DCAxy cut  

  TH2F *fM2vsP_DCAxyCut[2][18];        //! M^{2} vs. P with a DCAxy cut (2sigma TPC cut)
  
  TH1D *hDCAxy_pbin[2][18][nbin];           //! DCA distribution in 2 sigma TPC cut for each particle species, in p bins

  TH1D *hM2CutDCAxy_pbin[2][18][nbin];      //! M^{2} IN DCA cut (in 2 sigma TPC cut), in p bins

  TH1D *hDCAz_pbin[2][18][nbin];            //! DCAz distribution in 2 sigma TPC cut for each particle species, in p bins

  TH1D *hM2CutGroundDCAxy_pbin[2][18][nbin];//! M^{2} OUT DCA cut (in 2 sigma TPC cut), in p bins
   
  TH1D *hDCAxy_pTpcbin[2][18][nbin];           //! DCA distribution in 2 sigma TPC cut for each particle species, in pTPC bins

  TH1D *hM2CutDCAxy_pTpcbin[2][18][nbin];      //! M^{2} IN DCA cut (in 2 sigma TPC cut), in pTPC bins

  TH1D *hDCAz_pTpcbin[2][18][nbin];            //! DCAz distribution in 2 sigma TPC cut for each particle species, in pTPC bins

  TH1D *hM2CutGroundDCAxy_pTpcbin[2][18][nbin];//! M^{2} OUT DCA cut (in 2 sigma TPC cut), in pTPC bins

  TH1D *hM2BkgMism[2][3][nbin];                //! M2 from mismatch background in each momentum bin

  TH1F *hNminTPCcl[2];

  ClassDef(AliAnalysisNucleiMass, 1);
};

#endif
