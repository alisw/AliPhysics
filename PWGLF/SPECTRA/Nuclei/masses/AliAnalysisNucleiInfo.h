#ifndef ALIANALYSISNUCLEIINFO_H
#define ALIANALYSISNUCLEIINFO_H

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
class TH1I;
class TF1;
class TF2;
class TGraph;
class AliESDtrackCuts;
class TProfile;
class TFile;
class TObject;

class AliAnalysisNucleiInfo : public AliAnalysisTaskSE {
 public:
  AliAnalysisNucleiInfo();
  AliAnalysisNucleiInfo(const char *name);
  
  virtual ~AliAnalysisNucleiInfo();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);


  //Cuts on the tracks
  void SetFilterBit(Int_t TestFilterBit=16) {FilterBit=TestFilterBit;}
  //Geometrical cuts
  void SetEtaLimit(Double_t etaMin=-0.8, Double_t etaMax=0.8) {EtaLimit[0]=etaMin;EtaLimit[1]=etaMax;}
  void SetDCACut(Double_t DCAxyCUT=1000.0, Double_t DCAzCUT=1000.0) {DCAxyCut=DCAxyCUT; DCAzCut=DCAzCUT;}
  //Other cuts 
  void SetNsigmaTPCCut(Double_t nSigmaTpcCut=2) {NsigmaTpcCut=nSigmaTpcCut;}
  void SetStartTimeTofRes(Double_t startTimeTofRes=9999.9){StartTimeTofRes=startTimeTofRes;}

 private:
  AliAnalysisNucleiInfo(const AliAnalysisNucleiInfo &old); 
  AliAnalysisNucleiInfo& operator=(const AliAnalysisNucleiInfo &source);
    
  static const Int_t nBconf=2;                     // Number of Magnetic field configuration (B++ and B--)
  static const Int_t nPart=9;                      // Number of particle type: e,mu,pi,K...
  static const Int_t nSpec=18;                     // Number of particle species: particles: e+,e-,mu+,mu-,...
    
  //Variables settings with public methods:
  Int_t FilterBit;                                 // Filter Bit to be used
  Double_t EtaLimit[2];                            // Eta windows in analysis
  Double_t DCAxyCut;                               // Cut on DCA-xy
  Double_t DCAzCut;                                // Cut on DCA-z
  Double_t NsigmaTpcCut;                           // number of sigma Tpc Cut
  Double_t StartTimeTofRes;
  
  //other:
  Int_t iBconf;                                   //! If Magnetic Field configuration is down or up
  Bool_t kTOF;                                    //! kTOFout and kTIME required
  
  static const Int_t iTriggerSel=0;             // -99->no trigger required ; 0-> if kMB ; 16-> if kCentral ; 17-> if kSemiCentral ; -2 -> No MB, No Central and No SemiCentral  

  AliAODEvent* fAOD;                              //! AOD object
  AliESDEvent* fESD;                              //! ESD object
  AliVEvent* fEvent;                              //! general object
  AliPIDResponse *fPIDResponse;                   //! pointer to PID response
  TList *fList[nBconf];                           //! lists for slot
  
  TH1I *htriggerbits[nBconf][2];                  //! Trigger bits distribution
  TH1F *htemp[nBconf];                            //! Temp. plot: avoid a problem with the merge of the output when a TList is empty (of the opposite magnetic field configuration)
  TH1F *hZvertex[nBconf][2];                      //! z-vertex distribution before and after the cuts on the event
  
  TH1F *hEta[nBconf];                             //! Eta distribution of the tracks
  TH1F *hPhi[nBconf];                             //! Phi particle distribution
  TH1I *hNtrackAtTof[nBconf];                        //! Number of the tracks when kTOF is required
  
  //TPC info:
  TH2F *fdEdxVSp[nBconf][2];                      //! dedx vs pTpc
  TProfile *hDeDxExp[nBconf][9];                  //! TPC spline used
  TH2F *fNsigmaTpc[nBconf][18];                   //! NsigmaTPC vs. pTpc
  
  //TOF info:
  TH2F *fBetaTofVSp[nBconf][2];                   //! beta vs pVtx
  TProfile *hBetaExp[nBconf][9];                  //! TOF expected beta
  TH2F *fNsigmaTof[nBconf][2][18];                //! NsigmaTOF vs. pT
  TH2F *fTofMinusExp[nBconf][2][18];              //! tof-t_exp w/o tpc
  TH1F *hStartTimeRes[nBconf];                    //! start time resolution
 
  //ITS info:
  TH2F *h2DCAap[nBconf][18];                      //! DCAxy vs DCAz with NsigmaTpcCut for each particle species
  
  ClassDef(AliAnalysisNucleiInfo, 2);
};

#endif
