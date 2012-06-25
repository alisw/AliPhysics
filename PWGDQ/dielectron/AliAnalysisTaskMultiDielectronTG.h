#ifndef ALIANALYSISTASKMULTIDIELECTRONTG_H
#define ALIANALYSISTASKMULTIDIELECTRONTG_H
/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//#####################################################
//#                                                   # 
//#        Basic Analysis task for Dielectron         #
//#          single event analysis                    #
//#                                                   #
//#  by WooJin J. Park, GSI / W.J.Park@gsi.de         #
//#     Ionut C. Arsene, GSI / I.C.Arsene@gsi.de      #
//#     Magnus Mager, CERN / Magnus.Mager@cern.ch     #
//#     Jens Wiechula, Uni HD / Jens.Wiechula@cern.ch #
//#                                                   #
//#####################################################

#include "TList.h"

#include "AliAnalysisTaskSE.h"

#include <vector>
#include <deque>
#include <cstdlib>


// #include "AliDielectronPID.h"

class AliDielectron;
class TH1D;
class AliAnalysisCuts;
class AliTriggerAnalysis;
class AliESDtrackCuts;

class AliDielectronSingleTG : public TObject 
{
 public:
  AliDielectronSingleTG():
    fCharge(0),
    fCentrality(0),
    fXv(0),
    fYv(0),
    fZv(0),
    fPx(0),
    fPy(0),
    fPz(0),
    fPt(0),
    fEta(0),
    fPhi(0),
    fTheta(0),
    fConv(0),
    fGst(0),
    Obj(0x0)
      {;
      }
    
    AliDielectronSingleTG(Int_t charge, Double_t cent, 
			  Double_t xv, Double_t yv, Double_t zv, 
			  Double_t px, Double_t py, Double_t pz, Double_t pt,
			  Double_t eta, Double_t phi, Double_t theta,
			  Int_t conv, Int_t ghost, AliVTrack *trk)
      : 
      fCharge(charge), 
      fCentrality(cent), 
      fXv(xv),
      fYv(yv),
      fZv(zv),
      fPx(px),
      fPy(py),
      fPz(pz),
      fPt(pt),
      fEta(eta),
      fPhi(phi),
      fTheta(theta),
      fConv(conv),
      fGst(ghost), Obj(0x0)
	{
	  SetTrack(trk);
	  ;
	}

      ~AliDielectronSingleTG() {;}

      
      void SetTrack(AliVTrack *trk) { Obj = trk;}
      Int_t Charge(void){ return fCharge;}
      Double_t  Phi(void){ return fPhi;}
      Double_t  Theta(void){ return fTheta;}
      Double_t  Px(void){ return fPx;}
      Double_t  Py(void){ return fPy;}
      Double_t  Pz(void){ return fPz;}
      Double_t  Xv(void){ return fPx;}
      Double_t  Yv(void){ return fPy;}
      Double_t  Zv(void){ return fPz;}
      Double_t  Pt(void){ return fPt;}
      AliVTrack *GetTrack(void){ return Obj;}
      void SetConvFlag(Int_t val){ fConv = val;}
      void SetGstFlag(Int_t val){ fGst = val;}
      Int_t GetConvFlag(void){ return fConv;}
      Int_t GetGstFlag(void){ return fGst;}
      
 protected:
      Int_t fCharge; 
      Double_t fCentrality; 
      Double_t fXv; 
      Double_t fYv;
      Double_t fZv;
      Double_t fPx;
      Double_t fPy;
      Double_t fPz;
      Double_t fPt;
      Double_t fEta;
      Double_t fPhi;
      Double_t fTheta;
      Int_t fConv;
      Int_t fGst;
      AliVTrack *Obj;
      
      ClassDef(AliDielectronSingleTG, 2) // Event pool class              

};



class AliAnalysisTaskMultiDielectronTG : public AliAnalysisTaskSE {
  
public:
  AliAnalysisTaskMultiDielectronTG();
  AliAnalysisTaskMultiDielectronTG(const char *name);
  virtual ~AliAnalysisTaskMultiDielectronTG();

  enum ETriggerLogig {kAny, kExact};

  virtual void UserExec(Option_t *option);
  virtual void UserCreateOutputObjects();
  virtual void FinishTaskOutput();
  //temporary
//   virtual void NotifyRun(){AliDielectronPID::SetCorrVal((Double_t)fCurrentRunNumber);}
  
  void UsePhysicsSelection(Bool_t phy=kTRUE) {fSelectPhysics=phy;}
  void SetTriggerMask(ULong64_t mask) {fTriggerMask=mask;}
  UInt_t GetTriggerMask() const { return fTriggerMask; }
  void SetExcludeTriggerMask(ULong64_t mask) {fExcludeTriggerMask=mask;}
  UInt_t GetExcludeTriggerMask() const { return fExcludeTriggerMask; }
  void SetTriggerLogic(ETriggerLogig log) {fTriggerLogic=log;}
  ETriggerLogig GetTriggerLogic() const {return fTriggerLogic;}

  void SetEventFilter(AliAnalysisCuts * const filter) {fEventFilter=filter;}
  void SetTriggerOnV0AND(Bool_t v0and=kTRUE)    { fTriggerOnV0AND=v0and;    }
  void SetRejectPileup(Bool_t pileup=kTRUE)     { fRejectPileup=pileup;     }
  void AddDielectron(AliDielectron * const die) { fListDielectron.Add(die); }


  void reject_conversion(double val){d_conv_phiv = val;}
  void enable_v0mixing(Bool_t val){d_v0_mixing = val;}
  void check_ghost_pairs(vector<AliDielectronSingleTG*> e1);
  void fill_pair(AliDielectronSingleTG* e1, AliDielectronSingleTG* e2, int type, AliDielectron *die);
  bool PairTrackcut(double var1);
  void calc_vars(AliDielectronSingleTG* e1, AliDielectronSingleTG* e2, 
		 double &mass, double &phiv, double &px, double &py, double&pz,
		 double &pt, double &e, double &phi, double &eta, double &cos, double &psi);
  void calc_pair(vector<AliDielectronSingleTG*> e1, vector<AliDielectronSingleTG*> e2, AliDielectron *die, Int_t idie);
  void randomize_pool(vector<AliDielectronSingleTG*> e1, vector<AliDielectronSingleTG*> e2);        
  void reshuffle_buffer(vector<AliDielectronSingleTG*> ve, deque<AliDielectronSingleTG*> pool);

protected:
  enum {kAllEvents=0, kSelectedEvents, kV0andEvents, kFilteredEvents, kPileupEvents, kNbinsEvent};
  TList fListDielectron;             // List of dielectron framework instances
  TList fListHistos;                 //! List of histogram manager lists in the framework classes
  TList fListCF;                     //! List with CF Managers
  TList *tQAElectron;                     //! List with CF Managers
  

  Bool_t fSelectPhysics;             // Whether to use physics selection
  UInt_t fTriggerMask;               // Event trigger mask
  UInt_t fExcludeTriggerMask;        // Triggers to exclude from the analysis
  Bool_t fTriggerOnV0AND;            // if to trigger on V0and
  Bool_t fRejectPileup;              // pileup rejection wanted

  ETriggerLogig fTriggerLogic;       // trigger logic: any or all bits need to be matching
  
  AliTriggerAnalysis *fTriggerAnalysis; //! trigger analysis class

  AliAnalysisCuts *fEventFilter;     // event filter

  AliESDtrackCuts *fCutsMother;   

  TH1D *fEventStat;                  //! Histogram with event statistics
  TH1D *fEvent;
  TH2D *fdEdXvsPt;
  TH2D *fdEdXnSigmaElecvsPt;
  TH2D *fdEdXvsPtTOF;
  TH2D *fdEdXnSigmaElecvsPtTOF;
  TH2D *fTOFbetavsPt;
  TH2D *fTOFnSigmaElecvsPt;
  TH2F *hNCrossedRowsTPC; 
  TH2F *hChi2ClusTPC;
  TH2F *hRatioCrossClusTPC;

  Double_t fgValues[AliDielectronVarManager::kNMaxValues];
  std::vector<AliDielectronSingleTG*>  vem;
  std::vector<AliDielectronSingleTG*>  vep;
  std::vector<AliDielectronSingleTG*>  vem_tmp;
  std::vector<AliDielectronSingleTG*>  vep_tmp;
  Double_t d_conv_phiv; 
  Double_t bz;
  Bool_t d_v0_mixing;


  //Buffer for event mixing
  static const int NBUF=100; //depth of buffer
  static const int NMix=40; //# of events mixed (for +-)
  //static const int NMix=2; //# of events mixed (for +-)
  
  
  static const int NRPBIN=12;
  static const int NZBIN=10;
  static const int NCENT=10;
  static const int NDIE=10;
  int d_ibuf[NDIE][NZBIN][NCENT][NRPBIN];
  std::vector<AliDielectronSingleTG*> d_vep[NBUF][NDIE][NZBIN][NCENT][NRPBIN];
  std::vector<AliDielectronSingleTG*> d_vem[NBUF][NDIE][NZBIN][NCENT][NRPBIN];
  
  static const unsigned int MAXPOOL=500;
  //static const unsigned int MAXPOOL=50;
  static const int MAX_TRY=3;
  std::deque<AliDielectronSingleTG*> d_poolp[NDIE][NZBIN][NCENT][NRPBIN];
  std::deque<AliDielectronSingleTG*> d_poolm[NDIE][NZBIN][NCENT][NRPBIN]; 
  
  AliAnalysisTaskMultiDielectronTG(const AliAnalysisTaskMultiDielectronTG &c);
  AliAnalysisTaskMultiDielectronTG& operator= (const AliAnalysisTaskMultiDielectronTG &c);
  
  ClassDef(AliAnalysisTaskMultiDielectronTG, 2); //Analysis Task handling multiple instances of AliDielectron
};
#endif
