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
#include "AliDielectronVarManager.h"

using std::vector;
using std::deque;

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
    fObj(0x0)
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
      fGst(ghost), fObj(0x0)
	{
	  SetTrack(trk);
	  ;
	}

      AliDielectronSingleTG(const AliDielectronSingleTG&);
      AliDielectronSingleTG &operator=(const AliDielectronSingleTG&);
	


      ~AliDielectronSingleTG() {;}

      
      void SetTrack(AliVTrack * const trk) { fObj = trk;}
      virtual Int_t Charge(void) const { return fCharge;}
      Double_t  Phi(void) const { return fPhi;}
      Double_t  Eta(void) const { return fEta;}
      Double_t  Theta(void) const { return fTheta;}
      Double_t  Px(void) const { return fPx;}
      Double_t  Py(void) const { return fPy;}
      Double_t  Pz(void) const { return fPz;}
      Double_t  Xv(void) const { return fPx;}
      Double_t  Yv(void) const { return fPy;}
      Double_t  Zv(void) const { return fPz;}
      Double_t  Pt(void) const { return fPt;}
      AliVTrack *GetTrack(void) const { return fObj;}
      void SetConvFlag(Int_t val){ fConv = val;}
      void SetGstFlag(Int_t val){ fGst = val;}
      Int_t GetConvFlag(void) const { return fConv;}
      Int_t GetGstFlag(void) const { return fGst;}
      
 protected:
      Int_t fCharge;     ///charge of track
      Double_t fCentrality;  // centrality 
      Double_t fXv;  // vertex in X
      Double_t fYv;  // vertex in Y
      Double_t fZv;  // vertex in Z
      Double_t fPx;  // Momentum in X
      Double_t fPy;  // Momentum in Y
      Double_t fPz;  // Momentum in Z
      Double_t fPt;  // Momentum in Transverse
      Double_t fEta; // Particle Eta
      Double_t fPhi; // Particle Phi
      Double_t fTheta; //Particle Theta
      Int_t fConv; /// Conversion Flag
      Int_t fGst;  /// Ghost flag
      AliVTrack *fObj;  ///AliVTrack
      
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


  void RejectOP(double val){fdop = val;} ///To reject conversions
  void RejectConversion(double val, double mass){fdconvphiv = val; fdconvMee = mass;} ///To reject conversions
  void EnableV0mixing(Bool_t val){fdv0mixing = val;}   ///Enable V0 mixing  
  void CheckGhostPairs(std::vector<AliDielectronSingleTG*> e1); ///identify ghost pairs in like sign pais
  Bool_t CheckGhost(std::vector<AliDielectronSingleTG*> e1, std::vector<AliDielectronSingleTG*> e2); ///check ghost pairs for like sign and mixed like-sign pais
  void RejectPairs(std::vector<AliDielectronSingleTG*> e1, std::vector<AliDielectronSingleTG*> e2, Int_t idie); ///identify conversions for the rejection
  void FillPair(AliDielectronSingleTG* e1, AliDielectronSingleTG* e2, int type, AliDielectron *die, Int_t idie); /// Fill Pairs
  bool PairTrackcut(double var1, double var2, double var3, int idie); /// Pair cuts
  void CalcVars(AliDielectronSingleTG* e1, AliDielectronSingleTG* e2, 
		 double &mass, double &phiv, double &px, double &py, double&pz,
		double &pt, double &e, double &phi, double &eta, double &cos, double &psi); /// Calcualate kinematic variables
  void CalcPair(std::vector<AliDielectronSingleTG*> e1, std::vector<AliDielectronSingleTG*> e2, AliDielectron *die, Int_t idie);  ///Process Pairs
  void RandomizePool(std::vector<AliDielectronSingleTG*> e1, std::vector<AliDielectronSingleTG*> e2);         ///Randimize pairs
  void ReshuffleBuffer(std::vector<AliDielectronSingleTG*> ve, std::deque<AliDielectronSingleTG*> pool);  ///ReshuffleBuffer

  Double_t GetPhiv(AliDielectronSingleTG* e1, AliDielectronSingleTG* e2); /// calculate phiv
  Double_t GetOpeningAngle(AliDielectronSingleTG* e1, AliDielectronSingleTG* e2); /// calculate opening angle 
  Double_t GetMass(AliDielectronSingleTG* e1, AliDielectronSingleTG* e2); /// calculate Mass 
  void SetRejBGPairs(bool Val1, bool Val2){ fBGRejUnlike = Val1 ; fBGRejLike = Val2 ; } /// SetFlag to Enable Rejection of ghost BG pairs 
  void SetPairCuts(Int_t Types[20]){
    for(int i=0;i<20;i++){
      fRejectPairFlag[i] = Types[i];
    }
  }

protected:
  enum {kAllEvents=0, kSelectedEvents, kV0andEvents, kFilteredEvents, kPileupEvents, kNbinsEvent};
  TList fListDielectron;             // List of dielectron framework instances
  TList fListHistos;                 //! List of histogram manager lists in the framework classes
  TList fListCF;                     //! List with CF Managers
  TList *fQAElectron;                     //! List with CF Managers
  

  Bool_t fSelectPhysics;             // Whether to use physics selection
  UInt_t fTriggerMask;               // Event trigger mask
  UInt_t fExcludeTriggerMask;        // Triggers to exclude from the analysis
  Bool_t fTriggerOnV0AND;            // if to trigger on V0and
  Bool_t fRejectPileup;              // pileup rejection wanted

  enum PairRejType { NoRej=0, RejPairOp, RejPairPv, CutPairOp, CutPairPv};
  
  Int_t fRejectPairFlag[20];

  ETriggerLogig fTriggerLogic;       // trigger logic: any or all bits need to be matching
  
  AliTriggerAnalysis *fTriggerAnalysis; //! trigger analysis class

  AliAnalysisCuts *fEventFilter;     // event filter

  AliESDtrackCuts *fCutsMother;    /// Mother Cuts for QA

  TH1D *fEventStat;                  //! Histogram with event statistics
  TH1D *fEvent;                      // Centrality
  TH2D *fdEdXvsPt;                   // TPC dE/dx 
  TH2D *fdEdXnSigmaElecvsPt;         // TPC nSigmaEle vs. pt
  TH2D *fdEdXvsPtTOF;                // TPC dE/dx with TOF cut
  TH2D *fdEdXnSigmaElecvsPtTOF;      // TPC nSigmaEle vs. pt with TOF Cuts
  TH2D *fTOFbetavsPt;                // TOF beta vs. pT
  TH2D *fTOFnSigmaElecvsPt;          // TOF nSigma Electron vs. pT
  TH2F *fNCrossedRowsTPC;            // TPC NCrossedRows vs. pT
  TH2F *fChi2ClusTPC;                // TPC Chi2 Per Cluster 
  TH2F *fRatioCrossClusTPC;          // TPC Crossed rows per finable Clusters

  Double_t fgValues[AliDielectronVarManager::kNMaxValues];   /// Track/Pair information from AliDielectronVarManager
  std::vector<AliDielectronSingleTG*>  fVem;      /// Lists of electrons
  std::vector<AliDielectronSingleTG*>  fVep;      /// Lists of positions
  std::vector<AliDielectronSingleTG*>  fVemtmp;   /// template for electron lists
  std::vector<AliDielectronSingleTG*>  fVeptmp;   /// template for positron lists
  Double_t fdconvphiv;      /// PhiCut
  Double_t fdconvMee;      /// MassCut
  Double_t fdop;      /// Opening angle Cut
  Double_t fbz;            /// Magnetic field
  Bool_t fdv0mixing;       /// Mixing using V0 
  Bool_t fBGRejUnlike;     //// Ghost rejection flag for event mixing (unlike pairs)
  Bool_t fBGRejLike;     //// Ghost rejection flag for event mixing (like pairs)

  //Buffer for event mixing
  static const int fgkNBUF=100; //depth of buffer
  static const int fgkNMix=40; //# of events mixed (for +-)
  //static const int NMix=2; //# of events mixed (for +-)
  
  
  static const int fgkNRPBIN=12;    ///Number of RPbin for mixed event 
  static const int fgkNZBIN=10;     ///Number of zbin for mixed event 
  static const int fgkNCENT=10;     ///Number of centrality for mixed event 
  static const int fgkNDIE=20;      ///maximum number of cuts for AliDielectron
  int fibuf[fgkNDIE][fgkNZBIN][fgkNCENT][fgkNRPBIN];    ///buffer occupation for mixed event
  std::vector<AliDielectronSingleTG*> fvep[fgkNBUF][fgkNDIE][fgkNZBIN][fgkNCENT][fgkNRPBIN];   //// positron buffer for mixing
  std::vector<AliDielectronSingleTG*> fvem[fgkNBUF][fgkNDIE][fgkNZBIN][fgkNCENT][fgkNRPBIN];   //// electron buffer for mixing
    
  static const unsigned int fgkMAXPOOL=500;   ////maximum pool for mixing
  //static const unsigned int MAXPOOL=50;
  static const int fgkMAXTRY=3;    ///try to shuffle 
  std::deque<AliDielectronSingleTG*> fpoolp[fgkNDIE][fgkNZBIN][fgkNCENT][fgkNRPBIN]; ///pool for positrons
  std::deque<AliDielectronSingleTG*> fpoolm[fgkNDIE][fgkNZBIN][fgkNCENT][fgkNRPBIN]; ///pool for electrons
  
  AliAnalysisTaskMultiDielectronTG(const AliAnalysisTaskMultiDielectronTG &c);
  AliAnalysisTaskMultiDielectronTG& operator= (const AliAnalysisTaskMultiDielectronTG &c);
  
  ClassDef(AliAnalysisTaskMultiDielectronTG, 2); //Analysis Task handling multiple instances of AliDielectron
};
#endif
