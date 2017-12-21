#ifndef ALIANALYSISTASKSECHARMFRACTION_H
#define ALIANALYSISTASKSECHARMFRACTION_H

/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///*************************************************************************
/// \class Class AliAnalysisTaskSECharmFraction
/// \brief AliAnalysisTask for the extraction of the fraction of prompt charm
/// using the charm hadron impact parameter to the primary vertex
///
///
/// \author Author: Andrea Rossi andrea.rossi@pd.infn.it
//*************************************************************************

class TH1F;
class TH2F;
class TF1;
class AliAODDEvent;
class AliAODMCHeader;
class AliAODRecoDecayHF2Prong;
class AliAODRecoDecayHF;
class AliAODMCParticle;
class AliAnalysisVertexingHF;
class AliRDHFCutsD0toKpi;
class AliNormalizationCounter;
class AliVertexingHFUtils;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskSECharmFraction : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskSECharmFraction();
  AliAnalysisTaskSECharmFraction(const char *name);
  AliAnalysisTaskSECharmFraction(const char *name,AliRDHFCutsD0toKpi *cutsA,AliRDHFCutsD0toKpi *cutsB);

  virtual ~AliAnalysisTaskSECharmFraction(); 

  /// Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void Init();
  virtual void LocalInit() {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);  
  void SetReadMC(Bool_t readMC=kTRUE){fReadMC=readMC;}
  void SetSplitMassD0D0bar(Bool_t splitD0D0bar=kTRUE){fsplitMassD0D0bar=splitD0D0bar;}
  Bool_t GetIsSplitMassD0D0bar(){return fsplitMassD0D0bar;}
  void SetUsePID(Bool_t pid){fusePID=pid;}
  void SetAnalyzeLikeSign(Bool_t likesign=kFALSE){fLikeSign=likesign;}
  void SetNMaxTrForVtx(const Int_t ntrMaxforVtx){fNtrMaxforVtx=ntrMaxforVtx;}
  Int_t GetNMaxTrForVtx(){return fNtrMaxforVtx;}
  void SetPtBins(Int_t nbins,const Float_t *ptbins);
  void SetSignalInvMassCut(const Double_t signalInvMassCut=0.027){fsignalInvMassCut=signalInvMassCut;}
  void SetLargeInvMassCut(const Double_t largeInvMassCut=2.){flargeInvMassCut=largeInvMassCut;}
  void SetSideBandInvMassCut(const Double_t sidebandInvMassCut=0.054){// default value ~ 2x3 times inv mass resol.: a factor 2 is applied w.r.t. 3sigma, should be safe enough to exclude most of the reflections 
    fsidebandInvMassCut=sidebandInvMassCut;  
  }
  void SetSideBandInvMassWindow(const Double_t sidebandInvMassWindow=0.108){//~ 6 times inv. mass resol.
    fsidebandInvMassWindow=sidebandInvMassWindow;
  }  
  void SetAcceptanceCut(const Double_t eta=0.8,const Double_t nITSpoints=5.,const Double_t nSPDpoints=2.){fAcceptanceCuts[0]=eta;fAcceptanceCuts[1]=nITSpoints;fAcceptanceCuts[2]=nSPDpoints;}
  void SetStandardMassSelection();
  Int_t SetStandardCuts();
  void CheckInvMassD0(AliAODRecoDecayHF2Prong *d,Double_t &invMassD0,Double_t &invMassD0bar,Bool_t &isPeakD0,Bool_t &isPeakD0bar,Bool_t &isSideBandD0,Bool_t &isSideBandD0bar);
  void SetAnalysisLevel(Int_t level){fFastAnalysis=level;}
  void SetCheckBitD0flag(Bool_t checkfl){fcheckD0Bit=checkfl;}
  Bool_t GetCheckBitD0flag(){return fcheckD0Bit;}
  Int_t GetAnalysisLevel(){return fFastAnalysis;}
  Int_t CheckOrigin(const TClonesArray* arrayMC, const AliAODMCParticle *mcPartCandidate)const;
  AliAODRecoDecayHF *GetD0toKPiSignalType(const AliAODRecoDecayHF2Prong *d,TClonesArray *arrayMC,Int_t &signaltype,Double_t &massMumTrue,Double_t *primaryVtx,Int_t &isD0D0bar);
  AliAODRecoDecayHF *GetD0toKPiSignalTypeObsolete(const AliAODRecoDecayHF2Prong *d,TClonesArray *arrayMC,Int_t &signaltype,Double_t &massMumTrue,Double_t *primaryVtx);
  AliAODRecoDecayHF* ConstructFakeTrueSecVtx(const AliAODMCParticle *b1,const AliAODMCParticle *b2,const AliAODMCParticle *mum,Double_t *primaryVtxTrue);
  void SetUseMC(Bool_t useMC){fUseMC=useMC;}
  Bool_t SpecialSelD0(AliAODRecoDecayHF2Prong *d,Int_t &nusedforVtx);
  Bool_t FillAziList(AliAODEvent *aod,Double_t azilist[30000],Int_t trkIDlist[30000],Int_t &nprim)const;
  void FillAziHistos(AliAODRecoDecayHF2Prong *d,TList *&list,Int_t ptbin,Double_t azilist[30000],Int_t trkIDlist[30000],Int_t nprim,Int_t okD0,Int_t okD0bar,Bool_t isPeakD0,Bool_t isPeakD0bar,Bool_t isSideBandD0,Bool_t isSideBandD0bar)const;

  AliAODVertex* GetPrimaryVtxSkipped(AliAODEvent *aodev,AliAODRecoDecayHF2Prong *d);
  void SetRejecCandidateMCUpgrade(Bool_t selection){fselectForUpgrade=selection;}
  void SetSkipEventSelection(Bool_t skip){fskipEventSelection=skip;}
  void SetMaxZvtxForSkipEventSelection(Double_t zmax){fZvtxUpgr=zmax;}
  void SetPtWeightsFromDataPbPb276overLHC12a17a();
  void SetFillImpParTree(Bool_t fillimppar){fFillTree=fillimppar;}
  void SetAODMismatchProtection(Int_t opt=1) {fAODProtection=opt;}
  void SetLightOutput(Bool_t lightOutput){fLightOutput=lightOutput;}

  /* ######### THE FOLLOWING IS FOR FURTHER IMPLEMENATION ############
     Int_t GetPtBin(Double_t pt)const;
     void SetD0Cuts(Int_t ptbin,Double_t &*d0cutsLoose,Double_t &*d0cutsTight);
     
     //  void InvMassSelection();
     
     void SetCheckMC(Bool_t checkMC){fcheckMC=checkMC;}
     void SetCheckMC_D0(Bool_t check_D0){fcheckMCD0=check_D0;}
     void SetCheckMC_2prongs(Bool_t check2prongs){fcheckMC2prongs=check2prongs;}
     void SetCheckMC_prompt(Bool_t checkprompt){fcheckMCprompt=checkprompt;}
     void SetCheckMC_fromB(Bool_t checkfromB){fcheckMCfromB=checkfromB;}
     void SetCheckMC_fromDstar(Bool_t skipD0star){fSkipD0star=skipD0star;}
     void SetUseCuts(Bool_t usecuts){fD0usecuts=usecuts;}
     void SetSideBands(Double_t sideband){fSideBands=sideband;}
     void SetStudyPureBackground(Bool_t back){fStudyPureBackground=back;}
  */
  AliRDHFCutsD0toKpi* GetLooseCut(){
    return fCutsLoose;
  }
  AliRDHFCutsD0toKpi* GetTightCut(){
    return fCutsTight;
  }
  /* void SetCutFunction(Int_t (*setcuts)(AliAnalysisTaskSECharmFraction*,Double_t,Double_t)){
     fSetCuts=setcuts;
     fStandCuts=kFALSE;
     }
  */
  //  Int_t SetCuts(AliAnalysisTaskSECharmFraction *alchfr,Double_t pt,Double_t invMassCut);
  
 private:
  Bool_t FillHistos(AliAODRecoDecayHF2Prong *d,TList *&list,Int_t ptbin,Int_t okD0,Int_t okD0bar,Double_t invMassD0,Double_t invMassD0bar,Bool_t isPeakD0,Bool_t isPeakD0bar,Bool_t isSideBandD0,Bool_t isSideBandD0bar,Double_t massmumtrue,AliAODRecoDecayHF *aodDMC,Double_t *vtxTrue,Int_t isD0D0barMC,Double_t bField);
  void FillHistoMCproperties(TClonesArray *arrayMC);

  AliRDHFCutsD0toKpi *fCutsLoose;        /// Loose cuts object
  AliRDHFCutsD0toKpi *fCutsTight;      /// Vertexer heavy flavour
  Int_t fFastAnalysis;                  /// Level of analysis speed: default is 1, switch it to 2 to fill the THnSparse
  Bool_t  fLightOutput;                 ///Flag to fill just the output for the simpler analysis 
  Bool_t  fReadMC;                          /// Flag To switch on/off access to MC 
  Bool_t fcheckD0Bit;                       /// Flag to check the D0 bit flag
  Bool_t  fsplitMassD0D0bar;                /// Flag to use two shistos for D0 and D0bar invariant masses
  Bool_t  fLikeSign;                        /// Flag to analyse Like Sign array
  Bool_t  fusePID;                          /// Flag to use PID
  Double_t    fmD0PDG;                      ///  MC D0 mass
  Int_t        fnbins;                      /// Number of pt bins
  Float_t *fptbins;                        //[fnbins] ptbins 
  Int_t fNtrMaxforVtx;                      /// N Max acceptable tracks used for vertex (0,1,2)
  Double_t fptAll;                          //!<! Sum of pt of the reco tracks
  Double_t fptAllSq;                        //!<! Sum of the square of the pt of the reco tracks
  Double_t fptMax[3];                       //!<! Three largest track pt in the event
  Double_t fAcceptanceCuts[3];                /// array with acceptance cuts
  Double_t fsignalInvMassCut;               /// invariant mass cut to define signal region
  Double_t flargeInvMassCut;                /// invariant mass cut to accept all inv mass window
  Double_t fsidebandInvMassCut;             /// invariant mass cut to define side band region lower limit
  Double_t fsidebandInvMassWindow;          /// invariant mass cut to define side band region width
  Bool_t fUseMC;                            /// flag to use or not MC info
  Bool_t fCleanCandOwnVtx;                  /// flag to switch on/off cleaning of the candidate own vtx
  Int_t fAODProtection;  /// flag to activate protection against AOD-dAOD mismatch.
  TH1F *fNentries;                          //!<! histo for #AOD analysed, container 1
  TH1F *fSignalType;                        //!<! histo for the type of MC signal , container 2
  TH1F *fSignalTypeLsCuts;                 //!<! histo for the type of MC signal with loose cuts , container 3
  TH1F *fSignalTypeTghCuts;                //!<! histo for the type of MC signal with tight cuts, container 4
  AliNormalizationCounter *fCounter;        //!<! counter for the normalization 
  TList *flistMCproperties;               //!<! TLists for MC properties of D0 w.r.t. B mesons and c quarks cntainer 5
  TList *flistNoCutsSignal;               //!<! TList for signal (D prompt) with nocuts, container 6
  TList *flistNoCutsBack;               //!<! TList for background with nocuts, container 7
  TList *flistNoCutsFromB;               //!<! TList for D from B or D from Dstar from Bwith nocuts, container 8
  TList *flistNoCutsFromDstar;               //!<! TList for D from Dstar with nocuts, container 9
  TList *flistNoCutsOther;               //!<! TList for others with nocuts, container 10
  TList *flistLsCutsSignal;               //!<! TList for signal (D prompt) with loose cuts, container 11
  TList *flistLsCutsBack;               //!<! TList for background with loose cuts, container 12
  TList *flistLsCutsFromB;               //!<! TList for D from B or D from Dstar from B with loose cuts, container 13
  TList *flistLsCutsFromDstar;               //!<! TList for D from Dstar with loose cuts, container 14
  TList *flistLsCutsOther;               //!<! TList for others with loose cuts, container 15
  TList *flistTghCutsSignal;               //!<! TList for signal (D prompt) with tight cuts, container 16
  TList *flistTghCutsBack;               //!<! TList for backgrnd with tight cuts, container 17
  TList *flistTghCutsFromB;               //!<! TList for D from B or D from Dstar from Bwith tight cuts, container 18
  TList *flistTghCutsFromDstar;               //!<! TList for D from Dstar Dstar with tight cuts, container 19
  TList *flistTghCutsOther;               //!<! TList for others with tight cuts, container 20
  AliVertexingHFUtils *fVertUtil;         /// vertexing HF Util
  Bool_t fselectForUpgrade;               /// switch to reject candidates from HIJING and not Pythia for upgrade studies
  Bool_t fskipEventSelection;               /// switch to skip event selection (for upgrade studies)
  Double_t fZvtxUpgr;                       /// cut value on max zvtx used ONLY if fskipEventSelection is kTRUE
  TF1 *fWeightPt;                           /// function with pt weights used only for MC histos for reflections and signal mass shape
  /*  Bool_t       fD0usecuts;            /// Switch on the use of the cuts             TO BE IMPLEMENTED 
      Bool_t       fcheckMC;              ///  Switch on MC check: minimum check is same mother  TO BE IMPLEMENTED
      Bool_t       fcheckMCD0;           ///  check the mother is a D0                  TO BE IMPLEMENTED
      Bool_t       fcheckMC2prongs;         ///  check the decay is in two prongs       TO BE IMPLEMENTED  
      Bool_t       fcheckMCprompt;       ///  check the D0 comes from a c quark         TO BE IMPLEMENTED
      Bool_t       fcheckMCfromB;        ///  check the D0 comes from a b quark         TO BE IMPLEMENTED
      Bool_t       fSkipD0star;           /// skip if D0 comes from a D*                TO BE IMPLEMENTED
      Bool_t  fStudyd0fromBTrue;         /// Flag for analyze true impact par of D0 from B       TO BE IMPLEMENTED 
      Bool_t  fStudyPureBackground;      /// Flag to study the background (reverse the selection on the signal)     TO BE IMPLEMENTED 
      Double_t  fSideBands;                ///Side bands selection (see cxx)            TO BE IMPLEMENTED
  */ 
  Bool_t fFillTree;                        /// Switch to fill the trees of variabls for the impact parameter unbinned fit
  TTree *fTreeNCsign;                      //!<!
  TTree *fTreeNCback;                      //!<!
  TTree *fTreeNCfromB;                     //!<!
  TTree *fTreeRecoNCfromB;                 //!<!
  TTree *fTreeNCfromDstar;                 //!<!
  TTree *fTreeNCother;                     //!<!
  TTree *fTreeLSCsign;                     //!<!
  TTree *fTreeLSCback;                     //!<!
  TTree *fTreeLSCfromB;                    //!<!
  TTree *fTreeRecoLSCfromB;                //!<!
  TTree *fTreeLSCfromDstar;                //!<!
  TTree *fTreeLSCother;                    //!<!
  TTree *fTreeTGHCsign;                    //!<!
  TTree *fTreeTGHCback;                    //!<!
  TTree *fTreeTGHCfromB;                   //!<!
  TTree *fTreeRecoTGHCfromB;               //!<!
  TTree *fTreeTGHCfromDstar;               //!<!
  TTree *fTreeTGHCother;                   //!<!
  Double_t *fVariablesTreeNCsign;          //!<!
  Double_t *fVariablesTreeNCback;          //!<!
  Double_t *fVariablesTreeNCfromB;         //!<!
  Double_t *fVariablesTreeRecoNCfromB;     //!<!
  Double_t *fVariablesTreeNCfromDstar;     //!<!
  Double_t *fVariablesTreeNCother;         //!<! 
  Double_t *fVariablesTreeLSCsign;         //!<!
  Double_t *fVariablesTreeLSCback;         //!<!
  Double_t *fVariablesTreeLSCfromB;        //!<!
  Double_t *fVariablesTreeRecoLSCfromB;    //!<!
  Double_t *fVariablesTreeLSCfromDstar;    //!<!
  Double_t *fVariablesTreeLSCother;        //!<!
  Double_t *fVariablesTreeTGHCsign;        //!<!
  Double_t *fVariablesTreeTGHCback;        //!<!
  Double_t *fVariablesTreeTGHCfromB;       //!<!
  Double_t *fVariablesTreeRecoTGHCfromB;   //!<!
  Double_t *fVariablesTreeTGHCfromDstar;   //!<!
  Double_t *fVariablesTreeTGHCother;       //!<!

  AliAnalysisTaskSECharmFraction(const AliAnalysisTaskSECharmFraction&); // not implemented
  AliAnalysisTaskSECharmFraction& operator=(const AliAnalysisTaskSECharmFraction&); // not implemented
  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskSECharmFraction,8); /// analysis task for prompt charm fraction
  /// \endcond 
};

#endif
