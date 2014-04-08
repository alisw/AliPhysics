/*************************************************************************
 *                                                                       *
 * Task for Jet Chemistry Analysis in PWG-JE Jet Task Force Train        *
 *                                                                       *
 *                                                                       *
 *  contact:                                                             *
 *  Alice Zimmermann                                                     *
 *  zimmermann@physi.uni-heidelberg.de                                   *
 *                                                                       *
 *                                                                       *  
 *                                                                       *
 *                                                                       *
 *************************************************************************/
   

#ifndef ALIANALYSISTASKJETCHEM_H
#define ALIANALYSISTASKJETCHEM_H

/* Copyright(c) 1998-2013, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

class AliAODv0;
class AliAODVertex;
class AliAnalysisCentralitySelector;
class AliPIDResponse;
class TString;
class TList;
class AliAODMCParticle;
class AliAODTrack;

#include "AliAnalysisTaskFragmentationFunction.h"
#include "AliPID.h"
#include "AliAnalysisTaskSE.h"
#include "AliAODMCHeader.h"

class AliAnalysisTaskJetChem : public AliAnalysisTaskFragmentationFunction {

 public:
  
  //----------------------------------------
  class AliFragFuncHistosInvMass : public TObject
  {
    
    public:
    
    AliFragFuncHistosInvMass(const char* name = "FFIMhistos", 
			     Int_t nJetPt = 0, Float_t jetPtMin = 0, Float_t jetPtMax = 0,
			     Int_t nInvMass = 0, Float_t invMassMin=0, Float_t invMassMax=0,  
			     Int_t nPt = 0, Float_t ptMin = 0, Float_t ptMax = 0,
			     Int_t nXi = 0, Float_t xiMin = 0, Float_t xiMax = 0,
			     Int_t nZ  = 0, Float_t zMin  = 0, Float_t zMax  = 0);
    AliFragFuncHistosInvMass(const AliFragFuncHistosInvMass& copy);
    AliFragFuncHistosInvMass& operator=(const AliFragFuncHistosInvMass &o);
    virtual ~AliFragFuncHistosInvMass();
    
    virtual void DefineHistos();
    virtual void FillFF(Float_t trackPt, Float_t invM, Float_t jetPt,Bool_t incrementJetPt);
    virtual void AddToOutput(TList* list) const;

  private:

    Int_t   fNBinsJetPt;    // FF histos bins
    Float_t fJetPtMin;      // FF histos limits
    Float_t fJetPtMax;      // FF histos limits
    Int_t   fNBinsInvMass;  // FF histos bins
    Float_t fInvMassMin;    // FF histos limits
    Float_t fInvMassMax;    // FF histos limits
    Int_t   fNBinsPt;       // FF histos bins
    Float_t fPtMin;         // FF histos limits
    Float_t fPtMax;         // FF histos limits
    Int_t   fNBinsXi;       // FF histos bins
    Float_t fXiMin;         // FF histos limits
    Float_t fXiMax;         // FF histos limits
    Int_t   fNBinsZ;        // FF histos bins
    Float_t fZMin;          // FF histos limits
    Float_t fZMax;          // FF histos limits
  


    TH3F*   fh3TrackPt;     //! FF: track transverse momentum 
    TH3F*   fh3Xi;          //! FF: xi 
    TH3F*   fh3Z;           //! FF: z  
    TH1F*   fh1JetPt;       //! jet pt of all jets 

    TString fNameFF;        // histo names prefix
    
    ClassDef(AliFragFuncHistosInvMass, 1);
  };
  

 //----------------------------------------
  class AliFragFuncHistosPhiCorrInvMass : public TObject
  {
				   
    public:
    
    AliFragFuncHistosPhiCorrInvMass(const char* name = "FFPhiCorrIMhistos", 
				    Int_t nPt = 0, Float_t ptMin = 0, Float_t ptMax = 0,
				    Int_t nPhi = 0, Float_t phiMin = 0, Float_t phiMax = 0,
				    Int_t nInvMass = 0, Float_t invMassMin=0, Float_t invMassMax=0);

    AliFragFuncHistosPhiCorrInvMass(const AliFragFuncHistosPhiCorrInvMass& copy);
    AliFragFuncHistosPhiCorrInvMass& operator=(const AliFragFuncHistosPhiCorrInvMass &o);
    virtual ~AliFragFuncHistosPhiCorrInvMass();
    
    virtual void DefineHistos();
    virtual void FillPhiCorr(Float_t pt, Float_t phi, Float_t invM);
    virtual void AddToOutput(TList* list) const;

  private:

    Int_t   fNBinsPt;       // FF histos bins
    Float_t fPtMin;         // FF histos limits
    Float_t fPtMax;         // FF histos limits

    Int_t   fNBinsPhi;      // FF histos bins
    Float_t fPhiMin;        // FF histos limits
    Float_t fPhiMax;        // FF histos limits
    
    Int_t   fNBinsInvMass;  // FF histos bins
    Float_t fInvMassMin;    // FF histos limits
    Float_t fInvMassMax;    // FF histos limits
  
    TH3F*   fh3PhiCorr;     //! FF: phi correlation histo 

    TString fNamePhiCorr;   // histo names prefix
    
    ClassDef(AliFragFuncHistosPhiCorrInvMass, 1);
  };
  
  //----------------------------------------

  AliAnalysisTaskJetChem(); 
  AliAnalysisTaskJetChem(const char *name);
  AliAnalysisTaskJetChem(const  AliAnalysisTaskJetChem &copy);
  AliAnalysisTaskJetChem& operator=(const  AliAnalysisTaskJetChem &o);
  virtual ~AliAnalysisTaskJetChem();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);

  static  void   SetProperties(TH3F* h,const char* x, const char* y,const char* z);

  Bool_t IsAccepteddEdx(Double_t mom, Double_t signal, AliPID::EParticleType n, Double_t cutnSig) const;//not used anymore
  Bool_t IsK0InvMass(Double_t mass) const; 
  Int_t  GetListOfV0s(TList *list, Int_t type, Int_t particletype, AliAODVertex* primVertex, AliAODEvent* aod);
  Int_t  GetListOfParticles(TList *list, Int_t type, Int_t particletype, AliAODVertex* primVertex);
  Int_t  GetListOfMCParticles(TList *outputlist, Int_t particletype, AliAODEvent* mcaodevent);
  void   GetTracksInCone(TList* inputlist, TList* outputlist, const AliAODJet* jet, Double_t radius, Double_t& sumPt, Double_t minPt,  Double_t maxPt, Bool_t& isBadPt);
  void   GetTracksInPerpCone(TList* inputlist, TList* outputlist, const AliAODJet* jet, Double_t radius, Double_t& sumPerpPt);
  Bool_t MCLabelCheck(AliAODv0* v0, Int_t particletype, const AliAODTrack* trackNeg, const AliAODTrack* trackPos, TList *listmc, Int_t& negDaughterpdg, Int_t& posDaughterpdg, Int_t& motherType, Int_t& v0Label, Double_t& MCPt, Bool_t& fPhysicalPrimary, Int_t& MCv0PDGCode);
  Bool_t IsParticleMatching(const AliAODMCParticle* mcp0, Int_t v0Label);
  Bool_t DaughterTrackCheck(AliAODv0* v0, Int_t& nnum, Int_t& pnum);
  Int_t  IsTrackInjected(AliAODv0 *v0, AliAODMCHeader *header, TClonesArray *arrayMC);
  TString GetGenerator(Int_t label, AliAODMCHeader* header);
  Double_t SmearJetPt(Double_t jetPt, Int_t cl, Double_t jetRadius, Double_t ptmintrack, Double_t& jetPtSmear);
  AliAODJet* GetMedianCluster();
  
  virtual void SetK0Type(Int_t i){ fK0Type = i; }
  virtual void SetFilterMaskK0(UInt_t i) {fFilterMaskK0 = i;}

  Bool_t IsLaInvMass(Double_t mass) const; //La and ALa mass check
  virtual void SetLaType(Int_t i){ fLaType = i; }
  virtual void SetFilterMaskLa(UInt_t i) {fFilterMaskLa = i;}

  virtual void SetALaType(Int_t i){ fALaType = i; }
  virtual void SetFilterMaskALa(UInt_t i) {fFilterMaskALa = i;}

  virtual void SetSelectArmenteros(Bool_t b) {IsArmenterosSelected = b;}
  //virtual void SetEventSelectionMask(UInt_t i){fEvtSelectionMask = i;}  //already inherited by AliAnalysisFragmentationFunctionTask
  //virtual void UsePhysicsSelection(Bool_t b) {fUsePhysicsSelection = b;} //already inherited by AliAnalysisFragmentationFunctionTask

  void CalculateInvMass(AliAODv0* v0vtx, Int_t particletype, Double_t& invM, Double_t& trackPt);
  
  Bool_t AcceptBetheBloch(AliAODv0 *v0, AliPIDResponse *PIDResponse, Int_t particletype); //don't use this method for MC Analysis


  Double_t MyRapidity(Double_t rE, Double_t rPz) const;

  //-- K0s

  void   SetFFInvMassHistoBins(Int_t nJetPt = 39, Float_t jetPtMin = 5., Float_t jetPtMax = 200., //previous 19, 5.,100.
			       Int_t nInvM = 400, Float_t invMMin = 0.4,  Float_t invMMax = 0.6, //previous 0.4 to 0.6
			       Int_t nPt = 200, Float_t ptMin = 0., Float_t ptMax = 20.,         //previous 0. to 10.
			       Int_t nXi = 35, Float_t xiMin = 0., Float_t xiMax = 7.,
			       Int_t nZ = 11,  Float_t zMin = 0.,  Float_t zMax = 1.1)
  { fFFIMNBinsJetPt = nJetPt; fFFIMJetPtMin = jetPtMin; fFFIMJetPtMax = jetPtMax; 
    fFFIMNBinsInvM = nInvM; fFFIMInvMMin = invMMin; fFFIMInvMMax = invMMax; fFFIMNBinsPt = nPt; fFFIMPtMin = ptMin; fFFIMPtMax = ptMax; 
    fFFIMNBinsXi = nXi; fFFIMXiMin = xiMin; fFFIMXiMax = xiMax; fFFIMNBinsZ  = nZ;  fFFIMZMin  = zMin;  fFFIMZMax  = zMax; }


  void   SetPhiCorrInvMassHistoBins(Int_t nPt = 40, Float_t ptMin = 0., Float_t ptMax = 20., 
				    Int_t nPhi = 20, Float_t phiMin = 0., Float_t phiMax = 2*TMath::Pi(),
				    Int_t nInvM = 240, Float_t invMMin = 0.4,  Float_t invMMax = 0.6)
				    
  { fPhiCorrIMNBinsPt = nPt; fPhiCorrIMPtMin = ptMin; fPhiCorrIMPtMax = ptMax;
    fPhiCorrIMNBinsPhi = nPhi; fPhiCorrIMPhiMin = phiMin; fPhiCorrIMPhiMax = phiMax;
    fPhiCorrIMNBinsInvM = nInvM; fPhiCorrIMInvMMin = invMMin; fPhiCorrIMInvMMax = invMMax;
  }
  
  // --La and ALa

  void   SetFFInvMassLaHistoBins(Int_t nJetPt = 39, Float_t jetPtMin = 5., Float_t jetPtMax = 200., //La
				 //Int_t nInvM = 140, Float_t invMMin = 1.06,  Float_t invMMax = 1.2,//original inv. mass range, now I shifted to Vits slightly larger mass window
			       Int_t nInvM = 200, Float_t invMMin = 1.05,  Float_t invMMax = 1.25,
			       Int_t nPt = 200, Float_t ptMin = 0., Float_t ptMax = 20., 
			       Int_t nXi = 35, Float_t xiMin = 0., Float_t xiMax = 7.,
			       Int_t nZ = 11,  Float_t zMin = 0.,  Float_t zMax = 1.1)

  { fFFIMLaNBinsJetPt = nJetPt; fFFIMLaJetPtMin = jetPtMin; fFFIMLaJetPtMax = jetPtMax; 
    fFFIMLaNBinsInvM = nInvM; fFFIMLaInvMMin = invMMin; fFFIMLaInvMMax = invMMax; fFFIMLaNBinsPt = nPt; fFFIMLaPtMin = ptMin; fFFIMLaPtMax = ptMax; 
    fFFIMLaNBinsXi = nXi; fFFIMLaXiMin = xiMin; fFFIMLaXiMax = xiMax; fFFIMLaNBinsZ  = nZ;  fFFIMLaZMin  = zMin;  fFFIMLaZMax  = zMax; }


  void   SetPhiCorrInvMassLaHistoBins(Int_t nPt = 40, Float_t ptMin = 0., Float_t ptMax = 20.,    //La
				    Int_t nPhi = 20, Float_t phiMin = 0., Float_t phiMax = 2*TMath::Pi(),
				    Int_t nInvM = 50, Float_t invMMin = 1.06,  Float_t invMMax = 1.2)
				    
  { fPhiCorrIMLaNBinsPt = nPt; fPhiCorrIMLaPtMin = ptMin; fPhiCorrIMLaPtMax = ptMax;
    fPhiCorrIMLaNBinsPhi = nPhi; fPhiCorrIMLaPhiMin = phiMin; fPhiCorrIMLaPhiMax = phiMax;
    fPhiCorrIMLaNBinsInvM = nInvM; fPhiCorrIMLaInvMMin = invMMin; fPhiCorrIMLaInvMMax = invMMax;
  }


 
  // consts

  enum { kTrackUndef =0, kOnFly, kOnFlyPID, kOnFlydEdx, kOnFlyPrim, kOffl, kOfflPID, kOffldEdx, kOfflPrim };  
  enum { kK0, kLambda, kAntiLambda };  
 
  //--
  Bool_t   fAnalysisMC;
  Double_t fDeltaVertexZ;
  Double_t fCuttrackNegNcls;
  Double_t fCuttrackPosNcls; 
  Double_t fCutPostrackRap;
  Double_t fCutNegtrackRap;
  Double_t fCutRap;
  Double_t fCutPostrackEta;
  Double_t fCutNegtrackEta;
  Double_t fCutEta;
  Double_t fCutV0cosPointAngle;
  Double_t fCutChi2PosDaughter;
  Double_t fCutChi2NegDaughter;
  Bool_t   fKinkDaughters;
  Bool_t   fRequireTPCRefit;
  Double_t fCutArmenteros; 
  Double_t fCutV0DecayMin;
  Double_t fCutV0DecayMax;
  Double_t fCutV0totMom;
  Double_t fCutDcaV0Daughters;
  Double_t fCutDcaPosToPrimVertex;
  Double_t fCutDcaNegToPrimVertex;
  Double_t fCutV0RadiusMin;
  Double_t fCutV0RadiusMax;
  Double_t fCutBetheBloch;
  Double_t fCutRatio;

  // cuts
  void SetCuttrackPosNcls(Double_t posNcls){fCuttrackPosNcls=posNcls; Printf("AliAnalysisTaskJetChem:: SetCuttrackPosNcls %f",posNcls);}
  void SetCuttrackNegNcls(Double_t negNcls){fCuttrackNegNcls=negNcls; Printf("AliAnalysisTaskJetChem:: SetCuttrackNegNcls %f",negNcls);}
  void SetCuttrackPosRap(Double_t posRap){fCutPostrackRap=posRap; Printf("AliAnalysisTaskJetChem:: SetCuttrackPosRap %f",posRap);}
  void SetCuttrackNegRap(Double_t negRap){fCutNegtrackRap=negRap; Printf("AliAnalysisTaskJetChem:: SetCuttrackNegRap %f",negRap);}
  void SetCutV0Rap(Double_t v0Rap){fCutRap=v0Rap; Printf("AliAnalysisTaskJetChem:: SetCutV0Rap %f",v0Rap);}
  void SetCuttrackPosEta(Double_t posEta){fCutPostrackEta=posEta; Printf("AliAnalysisTaskJetChem:: SetCuttrackPosEta %f",posEta);}
  void SetCuttrackNegEta(Double_t negEta){fCutNegtrackEta=negEta; Printf("AliAnalysisTaskJetChem:: SetCuttrackNegEta %f",negEta);}
  void SetCutV0Eta(Double_t v0Eta){fCutEta=v0Eta; Printf("AliAnalysisTaskJetChem:: SetCutV0Eta %f",v0Eta);}
  void SetCosOfPointingAngle(Double_t cospointAng){fCutV0cosPointAngle=cospointAng; Printf("AliAnalysisTaskJetChem:: SetCosOfPointingAngle %f",cospointAng);}
  void SetChi2CutPosDaughter(Double_t chi2PosDaughter){fCutChi2PosDaughter=chi2PosDaughter; Printf("AliAnalysisTaskJetChem:: SetChi2CutPosDaughter %f",chi2PosDaughter);}
  void SetChi2CutNegDaughter(Double_t chi2NegDaughter){fCutChi2NegDaughter=chi2NegDaughter; Printf("AliAnalysisTaskJetChem:: SetChi2CutNegDaughter %f",chi2NegDaughter);}
  void SetAcceptKinkDaughters(Bool_t isKinkDaughtersAccepted){fKinkDaughters=isKinkDaughtersAccepted; Printf("AliAnalysisTaskJetChem:: SetAcceptKinkDaughters %i", isKinkDaughtersAccepted);}
  void SetRequireTPCRefit(Bool_t isTPCRefit){fRequireTPCRefit=isTPCRefit; Printf("AliAnalysisTaskJetChem:: SetRequireTPCRefit %i", isTPCRefit);}
  void SetCutArmenteros(Double_t armenteros){fCutArmenteros=armenteros; Printf("AliAnalysisTaskJetChem:: SetCutArmenteros %f", armenteros);}
  void SetCutV0DecayMin(Double_t decayMin){fCutV0DecayMin=decayMin; Printf("AliAnalysisTaskJetChem:: SetCutDecayMin %f", decayMin);}
  void SetCutV0DecayMax(Double_t decayMax){fCutV0DecayMax=decayMax; Printf("AliAnalysisTaskJetChem:: SetCutDecayMax %f", decayMax);}
  void SetCutV0totMom(Double_t v0totMom){fCutV0totMom=v0totMom; Printf("AliAnalysisTaskJetChem:: SetCutV0totMom %f", v0totMom);}
  void SetCutDcaV0Daughters(Double_t dcav0daughters){fCutDcaV0Daughters=dcav0daughters; Printf("AliAnalysisTaskJetChem:: SetCutDcaV0Daughters %f", dcav0daughters);}
  void SetCutDcaPosToPrimVertex(Double_t dcaPosToPrimVertex){fCutDcaPosToPrimVertex=dcaPosToPrimVertex; Printf("AliAnalysisTaskJetChem:: SetCutDcaPosToPrimVertex %f", dcaPosToPrimVertex);}
  void SetCutDcaNegToPrimVertex(Double_t dcaNegToPrimVertex){fCutDcaNegToPrimVertex=dcaNegToPrimVertex; Printf("AliAnalysisTaskJetChem:: SetCutDcaNegToPrimVertex %f", dcaNegToPrimVertex);}
  void SetCutV0RadiusMin(Double_t v0RadiusMin){fCutV0RadiusMin=v0RadiusMin; Printf("AliAnalysisTaskJetChem:: SetCutV0RadiusMin %f", v0RadiusMin);}
  void SetCutV0RadiusMax(Double_t v0RadiusMax){fCutV0RadiusMax=v0RadiusMax; Printf("AliAnalysisTaskJetChem:: SetCutV0RadiusMax %f", v0RadiusMax);}
  void SetCutBetheBloch(Double_t cutBetheBloch){fCutBetheBloch=cutBetheBloch; Printf("AliAnalysisTaskJetChem:: SetCutBetheBloch %f", cutBetheBloch);}
  void SetCutRatioTPC(Double_t cutRatioTPC){fCutRatio=cutRatioTPC; Printf("AliAnalysisTaskJetChem:: SetCutRatioTPC %f", cutRatioTPC);}
  void SetAnalysisMC(Bool_t analysisMC) {fAnalysisMC = analysisMC;}
  void SetDeltaZVertexCut(Float_t deltaVtxZ){fDeltaVertexZ = deltaVtxZ;}

 private:
  
  Int_t fK0Type;                                           //! K0 cuts
  UInt_t fFilterMaskK0;                                    //! K0 legs cuts
  TList* fListK0s;                                         //! K0 list 
  AliPIDResponse *fPIDResponse;	                           // PID

  AliFragFuncQATrackHistos*  fV0QAK0;                      //! track QA: V0s in K0 inv mass range
  AliFragFuncHistos*         fFFHistosRecCutsK0Evt;        //! inclusive FF for K0 evt
  AliFragFuncHistosInvMass*  fFFHistosIMK0AllEvt;          //! K0 pt spec for all events
  AliFragFuncHistosInvMass*  fFFHistosIMK0Jet;             //! K0 FF all dPhi   
  AliFragFuncHistosInvMass*  fFFHistosIMK0Cone;            //! K0 FF jet cone   
  AliFragFuncHistosPhiCorrInvMass*  fFFHistosPhiCorrIMK0;  //! K0 correlation to jet axis 
  
  Int_t fLaType;                                           //! La cuts
  UInt_t fFilterMaskLa;                                    //! La legs cuts
  TList* fListLa;                                          //! La list 
  
  AliFragFuncHistosInvMass*  fFFHistosIMLaAllEvt;          //! La pt spec for all events
  AliFragFuncHistosInvMass*  fFFHistosIMLaJet;             //! La FF all dPhi   
  AliFragFuncHistosInvMass*  fFFHistosIMLaCone;            //! La FF jet cone   
  AliFragFuncHistosPhiCorrInvMass*  fFFHistosPhiCorrIMLa;  //! La correlation to jet axis 
  
  Int_t fALaType;                                          //! ALa cuts

  UInt_t fFilterMaskALa;                                   //! ALa legs cuts
  TList* fListALa;                                         //! ALa list 
  TList* fListFeeddownLaCand;                              //! feeddown from Xi (-,0) 
  TList* fListFeeddownALaCand;                             //! feeddown from Xibar (+,0) 
  TList* jetConeFDLalist;                                  //! feeddown from Xi (-,0) in jet cone
  TList* jetConeFDALalist;                                 //! feeddown from Xibar (+,0) in jet cone
  TList* fListMCgenK0s;                                    //! MC generated K0s
  TList* fListMCgenLa;                                     //! MC generated La                 
  TList* fListMCgenALa;                                    //! MC generated ALa              
  TList* fListMCgenK0sCone;                                //! MC generated K0s in cone around jet axis, particles are from fragmentation but also from underlying event  
  TList* fListMCgenLaCone;                                 //! MC generated Lambdas in cone around jet axis, particles are from fragmentation but also from underlying event
  TList* fListMCgenALaCone;                                //! MC generated Antilambdas in cone around jet axis, particles are from fragmentation but also from underlying event

  Bool_t IsArmenterosSelected;                             //Armenteros-Podolanski Cut (is/isn't) applied  
 
  AliFragFuncHistosInvMass*  fFFHistosIMALaAllEvt;          //! ALa pt spec for all events
  AliFragFuncHistosInvMass*  fFFHistosIMALaJet;             //! ALa FF all dPhi   
  AliFragFuncHistosInvMass*  fFFHistosIMALaCone;            //! ALa FF jet cone   
  AliFragFuncHistosPhiCorrInvMass*  fFFHistosPhiCorrIMALa;  //! ALa corelation to jet axis 
  
  // histogram bins 
  


  //--K0s 
  
  Int_t   fFFIMNBinsJetPt;    // FF histos bins
  Float_t fFFIMJetPtMin;      // FF histos limits
  Float_t fFFIMJetPtMax;      // FF histos limits
  
  Int_t   fFFIMNBinsInvM;     // FF histos bins
  Float_t fFFIMInvMMin;       // FF histos bins
  Float_t fFFIMInvMMax;       // FF histos bins
  
  Int_t   fFFIMNBinsPt;       // FF histos bins
  Float_t fFFIMPtMin;         // FF histos limits
  Float_t fFFIMPtMax;         // FF histos limits
  
  Int_t   fFFIMNBinsXi;       // FF histos bins
  Float_t fFFIMXiMin;         // FF histos limits
  Float_t fFFIMXiMax;         // FF histos limits
  
  Int_t   fFFIMNBinsZ;        // FF histos bins
  Float_t fFFIMZMin;          // FF histos limits
  Float_t fFFIMZMax;          // FF histos limits
  
  //--La
  
  Int_t   fFFIMLaNBinsJetPt;    // FF histos bins
  Float_t fFFIMLaJetPtMin;      // FF histos limits
  Float_t fFFIMLaJetPtMax;      // FF histos limits
  
  Int_t   fFFIMLaNBinsInvM;     // FF histos bins
  Float_t fFFIMLaInvMMin;       // FF histos bins
  Float_t fFFIMLaInvMMax;       // FF histos bins
  
  Int_t   fFFIMLaNBinsPt;       // FF histos bins
  Float_t fFFIMLaPtMin;         // FF histos limits
  Float_t fFFIMLaPtMax;         // FF histos limits
  
  Int_t   fFFIMLaNBinsXi;       // FF histos bins
  Float_t fFFIMLaXiMin;         // FF histos limits
  Float_t fFFIMLaXiMax;         // FF histos limits
  
  Int_t   fFFIMLaNBinsZ;        // FF histos bins
  Float_t fFFIMLaZMin;          // FF histos limits
  Float_t fFFIMLaZMax;          // FF histos limits
  

  //--K0s
  
  Int_t fPhiCorrIMNBinsPt;    // FF histos bins
  Float_t fPhiCorrIMPtMin;    // FF histos limits
  Float_t fPhiCorrIMPtMax;    // FF histos limits
  
  Int_t fPhiCorrIMNBinsPhi;   // FF histos bins
  Float_t fPhiCorrIMPhiMin;   // FF histos limits
  Float_t fPhiCorrIMPhiMax;   // FF histos limits
  
  Int_t fPhiCorrIMNBinsInvM;  // FF histos bins
  Float_t fPhiCorrIMInvMMin;  // FF histos limits
  Float_t fPhiCorrIMInvMMax;  // FF histos limits
  
  //--La

  Int_t fPhiCorrIMLaNBinsPt;    // FF histos bins
  Float_t fPhiCorrIMLaPtMin;    // FF histos limits
  Float_t fPhiCorrIMLaPtMax;    // FF histos limits

  Int_t fPhiCorrIMLaNBinsPhi;   // FF histos bins
  Float_t fPhiCorrIMLaPhiMin;   // FF histos limits
  Float_t fPhiCorrIMLaPhiMax;   // FF histos limits
		
  Int_t   fPhiCorrIMLaNBinsInvM;  // FF histos bins
  Float_t fPhiCorrIMLaInvMMin;  // FF histos limits
  Float_t fPhiCorrIMLaInvMMax;  // FF histos limits


  // Histograms
  
  TH1F* fh1EvtAllCent; 
  TH1F* fh1Evt;                      
  TH1F* fh1K0Mult;                   
  TH1F* fh1dPhiJetK0;                
  TH1F* fh1LaMult;                   
  TH1F* fh1dPhiJetLa;                
  TH1F* fh1ALaMult;                  
  TH1F* fh1dPhiJetALa; 
  TH1F* fh1JetEta;        
  TH1F* fh1JetPhi;                   
  TH2F* fh2JetEtaPhi;  
  TH1F* fh1V0JetPt; 
  TH2F* fh2FFJetTrackEta; //charged jet track eta distribution                 
  TH1F* fh1trackPosNCls;             
  TH1F* fh1trackNegNCls; 
  TH1F* fh1trackPosRap;              
  TH1F* fh1trackNegRap;              
  TH1F* fh1V0Rap;              
  TH1F* fh1trackPosEta;              
  TH1F* fh1trackNegEta;              
  TH1F* fh1V0Eta;                    
  TH1F* fh1V0totMom;                 
  TH1F* fh1CosPointAngle;            
  TH1F* fh1Chi2Pos;                  
  TH1F* fh1Chi2Neg;                  
  TH1F* fh1DecayLengthV0;            
  TH2F* fh2ProperLifetimeK0sVsPtBeforeCut;
  TH2F* fh2ProperLifetimeK0sVsPtAfterCut;
  TH1F* fh1ProperLifetimeV0BeforeCut;
  TH1F* fh1ProperLifetimeV0AfterCut; 
  TH1F* fh1V0Radius;                 
  TH1F* fh1DcaV0Daughters;           
  TH1F* fh1DcaPosToPrimVertex;       
  TH1F* fh1DcaNegToPrimVertex;        
  TH2F* fh2ArmenterosBeforeCuts;     
  TH2F* fh2ArmenterosAfterCuts;      
  TH2F* fh2BB3SigProton;             
  TH2F* fh2BBLaPos;                  
  TH2F* fh2BBLaNeg;                  
  TH1F* fh1CrossedRowsOverFindableNeg;
  TH1F* fh1CrossedRowsOverFindablePos;
  TH1F* fh1PosDaughterCharge;
  TH1F* fh1NegDaughterCharge;
  TH1F* fh1PtMCK0s;
  TH1F* fh1PtMCLa;
  TH1F* fh1PtMCALa;
  TH1F* fh1EtaK0s;
  TH1F* fh1EtaLa;
  TH1F* fh1EtaALa;  
  TH3F* fh3InvMassEtaTrackPtK0s;
  TH3F* fh3InvMassEtaTrackPtLa;
  TH3F* fh3InvMassEtaTrackPtALa;
  TH1F* fh1noAssociatedK0s;
  TH1F* fh1TrackMultCone;
  TH2F* fh2TrackMultCone;
  TH2F* fh2MCgenK0Cone;
  TH2F* fh2MCgenLaCone;
  TH2F* fh2MCgenALaCone;
  TH2F* fh2MCEtagenK0Cone;
  TH2F* fh2MCEtagenLaCone;
  TH2F* fh2MCEtagenALaCone;
  TH1F* fh1FFIMK0ConeSmear;
  TH1F* fh1FFIMLaConeSmear;
  TH1F* fh1FFIMALaConeSmear;
  TH3F* fh3MCrecK0Cone;
  TH3F* fh3MCrecLaCone;
  TH3F* fh3MCrecALaCone;
  TH3F* fh3MCrecK0ConeSmear;
  TH3F* fh3MCrecLaConeSmear;
  TH3F* fh3MCrecALaConeSmear;
  TH3F* fh3SecContinCone;
  TH3F* fh3StrContinCone;
  TH3F* fh3IMK0PerpCone;
  TH3F* fh3IMLaPerpCone;
  TH3F* fh3IMALaPerpCone;
  TH3F* fh3IMK0MedianCone;
  TH3F* fh3IMLaMedianCone;
  TH3F* fh3IMALaMedianCone;
  TH1F* fh1MCMultiplicityPrimary;
  TH1F* fh1MCMultiplicityTracks;
  TH1F* fh1MCmotherLa;
  TH1F* fh1MCmotherALa;
  TH3F* fh3FeedDownLa;
  TH3F* fh3FeedDownALa;     
  TH1F* fh1MCProdRadiusK0s;
  TH1F* fh1MCProdRadiusLambda;
  TH1F* fh1MCProdRadiusAntiLambda;
  TH1F* fh1MCPtV0s;
  TH1F* fh1MCPtK0s; 
  TH1F* fh1MCPtLambda; 
  TH1F* fh1MCPtAntiLambda;
  TH1F* fh1MCXiPt;
  TH1F* fh1MCXibarPt;
  TH2F* fh2MCEtaVsPtK0s;
  TH2F* fh2MCEtaVsPtLa;
  TH2F* fh2MCEtaVsPtALa;
  TH1F* fh1MCRapK0s; 
  TH1F* fh1MCRapLambda;
  TH1F* fh1MCRapAntiLambda;
  TH1F* fh1MCEtaAllK0s; 
  TH1F* fh1MCEtaK0s; 
  TH1F* fh1MCEtaLambda;
  TH1F* fh1MCEtaAntiLambda;



  ClassDef(AliAnalysisTaskJetChem, 3);
};

#endif


