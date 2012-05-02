#ifndef ALISPDCLUSTTASK_H
#define ALISPDCLUSTTASK_H

///////////////////////////////////////////////////////////////////////////
// Class AliSPDClustTask                                            //
// Analysis task to produce data and MC histos needed for tracklets      //
// dNdEta extraction in multiple bins in one go                          //
// Author:  ruben.shahoyan@cern.ch                                       //
///////////////////////////////////////////////////////////////////////////

class TH1F; 
class TH2F;
class AliESDEvent;
class TList;

class AliMCParticle;
class AliITSMultReconstructor;
class AliESDTrackCuts;

#include "../ITS/AliITSsegmentationSPD.h"
#include "AliAnalysisTaskSE.h"
#include "AliTriggerAnalysis.h" 
#include <TMath.h>

class AliSPDClustTask : public AliAnalysisTaskSE {
 public:
  enum {kHTracklets=0,
	kHPtPion=50,kHPtKaon=51,kHPtProton=52,kHPtK0=53,kHPtLambda0=54,
	kHClusters=100,
	kClTypevsEta=0,kClZ=1,kClEta=2,kClTypevsEtaW=3,kClZW=4,kClEtaW=5,
	kClZPions=6, kClEtaPions=7, kClZPionsW=8, kClEtaPionsW=9,
	kClZKaons=10,kClEtaKaons=11,kClZKaonsW=12,kClEtaKaonsW=13,
	kClZProtons=14,kClEtaProtons=15,kClZProtonsW=16,kClEtaProtonsW=17,
	kClZK0s=18,kClEtaK0s=19,kClZK0sW=20,kClEtaK0sW=21,
	kClZLambda0s=22,kClEtaLambda0s=23,kClZLambda0sW=24,kClEtaLambda0sW=25}; // to facilitated access to histos, see BookHistos
  //
  AliSPDClustTask(const char *name = "AliSPDClustTask");
  virtual ~AliSPDClustTask(); 
  
  virtual void  UserCreateOutputObjects();
  virtual void  UserExec(Option_t *option);
  virtual void  Terminate(Option_t *);
  void       SetUseMC(Bool_t mc = kFALSE)              {fUseMC = mc;}
  TObjArray* BookHistos();
  void       FillHistos(AliStack *stack);
  // RS
  void       SetNStdDev(Float_t f=1.)           {fNStdDev = f<1e-5 ? 1e-5:f;}
  void       SetScaleDThetaBySin2T(Bool_t v=kFALSE) {fScaleDTBySin2T = v;}
  void       SetCutOnDThetaX(Bool_t v=kFALSE)   {fCutOnDThetaX = v;}
  void       SetPhiWindow(float w=0.08)         {fDPhiWindow   = w<1e-5 ? 1e-5:w;}
  void       SetThetaWindow(float w=0.025)      {if (w<0) fCutOnDThetaX=kTRUE; fDThetaWindow = TMath::Abs(w)<1e-5 ? 1e-5:TMath::Abs(w);}
  void       SetPhiShift(float w=0.0045)        {fDPhiShift = w;}
  void       SetPhiOverlapCut(float w=0.005)    {fPhiOverlapCut = w;}
  void       SetZetaOverlapCut(float w=0.05)    {fZetaOverlap = w;}
  void       SetRemoveOverlaps(Bool_t w=kFALSE) {fRemoveOverlaps = w;}
  //
  void       SetDPhiSCut(Float_t c=0.06)        {fDPhiSCut = c;}
  void       SetNStdCut(Float_t c=1.0)          {fNStdCut = c;}
  //
  void       SetEtaCut(Float_t etaCut)          {fEtaMax = TMath::Abs(etaCut); fEtaMin= -fEtaMax;}
  void       SetEtaMin(Float_t etaMin)          {fEtaMin = etaMin;}
  void       SetEtaMax(Float_t etaMax)          {fEtaMax = etaMax;}
  void       SetZVertexMin(Float_t z)           {fZVertexMin = z;}
  void       SetZVertexMax(Float_t z)           {fZVertexMax = z;}
  //
  void       SetInput(const char *filename);
  Int_t      FindMotherParticle(AliStack* stack, Int_t i);
  Double_t   PtWeight(Double_t pt, Int_t pdgCode);
  //
 protected:
  void       InitMultReco();
  void       FillClusterInfo();
  //
 protected:
  TList*       fOutput;                   // output list send on output slot 1 
  //
  Bool_t       fUseMC; 
  TObjArray*   fHistos;                   //! histos array
  Float_t      fVtx[3];                   //! event vertex
  //
  // Settings for the reconstruction
  // tracklet reco settings
  Float_t      fEtaMin;                    // histos filled only for this eta range
  Float_t      fEtaMax;                    // histos filled only for this eta range
  Float_t      fZVertexMin;                // min Z vtx to process
  Float_t      fZVertexMax;                // max Z vtx to process
  //
  Bool_t       fScaleDTBySin2T;            // request dTheta scaling by 1/sin^2(theta)
  Bool_t       fCutOnDThetaX;              // if true, apart from NStdDev cut apply also the cut on dThetaX
  Float_t      fNStdDev;                   // cut on weighted distance
  Float_t      fDPhiWindow;                // max dPhi
  Float_t      fDThetaWindow;              // max dTheta
  Float_t      fDPhiShift;                 // mean bend
  Float_t      fPhiOverlapCut;             // overlaps cut in phi
  Float_t      fZetaOverlap;               // overlaps cut in Z
  Bool_t       fRemoveOverlaps;            // request overlaps removal
  //
  Float_t      fDPhiSCut;                  // cut on signal dphiS
  Float_t      fNStdCut;                   // cut on signal weighted distance
  //
  AliITSMultReconstructor *fMultReco;              //! mult.reco object
  TTree*       fRPTree;                    //! tree of recpoints
  AliStack*    fStack;                     //! MC stack
  AliMCEvent*  fMCEvent;                   //! MC Event
  Float_t      fESDVtx[3];                 //  ESD vertex
  //
  Bool_t fDontMerge;                       // no merging requested
  //
  TH1F*  fhPtPionIn; // Input histogram containing the pion spectra weight
  TH1F*  fhPtKaonIn; // Input histogram containing the pion spectra weight
  TH1F*  fhPtProtonIn; // Input histogram containing the pion spectra weight
 private:    
  AliSPDClustTask(const AliSPDClustTask&); // not implemented
  AliSPDClustTask& operator=(const AliSPDClustTask&); // not implemented 
  
  ClassDef(AliSPDClustTask, 1);  
};


#endif
