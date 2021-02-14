//--- Task for investigation of the DoubleHyperHydrogen4 ---
//---     Author: Janik Ditzel; janik.ditzel@cern.ch     ---

#ifndef ALIANALYSISTASKDOUBLEHYPNUCTREE_H
#define ALIANALYSISTASKDOUBLEHYPNUCTREE_H

class TH1F;
class TH2F;
class AliESDEvent;
class AliESDpid;
class AliESDtrackCuts;
class AliESDv0;
class AliESDVertex;
class AliESDInputHandler;
class AliESDtrack;

#include "AliAnalysisTaskSE.h"
#include "AliStack.h"
#include "AliEventCuts.h"
#include "TObject.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include <TClonesArray.h>
#include "AliPID.h"
#include "AliVertexerTracks.h"
#include "AliMCParticle.h"


class AliAnalysisTaskDoubleHypNucTree : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskDoubleHypNucTree();
  AliAnalysisTaskDoubleHypNucTree(const char *name);
  virtual ~AliAnalysisTaskDoubleHypNucTree();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(const Option_t*);
  void SelectPIDcheckOnly(Bool_t pidch = kFALSE) {fPIDCheckOnly = pidch;};
  void SetTriggerMask(UInt_t triggerMask = AliVEvent::kAny) {fTriggerMask = triggerMask;};

enum PdgCodeType_t {
    kPDGPionPlus,
    kPDGPionMinus,
    kPDGProton,
    kPDGAntiProton,    
    kPDGDeuteron,
    kPDGAntiDeuteron,
    kPDGTriton,
    kPDGAntiTriton,
    kPDGHelium3,
    kPDGAntiHelium3,
    kPDGHelium4,
    kPDGAntiHelium4,
    kPDGHyperHydrogen3,
    kPDGAntiHyperHydrogen3,
    kPDGHyperHydrogen4,
    kPDGAntiHyperHydrogen4,
    kPDGHyperHelium4,
    kPDGAntiHyperHelium4,
    kPDGHyperHelium5,
    kPDGAntiHyperHelium5,
    kPDGDoubleHyperHydrogen4,
    kPDGAntiDoubleHyperHydrogen4,
    kPdgCodeMax  // Number of enum entries                                                                                                                                          
  };

  static const Int_t fgkPdgCode[];

 private:
  //--General--
  AliESDInputHandler    *fInputHandler;        
  AliESDpid             *fPID;                 
  AliESDEvent           *fESDevent;     
  AliMCEvent            *mcEvent;
  AliStack              *fStack;     
  AliEventCuts          fEventCuts; 
  AliESDtrackCuts       *trackCutsNuclei;
  AliESDtrackCuts       *trackCutsStrong;
  AliESDtrackCuts       *trackCutsSoft;
  TTree                 *fTree;
  TTree                 *fTreeGen;
  TList                 *fHistogramList; 
  TH2F                  *fHistdEdx;            
  TH1F                  *fHistNumEvents;       
  TH1F			*fHistTrigger;	 	 
  Double_t              fBetheParamsHe[6];  
  Double_t              fBetheParamsT[6];
  Bool_t                fPIDCheckOnly;  
  Int_t                 fMCtrue; 
  UInt_t		fTriggerMask; 
  Int_t                 fvariante;
  //--Trigger--
  Int_t                 MB, HMV0, HMSPD, HNU, HQU;
  //--Tracks-- 
  AliESDtrack           *track1;
  AliESDtrack           *track2;
  AliESDtrack           *track3;
  AliESDtrack           *track4;
  Double_t              ptot1, sign1, ptot2, sign2, ptot3, sign3, ptot4, sign4;
  Bool_t                He3Pos1, He3Neg1, pPos2, pNeg2, piPos3, piNeg3, piPos4, piNeg4;
  //--Vertex Reco--
  AliESDVertex          *secVertex;
  AliESDVertex          *tertVertex;
  AliVertexerTracks     *secvertexer;
  AliVertexerTracks     *tertvertexer;   
  AliExternalTrackParam *exTrack;
  TObjArray             *trkArray;
  TObjArray             *trkArray1;
  Double_t              PrimVertex[3];
  Double_t              SecVertex[3];
  Double_t              TertVertex[3];
  Double_t              cov[21];
  Double_t              cov0[21];
  Double_t              cov1[21];
  Double_t              cov2[21];
  UShort_t              id[2];
  Double_t              pxpypz[3];
  Double_t              xyz[3];
  Double_t              dd[3];  
  Short_t               sign; 
  //--MC--
  Int_t                 noCombinatoricsMC;
  Int_t                 lessCombinatoricsMC;
  Int_t                 stackN;
  Int_t                 label1;
  Int_t                 label2;
  Int_t                 label3;
  Int_t                 label4;   
  Int_t                 labelMother1;
  Int_t                 labelMother2;
  Int_t                 labelMother3;
  Int_t                 labelMother4;  
  Int_t                 labelGrandMother1;
  Int_t                 labelGrandMother2;
  Int_t                 labelGrandMother3;
  Int_t                 labelGrandMother4;
  Long_t                PDGCodeMother;
  AliMCParticle         *FirstDaughter;
  AliMCParticle         *SecondDaughter;
  AliMCParticle         *ThirdDaughter;
  AliMCParticle         *FourthDaughter;
  AliMCParticle         *ParticleMother;
  AliMCParticle         *ParticleMother1;
  AliMCParticle         *ParticleMother2;
  AliMCParticle         *ParticleMother3;
  AliMCParticle         *ParticleMother4;
  AliMCParticle         *ParticleGrandMother1;
  AliMCParticle         *ParticleGrandMother2;
  AliMCParticle         *ParticleGrandMother3;
  AliMCParticle         *ParticleGrandMother4;
  //--Particle--
  TLorentzVector        *sublorentzsum;
  TLorentzVector        *particle1;
  TLorentzVector        *particle2;
  TLorentzVector        *particle3;
  TLorentzVector        *particle4;
  TLorentzVector        *lorentzsum;
  TVector3              *h;
  Double_t              xthiss; 
  Double_t              xpp; 
  //saved in Tree
  Int_t		        fPeriod; 
  Int_t                 fMagneticField;
  Int_t                 fCentrality;
  Int_t                 frunnumber;
  Int_t                 mctruth;
  Int_t                 feventclass;
  Int_t                 fTrigMB, fTrigHMV0, fTrigHMSPD, fTrigHNU, fTrigHQU;
  Int_t                 fPDGMother, fChargeMother;
  Int_t                 fNclsDaughter, fNclsITSDaughter, fNclsDaughter1, fNclsITSDaughter1, fNclsDaughter2, fNclsITSDaughter2, fNclsDaughter3, fNclsITSDaughter3;
  Float_t               fEDaughter, fpDaughter, fptDaughter,  fpxDaughter,  fpyDaughter,  fpzDaughter, fyDaughter, fdEdxDaughter, fdEdxSigmaDaughter, fDcaDaughter, fDcazDaughter, fDcaSecDaughter,  fChi2Daughter,  fEtaDaughter, fPhiDaughter, fGeoLengthDaughter, fTOFSignalDaughter;
  Float_t               fEDaughter1, fpDaughter1, fptDaughter1, fpxDaughter1, fpyDaughter1, fpzDaughter1, fyDaughter1, fdEdxDaughter1, fdEdxSigmaDaughter1, fDcaDaughter1, fDcaSecDaughter1, fDcazDaughter1, fChi2Daughter1, fEtaDaughter1, fPhiDaughter1, fGeoLengthDaughter1, fTOFSignalDaughter1;
  Float_t               fEDaughter2, fpDaughter2, fptDaughter2, fpxDaughter2, fpyDaughter2, fpzDaughter2, fyDaughter2, fdEdxDaughter2, fdEdxSigmaDaughter2, fDcaDaughter2, fDcaSecDaughter2, fDcazDaughter2, fChi2Daughter2, fEtaDaughter2, fPhiDaughter2, fGeoLengthDaughter2, fTOFSignalDaughter2;
  Float_t               fEDaughter3, fpDaughter3, fptDaughter3, fpxDaughter3, fpyDaughter3, fpzDaughter3, fyDaughter3, fdEdxDaughter3, fdEdxSigmaDaughter3, fDcaDaughter3, fDcaSecDaughter3, fDcazDaughter3, fChi2Daughter3, fEtaDaughter3, fPhiDaughter3, fGeoLengthDaughter3, fTOFSignalDaughter3;
  Float_t               fEDaughter4, fpDaughter4, fptDaughter4, fpxDaughter4, fpyDaughter4, fpzDaughter4, fyDaughter4, fDcaDaughter4, fDcazDaughter4, fDcaSecDaughter4;
  Float_t               fSigmaYXDaughter, fSigmaXYZDaughter, fSigmaZDaughter, fSigmaYXDaughter1, fSigmaXYZDaughter1, fSigmaZDaughter1, fSigmaYXDaughter2, fSigmaXYZDaughter2, fSigmaZDaughter2, fSigmaYXDaughter3, fSigmaXYZDaughter3, fSigmaZDaughter3, fSigmaYXDaughter4, fSigmaXYZDaughter4, fSigmaZDaughter4, fPtUncertDaughter, fPtUncertDaughter1, fPtUncertDaughter2, fPtUncertDaughter3, fPtUncertDaughter4;
  Float_t               fDCA2B, fDCA3B1, fDCA3B2, fDCA3B3, fPA, fSubPA, fSubPA2, fDecAngle, farmalpha, farmpt;
  Float_t               fPrimVertexX, fPrimVertexY, fPrimVertexZ,fSecVertexX, fSecVertexY, fSecVertexZ, fTertVertexX, fTertVertexY, fTertVertexZ;
  Float_t               fmMother, fEMother, fpxMother, fpyMother, fpzMother, fptMother, fpMother, fyMother, fctMother;
  Float_t               fmSubMother, fESubMother, fpxSubMother, fpySubMother, fpzSubMother, fptSubMother, fpSubMother, fySubMother, fctSubMother;
  Float_t               fthetaP, fthetaN;
  

  Bool_t TriggerSelection();

  Float_t GetInvPtDevFromBC(Int_t b, Int_t c);
  
  Double_t Bethe(const AliESDtrack& track, Double_t mass, Int_t charge, Double_t* params);
  
  Double_t GeoLength(const AliESDtrack& track);
  
  Double_t GetTOFSignal(const AliESDtrack& track);
  
  void dEdxCheck();
  
  void AnalysisLoop();
  void CreateSecVertex();
  void SetDaughterInformation();

  void SetBetheBlochParams(Int_t runNumber);

  void MCGenerated();                                      
  
  void MCFourBodyDecay (Int_t stackN, AliMCParticle *tparticleMother, Long_t PDGMother, Long_t PDGFirstDaughter, Long_t PDGSecondDaughter, Long_t PDGThirdDaughter, Long_t PDGFourthDaughter, Double_t massFirstDaughter, Double_t massSecondDaughter, Double_t massThirdDaughter, Double_t massFourthDaughter); 
  
  Int_t GetLabel(Int_t labelFirstMother, Int_t particlePdgCode, Int_t motherparticlePdgCode);
  Int_t GetLabel(Int_t labelFirstMother, Int_t particlePdgCode);

  void ResetVals(TString mode);


  AliAnalysisTaskDoubleHypNucTree(const AliAnalysisTaskDoubleHypNucTree&);
  AliAnalysisTaskDoubleHypNucTree &operator=(const AliAnalysisTaskDoubleHypNucTree&);

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskDoubleHypNucTree, 3);
  /// \endcond
};
#endif
