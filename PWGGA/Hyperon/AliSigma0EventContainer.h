#ifndef AliSigma0EventContainer_H
#define AliSigma0EventContainer_H

#include "AliSigma0ParticleBase.h"
#include "AliSigma0ParticlePhotonMother.h"
#include "AliSigma0ParticleV0.h"
#include "Riostream.h"
#include "TObject.h"

#include "TH1.h"
#include "TH2.h"
#include "TList.h"

using std::vector;
#include <deque>
using std::deque;

class AliPIDResponse;

class AliSigma0EventContainer : public TObject {
 public:
  AliSigma0EventContainer();
  AliSigma0EventContainer(const AliSigma0EventContainer &);
  AliSigma0EventContainer &operator=(const AliSigma0EventContainer &);
  virtual ~AliSigma0EventContainer();

  void ProcessEvent(AliVEvent *inputEvent, AliMCEvent *mcEvent,
                    vector<AliSigma0ParticleBase> &ParticleContainer,
                    vector<AliSigma0ParticleBase> &AntiParticleContainer,
                    vector<AliSigma0ParticleV0> &V0Container,
                    vector<AliSigma0ParticleV0> &AntiV0Container,
                    vector<AliSigma0ParticlePhotonMother> &SigmaContainer,
                    vector<AliSigma0ParticlePhotonMother> &AntiSigmaContainer,
                    vector<AliAODConversionPhoton> &PhotonContainer);

  int GetMultiplicityITS(AliVEvent *inputEvent);

  void CorrelationFunction();
  void PhiStar();
  void EventMixing();
  void FillEventBuffer();

  void TrackCleaner(vector<AliSigma0ParticleBase> &ParticleContainer,
                    vector<AliSigma0ParticleBase> &AntiParticleContainer,
                    vector<AliSigma0ParticleV0> &V0Container,
                    vector<AliSigma0ParticleV0> &AntiV0Container);

  int GetZvertexBins(float zVertex) const;
  int GetMultiplicityBins(int multiplicity) const;

  void SetIsMC(bool isMC) { fIsMC = isMC; }
  void SetExtendedQA(bool isExtended) { fIsExtendedQA = isExtended; }
  void SetUseEmptyEvents(bool use) { fUseEmptyEvents = use; }
  void SetPhotonMixingDepth(int mix) { fMixingDepthPhoton = mix; }
  void SetProtonMixingDepth(int mixP) { fMixingDepthProton = mixP; }
  void SetLambdaMixingDepth(int mixL) { fMixingDepthLambda = mixL; }
  void SetSigmaMixingDepth(int mixS) { fMixingDepthSigma = mixS; }
  void SetZvertexBins(float max, float min, float step);

  void SetSigmaMass(float mass) { fMassSigma = mass; }
  void SetSigmaMassCut(float cut) { fSigmaMassCut = cut; }
  void SetSigmaSideband(float low, float up) {
    fSigmaSidebandLow = low;
    fSigmaSidebandUp = up;
  }

  void InitCutHistograms();
  TList *GetCutHistograms() const { return fHistograms; }

 protected:
  TList *fHistograms;
  TList *fHistogramsMC;
  TList *fHistogramsExtended;

  bool fIsMC;
  bool fIsExtendedQA;
  bool fUseEmptyEvents;

  AliVEvent *fInputEvent;  //!
  AliMCEvent *fMCEvent;    //!

  int fMixingDepthPhoton;
  int fMixingDepthProton;
  int fMixingDepthLambda;
  int fMixingDepthSigma;

  int fZVertex;
  int fMultiplicity;

  float fMassSigma;
  float fSigmaMassCut;
  float fSigmaSidebandLow;
  float fSigmaSidebandUp;

  //    int                 fNZVertexBins;
  //    vector<int>         fZVertexBins; //

  vector<AliSigma0ParticleBase> fProton;             //!
  vector<AliSigma0ParticleBase> fAntiProton;         //!
  vector<AliSigma0ParticleV0> fLambda;               //!
  vector<AliSigma0ParticleV0> fAntiLambda;           //!
  vector<AliSigma0ParticlePhotonMother> fSigma;      //!
  vector<AliSigma0ParticlePhotonMother> fAntiSigma;  //!
  vector<AliAODConversionPhoton> fPhoton;            //!

  deque<vector<AliSigma0ParticleBase> > fProtonMixed[10][13];      //!
  deque<vector<AliSigma0ParticleBase> > fAntiProtonMixed[10][13];  //!

  deque<vector<AliSigma0ParticleV0> > fLambdaMixed[10][13];      //!
  deque<vector<AliSigma0ParticleV0> > fAntiLambdaMixed[10][13];  //!

  deque<vector<AliSigma0ParticlePhotonMother> > fSigmaMixed[10];      //!
  deque<vector<AliSigma0ParticlePhotonMother> > fAntiSigmaMixed[10];  //!

  deque<vector<AliAODConversionPhoton> > fPhotonMixed[10][13];  //!

  // Histograms
  // =====================================================================
  TH1F *fHistCuts;                                 //
  TH1F *fHistKProtonProton;                        //
  TH1F *fHistKProtonProtonMixed;                   //
  TH1F *fHistKProtonProtonMulti[13];               //
  TH1F *fHistKProtonProtonMixedMulti[13];          //
  TH1F *fHistKAntiProtonAntiProton;                //
  TH1F *fHistKAntiProtonAntiProtonMixed;           //
  TH1F *fHistKAntiProtonAntiProtonMulti[13];       //
  TH1F *fHistKAntiProtonAntiProtonMixedMulti[13];  //
  TH1F *fHistKProtonLambda;                        //
  TH1F *fHistKProtonLambdaMixed;                   //
  TH1F *fHistKProtonLambdaMulti[13];               //
  TH1F *fHistKProtonLambdaMixedMulti[13];          //
  TH1F *fHistKAntiProtonAntiLambda;                //
  TH1F *fHistKAntiProtonAntiLambdaMixed;           //
  TH1F *fHistKAntiProtonAntiLambdaMulti[13];       //
  TH1F *fHistKAntiProtonAntiLambdaMixedMulti[13];  //
  TH1F *fHistKLambdaLambda;                        //
  TH1F *fHistKLambdaLambdaMixed;                   //
  TH1F *fHistKLambdaLambdaMulti[13];               //
  TH1F *fHistKLambdaLambdaMixedMulti[13];          //
  TH1F *fHistKAntiLambdaAntiLambda;                //
  TH1F *fHistKAntiLambdaAntiLambdaMixed;           //
  TH1F *fHistKAntiLambdaAntiLambdaMulti[13];       //
  TH1F *fHistKAntiLambdaAntiLambdaMixedMulti[13];  //
  TH1F *fHistKProtonSigma;                         //
  TH1F *fHistKProtonSigmaMixed;                    //
  TH1F *fHistKProtonSigmaMulti[13];                //
  TH1F *fHistKProtonSigmaMixedMulti[13];           //
  TH1F *fHistKAntiProtonAntiSigma;                 //
  TH1F *fHistKAntiProtonAntiSigmaMixed;            //
  TH1F *fHistKAntiProtonAntiSigmaMulti[13];        //
  TH1F *fHistKAntiProtonAntiSigmaMixedMulti[13];   //

  TH1F *fHistKProtonSigmaSidebandLower;               //
  TH1F *fHistKProtonSigmaMixedSidebandLower;          //
  TH1F *fHistKAntiProtonAntiSigmaSidebandLower;       //
  TH1F *fHistKAntiProtonAntiSigmaMixedSidebandLower;  //
  TH1F *fHistKProtonSigmaSidebandUpper;               //
  TH1F *fHistKProtonSigmaMixedSidebandUpper;          //
  TH1F *fHistKAntiProtonAntiSigmaSidebandUpper;       //
  TH1F *fHistKAntiProtonAntiSigmaMixedSidebandUpper;  //

  TH2F *fHistKProtonSigmaVsLambda;          //
  TH2F *fHistKAntiProtonAntiSigmaVsLambda;  //

  TH1F *fHistKProtonSigmaLambda;               //
  TH1F *fHistKProtonSigmaMixedLambda;          //
  TH1F *fHistKAntiProtonAntiSigmaLambda;       //
  TH1F *fHistKAntiProtonAntiSigmaMixedLambda;  //

  TH2F *fHistMCProtonProtonMomentum;          //
  TH2F *fHistMCAntiProtonAntiProtonMomentum;  //
  TH2F *fHistMCProtonLambdaMomentum;          //
  TH2F *fHistMCAntiProtonAntiLambdaMomentum;  //
  TH2F *fHistMCLambdaLambdaMomentum;          //
  TH2F *fHistMCAntiLambdaAntiLambdaMomentum;  //
  TH2F *fHistMCProtonSigmaMomentum;           //
  TH2F *fHistMCAntiProtonAntiSigmaMomentum;   //

  TH2F *fHistProtonProtonPhiStar[9];                     //
  TH2F *fHistAntiProtonAntiProtonPhiStar[9];             //
  TH2F *fHistProtonLambdaProtonPhiStar[9];               //
  TH2F *fHistAntiProtonAntiLambdaProtonPhiStar[9];       //
  TH2F *fHistProtonLambdaPionPhiStar[9];                 //
  TH2F *fHistAntiProtonAntiLambdaPionPhiStar[9];         //
  TH2F *fHistProtonProtonPhiStarMixed[9];                //
  TH2F *fHistAntiProtonAntiProtonPhiStarMixed[9];        //
  TH2F *fHistProtonLambdaProtonPhiStarMixed[9];          //
  TH2F *fHistAntiProtonAntiLambdaProtonPhiStarMixed[9];  //
  TH2F *fHistProtonLambdaPionPhiStarMixed[9];            //
  TH2F *fHistAntiProtonAntiLambdaPionPhiStarMixed[9];    //

  TH1F *fHistProtonLambdaNShared;          //
  TH1F *fHistProtonAntiLambdaNShared;      //
  TH1F *fHistAntiProtonLambdaNShared;      //
  TH1F *fHistAntiProtonAntiLambdaNShared;  //
  TH1F *fHistLambdaNShared;                //
  TH1F *fHistAntiLambdaNShared;            //

  TH1F *fHistSigmaMixedPt;              //
  TH1F *fHistSigmaMixedInvMass;         //
  TH2F *fHistSigmaMixedInvMassPt;       //
  TH2F *fHistSigmaMixedInvMassEta;      //
  TH1F *fHistAntiSigmaMixedPt;          //
  TH1F *fHistAntiSigmaMixedInvMass;     //
  TH2F *fHistAntiSigmaMixedInvMassPt;   //
  TH2F *fHistAntiSigmaMixedInvMassEta;  //

  TH1F *fHistProtonPt;              //
  TH1F *fHistAntiProtonPt;          //
  TH1F *fHistLambdaPt;              //
  TH1F *fHistLambdaInvMass;         //
  TH1F *fHistLambdaInvMassRec;      //
  TH1F *fHistAntiLambdaPt;          //
  TH1F *fHistAntiLambdaInvMass;     //
  TH1F *fHistAntiLambdaInvMassRec;  //
  TH1F *fHistSigmaPt;               //
  TH1F *fHistSigmaInvMass;          //
  TH1F *fHistSigmaInvMassRec;       //
  TH1F *fHistAntiSigmaPt;           //
  TH1F *fHistAntiSigmaInvMass;      //
  TH1F *fHistAntiSigmaInvMassRec;   //

 private:
  ClassDef(AliSigma0EventContainer, 2)
};

#endif
