#ifndef ALIANALYSISTASKJETCHARGEFLAVOURPB_H
#define ALIANALYSISTASKJETCHARGEFLAVOURPB_H


class TTree;
class TChain;
class AliAnalysisManager;

#include "AliAnalysisTaskEmcalJet.h"


class AliAnalysisTaskJetChargeFlavourPb : public AliAnalysisTaskEmcalJet {
 public:


  AliAnalysisTaskJetChargeFlavourPb();
  AliAnalysisTaskJetChargeFlavourPb(const char *name);
  virtual ~AliAnalysisTaskJetChargeFlavourPb();

  void                                UserCreateOutputObjects();
  void                                Terminate(Option_t *option);

  //Setters
  void SetJetContainer(Int_t c)            { fContainer = c; }
  void SetJetPtMinThreshold(Float_t f)     { fPtThreshold = f; }
  void SetCentralitySelectionOn(Bool_t t)  { fCentSelectOn = t; }
  void SetMinCentrality(Float_t t)         { fCentMin = t; }
  void SetMaxCentrality(Float_t t)         { fCentMax = t; }
  void SetJetRadius(Float_t t)            { fJetRadius = t; }
  void SetK(Float_t t)                    { JetChargeK = t; }

 protected:
  // Obligatory Functions
  Bool_t                              RetrieveEventObjects();
  Bool_t                              Run();
  Bool_t                              FillHistograms();
  // Essential variables
  // Jet container to be analysed: 0 for raw, 1 for subtracted
  Int_t                               fContainer;                           //
  // Jet lower Pt threshold
  Float_t                             fPtThreshold;                         //
  // Switch on/off centrality selection
  Bool_t                              fCentSelectOn;                        //
  // Min centrality value
  Float_t                             fCentMin;                             //
  // Max centrality value
  Float_t                             fCentMax;                             //
  // Jet radius
  Float_t                             fJetRadius;
  // Value of K the scaling factor on jet charge default k = 0.5
  Float_t                             JetChargeK = 0.5;                     //

  // Here is the TTree
  TTree                               *fTreeJets;                                                  //!
  // These are the branch variables; there are nBranches of them

  Float_t Pt;                     //
  Float_t Phi;                    //
  Float_t Eta;                    //
  Float_t JetCharge;              //
  Float_t JetChargeMid;              //
  Float_t JetChargeHigh;              //
  Float_t ParticlePt;             //
  Float_t ParticlePhi;            //
  Float_t ParticleEta;            //
  Float_t ParticleJetCharge;      //
  Float_t ParticleJetChargeMid;      //
  Float_t ParticleJetChargeHigh;      //
  Float_t LeadingTrackPt;         //
  Int_t PdgCode;                  //
  Int_t PtMatchedPdgCode;         //
  Int_t GeoMatchedPdgCode;        //
  Float_t ProgenetorFraction;     //

  Int_t nParticles;        //
  Int_t nParticlesParticle;        //

  
  TChain                              *pChain;                                                     //



  ClassDef(AliAnalysisTaskJetChargeFlavourPb, 11)
};
#endif
