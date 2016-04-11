#ifndef ALIANALYSISTASKCHARGEDJETSHADRONTOY_H
#define ALIANALYSISTASKCHARGEDJETSHADRONTOY_H
/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

class TClonesArray;
class TString;
class TRandom3;

#include "AliAnalysisTaskSE.h"

// Toy model to create an event for charged jet-hadron correlations

class AliAnalysisTaskChargedJetsHadronToy : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskChargedJetsHadronToy();
  virtual ~AliAnalysisTaskChargedJetsHadronToy();

  virtual void  UserCreateOutputObjects();
  virtual void  UserExec(Option_t *);
  virtual void  Terminate(Option_t *) {}
  void          ExecOnce();

  // ### SETTERS/GETTERS
  void                        SetCreateUE(Bool_t val)                 {fCreateUE = val;}
  void                        SetCreateJets(Bool_t val)               {fCreateJets = val;}
  void                        SetUEMultiplicityDistribution(TH1* val) {fUEMultDistribution = val;}
  void                        SetUEDistribution(TF1* val)             {fUEDistribution = val;}
  void                        SetUEMultiplicity(Int_t val)            {fUEMultiplicity = val;}
  void                        SetGeneratedJetParticleDistribution(TF1* val)   {fGeneratedJetParticleDistribution = val;}
  void                        SetGeneratedJetPtDistribution(TH1* val) {fGeneratedJetPtDistribution = val;}
  void                        SetGeneratedJetCount(Int_t val)         {fGeneratedJetCount = val;}
  void                        SetGeneratedJetPtRange(Double_t min, Double_t max) {fGeneratedJetPtMin = min; fGeneratedJetPtMax = max;}
  void                        SetGeneratedJetWidthPhi(Double_t val)   {fGeneratedJetWidthPhi = val;}
  void                        SetGeneratedJetWidthEta(Double_t val)   {fGeneratedJetWidthEta = val;}
  void                        SetGeneratedJetMinEta(Double_t val)     {fGeneratedJetMinEta = val;}
  void                        SetGeneratedJetMaxEta(Double_t val)     {fGeneratedJetMaxEta = val;}

  void                        SetInputTracksName(const char* val)     {fInputArrTracksName = val;}
  void                        SetOutputTracksName(const char* val)    {fOutputArrTracksName = val;}
  void                        SetGeneratedJetsName(const char* val)   {fGeneratedJetsArrName = val;}

private:
  // ### Settings
  Bool_t                      fCreateUE;                          // create UE in toymodel
  Bool_t                      fCreateJets;                        // create jets in toymodel
  TH1*                        fUEMultDistribution;                // histogram for multiplicity distribution
  TF1*                        fUEDistribution;                    // function for particle pt distribution
  Int_t                       fUEMultiplicity;                    // multiplicity in UE
  TF1*                        fGeneratedJetParticleDistribution;  // function for particle pt distribution in jets
  TH1*                        fGeneratedJetPtDistribution;        // pt distribution used to produce jets
  Int_t                       fGeneratedJetCount;                 // count of generated jets

  Double_t                    fGeneratedJetPtMin;                 // generated jets pT (min)
  Double_t                    fGeneratedJetPtMax;                 // generated jets pT (max)
  Double_t                    fGeneratedJetWidthPhi;              // width of generated jets in phi
  Double_t                    fGeneratedJetWidthEta;              // width of generated jets in eta
  Double_t                    fGeneratedJetMinEta;                // generated jets min eta
  Double_t                    fGeneratedJetMaxEta;                // generated jets max eta

  // ### Input/output settings+arrays
  TClonesArray*               fInputArrTracks;                    //! input array containing tracks from events
  TString                     fInputArrTracksName;                // Name of the TClonesArray that will be loaded
  TClonesArray*               fOutputArrTracks;                   //! array holding tracks from toy model
  TString                     fOutputArrTracksName;               // Name of the destination TClonesArray
  TClonesArray*               fGeneratedJetsArr;                  //! array holding generated jets from toy model
  TString                     fGeneratedJetsArrName;              // Name of the destination TClonesArray

  // ### Misc
  TF1*                        fDistEtaGaussian;                   //! function for gaussian distribution in toy
  TF1*                        fDistPhiGaussian;                   //! function for gaussian distribution in toy
  TRandom3*                   fRandom;                            //! random number generator
  Bool_t                      fInitialized;                       // internal state when ExecOnce has been executed

  void                        AssembleEvent();

  AliAnalysisTaskChargedJetsHadronToy(const AliAnalysisTaskChargedJetsHadronToy&);            // not implemented
  AliAnalysisTaskChargedJetsHadronToy &operator=(const AliAnalysisTaskChargedJetsHadronToy&); // not implemented

  ClassDef(AliAnalysisTaskChargedJetsHadronToy, 1);
};

#endif
