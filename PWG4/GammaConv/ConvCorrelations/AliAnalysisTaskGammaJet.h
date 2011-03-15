/* This file is property of and copyright                                 *
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/// @file   AliAnalysisTaskGammaJet.h
/// @author Svein Lindal
/// @brief  Class used to run isolation studies of conversion gamma / pions
 


#ifndef ALIANALYSISTASKGAMMAJET_H
#define ALIANALYSISTASKGAMMAJET_H



class TH1F;
class AliESDEvent;
class AliGammaConversionAODObject;
class AliAODConversionParticle;
class AliAODPWG4ParticleCorrelation;
class AliAODPWG4Particle;
class TObjArray;
class TString;
class TClonesArray;


#include "AliAnaConvIsolation.h"
#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskGammaJet : public AliAnalysisTaskSE {

public:
 
  AliAnalysisTaskGammaJet(); 
  AliAnalysisTaskGammaJet(const char *name);
  virtual ~AliAnalysisTaskGammaJet();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

  void SetDeltaAODFileName(TString string) { fDeltaAODFileName = string;}
 
  AliAnaConvIsolation * GetIsolation() const {return fAnaIsolation;}
  void SetIsolation( AliAnaConvIsolation * isolation) { fAnaIsolation = isolation; }

  void SetConversionCutId(TString cut) { fConversionCutString = Form("GammaConv_%s", cut.Data());}

  void SetMinPt(Float_t minpt) { fMinPt = minpt;}
  void SetMinNTracks(Int_t nTracks) { fMinNTracks = nTracks; }

  void AddIsolationAna(TObject * isoAna) { fAnaIsolationArray->Add(isoAna);}
  void AddPhotonHadronAna(TObject * ana) { fAnaPhotonArray->Add(ana);}
  void AddPhotonJetAna(TObject * ana) { fAnaPhotonJetArray->Add(ana);}
  void AddPionHadronAna(TObject * ana) { fAnaPionArray->Add(ana);}

  void GetPionGrandChildren(const AliAODConversionParticle * const pion, const TClonesArray * photons, Int_t* trackLabels);

 private:

  void NotifyRun();
  Bool_t UserNotify();


  //Get the AOD event from whereever it might be accessible
  AliAODEvent * GetAODEvent();

  //Get Conversion gammas branch
  TClonesArray * GetConversionGammas(const AliAODEvent * aodEvent);
  TClonesArray * GetConversionPions(const AliAODEvent * aodEvent);

  //Process conv gamma
  void ProcessConvGamma( const TClonesArray * const convGamma, const TClonesArray * const pions, const TClonesArray * const tracks);
  void ProcessPions( const TClonesArray * const pions, const TClonesArray * const photons, const TClonesArray * const tracks);

  //Process calorimeters
  void ProcessCalorimeters( const AliAODEvent * const aodEvent );
  

  //Does any pions have given photon (iPhot) index listed as daughter
  Bool_t IsDecayPion(Int_t iPhot, const TClonesArray * const pions); //see above


  TList       *fOutputList;       //! Output list
  TString     fDeltaAODFileName;  //! File where Gamma Conv AOD is located, if not in default AOD
  TString fConversionCutString;   //! The cut string of the conversion analysis used to produce input AOD


  AliAnaConvIsolation * fAnaIsolation;

  TObjArray * fAnaIsolationArray; //!Array of isolation analysis objects
  TObjArray * fAnaPionArray;      //!Array of pion - hadron ana objects
  TObjArray * fAnaPhotonArray;    //!Array of photon - hadron ana objects
  TObjArray * fAnaPhotonJetArray; //!Array of photon - jet ana objects

  Float_t fMinPt; //Minimum pt for leading particles
  Int_t fMinNTracks; //Minimum number of tracks in event
  
  AliAnalysisTaskGammaJet(const AliAnalysisTaskGammaJet&); // not implemented
  AliAnalysisTaskGammaJet& operator=(const AliAnalysisTaskGammaJet&); // not implemented
  
  ClassDef(AliAnalysisTaskGammaJet, 4); // example of analysis
};

#endif
