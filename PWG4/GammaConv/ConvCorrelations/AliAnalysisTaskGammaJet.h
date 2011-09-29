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
class AliAODConversionPhoton;
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

  void SetGammaCutId(TString cut) { fGammaCutString = Form("GammaConv_%s", cut.Data());}
  void SetPionCutId(TString cut) { fPionCutString = Form("GammaConv_%s", cut.Data());}

  void SetMinPt(Float_t minpt) { fMinPt = minpt;}
  void SetMinNTracks(Int_t nTracks) { fMinNTracks = nTracks; }

  void AddIsolationAna(TObject * isoAna) { fAnaIsolationArray->Add(isoAna);}
  void AddPhotonHadronAna(TObject * ana) { fAnaPhotonArray->Add(ana);}
  void AddPhotonJetAna(TObject * ana) { fAnaPhotonJetArray->Add(ana);}
  void AddPionJetAna(TObject * ana) { fAnaPionJetArray->Add(ana);}
  void AddPionHadronAna(TObject * ana) { fAnaPionArray->Add(ana);}

  void GetPionGrandChildren(const AliAODConversionPhoton * const pion, const TClonesArray * photons, Int_t* trackLabels);
  void SetEtaLimits(Float_t eta) { fEtaLimit = eta; }


 private:

  void NotifyRun();
  Bool_t UserNotify();
  Bool_t EventIsSynced(const TClonesArray * const tracks, const TClonesArray * const convGamma, const TClonesArray * const pions);
  Bool_t BothTracksPresent(const AliAODConversionPhoton * const photon, const TClonesArray * const tracks) const;
  Bool_t BothGammaPresent(const AliAODConversionPhoton * const pion, const TClonesArray * const photons, const TClonesArray * const tracks) const;
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
  Float_t     fEtaLimit;
  TString     fDeltaAODFileName;  //! File where Gamma Conv AOD is located, if not in default AOD
  TString     fGammaCutString;   //! The cut string of the conversion analysis used to produce input AOD
  TString     fPionCutString;   //! The cut string of the conversion analysis used to produce input AOD


  AliAnaConvIsolation * fAnaIsolation;

  TObjArray * fAnaIsolationArray; //!Array of isolation analysis objects
  TObjArray * fAnaPionArray;      //!Array of pion - hadron ana objects
  TObjArray * fAnaPhotonArray;    //!Array of photon - hadron ana objects
  TObjArray * fAnaPhotonJetArray; //!Array of photon - jet ana objects
  TObjArray * fAnaPionJetArray; //!Array of photon - jet ana objects

  Float_t fMinPt; //Minimum pt for leading particles
  Int_t fMinNTracks; //Minimum number of tracks in event
  
  TH2F * fhTracksMissingPt[2];
 


  
  AliAnalysisTaskGammaJet(const AliAnalysisTaskGammaJet&); // not implemented
  AliAnalysisTaskGammaJet& operator=(const AliAnalysisTaskGammaJet&); // not implemented
  
  ClassDef(AliAnalysisTaskGammaJet, 4); // example of analysis
};

#endif
