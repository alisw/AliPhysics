#ifndef ALIANABTAG_H
#define ALIANABTAG_H
/* Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id:  $ */

//_________________________________________________________________________
//
// Class for the electron identification and B-tagging.
// Clusters from EMCAL matched to tracks are selected 
// and kept in the AOD. Few histograms produced.
//

//-- Author T.R.P-R.Aronsson

// --- ROOT system ---
class TH2F ;
class TString ;
class TNtuple ;
class TH3F;

// --- ANALYSIS system ---
#include "AliAnaPartCorrBaseClass.h"

class AliAODMCParticle;
class AliCaloTrackReader;
class AliAODTrack;
class TList ;

class AliAnaBtag : public AliAnaPartCorrBaseClass {

public: 
  AliAnaBtag() ; // default ctor
  virtual ~AliAnaBtag() ; //virtual dtor
private:
	AliAnaBtag(const AliAnaBtag & g) ; // cpy ctor
	AliAnaBtag & operator = (const AliAnaBtag & g) ;//cpy assignment
	
public:
	
  TList *  GetCreateOutputObjects();

  //Main functions
  void Init();
  void MakeAnalysisFillAOD()  ;
  void MakeAnalysisFillHistograms() ; 
  
  //B-tagging
  Int_t GetDVMBtag(AliAODTrack * tr, Int_t &pairs, Int_t &start, Int_t &stop);//Main tagger

  //Temporary local method to get DCA
  Bool_t GetDCA(const AliAODTrack* tr,Double_t imp[2], Double_t cov[3]);
  Bool_t PhotonicV0(Int_t trackId); //check with V0 list

  //used with the jet tracks to tag bjets
  Bool_t CheckIfBjet(const AliAODTrack* track);
  Bool_t IsMcBJet(Double_t x, Double_t y);
  Bool_t IsMcDJet(Double_t x, Double_t y);



  Double_t GetDrCut() const { return fDrCut; }
  Double_t GetPairDcaCut() const { return fPairDcaCut; }
  Double_t GetDecayLenCut() const { return fDecayLenCut; }
  Double_t GetImpactCut() const { return fImpactCut; }
  Double_t GetAssocPtCut() const { return fAssocPtCut; }
  Double_t GetMassCut() const { return fMassCut; }
  Double_t GetSdcaCut() const { return fSdcaCut; }
  Int_t    GetITSCut() const { return fITSCut; }
  Int_t    GetNTagTrackCut() const { return fNTagTrkCut; }
  Double_t GetIPSigCut() const { return fIPSigCut; }


  void SetWriteNtuple(Int_t w) { fWriteNtuple = w; }
  void SetDrCut(Double_t dr)  { fDrCut = dr; }
  void SetPairDcaCut(Double_t pdca) { fPairDcaCut = pdca; }
  void SetDecayLenCut(Double_t dlen) { fDecayLenCut = dlen; }
  void SetImpactCut(Double_t imp) { fImpactCut = imp; }
  void SetAssocPtCut(Double_t pt) { fAssocPtCut = pt; }
  void SetMassCut(Double_t mass) { fMassCut = mass; }
  void SetSdcaCut(Double_t sdca) { fSdcaCut = sdca; }
  void SetITSCut(Int_t its) { fITSCut = its; }
  void SetNTagTrackCut(Int_t ntr) { fNTagTrkCut = ntr; }
  void SetIPSigCut(Double_t ips) { fIPSigCut = ips; }


  void InitParameters();
  void Terminate(TList * outputList);
         
  private:
  //For DVM B-tag method
  Double_t ComputeSignDca(AliAODTrack *tr, AliAODTrack *tr2 , Double_t &masscut, Double_t &pdcacut, Double_t &massphoton, Double_t &decay);

  AliAODMCParticle* GetMCParticle(Int_t part);


  private:
  //NTuples!
  Int_t fWriteNtuple;    //Will always be no, but might be set to yes in config file.
  TNtuple * electrons;   //Electrons
  TNtuple * pairs;       //Pairs to the electrons
  TNtuple * events;
  Int_t fEventNumber;    // For Ntuple to label events (starts at 0)
  Int_t fNElec;
  Int_t fNElecEv;
  Int_t fNPair;

  //DVM B-tagging
  Double_t fDrCut;       //max dR
  Double_t fPairDcaCut;  //max pair-DCA
  Double_t fDecayLenCut; //max 3d-decaylength
  Double_t fImpactCut;   //max track impact param
  Double_t fAssocPtCut;  //min associated pt
  Double_t fMassCut;     //min Minv cut
  Double_t fSdcaCut;     //min signDca
  Int_t   fITSCut;       //min ITS hits (both)
  //IP Sig B-tagging
  Int_t    fNTagTrkCut;  //min number of tracks required for IP sig tag
  Double_t fIPSigCut;    //min IP significance cut

  Double_t fJetEtaCut;   //max eta for jets
  Double_t fJetPhiMin;   //min phi for jets
  Double_t fJetPhiMax;   //max phi for jets

  
  
  //Histograms
  TH1F * fhEmcalElectrons;              //All electrons, as id:ed by EMCAL
  TH1F * fhTRDElectrons;                //Electrons from TRD
  TH1F * fhTPCElectrons;                //Electrons from TPC
  TH1F * fhEmcalMCE;                    //All electrons, as id:ed by EMCAL MC
  TH1F * fhTRDMCE;                      //Electrons from TRD MC
  TH1F * fhTPCMCE;                      //Electrons from TPC MC
  TH1F * fhEmcalMCEFake;                //All electrons, as id:ed by EMCAL MC, fake
  TH1F * fhTRDMCEFake;                  //Electrons from TRD MC, fake
  TH1F * fhTPCMCEFake;                  //Electrons from TPC MC, fake
  TH1F * fhEmcalMCP;                    //Pions, as id:ed by EMCAL
  TH1F * fhTRDMCP;                      //Pions from TRD
  TH1F * fhTPCMCP;                      //Pions from TPC
  TH1F * fhEmcalMCK;                    //Kaons from EMCAL
  TH1F * fhTRDMCK;                      //Kaons from TRD
  TH1F * fhTPCMCK;                      //Kaons from TPC
  TH1F * fhSpecies;                     //PDG of id:ed electrons
  TH1F * fhDVM1;                        //initial b-tag dvm1
  TH1F * fhDVM2;                        //initial b-tag dvm2 
  TH1F * fhNVTX;                        //NVtx of all btags 
  TH1F * fhNVTXMC;                        //NVtx of MC btags 
  TH2F * fhJets;                        //All jets 2d
  TH2F * fhJetsAllEtaPhi;               //Eta phi for all jets
  TH2F * fhJetsLeadingBElectronEtaPhi;  //deta dphi btw leading jet and bele
  TH2F * fhDVM1EtaPhi;                  //eta phi for b-electrons
  TH2F * fhBJetElectronDetaDphi;        //eta phi for jet with b-electrons
  TH2F * fhClusterMap;                  //2D eta-phi of EMCAL clusters
  TH1F * fhClusterEnergy;               //cluster E of EMCAL
  TH1F * fhTestalle;					//Temp histo for EMCAL cluster energy
  TH1F * fhResidual;                    //Residuals from trackmatching
  TH1F * fhPairPt;                      //Pairs


  //Analysis of electrons
  TH2F * fhElectrons;
  //Analysis for tracks from esd/aod
  TH2F * fhTracks;


  ClassDef(AliAnaBtag,2)

} ;
 

#endif//ALIANABTAG_H



