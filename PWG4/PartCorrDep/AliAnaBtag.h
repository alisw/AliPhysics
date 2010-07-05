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

//-- Author T.R.P-R.Aronsson, J.Klay

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
  Int_t GetDVMBtag(AliAODTrack * tr); //returns # tracks from secvtx

  //Temporary local method to get DCA because AliAODTrack is stupid
  Bool_t GetDCA(const AliAODTrack* tr,Double_t imp[2], Double_t cov[3]);
  Bool_t PhotonicV0(Int_t trackId); //check with V0 list

  //check if track has been flagged as a non-photonic or DVM electron
  //used with the jet tracks to tag bjets
  Bool_t CheckIfBjet(const AliAODTrack* track);
  Bool_t IsMcBJet(Double_t x, Double_t y);
  Bool_t IsMcDJet(Double_t x, Double_t y);

  void Print(const Option_t * opt)const;
  
  TString GetCalorimeter()   const {return fCalorimeter ; }
  Double_t GetpOverEmin()   const {return fpOverEmin ; }
  Double_t GetpOverEmax()   const {return fpOverEmax ; }
  Bool_t GetWriteNtuple()   const {return fWriteNtuple ; }

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
  Double_t GetMinClusEne() const { return fMinClusEne; }

  void SetCalorimeter(TString det)    {fCalorimeter = det ; }
  void SetpOverEmin(Double_t min)     {fpOverEmin = min ; }
  void SetpOverEmax(Double_t max)     {fpOverEmax = max ; }
  void SetResidualCut(Double_t cut)     {fResidualCut = cut ; }
  void SetWriteNtuple(Bool_t val)     {fWriteNtuple = val ; }

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
  void SetMinClusEne(Double_t ene) { fMinClusEne = ene; }

  void InitParameters();
  void Terminate(TList * outputList);
         
  private:
  //For DVM B-tag method
  Double_t ComputeSignDca(AliAODTrack *track, AliAODTrack *track2, float cut1, double pdcacut);

  //Int_t GetMCSource(Int_t mctag);
  AliAODMCParticle* GetMCParticle(Int_t part);


  private:
  TString  fCalorimeter;  //! Which detector? EMCAL or PHOS
  Double_t fpOverEmin;    //! Minimum p/E value for Electrons
  Double_t fpOverEmax;    //! Maximum p/E value for Electrons
  Double_t fResidualCut;  //! Track-cluster matching distance
  Double_t fMinClusEne;   //! Min clus energy for matching

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

  Bool_t  fWriteNtuple; //flag for filling ntuple or not
  

  ///////////////////////////////////////
  //Output histograms and Ntuples

  TNtuple * fEleNtuple;                 //Ntuple for electrons if neede


  TH1F * fhEmcalElectrons;              //All electrons, as id:ed by EMCAL
  TH1F * fhTRDElectrons;                //Electrons from TRD
  TH1F * fhTPCElectrons;                //Electrons from TPC
  TH1F * fhDVM1;                        //initial b-tag dvm1
  TH1F * fhDVM2;                        //initial b-tag dvm2 
  TH2F * fhJets;                        //All jets 2d
  TH2F * fhJetsAllEtaPhi;               //Eta phi for all jets
  TH2F * fhJetsLeadingBElectronEtaPhi;  //deta dphi btw leading jet and bele
  TH2F * fhDVM1EtaPhi;                  //eta phi for b-electrons
  TH2F * fhBJetElectronDetaDphi;        //eta phi for jet with b-electrons
  TH1F * fhClusterEnergy;               //cluster E of EMCAL
  TH1F * fhTestalle;
  TH1F * fhResidual;                    //Residuals from trackmatching

  //Analysis of electrons
  TH2F * fhElectrons;
  //Analysis for tracks from esd/aod
  TH2F * fhTracks;


  ClassDef(AliAnaBtag,12)

} ;
 

#endif//ALIANABTAG_H



