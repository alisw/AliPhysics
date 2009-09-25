#ifndef ALIANAELECTRON_H
#define ALIANAELECTRON_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id:  $ */

//_________________________________________________________________________
//
// Class for the electron identification.
// Clusters from EMCAL matched to tracks are selected 
// and kept in the AOD. Few histograms produced.
//

//-- Author: J.L. Klay (Cal Poly)

// --- ROOT system ---
class TH2F ;
class TString ;
class TNtuple ;
class TH3F;

// --- ANALYSIS system ---
#include "AliAnaPartCorrBaseClass.h"

class AliCaloTrackReader;
class AliAODTrack;
class TList ;

class AliAnaElectron : public AliAnaPartCorrBaseClass {

public: 

  AliAnaElectron() ; // default ctor
  AliAnaElectron(const AliAnaElectron & g) ; // cpy ctor
  AliAnaElectron & operator = (const AliAnaElectron & g) ;//cpy assignment
  virtual ~AliAnaElectron() ; //virtual dtor
  
  TList *  GetCreateOutputObjects();

  void Init();

  void MakeAnalysisFillAOD()  ;
  
  void MakeAnalysisFillHistograms() ; 
  
  //B-tagging
  Int_t GetDVMBtag(AliAODTrack * tr); //returns # tracks from secvtx

  //Temporary local method to get DCA because AliAODTrack is stupid
  Bool_t GetDCA(const AliAODTrack* tr,Double_t imp[2], Double_t cov[3]);

  Bool_t IsItPhotonic(const AliAODPWG4Particle* part); //check with track list
  Bool_t IsItPhotonic2(const AliAODPWG4Particle* part); //check with V0 list

  //check if track has been flagged as a non-photonic or DVM electron
  //used with the jet tracks to tag bjets
  Bool_t CheckTrack(const AliAODTrack* track,const char* type);  

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

  void InitParameters();

  void Terminate(TList * outputList);
  void ReadHistograms(TList * outputList); //Fill histograms with
					   //histograms in ouput list,
					   //needed in Terminate.            
  private:
  //For DVM B-tag method
  Double_t ComputeSignDca(AliAODTrack *track, AliAODTrack *track2 , float cut1);
  //the 2 following functions are internal methods of the b-tagging
  //based on transverse impact parameter
  Double_t GetIPSignificance(AliAODTrack *tr, Double_t jetPhi);
  void GetImpactParamVect(Double_t Pxy[2], Double_t t[2], Double_t Vxy[2], Double_t ip[2]);

  private:
  TString  fCalorimeter;  //! Which detector? EMCAL or PHOS
  Double_t fpOverEmin;    //! Minimum p/E value for Electrons
  Double_t fpOverEmax;    //! Maximum p/E value for Electrons
  Double_t fResidualCut;  //! Track-cluster matching distance

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

  Bool_t  fWriteNtuple; //flag for filling ntuple or not

  TNtuple* fEleNtuple; //! testing ntuple

  //matching checks   
  TH1F *fh1pOverE;     //! p/E for track-cluster matches
  TH1F *fh1dR;         //! distance between projected track and cluster
  TH2F *fh2EledEdx;    //! dE/dx vs. momentum for electron candidates
  TH2F *fh2MatchdEdx;  //! dE/dx vs. momentum for all matches
  TH2F *fh2dEtadPhi;   //! DeltaEta vs. DeltaPhi of all track/cluster
		       //! pairs
  TH2F *fh2dEtadPhiMatched;   //! DeltaEta vs. DeltaPhi of matched
				//! track/cluster pairs
  TH2F *fh2dEtadPhiUnmatched;   //! DeltaEta vs. DeltaPhi of unmatched track/cluster pairs

  TH2F* fh2TrackPVsClusterE;     //!track momentum vs. cluster energy
  TH2F* fh2TrackPtVsClusterE;    //!track pt vs. cluster energy
  TH2F* fh2TrackPhiVsClusterPhi; //!track phi vs. cluster phi
  TH2F* fh2TrackEtaVsClusterEta; //!track eta vs. cluster eta

  //Photonic Electron checks
  TH1F* fh1OpeningAngle; //!opening angle between pairs of photon candidates
  TH1F* fh1MinvPhoton;   //!invariant mass distribution of electron pairs

  //Reconstructed
  TH1F * fhPtElectron;  //! Number of identified electron vs transverse momentum 
  TH2F * fhPhiElectron; //! Azimuthal angle of identified  electron vs transverse momentum 
  TH2F * fhEtaElectron; //! Pseudorapidity of identified  electron vs tranvserse momentum 

  TH1F * fhPtNPE;  //! Number of non-photonic electron vs transverse momentum 
  TH2F * fhPhiNPE; //! Azimuthal angle of non-photonic electron vs transverse momentum 
  TH2F * fhEtaNPE; //! Pseudorapidity of non-photonic electron vs tranvserse momentum 

  TH1F * fhPtPE;  //! Number of photonic electron vs transverse momentum 
  TH2F * fhPhiPE; //! Azimuthal angle of photonic electron vs transverse momentum 
  TH2F * fhEtaPE; //! Pseudorapidity of photonic electron vs tranvserse momentum 

  TH1F * fhPtConversion;  //! Number of conversion electron vs transverse momentum 
  TH2F * fhPhiConversion; //! Azimuthal angle of conversion  electron vs transverse momentum 
  TH2F * fhEtaConversion; //! Pseudorapidity of conversion electron vs tranvserse momentum 

  TH1F * fhPtBottom;  //! Number of bottom electron vs transverse momentum 
  TH2F * fhPhiBottom; //! Azimuthal angle of bottom  electron vs transverse momentum 
  TH2F * fhEtaBottom; //! Pseudorapidity of bottom electron vs tranvserse momentum 

  TH1F * fhPtCharm;  //! Number of charm electron vs transverse momentum 
  TH2F * fhPhiCharm; //! Azimuthal angle of charm  electron vs transverse momentum 
  TH2F * fhEtaCharm; //! Pseudorapidity of charm electron vs tranvserse momentum 

  TH1F * fhPtCFromB;  //! Number of charm from bottom electron vs transverse momentum 
  TH2F * fhPhiCFromB; //! Azimuthal angle of charm from bottom electron vs transverse momentum 
  TH2F * fhEtaCFromB; //! Pseudorapidity of charm from bottom electron vs tranvserse momentum 

  TH1F * fhPtDalitz;  //! Number of dalitz electron vs transverse momentum 
  TH2F * fhPhiDalitz; //! Azimuthal angle of dalitz  electron vs transverse momentum 
  TH2F * fhEtaDalitz; //! Pseudorapidity of dalitz electron vs tranvserse momentum 

  TH1F * fhPtWDecay;  //! Number of W-boson electron vs transverse momentum 
  TH2F * fhPhiWDecay; //! Azimuthal angle of W-boson  electron vs transverse momentum 
  TH2F * fhEtaWDecay; //! Pseudorapidity of W-boson electron vs tranvserse momentum 
  		
  TH1F * fhPtZDecay;  //! Number of Z-boson electron vs transverse momentum 
  TH2F * fhPhiZDecay; //! Azimuthal angle of Z-boson  electron vs transverse momentum 
  TH2F * fhEtaZDecay; //! Pseudorapidity of Z-boson electron vs tranvserse momentum 

  TH1F * fhPtAll;  //! Number of all electron vs transverse momentum 
  TH2F * fhPhiAll; //! Azimuthal angle of all  electron vs transverse momentum 
  TH2F * fhEtaAll; //! Pseudorapidity of all electron vs tranvserse momentum 

  TH1F * fhPtUnknown;  //! Number of unknown electron vs transverse momentum
  TH2F * fhPhiUnknown; //! Azimuthal angle of unknown  electron vs transverse momentum
  TH2F * fhEtaUnknown; //! Pseudorapidity of unknown electron vs tranvserse momentum

  TH1F * fhPtMisidentified;  //! Number of misidentified electron vs transverse momentum
  TH2F * fhPhiMisidentified; //! Azimuthal angle of misidentified electron vs transverse momentum
  TH2F * fhEtaMisidentified; //! Pseudorapidity of misidentified electron vs tranvserse momentum 

  TH1F* fhPtHadron;     //!Pt distribution of reco charged hadrons
			//!(pi,k,p) in EMCAL acceptance
  TH1F* fhPtEleTrkDet;  //!Pt distribution of reco electrons using
			//!pid info from tracking detectors only in
			//!EMCAL acceptance
  //event QA
  TH1F * fhImpactXY;    //! impact parameter of all tracks to primary vertex
  TH1F * fhRefMult;     //! refmult (sep 14)
  TH1F * fhRefMult2;    //! refmult2 (sep 14)

  //DVM B-tagging
  TH2F * fhDVMBtagCut1; //! DVM B-tagging result for cut1 (minv>1.0)
  TH2F * fhDVMBtagCut2; //! DVM B-tagging result for cut2 (minv>1.5)
  TH2F * fhDVMBtagCut3; //! DVM B-tagging result for cut3 (minv>1.8)
  TH2F * fhDVMBtagQA1;  //! DVM B-tagging : QA of pairDca vs decaylength
  TH2F * fhDVMBtagQA2;  //! DVM B-tagging : QA of signDca vs mass
  TH1F * fhDVMBtagQA3;  //! DVM B-tagging : QA (sep 14)
  TH1F * fhDVMBtagQA4;  //! DVM B-tagging : QA (sep 14)
  TH1F * fhDVMBtagQA5;  //! DVM B-tagging : QA (sep 14)
  //IPSig B-tagging
  TH1F * fhIPSigBtagQA1; //! IPSig B-tagging : QA of # tag tracks
  TH1F * fhIPSigBtagQA2; //! IPSig B-tagging : QA of IP sig

  //B-Jet histograms
  TH2F* fhJetType;       //! How many of each tag were found vs jet pt
  TH2F* fhBJetXsiFF;     //! B-tagged jet FF with xsi = log(pt_Jet/pt_Track)
  TH2F* fhBJetPtFF;      //! B-tagged jet FF with pt_Track
  TH2F* fhBJetEtaPhi;    //! B-tagged jet eta-phi distribution
  TH2F* fhNonBJetXsiFF;  //! Non b-tagged jet FF with xsi = log(pt_Jet/pt_Track)
  TH2F* fhNonBJetPtFF;   //! Non b-tagged jet FF with pt_Track
  TH2F* fhNonBJetEtaPhi; //! Non b-tagged jet eta-phi distribution

  //MC
  TNtuple *fMCEleNtuple; //! Ntuple of MC electrons
  TH1F* fhPtMCHadron;    //! Pt distribution of MC charged hadrons (pi,k,p) in EMCAL acceptance
  TH1F* fhPtMCBottom;    //! Pt distribution of MC bottom electrons in EMCAL
  TH1F* fhPtMCCharm;     //! Pt distribution of MC charm electrons in EMCAL 
  TH1F* fhPtMCCFromB;    //! Pt distribution of MC charm from bottom ele in EMCAL
  TH1F* fhPtMCConversion;//! Pt distribution of MC conversion electrons in EMCAL 
  TH1F* fhPtMCDalitz;    //! Pt distribution of MC Dalitz electrons in EMCAL
  TH1F* fhPtMCWDecay;    //! Pt distribution of MC W decay electrons in EMCAL
  TH1F* fhPtMCZDecay;    //! Pt distribution of MC Z decay electrons in EMCAL
  TH1F* fhPtMCUnknown;   //! Pt distribution of MC unknown electrons in EMCAL

  ClassDef(AliAnaElectron,6)

} ;
 

#endif//ALIANAELECTRON_H



