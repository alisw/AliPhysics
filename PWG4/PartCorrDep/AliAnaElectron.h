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
class TObjString;

// --- ANALYSIS system ---
#include "AliAnaPartCorrBaseClass.h"

class AliAODMCParticle;
class AliCaloTrackReader;
class AliAODTrack;
class TList ;

class AliAnaElectron : public AliAnaPartCorrBaseClass {

 public: 
  AliAnaElectron() ; // default ctor
  virtual ~AliAnaElectron() ; //virtual dtor
 private:
  AliAnaElectron(const AliAnaElectron & g) ; // cpy ctor
  AliAnaElectron & operator = (const AliAnaElectron & g) ;//cpy assignment
  
 public:
	
  TObjString * GetAnalysisCuts();
  TList      * GetCreateOutputObjects();

  void Init();

  void MakeAnalysisFillAOD()  ;
  
  void MakeAnalysisFillHistograms() ; 
  
  //B-tagging
  Int_t GetDVMBtag(AliAODTrack * tr); //returns # tracks from secvtx

  //Temporary local method to get DCA because AliAODTrack is stupid
  Bool_t GetDCA(const AliAODTrack* tr,Double_t imp[2], Double_t cov[3]);

  Bool_t PhotonicPrim(const AliAODPWG4Particle* part); //check with track list
  Bool_t PhotonicV0(Int_t trackId); //check with V0 list

  //check if track has been flagged as a non-photonic or DVM electron
  //used with the jet tracks to tag bjets
  Bool_t CheckTrack(const AliAODTrack* track,const char* type);  
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
  //For determining origin of electron
  Int_t GetMCSource(Int_t mctag);

  //Need a clean way to get the MC info.  An AliAODMCParticle object
  //is returned from whichever source we are operating on
  AliAODMCParticle* GetMCParticle(Int_t part);
  //Get MC B Parent pt
  Double_t GetBParentPt(Int_t label);
  //Get Number of particles in AliAODMCParticle array, if it exists
  Int_t GetNumAODMCParticles();

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

  ///////////////////////////////////////
  //RC = RECO only - these histos will be filled using only reco
  //information

  //event QA
  TH1F * fhImpactXY;    //! XY impact parameter of all tracks to primary vertex
  TH1F * fhRefMult;     //! refmult (tracks with |eta| < 0.5)
  TH1F * fhRefMult2;    //! refmult2 (tracks with |eta| < 0.5 & impXY,impZ < 1.0)

  //matching checks   
  TH3F *fh3pOverE;     //! p/E for track-cluster matches vs pt vs mult
  TH3F *fh3EOverp;     //! E/p for track-cluster matches vs pt vs mult
  TH3F *fh3pOverE2;     //! p/E for track-cluster matches vs pt vs mult
  TH3F *fh3EOverp2;     //! E/p for track-cluster matches vs pt vs mult
  TH3F *fh3pOverE3;     //! p/E for track-cluster matches vs pt vs mult
  TH3F *fh3EOverp3;     //! E/p for track-cluster matches vs pt vs mult

  //JLK
  TH2F *fh2pOverE;      //! p/E for track-cluster matches vs pt vs mult         
  TH2F *fh2EOverp;      //! E/p for track-cluster matches vs pt vs mult         
  TH2F *fh2pOverE2;     //! p/E for track-cluster matches vs pt vs mult         
  TH2F *fh2EOverp2;     //! E/p for track-cluster matches vs pt vs mult         
  //JLK

  TH1F *fh1dR;         //! distance between projected track and cluster
  TH2F *fh2EledEdx;    //! dE/dx vs. momentum for electron candidates
  TH2F *fh2MatchdEdx;  //! dE/dx vs. momentum for all matches
  TH2F *fh2dEtadPhi;   //! DeltaEta vs. DeltaPhi of all track/cluster pairs
  TH2F *fh2dEtadPhiMatched;   //! DeltaEta vs. DeltaPhi of matched track/cluster pairs
  TH2F *fh2dEtadPhiUnmatched;   //! DeltaEta vs. DeltaPhi of unmatched track/cluster pairs

  TH2F* fh2TrackPVsClusterE;     //!track momentum vs. cluster energy
  TH2F* fh2TrackPtVsClusterE;    //!track pt vs. cluster energy
  TH2F* fh2TrackPhiVsClusterPhi; //!track phi vs. cluster phi
  TH2F* fh2TrackEtaVsClusterEta; //!track eta vs. cluster eta

  //Photonic Electron checks
  TH1F* fh1OpeningAngle; //!opening angle between pairs of photon candidates
  TH1F* fh1MinvPhoton;   //!invariant mass distribution of electron pairs

  //Reconstructed electrons
  TH1F * fhPtElectron;  //! Number of identified electron vs transverse momentum 
  TH2F * fhPhiElectron; //! Azimuthal angle of identified  electron vs transverse momentum 
  TH2F * fhEtaElectron; //! Pseudorapidity of identified  electron vs tranvserse momentum 

  TH1F * fhPtNPE;  //! Number of non-photonic electron vs transverse momentum 
  TH2F * fhPhiNPE; //! Azimuthal angle of non-photonic electron vs transverse momentum 
  TH2F * fhEtaNPE; //! Pseudorapidity of non-photonic electron vs tranvserse momentum 

  TH1F * fhPtPE;  //! Number of photonic electron vs transverse momentum 
  TH2F * fhPhiPE; //! Azimuthal angle of photonic electron vs transverse momentum 
  TH2F * fhEtaPE; //! Pseudorapidity of photonic electron vs tranvserse momentum 

  //These next set do use some MC info.  The first bin of the second
  //dimension is filled for both REAL and MC data, other bins filled
  //only if MC
  //Histograms for comparison to tracking detectors
  TH2F* fhPtHadron;        //!Pt distribution of reco charged hadrons
                           //!(pi,k,p) in EMCAL acceptance
  TH2F* fhPtNPEleTPC;      //!Pt distribution of non-photonic reco electrons using
			   //!just TPC dEdx info in EMCAL acceptance
  TH2F* fhPtNPEleTPCTRD;   //!Pt distribution of non-photonic reco electrons using
			   //!pid info from tracking detectors only in EMCAL acceptance
  TH2F* fhPtNPEleTTE;      //!Pt distribution of non-photonic reco
			   //!electrons using pid info from TPC+TRD+EMCAL
			   //!in EMCAL acceptance
  TH2F* fhPtNPEleEMCAL;    //!Pt distribution of non-photonic reco
			   //!electrons using EMCAL only
			   //!in EMCAL acceptance

  //DVM B-tagging
  TH2F * fhDVMBtagCut1; //! DVM B-tagging result for cut1 (minv>1.0)
  TH2F * fhDVMBtagCut2; //! DVM B-tagging result for cut2 (minv>1.5)
  TH2F * fhDVMBtagCut3; //! DVM B-tagging result for cut3 (minv>1.8)
  TH2F * fhDVMBtagQA1;  //! DVM B-tagging : QA of pairDca vs decaylength
  TH2F * fhDVMBtagQA2;  //! DVM B-tagging : QA of signDca vs mass
  TH1F * fhDVMBtagQA3;  //! DVM B-tagging : QA number of ITS clusters
  TH1F * fhDVMBtagQA4;  //! DVM B-tagging : QA prim vtx impXY
  TH1F * fhDVMBtagQA5;  //! DVM B-tagging : QA prim vtx impZ
  //IPSig B-tagging
  TH1F * fhIPSigBtagQA1; //! IPSig B-tagging : QA of # tag tracks
  TH1F * fhIPSigBtagQA2; //! IPSig B-tagging : QA of IP sig
  TH1F * fhTagJetPt1x4;  //! IPSig B-tagging : result for (1 track, ipSignif>4)
  TH1F * fhTagJetPt2x3;  //! IPSig B-tagging : result for (2 track, ipSignif>3)
  TH1F * fhTagJetPt3x2;  //! IPSig B-tagging : result for (3 track, ipSignif>2)
  TH1F * fhePlusTagJetPt1x4;  //! IPSig B-tagging : eJet + result for (1 track, ipSignif>4)
  TH1F * fhePlusTagJetPt2x3;  //! IPSig B-tagging : eJet + result for (2 track, ipSignif>3)
  TH1F * fhePlusTagJetPt3x2;  //! IPSig B-tagging : eJet + result for (3 track, ipSignif>2)

  //B-Jet histograms
  TH2F* fhJetType;       //! How many of each tag were found vs jet pt
  TH2F* fhLeadJetType;   //! How many leading of each tag were found vs jet pt
  TH2F* fhBJetXsiFF;     //! B-tagged jet FF with xsi = log(pt_Jet/pt_Track)
  TH2F* fhBJetPtFF;      //! B-tagged jet FF with pt_Track
  TH2F* fhBJetEtaPhi;    //! B-tagged jet eta-phi distribution
  TH2F* fhNonBJetXsiFF;  //! Non b-tagged jet FF with xsi = log(pt_Jet/pt_Track)
  TH2F* fhNonBJetPtFF;   //! Non b-tagged jet FF with pt_Track
  TH2F* fhNonBJetEtaPhi; //! Non b-tagged jet eta-phi distribution

  ///////////////////////////////////////////////////////////////////
  //MC = From here down, the histograms use MC information, so they will
  //only be filled in simulations
  TNtuple* fEleNtuple; //! testing ntuple

  TH2F * fhPhiConversion; //! Azimuthal angle of conversion  electron vs transverse momentum 
  TH2F * fhEtaConversion; //! Pseudorapidity of conversion electron vs tranvserse momentum 

  //Histograms for comparison to tracking detectors
  TH2F* fhPtTrack;         //!Pt distribution of reco tracks with MC-ID

  TH2F* fhPtNPEBHadron;    //!correlate our best reconstructed
			   //b-electrons with the b-hadron momentum

  //For computing efficiency of IPSIG tag
  //these require that an MC b-Ancestor is present in the jet
  TH1F * fhBJetPt1x4;    //! IPSig B-tagging : result for (1 track, ipSignif>4)
  TH1F * fhBJetPt2x3;    //! IPSig B-tagging : result for (2 track, ipSignif>3)
  TH1F * fhBJetPt3x2;    //! IPSig B-tagging : result for (3 track, ipSignif>2)

  TH1F * fhFakeJetPt1x4;    //! IPSig B-tagging : fake result for (1 track, ipSignif>4)
  TH1F * fhFakeJetPt2x3;    //! IPSig B-tagging : fake result for (2 track, ipSignif>3)
  TH1F * fhFakeJetPt3x2;    //! IPSig B-tagging : fake result for (3 track, ipSignif>2)

  TH2F* fhDVMJet;        //! DVM jet algo check

  ////////////////////////////
  //MC Only Rate histograms

  TNtuple *fMCEleNtuple;  //! Ntuple of MC electrons

  TH2F* fhMCBJetElePt;    //! Pt of B-Jet vs pt of electron
  TH2F* fhMCBHadronElePt; //! Pt of B-hadrons vs pt of electron
  TH1F* fhPtMCHadron;     //! Pt distribution of MC charged hadrons (pi,k,p) in EMCAL acceptance
  TH2F* fhPtMCElectron;   //! Pt distribution of MC electrons from various sources in EMCAL
  TH2F* fhMCXYConversion; //! XY distribution of conversion electrons
  TH2F* fhMCRadPtConversion; //! Radius vs. pT distribution of conversion electrons

  ClassDef(AliAnaElectron,12)

} ;
 

#endif//ALIANAELECTRON_H



