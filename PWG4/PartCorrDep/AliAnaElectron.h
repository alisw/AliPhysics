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
  Double_t ComputeSignDca(AliAODTrack *track, AliAODTrack *track2 , float cut1);
  Int_t GetBtag(AliAODTrack * tr);

  void Print(const Option_t * opt)const;
  
  TString GetCalorimeter()   const {return fCalorimeter ; }
  Double_t GetpOverEmin()   const {return fpOverEmin ; }
  Double_t GetpOverEmax()   const {return fpOverEmax ; }
  Bool_t GetWriteNtuple()   const {return fWriteNtuple ; }

  void SetCalorimeter(TString det)    {fCalorimeter = det ; }
  void SetpOverEmin(Double_t min)     {fpOverEmin = min ; }
  void SetpOverEmax(Double_t max)     {fpOverEmax = max ; }
  void SetResidualCut(Double_t cut)     {fResidualCut = cut ; }
  void SetWriteNtuple(Bool_t val)     {fWriteNtuple = val ; }

  void InitParameters();

  void Terminate(TList * outputList);
  void ReadHistograms(TList * outputList); //Fill histograms with
					   //histograms in ouput list,
					   //needed in Terminate.            

  private:
  TString  fCalorimeter;  //! Which detector? EMCAL or PHOS
  Double_t fpOverEmin;    //! Minimum p/E value for Electrons
  Double_t fpOverEmax;    //! Maximum p/E value for Electrons
  Double_t fResidualCut;  //! Track-cluster matching distance

  //B-tagging
  Float_t fDrCut;       //max dR
  Float_t fPairDcaCut;  //max pair-DCA
  Float_t fDecayLenCut; //max 3d-decaylength
  Float_t fImpactCut;   //max track impact param
  Float_t fAssocPtCut;  //min associated pt
  Float_t fMassCut;     //min Minv cut
  Float_t fSdcaCut;     //min signDca
  Int_t   fITSCut;      //min ITS hits (both)

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

  //Reconstructed
  TH1F * fhPtElectron;  //! Number of identified electron vs transverse momentum 
  TH2F * fhPhiElectron; //! Azimuthal angle of identified  electron vs transverse momentum 
  TH2F * fhEtaElectron; //! Pseudorapidity of identified  electron vs tranvserse momentum 

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

  TH1F * fhPtPrompt;  //! Number of prompt electron vs transverse momentum 
  TH2F * fhPhiPrompt; //! Azimuthal angle of prompt  electron vs transverse momentum 
  TH2F * fhEtaPrompt; //! Pseudorapidity of prompt electron vs tranvserse momentum 

  TH1F * fhPtUnknown;  //! Number of unknown electron vs transverse momentum 
  TH2F * fhPhiUnknown; //! Azimuthal angle of unknown  electron vs transverse momentum 
  TH2F * fhEtaUnknown; //! Pseudorapidity of unknown electron vs tranvserse momentum 

  //B-tagging
  TH2F * fhBtagCut1; //! B-tagging result for cut1 (minv>1.0)
  TH2F * fhBtagCut2; //! B-tagging result for cut2 (minv>1.5)
  TH2F * fhBtagCut3; //! B-tagging result for cut3 (minv>1.8)

  //MC
  TNtuple *fMCEleNtuple; //! Ntuple of MC electrons

  ClassDef(AliAnaElectron,2)

} ;
 

#endif//ALIANAELECTRON_H



