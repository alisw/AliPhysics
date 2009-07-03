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

// --- ANALYSIS system ---
#include "AliAnaPartCorrBaseClass.h"

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
  
  void Print(const Option_t * opt)const;
  
  TString GetCalorimeter()   const {return fCalorimeter ; }
  Double_t GetpOverEmin()   const {return fpOverEmin ; }
  Double_t GetpOverEmax()   const {return fpOverEmax ; }

  void SetCalorimeter(TString det)    {fCalorimeter = det ; }
  void SetpOverEmin(Double_t min)     {fpOverEmin = min ; }
  void SetpOverEmax(Double_t max)     {fpOverEmax = max ; }
  void SetResidualCut(Double_t cut)     {fResidualCut = cut ; }

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
  TH2F *fh2OuterPtVsExtrapPt;
  TH2F *fh2OuterPhiVsExtrapPhi;
  TH2F *fh2OuterEtaVsExtrapEta;

  TH2F* fh2TrackPVsClusterE;
  TH2F* fh2TrackPtVsClusterE;
  TH2F* fh2TrackPhiVsClusterPhi;
  TH2F* fh2TrackEtaVsClusterEta;

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

  /*
  //MC
  TH1F * fhMCPtElectron;   //! pT of MC electrons 
  TH2F * fhMCPhiElectron;  //! Phi of MC electrons
  TH2F * fhMCEtaElectron;  //! eta of MC electrons

  TH1F * fhMCPtConversion;  //! Number of TRACKABLE MC conversion electron vs transverse momentum 
  TH2F * fhMCPhiConversion; //! Azimuthal angle of TRACKABLE MC conversion  electron vs transverse momentum 
  TH2F * fhMCEtaConversion; //! Pseudorapidity of TRACKABLE MC conversion electron vs tranvserse momentum 

  TH1F * fhMCPtBottom;  //! Number of MC bottom electron vs transverse momentum 
  TH2F * fhMCPhiBottom; //! Azimuthal angle of MC bottom  electron vs transverse momentum 
  TH2F * fhMCEtaBottom; //! Pseudorapidity of MC bottom electron vs tranvserse momentum 

  TH1F * fhMCPtCharm;  //! Number of MC charm electron vs transverse momentum 
  TH2F * fhMCPhiCharm; //! Azimuthal angle of MC charm  electron vs transverse momentum 
  TH2F * fhMCEtaCharm; //! Pseudorapidity of MC charm electron vs tranvserse momentum 

  TH1F * fhMCPtCFromB;  //! Number of MC charm from bottom electron vs transverse momentum 
  TH2F * fhMCPhiCFromB; //! Azimuthal angle of MC charm from bottom electron vs transverse momentum
  TH2F * fhMCEtaCFromB; //! Pseudorapidity of MC charm from bottom electron vs tranvserse momentum 

  TH1F * fhMCPtDalitz;  //! Number of MC dalitz electron vs transverse momentum 
  TH2F * fhMCPhiDalitz; //! Azimuthal angle of MC dalitz  electron vs transverse momentum 
  TH2F * fhMCEtaDalitz; //! Pseudorapidity of MC dalitz electron vs tranvserse momentum 

  TH1F * fhMCPtWDecay;  //! Number of MC W-boson electron vs transverse momentum 
  TH2F * fhMCPhiWDecay; //! Azimuthal angle of MC W-boson  electron vs transverse momentum 
  TH2F * fhMCEtaWDecay; //! Pseudorapidity of MC W-boson electron vs tranvserse momentum 
  		
  TH1F * fhMCPtZDecay;  //! Number of MC Z-boson electron vs transverse momentum 
  TH2F * fhMCPhiZDecay; //! Azimuthal angle of MC Z-boson  electron vs transverse momentum 
  TH2F * fhMCEtaZDecay; //! Pseudorapidity of MC Z-boson electron vs tranvserse momentum 

  TH1F * fhMCPtPrompt;  //! Number of prompt MC electron vs transverse momentum 
  TH2F * fhMCPhiPrompt; //! Azimuthal angle of prompt MC electron vs transverse momentum 
  TH2F * fhMCEtaPrompt; //! Pseudorapidity of prompt MC electron vs tranvserse momentum 

  TH1F * fhMCPtUnknown;  //! Number of unknown MC electron vs transverse momentum 
  TH2F * fhMCPhiUnknown; //! Azimuthal angle of unknown MC electron vs transverse momentum 
  TH2F * fhMCEtaUnknown; //! Pseudorapidity of unknown MC electron vs tranvserse momentum 
  */

  ClassDef(AliAnaElectron,1)

} ;
 

#endif//ALIANAELECTRON_H



