/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//-----------------------------------------------------------------------
// Author : R. Vernet, Consorzio Cometa - Catania (it)
//-----------------------------------------------------------------------
// Modification done by X. Lopez - LPC Clermont (fr)
//-----------------------------------------------------------------------
// Modification done by S. Ahn - LPC Clermont (fr), Konkuk University (kr)
//-----------------------------------------------------------------------

#ifndef ALICFMUONRESUPSILON_H
#define ALICFMUONRESUPSILON_H

#include "AliAnalysisTaskSE.h"

class AliCFContainer;
class AliCFManager;

class AliCFMuonResUpsilon : public AliAnalysisTaskSE {
  public:

  AliCFMuonResUpsilon();
  AliCFMuonResUpsilon(const Char_t* name);
  AliCFMuonResUpsilon& operator= (const AliCFMuonResUpsilon& c);
  AliCFMuonResUpsilon(const AliCFMuonResUpsilon& c);
  virtual ~AliCFMuonResUpsilon();

  // ANALYSIS FRAMEWORK STUFF to loop on data and fill output objects
  void     UserExec(Option_t *option);
  void     Terminate(Option_t *);
	void		 UserCreateOutputObjects();
   
  // CORRECTION FRAMEWORK RELATED FUNCTIONS
  void           SetCFManager(AliCFManager* const io)	  {fCFManager = io;}   // global correction manager
  AliCFManager * GetCFManager() const 									{return fCFManager;}           // get corr manager
  void           SetQAList(TList* const list) 					{fQAHistList = list;}

  // Setters
  Bool_t IsReadAODData() const 													{return fReadAODData;}
  void   SetReadAODData(Bool_t flag=kTRUE) 							{fReadAODData=flag;}
  void   SetReadMCInfo(Bool_t flag=kTRUE) 							{fReadMCInfo=flag;}
	void 	 SetRequirePdgCode(const Int_t PDG)							{fPDG=PDG;}
	void	 SetPtRange(const Double_t ptmin, const Double_t ptmax) {fPtMin=ptmin; fPtMax=ptmax; }
	void	 SetRapidityRange(const Double_t ymin, const Double_t ymax) {fYMin=ymin; fYMax=ymax; }
  void   SetTrigClassMuonName(TString name = "CMUS")		  {fTrigClassMuon=name;}
  void   SetTrigClassInteracName(TString name = "CINT") {fTrigClassInteraction=name;}
  void   SetTrigClassSideName(TString name[3])	  	{for(Int_t i=0;i<3;i++) fTrigClassSide[i]=name[i];}

 protected:
  
  Bool_t          fReadAODData;   // flag for AOD/ESD input files
  Bool_t          fReadMCInfo;   	// flag for reading MC info (ESD->Kinematics, AOD->MCbranch)
  AliCFManager   *fCFManager;   	// pointer to the CF manager
 	TH1D					 *fnevts;					// TH1 for event statistics
	Bool_t					fIsPhysSelMB;		// flag for the physics selection : MB
	Bool_t					fIsPhysSelMUON;	// flag for the physics selection : MUON
  TList          *fQAHistList;   	// list of QA histograms
	Int_t						fPDG;						// PDG code of resonance
	Double_t				fPtMin;					// min Pt of resonance
	Double_t				fPtMax;					// max Pt of resonance
	Double_t				fYMin;					// min rapidity of resonance
	Double_t				fYMax;					// max rapidity of resonance

  // CUTS ON THE FIRED TRIGGER CLASS
  TString					fTrigClassMuon		;    // name of the muon trigger class (CMU by default)
  TString					fTrigClassInteraction	;    // name of the interaction trigger class (CINT by default)
  TString					fTrigClassSide[4]	;    // name of the muon trigger classes containing the side
  Bool_t					fDistinguishTrigClass	;    // flag to activate the cut on the fired trigger class


  Double_t Imass(Float_t e1, Float_t px1, Float_t py1, Float_t pz1,	Float_t e2, Float_t px2, Float_t py2, Float_t p2) const;
 	Double_t Rap(Float_t e, Float_t pz) const;
  
  Double_t CostCS(Double_t px1, Double_t py1, Double_t pz1, Double_t e1, Double_t charge1,
                        Double_t px2, Double_t py2, Double_t pz2, Double_t e2, Double_t energy);
  Double_t CostHE(Double_t px1, Double_t py1, Double_t pz1, Double_t e1, Double_t charge1,
                        Double_t px2, Double_t py2, Double_t pz2, Double_t e2);
  Double_t PhiCS(Double_t px1, Double_t py1, Double_t pz1, Double_t e1, Double_t charge1,
                       Double_t px2, Double_t py2, Double_t pz2, Double_t e2, Double_t energy);
  Double_t PhiHE(Double_t px1, Double_t py1, Double_t pz1, Double_t e1, Double_t charge1,
                       Double_t px2, Double_t py2, Double_t pz2, Double_t e2, Double_t energy);
  
  ClassDef(AliCFMuonResUpsilon,1);
};

#endif
