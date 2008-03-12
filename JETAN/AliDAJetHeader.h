#ifndef ALIDAJETHEADER_H
#define ALIDAJETHEADER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
//---------------------------------------------------------------------
// Jet header class for Deterministic Annealing
// Stores the parameters of the DA jet algorithm
// Author: Davide Perrino (davide.perrino@ba.infn.it)
//---------------------------------------------------------------------

#include <AliJetHeader.h>

class AliDAJetHeader : public AliJetHeader
{
  public:

	AliDAJetHeader();
	virtual ~AliDAJetHeader() {}

	void SetDirectory	(Char_t *dir    ) { fDirectory = dir; }
	void SetFileOut 	(Char_t *fout   ) { fFileOut=fout;  }
	void SetPythiaOnly	(Bool_t pyt     ) { fPytOnly=pyt;   }
	void SetPtCut		(Double_t ptcut ) { fPtCut =ptcut;  }
	void SetEtaCut		(Double_t etacut) { fEtaCut=etacut; }
	void ChargedOnly	(Bool_t charged ) { fChgOnly=charged; }
	void SelectJets		(Bool_t seljets ) { fSelectJets=seljets; }
	void SetNclust		(Int_t ncl      ) { fNclustMax=ncl ; fFixedCl=kTRUE; }
	void SetEtMin		(Float_t etmin  ) { fEtMin =etmin;  }

	Char_t*  GetDirectory()	const { return fDirectory; }
	Char_t*  GetFileOut()	const { return fFileOut; }
	Bool_t   GetPythiaOnly()const { return fPytOnly; }
	Double_t GetPtCut()		const { return fPtCut;   }
	Double_t GetEtaCut()	const { return fEtaCut;  }
	Bool_t   GetChgOnly()	const { return fChgOnly; }
	Bool_t   GetSelJets()	const { return fSelectJets; }
	Int_t    GetNclustMax() const { return fNclustMax; }
	Bool_t   GetFixedCl()	const { return fFixedCl; }
	Float_t  GetEtMin()		const { return fEtMin;   }

  protected:
	AliDAJetHeader(const AliDAJetHeader &jh);
	AliDAJetHeader& operator=(const AliDAJetHeader &jh);
	Char_t	   *fDirectory;                                 //directory name 
	Char_t	   *fFileOut;                                   //output file name 
	Bool_t		fPytOnly;						//use only data from PYTHIA
	Double_t	fPtCut;							//cut on transverse momentum
	Double_t	fEtaCut;						//cut on absolute eta
	Bool_t		fChgOnly;						//flag on charged particles
	Bool_t		fSelectJets;					//select jets among clusters
	Int_t		fNclustMax;						//number of clusters when to stop annealing
	Bool_t		fFixedCl;						//use a fixed fNclustMax
	Float_t		fEtMin;							//minimum energy for found jets

	ClassDef(AliDAJetHeader,1)
};

#endif
