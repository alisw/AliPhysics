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

	void SelectJets		(Bool_t seljets ) { fSelectJets=seljets; }
	void SetNclust		(Int_t ncl      ) { fNclustMax=ncl ; fFixedCl=kTRUE; }
	void SetEtMin		(Float_t etmin  ) { fEtMin =etmin;  }

	Bool_t   GetSelJets()	const { return fSelectJets; }
	Int_t    GetNclustMax() const { return fNclustMax; }
	Bool_t   GetFixedCl()	const { return fFixedCl; }
	Float_t  GetEtMin()		const { return fEtMin;   }

  protected:
	AliDAJetHeader(const AliDAJetHeader &jh);
	AliDAJetHeader& operator=(const AliDAJetHeader &jh);
	Bool_t		fSelectJets;					//select jets among clusters
	Int_t		fNclustMax;						//number of clusters when to stop annealing
	Bool_t		fFixedCl;						//use a fixed fNclustMax
	Float_t		fEtMin;							//minimum energy for found jets

	ClassDef(AliDAJetHeader,2)
};

#endif
