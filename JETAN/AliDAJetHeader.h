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
	void SetRadius		(Float_t radius );
	void SetNclust		(Int_t ncl      ) { fNclustMax=ncl ; fFixedCl=kTRUE; }
	void SetEtMin		(Float_t etmin  ) { fEtMin =etmin;  }
	void SetNeff		(Int_t n        ) { fNeff = n; }
	void SetEtaEff		(Float_t eta    ) { fEtaEff = eta;  }

	Bool_t   GetSelJets()	const { return fSelectJets; }
	Float_t  GetRadius()	const { return fRadius;    }
	Int_t    GetNclustMax() const { return fNclustMax; }
	Bool_t   GetFixedCl()	const { return fFixedCl; }
	Float_t  GetEtMin()		const { return fEtMin;   }
	Int_t    GetNeff()		const { return fNeff;    }
	Float_t  GetEtaEff()	const { return fEtaEff;  }

  protected:
	AliDAJetHeader(const AliDAJetHeader &jh);
	AliDAJetHeader& operator=(const AliDAJetHeader &jh);
	Bool_t		fSelectJets;					// select jets among clusters
    Float_t 	fRadius;						// jet "radius"
	Int_t		fNclustMax;						// number of clusters when to stop annealing
	Bool_t		fFixedCl;						// use a fixed fNclustMax
	Float_t		fEtMin;							// minimum energy for found jets
	Int_t		fNeff;							// number of total input data, including fakes
	Float_t		fEtaEff;						// eta range in which fake tracks are generated

	ClassDef(AliDAJetHeader,3)
};

#endif
