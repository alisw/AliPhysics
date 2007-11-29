#ifndef ALITRDSEEDV1_H
#define ALITRDSEEDV1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  The TRD track seed                                                    //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#ifndef ALITRDSEED_H
#include "AliTRDseed.h"
#endif

#ifndef ALIRIEMAN_H
#include "AliRieman.h"
#endif

class TTreeSRedirector;

class AliRieman;

class AliTRDstackLayer;
class AliTRDcluster;
class AliTRDrecoParam;

class AliTRDseedV1 : public AliTRDseed
{

 public:

	AliTRDseedV1(Int_t layer = -1, AliTRDrecoParam *p=0x0);
	~AliTRDseedV1();
	AliTRDseedV1(const AliTRDseedV1 &ref, Bool_t owner=kFALSE);
	AliTRDseedV1& operator=(const AliTRDseedV1 &);

	Bool_t	AttachClustersIter(AliTRDstackLayer *layer, Float_t quality, Bool_t kZcorr = kFALSE, AliTRDcluster *c=0x0);
	Bool_t	AttachClustersProj(AliTRDstackLayer *layer, Float_t quality, Bool_t kZcorr = kFALSE, AliTRDcluster *c=0x0);
	static	Float_t FitRiemanTilt(AliTRDseedV1 * cseed, Bool_t terror);
	Bool_t  FitTracklet();
	
//	        Bool_t  AttachClusters(Double_t *dx, Float_t quality, Bool_t kZcorr=kFALSE, AliTRDcluster *c=0x0);
	inline Float_t GetChi2Z(const Float_t z = 0.) const;
	inline Float_t GetChi2Y(const Float_t y = 0.) const;
	       Float_t GetQuality(Bool_t kZcorr) const;
	       Int_t   GetLayer() const                       { return fLayer;    }

	inline void    Update(const AliRieman *rieman);
               void    Print(Option_t * /*o*/) const          { }
	       void    Print();

	       void    SetLayer(Int_t l)                      { fLayer     = l;   }
	       void    SetNTimeBins(Int_t nTB)                { fTimeBins  = nTB; }
	       void    SetRecoParam(AliTRDrecoParam *p)       { fRecoParam = p;   }

 protected:

	void Copy(TObject &ref) const;

 private:

	Int_t            fLayer;     //  layer for this seed
	Int_t            fTimeBins;  //  local copy of the DB info
	Bool_t           fOwner;     //  owner of the clusters
	AliTRDrecoParam *fRecoParam; //! local copy of the reco params 

	ClassDef(AliTRDseedV1, 1)    //  New TRD seed 

};

//____________________________________________________________
inline Float_t AliTRDseedV1::GetChi2Z(const Float_t z) const
{
	Float_t z1  = (z == 0.) ? fMeanz : z;
	Float_t chi = fZref[0] - z1;
	return chi*chi;
}

//____________________________________________________________
inline Float_t AliTRDseedV1::GetChi2Y(const Float_t y) const
{
	Float_t y1  = (y == 0.) ? fYfitR[0] : y;
	Float_t chi = fYref[0] - y1;
	return chi*chi;
}

//____________________________________________________________
inline void AliTRDseedV1::Update(const AliRieman *rieman)
{
	fZref[0] = rieman->GetZat(fX0);
	fZref[1] = rieman->GetDZat(fX0);
	fYref[0] = rieman->GetYat(fX0);
	fYref[1] = rieman->GetDYat(fX0);
}

#endif

