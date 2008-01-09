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

#ifndef ALIPID_H
#include "AliPID.h"
#endif

#ifndef ALIRIEMAN_H
#include "AliRieman.h"
#endif

class TTreeSRedirector;

class AliRieman;

class AliTRDstackLayer;
class AliTRDcluster;
class AliTRDrecoParam;
class AliTRDtrack;

class AliTRDseedV1 : public AliTRDseed
{

  public:

	enum {
	  knSlices = 10
	};

	AliTRDseedV1(Int_t layer = -1, AliTRDrecoParam *p=0x0);
	~AliTRDseedV1();
	AliTRDseedV1(const AliTRDseedV1 &ref);
	AliTRDseedV1& operator=(const AliTRDseedV1 &ref);

	Bool_t	AttachClustersIter(AliTRDstackLayer *layer, Float_t quality, Bool_t kZcorr = kFALSE
                                 , AliTRDcluster *c=0x0);
	Bool_t	AttachClusters(AliTRDstackLayer *layer, Bool_t kZcorr = kFALSE);
	void    CookdEdx(Int_t nslices);
	static	Float_t FitRiemanTilt(AliTRDseedV1 * cseed, Bool_t terror);
	Bool_t  Fit();

	       void      Init(AliTRDtrack *track);
	inline void      Init(const AliRieman *fit);
	
	inline Float_t   GetChi2Z(const Float_t z = 0.) const;
	inline Float_t   GetChi2Y(const Float_t y = 0.) const;
	       void      GetCovAt(Double_t x, Double_t *cov) const;
	       Float_t*  GetdEdx() {return &fdEdx[0];}
	       Double_t* GetdQdl() {return &fdQdl[0];}
	       Double_t  GetMomentum() const {return fMom;}
	       Int_t     GetN() const {return fN2;}
	       Float_t   GetQuality(Bool_t kZcorr) const;
	       Int_t     GetPlane() const                       { return fPlane;    }
	       Double_t* GetProbability();
	       Double_t  GetYat(Double_t x) const {return fYfitR[0] + fYfitR[1] * (x - fX0);}
	       Double_t  GetZat(Double_t x) const {return fZfitR[0] + fZfitR[1] * (x - fX0);}
				 
	       Bool_t    IsOwner() const {return fOwner;}
         void      Print(Option_t * /*o*/) const          { }
	       void      Print();
	       
	       void      SetdQdl(Double_t length);
	       void      SetMomentum(Double_t mom) {fMom = mom;}
	       void      SetOwner(Bool_t own = kTRUE);
	       void      SetPlane(Int_t p)                      { fPlane     = p;   }
	       void      SetRecoParam(AliTRDrecoParam *p)       { fRecoParam = p;   }

 protected:

	void Copy(TObject &ref) const;

 private:

	Int_t            fPlane;                  //  TRD plane
	Bool_t           fOwner;                  //  Toggle ownership of clusters
	Float_t          fMom;                    //  Momentum estimate for tracklet [GeV/c]
	Float_t          fdEdx[knSlices];         //  dE/dx measurements for tracklet
 	Double_t         fdQdl[knTimebins];       //  dQ/dl for all clusters attached to tracklet 
 	Double_t         fdQ[knTimebins];         //! dQ for all clusters attached to tracklet TODO migrate to AliTRDcluster
        Double_t         fProb[AliPID::kSPECIES]; //  PID probabilities
 	AliTRDrecoParam *fRecoParam;              //! Local copy of the reco params 

	ClassDef(AliTRDseedV1, 1)                 //  New TRD seed 

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
inline void AliTRDseedV1::Init(const AliRieman *rieman)
{
	fZref[0] = rieman->GetZat(fX0);
	fZref[1] = rieman->GetDZat(fX0);
	fYref[0] = rieman->GetYat(fX0);
	fYref[1] = rieman->GetDYat(fX0);
}

#endif

