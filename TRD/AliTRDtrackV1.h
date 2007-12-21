#ifndef ALITRDTRACKV1_H
#define ALITRDTRACKV1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  The TRD track                                                         //
//                                                                        //
//  Authors:                                                              //
//    Alex Bercuci <A.Bercuci@gsi.de>                                     //
//    Markus Fasel <M.Fasel@gsi.de>                                       //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#ifndef ALITRDTRACK_H
#include "AliTRDtrack.h"
#endif

class AliTRDseedV1;
class AliESDtrack;
class AliCluster;

class AliTRDtrackV1 : public AliTRDtrack
{
public:
	AliTRDtrackV1();
	AliTRDtrackV1(AliTRDseedV1 *trklts, const Double_t p[5], const Double_t cov[15], Double_t x, Double_t alpha);
	AliTRDtrackV1(const AliESDtrack &ref);
	AliTRDtrackV1(const AliTRDtrackV1 & /*ref*/);
	AliTRDtrackV1& operator=(const AliTRDtrackV1 &/*ref*/){ return *this; }

	Bool_t        CookPID();
	inline Int_t  GetNumberOfClusters() const;
	Double_t      GetPredictedChi2(const AliTRDseedV1 *tracklet) const;
        Double_t      GetPredictedChi2(const AliCluster* /*c*/) const                   { return 0.0; }
	const AliTRDseedV1* GetTracklet(Int_t plane) const {return plane >=0 && plane <6 ? &fTracklet[plane] : 0x0;}
	Int_t*        GetTrackletIndexes() {return &fTrackletIndex[0];}
	void          SetTracklet(AliTRDseedV1 *trklt, Int_t plane, Int_t index);
	Bool_t        Update(const  AliTRDseedV1 *tracklet, Double_t chi2);
virtual Bool_t        Update(const AliCluster*, Double_t, Int_t) { return kTRUE; }
	void          UpdateESDdEdx(AliESDtrack *t);

protected:
	AliTRDrecoParam *fRecoParam;       // reconstruction parameters
	ClassDef(AliTRDtrackV1, 1)         // development TRD track
};

//___________________________________________________________
inline Int_t AliTRDtrackV1::GetNumberOfClusters() const
{
/*	Int_t ntrklts = GetNumberOfTracklets();
	printf("AliTRDtrackV1::GetNumberOfClusters() %d\n", ntrklts);
	return ntrklts;*/
	Int_t ncls = 0;
	for(int ip=0; ip<6; ip++)
		if(fTrackletIndex[ip] >= 0) ncls+=fTracklet[ip].GetN();
		
	return ncls;
}

// //___________________________________________________________
// inline Int_t AliTRDtrackV1::GetNumberOfTracklets() const
// {
// 	Int_t ntrklt = 0;
// 	for(int ip=0; ip<6; ip++) if(fTrackletIndex[ip] >= 0) ntrklt++;
// 	return ntrklt;
// }

#endif


