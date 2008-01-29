#ifndef ALITRDSTACKLAYER_H
#define ALITRDSTACKLAYER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  A TRD layer in a single stack                                         //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#ifndef ALITRDPROPAGATIONLAYER_H
#include "AliTRDpropagationLayer.h"
#endif

#ifndef ALITRDCLUSTER_H
#include "AliTRDcluster.h"
#endif

class AliTRDrecoParam;

class AliTRDstackLayer : public AliTRDpropagationLayer
{
  public:
	enum{
		kMaxClustersLayer = 150,
		kMaxRows = 16
	};

	AliTRDstackLayer(Double_t z0 = 0., Double_t zLength = 0., UChar_t stackNr = 0
                       , AliTRDrecoParam *p=0x0);
	AliTRDstackLayer(const AliTRDpropagationLayer &layer, Double_t z0
                       , Double_t zLength, UChar_t stackNr, AliTRDrecoParam *p = 0x0);
	AliTRDstackLayer(const AliTRDpropagationLayer &layer);
	AliTRDstackLayer(const AliTRDstackLayer &layer);
	~AliTRDstackLayer();
	AliTRDstackLayer   &operator=(const AliTRDpropagationLayer &myLayer);
	AliTRDstackLayer   &operator=(const AliTRDstackLayer &myLayer);
	AliTRDcluster      *operator[](const Int_t i) const {
		return ((i < fN) && (i >= 0)) ? fClusters[i] : 0x0;
	}

	void           BuildIndices(Int_t iter = 0);
	void           BuildCond(AliTRDcluster *cl, Double_t *cond, UChar_t Layer, Double_t theta=0., Double_t phi=0.);
	AliTRDcluster* GetCluster(Int_t index) const {return index < fN ? fClusters[index] : 0x0;}
	Int_t          GetGlobalIndex(const Int_t index) const {return ((index < fN) && (index >= 0)) ? fIndex[index] : 0; }
	void           GetClusters(Double_t *cond, Int_t *index, Int_t& ncl, Int_t BufferSize = kMaxClustersLayer);
	AliTRDcluster* GetNearestCluster(Double_t *cond);

	Double_t       GetZ0()                           const { return fZ0; }
	Double_t       GetDZ0()                          const { return fZLength; }
	Int_t          GetNClusters()                    const { return fN; }
	UInt_t         GetStackNr()                      const { return fStackNr; }
	
	Bool_t         IsT0() const {return TestBit(1);}
	
	void           PrintClusters()                   const;
	Int_t          SearchNearestCluster(const Double_t y, const Double_t z, const Double_t Roady, const Double_t Roadz) const;
	void           SetRange(Float_t z0, Float_t zLength);
	void           SetNRows(const Int_t nRows){ fNRows = nRows; }
	void           SetStackNr(const UInt_t stackNr){ fStackNr = stackNr; }
	void           SetOwner(Bool_t own = kTRUE) {fOwner = own;}
	void           SetClustersArray(AliTRDcluster **cl, Int_t nClusters){fClusters = cl; fN = nClusters;}
	void           SetIndexArray(UInt_t *indexArray){fIndex = indexArray;}
	void           SetDebugStream(TTreeSRedirector *debug) {fDebugStream = debug;}
	void           SetRecoParam(AliTRDrecoParam *p) {fRecoParam = p;}
	void           SetT0(Bool_t set=kTRUE) {SetBit(1, set);}

private:
	void           Copy(TObject &o) const;
	Int_t          FindYPosition(Double_t y, UChar_t z, Int_t nClusters) const;
	Int_t          FindNearestYCluster(Double_t y, UChar_t z) const;

private:
	Bool_t            fOwner;               //  owner of the clusters
	UChar_t           fStackNr;             //  stack number in supermodule
	UChar_t           fNRows;               //  number of pad rows in the chamber
	UChar_t           fPositions[kMaxRows]; //  starting index of clusters in pad row 
	Double_t          fZ0;                  //  starting position of the layer in Z direction
	Double_t          fZLength;             //  length of the layer in Z direction
	AliTRDrecoParam  *fRecoParam;           //! reconstruction parameters
	TTreeSRedirector *fDebugStream;         //! debug streamer
	
	ClassDef(AliTRDstackLayer, 1)           //  stack propagation layer

};
#endif	// ALITRDSTACKLAYER_H_

