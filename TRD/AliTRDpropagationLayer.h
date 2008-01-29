#ifndef ALITRDPROPAGATIONLAYER_H
#define ALITRDPROPAGATIONLAYER_H   

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */ 

/* $Id$ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  The TRD propagation layer                                             //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#ifndef ROOT_TObject
#include "TObject.h"
#endif

class AliTRDcluster;

class AliTRDpropagationLayer : public TObject
{	

 public:	

  enum { kMaxClusterPerTimeBin = 2300
       , kZones                = 5    };
	
  AliTRDpropagationLayer();
  AliTRDpropagationLayer(Double_t x, Double_t dx, Double_t rho
		       , Double_t x0, Int_t tbIndex, Int_t plane);
  AliTRDpropagationLayer(const AliTRDpropagationLayer &p);
  virtual ~AliTRDpropagationLayer();
  AliTRDpropagationLayer &operator=(const AliTRDpropagationLayer &/*p*/) { return *this; }

  operator Int_t() const                                        { return fN;                    }
  AliTRDcluster *operator[](Int_t i)                            { return fClusters[i];          }
  virtual void Copy(TObject &o) const;
	
  void         SetZmax(Int_t cham, Double_t center, Double_t w) { fZc[cham]      = center;
                                                                  fZmax[cham]    = w;           }
  void         SetYmax(Double_t w, Double_t wsensitive)         { fYmax          = w;
                                                                  fYmaxSensitive = wsensitive;  }

  void         SetZ(Double_t* center, Double_t *w, Double_t *wsensitive);
  void         SetHoles(Bool_t* holes);
  //void         SetHole(Double_t Zmax, Double_t Ymax
  //		       , Double_t rho = 1.29e-3, Double_t x0 = 36.66
  //                   , Double_t Yc = 0.0, Double_t Zc = 0.0);
  void         SetSector(const Int_t sec)                       { fSec           = sec;         }
  void         SetX(Double_t x)                                 { fX             = x;           }
	void         SetT0(Bool_t set = kTRUE) {SetBit(1, set);}
	
  Double_t     GetYmax() const                                  { return fYmax;                 }
  Double_t     GetZmax(Int_t c) const                           { return fZmax[c];              }
  Double_t     GetZc(Int_t c) const                             { return fZc[c];                }
  UInt_t       GetIndex(Int_t i) const                          { return fIndex[i];             }
  Double_t     GetX() const                                     { return fX;                    }
  Double_t     GetdX() const                                    { return fdX;                   }
  Int_t        GetTimeBinIndex() const                          { return fTimeBinIndex;         }
  Int_t        GetPlane() const                                 { return fPlane;                }
  Int_t        GetSector() const                                { return fSec;                  }

  Bool_t       IsHole(Int_t zone) const                         { return fIsHole[zone];         }
  Bool_t       IsSensitive() const                              { return (fTimeBinIndex >= 0) ? kTRUE : kFALSE;}
	Bool_t       IsT0() const {return TestBit(1);}

  void         Clear(const Option_t * /*o*/)                    { ; } 
  void         Clear()                                          { for (Int_t i = 0; i < fN; i++) 
                                                                    fClusters[i] = NULL;
		                                                  fN = 0;                       }
	
  void         InsertCluster(AliTRDcluster *c, UInt_t index);
  Int_t        Find(Float_t y) const;
  Int_t        FindNearestCluster(Float_t y, Float_t z, Float_t maxroad, Float_t maxroadz) const;

 protected:

	Int_t         fN;                     // Number of clusters
	Int_t         fSec;                   // Sector mumber
	AliTRDcluster **fClusters;            //[fN] Array of pointers to clusters
	UInt_t        *fIndex;                //[fN] Array of cluster indexes
	Double_t      fX;                     // X coordinate of the middle plane
	Double_t      fdX;                    // Radial thickness of the time bin
	Double_t      fRho;                   // Default density of the material
	Double_t      fX0;                    // Default radiation length
	Int_t         fTimeBinIndex;          // Plane * F(local_tb)
	Int_t         fPlane;                 // Plane number
	Double_t      fZc[kZones];            // Z position of the center for 5 active areas
	Double_t      fZmax[kZones];          // Half of active area length in Z
	Double_t      fZmaxSensitive[kZones]; // Sensitive area for detection Z
	Bool_t        fIsHole[kZones];        // Is hole in given sector
	Double_t      fYmax;                  // Half of active area length in Y
	Double_t      fYmaxSensitive;         // Half of active area length in Y
	Bool_t        fHole;                  // kTRUE if there is a hole in the layer
	Double_t      fHoleZc;                // Z of the center of the hole
	Double_t      fHoleZmax;              // Half of the hole length in Z
	Double_t      fHoleYc;                // Y of the center of the hole
	Double_t      fHoleYmax;              // Half of the hole length in Y
	Double_t      fHoleRho;               // Density of the gas in the hole
	Double_t      fHoleX0;                // Radiation length of the gas in the hole

	ClassDef(AliTRDpropagationLayer, 1)   // Propagation layer for TRD tracker

};

#endif 

