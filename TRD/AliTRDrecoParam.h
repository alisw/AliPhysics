#ifndef ALITRDRECOPARAM_H
#define ALITRDRECOPARAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Parameter class for the TRD reconstruction                            //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#ifndef ALIDETECTORRECOPARAM_H
#include "AliDetectorRecoParam.h"
#endif

class AliTRDrecoParam : public AliDetectorRecoParam
{

  public:

	AliTRDrecoParam();
	~AliTRDrecoParam() { }

	Double_t GetChi2Y() const                 { return fkChi2Y;    }
	Double_t GetChi2Z() const                 { return fkChi2Z;    }
	Double_t GetFindableClusters() const      { return fkFindable; }
	Double_t GetMaxTheta() const              { return fkMaxTheta; }
	Double_t GetMaxPhi() const                { return fkMaxPhi;   }

	Double_t GetRoad0y() const                { return fkRoad0y;   }
	Double_t GetRoad0z() const                { return fkRoad0z;   }

	Double_t GetRoad1y() const                { return fkRoad1y;   }
	Double_t GetRoad1z() const                { return fkRoad1z;   }

	Double_t GetRoad2y() const                { return fkRoad2y;   }
	Double_t GetRoad2z() const                { return fkRoad2z;   }

	Double_t GetPlaneQualityThreshold() const { return fkPlaneQualityThreshold; }

	Double_t GetTrackLikelihood() const       { return fkTrackLikelihood;       }
	
	static   AliTRDrecoParam *GetLowFluxParam();
        static   AliTRDrecoParam *GetHighFluxParam();

 private:

	Double_t fkMaxTheta;              // Maximum theta
	Double_t fkMaxPhi;                // Maximum phi

	Double_t fkRoad0y;                // Road for middle cluster
	Double_t fkRoad0z;                // Road for middle cluster

	Double_t fkRoad1y;                // Road in y for seeded cluster
	Double_t fkRoad1z;                // Road in z for seeded cluster

	Double_t fkRoad2y;                // Road in y for extrapolated cluster
	Double_t fkRoad2z;                // Road in z for extrapolated cluster

	Double_t fkPlaneQualityThreshold; // Quality threshold
	Double_t fkFindable;              // Ratio of clusters from a track in one chamber which are at minimum supposed to be found.
	Double_t fkChi2Z;                 // Max chi2 on the z direction for seeding clusters fit
	Double_t fkChi2Y;                 // Max chi2 on the y direction for seeding clusters Rieman fit
	Double_t fkTrackLikelihood;       // Track likelihood for tracklets Rieman fit

	ClassDef(AliTRDrecoParam, 1)      // Reconstruction parameters for TRD detector

};
#endif
