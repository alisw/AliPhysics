/////////////////////////////////////////////////////////////////////////
//
// - AliEVE implementation -
// Containers for visualisation of TRD data structures 
//    - TRDHits - visualisation of MC Hits, Clusters (RecPoints)
//    - TRDDigits - visualisation of TRD digits
//
// by A.Bercuci (A.Bercuci@gsi.de)   Fri Oct 27 2006
///////////////////////////////////////////////////////////////////////

#ifndef ALIEVE_TRDData_H
#define ALIEVE_TRDData_H

#ifndef REVE_QuadSet_H
#include <Reve/QuadSet.h>
#endif

#ifndef REVE_BoxSet_H
#include <Reve/BoxSet.h>
#endif

#ifndef REVE_PointSet_H
#include <Reve/PointSet.h>
#endif

#include "AliTRDdataArrayI.h"

namespace Alieve {
	class TRDChamber;
	class TRDHits : public Reve::PointSet
	{
	public:
		TRDHits(const Text_t* name, Int_t n_points = 0);

		void PointSelected(Int_t n);

	ClassDef(TRDHits,1) // Base class for TRD hits visualisation
	};

		
	class TRDDigits : public Reve::OldQuadSet, public Reve::RenderElement
	{
	public:
		TRDDigits(TRDChamber *p);

		void			ComputeRepresentation();
		void			Paint(Option_t *opt="");
		void			Reset();
		void			SetShow(Bool_t log=kFALSE, Bool_t box=kFALSE){kLog = log; kBox=box;}
		void			SetData(AliTRDdataArrayI *digits);
		void			SetThreshold(Int_t t){fThreshold = t;}

	protected:
		TRDChamber *fChamber;
	
	private:
		Bool_t						kLog, kBox;
		UShort_t					fThreshold;
		Reve::BoxSet			fBoxes;
		AliTRDdataArrayI	fData;
		
		ClassDef(TRDDigits,1) // Digits visualisation for TRD
	};
}

#endif
