#ifndef ALIEVE_TRDModuleImp_H
#define ALIEVE_TRDModuleImp_H

/////////////////////////////////////////////////////////////////////////
//
// Implementation of TRDModule:
//    - TRDChamber - Data holder
//    - TRDNode    - Node structure 
//
// by A.Bercuci (A.Bercuci@gsi.de)   Fri Oct 27 2006
/////////////////////////////////////////////////////////////////////////

#ifndef REVE_RenderElement_H
#include <Reve/RenderElement.h>
#endif

#ifndef ALIEVE_TRDModule_H
#include "TRDModule.h"
#endif

class AliTRDpadPlane;
class AliTRDgeometry;
class AliTRDhit;
class AliTRDdataArrayI;
class AliTRDdigitsManager;
class TObjArray;

namespace Reve {
	class Track;
}
namespace Alieve {
	class TRDHits;
	class TRDDigits;

	class TRDChamber : public Reve::RenderElement, public TRDModule
	{
	friend class TRDDigits;
	public:
	
		TRDChamber(const Int_t det=0);
		virtual ~TRDChamber() {}
	
		TRDChamber(const TRDChamber&);
		TRDChamber& operator=(const TRDChamber&);
		
		void	AddHit(AliTRDhit *hit);
		void	LoadClusters(TObjArray *cs);
		void	LoadDigits(AliTRDdigitsManager *digits);
		void	LoadTracklets(TObjArray *ts);
		void	Paint(Option_t* option="");
		void	Reset();
		void	SetGeometry(AliTRDgeometry *geo);
		
	protected:
		TRDDigits	*fDigits;   // digits representation
		TRDHits		*fHits;     // hits representation
		TRDHits		*fRecPoints;// cluster representation
		std::vector<Reve::Track*> *fTracklets; // mcm tracklets

		// data representation section
		Int_t		rowMax; // number of rows for this pad plane
  	Int_t		colMax; // number of columns for this pad plane
  	Int_t		timeMax; // number of timebins
		Float_t	samplingFrequency; // sampling frequency
		Float_t	fX0; // radial distance from vertex to the chamber
		Int_t		fPla; // detector plane
		AliTRDpadPlane *fPadPlane; // pad plane object
		AliTRDgeometry *fGeo; // TRD geometry
	
	ClassDef(TRDChamber,1) // Holder for TRD chamber data
	};

	
	class TRDNode : public Reve::RenderElementListBase, public TRDModule
	{
	public:
		TRDNode(const char *typ, const Int_t det=0);
		void	Paint(Option_t* option="");
		void	Reset();

		void	Collapse(); // *MENU*
		void	Expand(); // *MENU*
		void	EnableListElements(); // *MENU*
		void	DisableListElements(); // *MENU*
		void	UpdateLeaves();
		void	UpdateNode();
		
		List_i begin(){return fChildren.begin();}
		List_i end(){return fChildren.end();}
	
	ClassDef(TRDNode, 1)
	};
}

#endif
