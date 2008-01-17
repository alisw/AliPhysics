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
#include <TEveQuadSet.h>
#endif

#ifndef REVE_BoxSet_H
#include <TEveBoxSet.h>
#endif

#ifndef REVE_PointSet_H
#include <TEvePointSet.h>
#endif

#ifndef ROOT_TGedFrame
#include <TGedFrame.h>
#endif

#include "AliTRDdataArrayI.h"

class AliTRDdigitsManager;
namespace Alieve {
	class TRDChamber;
	class TRDHits : public TEvePointSet
	{
	public:
		TRDHits(TRDChamber *p);

		void PointSelected(Int_t n);

	protected:
		TRDChamber *fParent;
	
	ClassDef(TRDHits,1) // Base class for TRD hits visualisation
	};

	class TRDHitsEditor : public TGedFrame
	{
	public:
		TRDHitsEditor(const TGWindow* p=0, Int_t width = 170, Int_t height = 30, UInt_t options = kChildFrame, Pixel_t back = GetDefaultFrameBackground());
		~TRDHitsEditor();

		virtual void SetModel(TObject* obj);

	protected:
		TRDHits* fM;

	ClassDef(TRDHitsEditor,1) // Editor for TRDHits
	};	


	class TRDDigits : public TEveQuadSet
	{
	friend class TRDDigitsEditor;
	public:
		TRDDigits(TRDChamber *p);

		void			ComputeRepresentation();
		void			Paint(Option_t *opt="");
		void			Reset();
		void			SetData(AliTRDdigitsManager *digits);

	protected:
		TRDChamber *fParent;
	
	private:
		TEveBoxSet			fBoxes;
		AliTRDdataArrayI	fData;
		
		ClassDef(TRDDigits,1) // Digits visualisation for TRD
	};
	
	class TRDDigitsEditor : public TGedFrame
	{
	public:
		TRDDigitsEditor(const TGWindow* p=0, Int_t width = 170, Int_t height = 30, UInt_t options = kChildFrame, Pixel_t back = GetDefaultFrameBackground());
		~TRDDigitsEditor();

		virtual void SetModel(TObject* obj);

	protected:
		TRDDigits* fM;

	ClassDef(TRDDigitsEditor,1) // Editor for TRDDigits
	};


	class TRDClusters : public TRDHits
	{
	public:
		TRDClusters(TRDChamber *p);

		void PointSelected(Int_t n);
	
	ClassDef(TRDClusters,1) // Base class for TRD clusters visualisation
	};

}

#endif
