#ifndef ALIEVE_TRDModule_H
#define ALIEVE_TRDModule_H

/////////////////////////////////////////////////////////////////////////
//
// - AliEVE implementation -
// The common structure of a TRD module (SM, Stack or Chamber)
//    - TRDModule - structure of TRD module for visualisation
//    - TRDModuleEditor - UI
//
// by A.Bercuci (A.Bercuci@gsi.de)   Fri Oct 27 2006
///////////////////////////////////////////////////////////////////////

#ifndef ROOT_TNamed
#include <TNamed.h>
#endif

#ifndef ROOT_TGedFrame
#include <TGedFrame.h>
#endif

class TObject;
class TGWindow;
class TGCheckButton;
class TGNumberEntry;
class TGColorSelect;

namespace Alieve {

	class TRDModule : public TNamed
	{
	friend class TRDModuleEditor;
	friend class TRDNode;
	friend class TRDChamber;
	public:
		TRDModule(const char *typ="XXX", const Int_t id=0);
		virtual ~TRDModule() {}
	
		virtual Bool_t GetDigitsBox(){return fDigitsBox;}
		virtual Bool_t GetDigitsLog(){return fDigitsLog;}
		virtual UShort_t GetDigitsThreshold(){return fDigitsThreshold;}
		virtual Int_t	GetID(){return fDet;}
		virtual void	Paint(Option_t* option="")=0;
		virtual void	Reset()=0;

	protected:
		// UI section
		Bool_t	fLoadHits, fRnrHits;
		Bool_t	fLoadDigits, fRnrDigits, fDigitsLog, fDigitsBox;
		Bool_t	kDigitsNeedRecompute;

		Bool_t	fLoadRecPoints, fRnrRecPoints;
		Bool_t	fLoadTracklets, fRnrTracklets;
		
		Int_t fDet; // detector number
		UShort_t	fDigitsThreshold; // digits threshold
	ClassDef(TRDModule,1) // Structure holder for TRD chamber
	};
	

	class TRDModuleEditor : public TGedFrame
	{
	public:
		TRDModuleEditor(const TGWindow* p=0, Int_t width = 170, Int_t height = 30, UInt_t options = kChildFrame, Pixel_t back = GetDefaultFrameBackground());
		~TRDModuleEditor();

		virtual void SetModel(TObject* obj);

		void	ModifyDigitsView();
		void	SetThreshold(Long_t thres);
		void	UpdateChamber();
		void	UpdateClusters(Pixel_t);
		void	UpdateHits(Pixel_t);
	
	protected:
		TRDModule* fM;

	private:
		TGCheckButton *fDisplayHits;
		TGColorSelect *fHitsColor;
		TGCheckButton *fDisplayDigits, *fToggleLog, *fToggleBox, *fThreshold;
		TGNumberEntry	*fThresValue;
		TGCheckButton *fDisplayClusters;
		TGColorSelect *fClustersColor;
		TGCheckButton *fDisplayTracks;
		
	ClassDef(TRDModuleEditor,1) // Editor for TRDModule
	};
}

#endif
