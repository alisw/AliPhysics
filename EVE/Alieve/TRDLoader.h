#ifndef ALIEVE_TRDLoader_H
#define ALIEVE_TRDLoader_H

////////////////////////////////////////////////////////////////////////
//                                                                      // - ALIEVE implementation -
// Loader for the TRD detector - base class 
//    - TRDLoader - loader of TRD data (simulation + measured)
//    - TRDLoaderEditor - UI
//
// by A.Bercuci (A.Bercuci@gsi.de)   Fri Oct 27 2006
////////////////////////////////////////////////////////////////////////

#include <TEveElement.h>
#include <TGedFrame.h>
#include <TString.h>

class AliTRDv1;
class AliTRDgeometry;

class TGNumberEntry;
class TGColorSelect;
class TGTextEntry;
class TTree;

class TEveGValuator;

namespace Alieve {
	class TRDChamber;
	class TRDLoaderManager;
	enum TRDDataTypes{
		kHits = 0,
		kDigits = 1,
		kClusters = 2,
		kTracks = 3,
		kRawRoot = 4,
		kRawData = 5
	};
	class TRDLoader : public TEveElementList
	{
	friend class TRDLoaderEditor;
	public:
		TRDLoader(const Text_t* n="TRDLoader", const Text_t* t=0x0);
		~TRDLoader();
		virtual void 		Paint(Option_t *option="");
		virtual void		SetDataType(TRDDataTypes type);
	protected:
		virtual void		AddChambers(int sm=-1, int stk=-1, int ly=-1);
		virtual TRDChamber*	GetChamber(int d);
		virtual Bool_t	GoToEvent(int ev);
		virtual Bool_t	LoadClusters(TTree *tC);
		virtual Bool_t	LoadDigits(TTree *tD);
		virtual Bool_t	LoadTracklets(TTree *tT);
		virtual Bool_t	Open(const char *file, const char *dir = ".");
		virtual void		Unload();

	protected:
		Bool_t	kLoadHits, kLoadDigits, kLoadClusters, kLoadTracks;
		Int_t		fSM, fStack, fLy; // supermodule, stack, layer
		TString	fFilename; // name of data file 
		TString	fDir; // data directory
		Int_t		fEvent; // current event to be displayed
			

		AliTRDv1					*fTRD; // the TRD detector
		AliTRDgeometry		*fGeo; // the TRD geometry
		
		ClassDef(TRDLoader, 1) // Alieve Loader class for the TRD detector
	};


	
	class TRDLoaderEditor : public TGedFrame
	{
	public:
		TRDLoaderEditor(const TGWindow* p=0, Int_t width = 170, Int_t height = 30, UInt_t options = kChildFrame, Pixel_t back = GetDefaultFrameBackground());
		~TRDLoaderEditor();

		virtual void	AddChambers();
		virtual void	FileOpen();
		virtual void	Load();
		virtual void	SetEvent(Double_t ev){fM->fEvent = (Int_t)ev;}
		virtual void	SetModel(TObject* obj);
		
	protected:
		TRDLoader		*fM;
		TGTextEntry		*fFile;
		TEveGValuator		*fEvent;
		TEveGValuator		*fSMNumber, *fStackNumber, *fPlaneNumber;
		
		ClassDef(TRDLoaderEditor,1) // Editor for TRDLoader
	};
	
}
#endif
