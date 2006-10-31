#ifndef ALIEVE_TRDLoader_H
#define ALIEVE_TRDLoader_H

////////////////////////////////////////////////////////////////////////
//                                                                      // - ALIEVE implementation -
// Loader for the TRD detector
//    - TRDLoader - loader of TRD data (simulation + measured)
//    - TRDLoaderEditor - UI
//
// by A.Bercuci (A.Bercuci@gsi.de)   Fri Oct 27 2006
////////////////////////////////////////////////////////////////////////

#ifndef REVE_RenderElement_H
#include <Reve/RenderElement.h>
#endif

#ifndef ROOT_TNamed
#include <TNamed.h>
#endif

#ifndef ROOT_TString
#include <TString.h>
#endif

#ifndef ROOT_TGedFrame
#include <TGedFrame.h>
#endif

class AliRunLoader;
class AliTRDv1;
class AliTRDgeometry;

class TGCheckButton;
class TGNumberEntry;
class TGColorSelect;
class TGTextEntry;
class TTree;
namespace Reve {
	class RGValuator;
}
namespace Alieve {
	class TRDChamber;
	enum DataTypes{
		kHits = 0,
		kDigits = 1,
		kClusters = 2,
		kESDs = 3 
	};
	class TRDLoader : public Reve::RenderElementListBase, public TNamed
	{
	friend class TRDLoaderEditor;
	public:
		TRDLoader(const Text_t* n="TRDLoader", const Text_t* t=0);
		~TRDLoader();

	protected:
		virtual void		AddChambers(const int sm=-1, const int stk=-1, const int ly=-1);
		virtual TRDChamber*	GetChamber(const int d);
		virtual Bool_t	GoToEvent(const int ev);
		virtual Bool_t	LoadClusters(TTree *tC);
		virtual Bool_t	LoadDigits(TTree *tD);
		virtual Bool_t	LoadHits(TTree *tH);
		virtual Bool_t	LoadTracklets(TTree *tT);
		virtual Bool_t	Open(const char *file, const char *dir = ".");
		virtual void 		Paint(Option_t *option="");
		virtual void		Unload();
		
	protected:
		Bool_t	kLoadHits, kLoadDigits, kLoadClusters, kLoadTracks;
		Int_t		fSM, fStack, fLy; // supermodule, stack, layer
		TString	fFilename; // name of data file 
		TString	fDir; // data directory
		Int_t		fEvent; // current event to be displayed
			
	private:
		AliTRDv1				*fTRD; // the TRD detector
		AliTRDgeometry	*fGeo; // the TRD geometry
		AliRunLoader		*fRunLoader; // Run Loader
		
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
		virtual void	SetModel(TObject* obj);
		
		
	protected:
		TRDLoader* fM;
  	TGTextEntry *fFile;
		Reve::RGValuator *fEvent;

		Reve::RGValuator *fSMNumber, *fStackNumber, *fPlaneNumber;
		TGCheckButton *fLoadHits, *fLoadDigits, *fLoadClusters, *fLoadTracks;
		
		ClassDef(TRDLoaderEditor,1) // Editor for TRDLoader
	};
}
#endif
