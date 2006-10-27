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

#ifndef ROOT_TGedFrame
#include <TGedFrame.h>
#endif

//#include <map>

class AliRunLoader;
class AliLoader;
class AliTRDhit;
class AliTRDv1;

class TGTextButton;
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
		virtual void	AddChambers(const int sm=-1, const int stk=-1, const int ly=-1);

		TRDChamber*	GetChamber(const int d);
		Bool_t			GoToEvent(const int ev);
		Bool_t			LoadClusters(TTree *tC);
		Bool_t			LoadDigits(TTree *tD);
		Bool_t			LoadHits(TTree *tH);
		Bool_t			LoadTracklets(TTree *tT);
		Bool_t			Open(const char *filename);
		virtual void Paint(Option_t *option="");
		void				Unload();
		
	protected:
		Bool_t kLoadHits, kLoadDigits, kLoadClusters, kLoadTracks;
	
	private:
		AliRunLoader	*fRunLoader;
		AliTRDv1			*fTRD;
		
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
		TGTextButton *fOpenFile;
		Reve::RGValuator *fEvent;

		TGNumberEntry *fSMNumber, *fStackNumber, *fPlaneNumber;
		TGCheckButton *fLoadHits, *fLoadDigits, *fLoadClusters, *fLoadESDs;
		
		ClassDef(TRDLoaderEditor,1) // Editor for TRDLoader
	};
}
#endif
