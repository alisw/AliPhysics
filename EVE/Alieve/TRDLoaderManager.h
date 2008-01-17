#ifndef ALIEVE_TRDLoaderManager_H
#define ALIEVE_TRDLoaderManager_H

////////////////////////////////////////////////////////////////////////
//                                                                      // - ALIEVE implementation -
// Loader manager for the TRD detector
//    - TRDLoaderManager - manager of TRD data loaders (simulation + measured)
//    - TRDLoaderManagerEditor - UI
//
// by A.Bercuci (A.Bercuci@gsi.de)   Mon Feb 26 2007
////////////////////////////////////////////////////////////////////////

#include <TEveElement.h>

#ifndef ROOT_TNamed
#include <TNamed.h>
#endif

#ifndef ROOT_TGedFrame
#include <TGedFrame.h>
#endif

class TGComboBox;
class TGTextButton;
class TClonesArray;

namespace Alieve {

	class TRDLoaderManager : public TEveElementList
	{
	friend class TRDLoaderManagerEditor;
	public:
		TRDLoaderManager(const Text_t* name="TRDLoader", const Text_t* title=0x0);
		~TRDLoaderManager();
		void 	Paint(Option_t *option);

	protected:
		void	Add(Int_t type, const Text_t *name, const Text_t *title=0x0);
		void	Remove(Int_t entry);
		
		ClassDef(TRDLoaderManager, 1) // Alieve loaders manager for TRD
	};

	class TRDLoaderManagerEditor : public TGedFrame
	{
	public:
		TRDLoaderManagerEditor(const TGWindow* p=0, Int_t width = 170, Int_t height = 30, UInt_t options = kChildFrame, Pixel_t back = GetDefaultFrameBackground());
		~TRDLoaderManagerEditor();
		virtual void	Add();
		virtual void	Remove(Int_t entry);
		virtual void	SetModel(TObject* obj);
		
	protected:
		TRDLoaderManager* fM;
	
	private:
		ULong_t bg;        // background color
		TGComboBox		*fSelector;
		TGTextButton	*fAdd, *fRemoveButton;
		TGGroupFrame 	*fGroupFrame;
		TClonesArray	*fRemove;

		ClassDef(TRDLoaderManagerEditor, 1)// Editor for TRDLoaderManager
	};
}

#endif

