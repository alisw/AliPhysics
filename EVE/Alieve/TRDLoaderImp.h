#ifndef ALIEVE_TRDLoaderImp_H
#define ALIEVE_TRDLoaderImp_H

////////////////////////////////////////////////////////////////////////
//                                                                     // - ALIEVE implementation -
// Single event loader for the TRD detector
//    - TRDLoaderSim - loader for simulations based on gAlice
//    - TRDLoaderRaw - loader for raw data
//
// by A.Bercuci (A.Bercuci@gsi.de)   Fri Oct 27 2006
////////////////////////////////////////////////////////////////////////

#ifndef ALIEVE_TRDLoader_H
#include "TRDLoader.h"
#endif


class AliRunLoader;
class AliTRDrawData;
class AliRawReaderDate;
class AliRawReaderRoot;

class TGCheckButton;
namespace Alieve {
	class TRDLoaderSim : public TRDLoader
	{
	friend class TRDLoaderSimEditor;
	public:
		TRDLoaderSim(const Text_t* n="TRDLoaderSim", const Text_t* t=0);
		~TRDLoaderSim();
		
		Bool_t			GoToEvent(int ev);
		Bool_t			LoadHits(TTree *tH);
		Bool_t			Open(const char *file, const char *dir=".");
			
	private:
		AliRunLoader			*fRunLoader; // Run Loader
	
		ClassDef(TRDLoaderSim, 1) // Alieve loader for the TRD detector (gAlice)
	};
	

	class TRDLoaderRaw : public TRDLoader
	{
	public:
		TRDLoaderRaw(const Text_t* n="TRDLoaderRaw", const Text_t* t=0);
		~TRDLoaderRaw();
		
		Bool_t			GoToEvent(int ev);
		Bool_t 			LoadEvent();
		Bool_t			Open(const char *file, const char *dir=".");
		void				SetDataType(TRDDataTypes type);

		
	private:
		void NextEvent(Bool_t rewindOnEnd=kTRUE);
	
	private:
		AliRawReaderDate	*fRawDateReader;
		AliRawReaderRoot	*fRawRootReader;
		AliTRDrawData			*fRaw;
		Bool_t						fDataRoot;
		Int_t							fEventOld;
		
		ClassDef(TRDLoaderRaw, 1) // Alieve loader for the TRD detector (raw)
	};

	class TRDLoaderSimEditor : public TGedFrame
	{
	public:
		TRDLoaderSimEditor(const TGWindow* p=0, Int_t width = 170, Int_t height = 30, UInt_t options = kChildFrame, Pixel_t back = GetDefaultFrameBackground());
		~TRDLoaderSimEditor();

		virtual void	SetModel(TObject* obj);
		virtual void	Toggle(Int_t id);
	
	protected:
		TRDLoaderSim* fM;
		TGCheckButton *fLoadHits, *fLoadDigits, *fLoadClusters, *fLoadTracks;
		
		ClassDef(TRDLoaderSimEditor,1) // Editor for TRDLoaderSim
	};
}
#endif
