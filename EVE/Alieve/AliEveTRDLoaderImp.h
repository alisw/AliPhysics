// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/
#ifndef ALIEVE_TRDLoaderImp_H
#define ALIEVE_TRDLoaderImp_H

////////////////////////////////////////////////////////////////////////
//                                                                     // - ALIEVE implementation -
// Single event loader for the TRD detector
//    - AliEveTRDLoaderSim - loader for simulations based on gAlice
//    - AliEveTRDLoaderRaw - loader for raw data
//
// by A.Bercuci (A.Bercuci@gsi.de)   Fri Oct 27 2006
////////////////////////////////////////////////////////////////////////

#ifndef ALIEVE_TRDLoader_H
#include "AliEveTRDLoader.h"
#endif


class AliRunLoader;
class AliTRDrawData;
class AliRawReaderDate;
class AliRawReaderRoot;

class TGCheckButton;

	class AliEveTRDLoaderSim : public AliEveTRDLoader
	{
	friend class AliEveTRDLoaderSimEditor;
	public:
		AliEveTRDLoaderSim(const Text_t* n="AliEveTRDLoaderSim", const Text_t* t=0);
		~AliEveTRDLoaderSim();

		Bool_t			GoToEvent(int ev);
		Bool_t			LoadHits(TTree *tH);
		Bool_t			Open(const char *file, const char *dir=".");

	private:
		AliRunLoader			*fRunLoader; // Run Loader

		ClassDef(AliEveTRDLoaderSim, 1) // Alieve loader for the TRD detector (gAlice)
	};


	class AliEveTRDLoaderRaw : public AliEveTRDLoader
	{
	public:
		AliEveTRDLoaderRaw(const Text_t* n="AliEveTRDLoaderRaw", const Text_t* t=0);
		~AliEveTRDLoaderRaw();

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

		ClassDef(AliEveTRDLoaderRaw, 1) // Alieve loader for the TRD detector (raw)
	};

	class AliEveTRDLoaderSimEditor : public TGedFrame
	{
	public:
		AliEveTRDLoaderSimEditor(const TGWindow* p=0, Int_t width = 170, Int_t height = 30, UInt_t options = kChildFrame, Pixel_t back = GetDefaultFrameBackground());
		~AliEveTRDLoaderSimEditor();

		virtual void	SetModel(TObject* obj);
		virtual void	Toggle(Int_t id);

	protected:
		AliEveTRDLoaderSim* fM;
		TGCheckButton *fLoadHits, *fLoadDigits, *fLoadClusters, *fLoadTracks;

		ClassDef(AliEveTRDLoaderSimEditor,1) // Editor for AliEveTRDLoaderSim
	};
#endif
