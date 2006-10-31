#ifndef ALIEVE_TRDLoaderSingle_H
#define ALIEVE_TRDLoaderSingle_H

////////////////////////////////////////////////////////////////////////
//                                                                     // - ALIEVE implementation -
// Single event loader for the TRD detector
//    - TRDLoader - loader of TRD data for one event
//
// by A.Bercuci (A.Bercuci@gsi.de)   Fri Oct 27 2006
////////////////////////////////////////////////////////////////////////

#ifndef ALIEVE_TRDLoade_H
#include "TRDLoader.h"
#endif

namespace Alieve {
	class TRDLoaderSingle : public TRDLoader
	{
	public:
		TRDLoaderSingle(const Text_t* n="TRDLoaderSingle", const Text_t* t=0);

//	protected:
		Bool_t			GoToEvent(const int ev);
		Bool_t			Open(const char *file, const char *dir=".");

		ClassDef(TRDLoaderSingle, 1) // Alieve sigle event loader class for the TRD detector
	};
}
#endif
