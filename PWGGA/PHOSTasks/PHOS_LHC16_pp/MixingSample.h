#ifndef MIXINGSAMPLE_H
#define MIXINGSAMPLE_H

// --- Root header files ---
#include <TList.h>
#include <TObject.h>
#include <TObjArray.h>

// --- Custom libraries ---
#include "PhotonSelection.h"

class MixingSample : public TObject
{
public:

	MixingSample();
	MixingSample(Int_t psize);
	virtual ~MixingSample();
	virtual TList * GetPool(EventFlags & e);
	virtual void UpdatePool(const TObjArray & clusters, EventFlags & e); 

protected:
	TList * fPool[10];
	Int_t fPoolSize;

private:
	MixingSample(const MixingSample &); // Not implemented
	MixingSample & operator = (const MixingSample &); // Not implemented

	ClassDef(MixingSample, 1)
};

#endif