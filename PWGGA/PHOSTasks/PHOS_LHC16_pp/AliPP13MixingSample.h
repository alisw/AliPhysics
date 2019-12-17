#ifndef ALIPP13MIXINGSAMPLE_H
#define ALIPP13MIXINGSAMPLE_H

// --- Root header files ---
#include <TList.h>
#include <TObject.h>
#include <TObjArray.h>

// --- Custom libraries ---
#include "AliPP13PhysicsSelection.h"

class AliPP13MixingSample : public TObject
{
public:

	AliPP13MixingSample();
	AliPP13MixingSample(Int_t psize);
	virtual ~AliPP13MixingSample();
	virtual TList * GetPool(EventFlags & e);
	virtual void UpdatePool(const TObjArray & clusters, EventFlags & e); 

protected:
	TList * fPool[10];
	Int_t fPoolSize;

private:
	AliPP13MixingSample(const AliPP13MixingSample &); // Not implemented
	AliPP13MixingSample & operator = (const AliPP13MixingSample &); // Not implemented

	ClassDef(AliPP13MixingSample, 1)
};

#endif
