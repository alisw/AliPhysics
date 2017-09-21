// --- Custom libraries ---
#include "MixingSample.h"


ClassImp(MixingSample)

MixingSample::MixingSample():
	TObject(),
	fPool(),
	fPoolSize(0)
{}

MixingSample::MixingSample(Int_t psize):
	TObject(),
	fPool(),
	fPoolSize(psize)
{
	for (Int_t i = 0; i < 10; ++i)
	{
		fPool[i] = new TList();
		fPool[i]->SetOwner(kTRUE);
	}
}

MixingSample::~MixingSample()
{
	for (Int_t i = 0; i < 10; ++i)
	{
		if (fPool[i]) delete fPool[i];
	}
}

TList * MixingSample::GetPool(EventFlags & e)
{
	Int_t zbin = Int_t((e.vtxBest[2] + 10.) / 2.);
	if (zbin < 0) zbin = 0;
	if (zbin > 9) zbin = 9;

	return fPool[zbin];
}

void MixingSample::UpdatePool(const TObjArray & clusters, EventFlags & e)
{
	TList * pool = GetPool(e);

	if (clusters.GetEntries() > 0)
		pool->AddFirst(clusters.Clone());

	if (pool->GetEntries() > fPoolSize)
	{
		TObject * tmp = pool->Last();
		pool->RemoveLast();
		delete tmp;
	}
}

