#ifndef ALICALOCELLSPHYSQA_H
#define ALICALOCELLSPHYSQA_H

#include <AliCaloCellsQA.h>

class AliCaloCellsPhysQA: public AliCaloCellsQA
{
public:
	AliCaloCellsPhysQA(): // Required for IO
		AliCaloCellsQA(),
		fEmin(0),
		fNCells(0),
		fTimingCut(0)
	{
	}

	AliCaloCellsPhysQA(Int_t nmods, Int_t det = kPHOS, Int_t startRunNumber = 100000, Int_t endRunNumber = 300000):
		AliCaloCellsQA(nmods, det, startRunNumber, endRunNumber),
		fEmin(0.3),
		fNCells(3),
		fTimingCut(999)
	{
	}

	virtual void SetPhysicsClusterCut(Double_t emin = 0.3, Int_t ncells = 3, Double_t time = 12.5e-9) 
	{
		fEmin = emin;
		fNCells = ncells;
		fTimingCut = time;
	}

protected:
	virtual Int_t CheckClusterGetSM(AliVCluster * clus);

	Double_t fEmin;
	Int_t fNCells;
	Double_t fTimingCut;

private:
	ClassDef(AliCaloCellsPhysQA, 1);
};

#endif