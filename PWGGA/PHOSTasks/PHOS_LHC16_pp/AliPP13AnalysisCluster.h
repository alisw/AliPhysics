#ifndef ALIPP13ANALYSISCLUSTER_H
#define ALIPP13ANALYSISCLUSTER_H

#include <TObject.h>

// --- Custom header files ---

// --- ROOT system ---
#include <TObject.h>

// --- AliRoot header files ---
#include <AliAODCaloCluster.h>

class AliPP13AnalysisCluster: public AliAODCaloCluster
{
public:
	AliPP13AnalysisCluster();
	AliPP13AnalysisCluster(const AliAODCaloCluster & c);

	void SetModule(Int_t t) { fModule = t; }
	void SetTrigger(Bool_t t) { fTrigger = t; }
	void SetTRU(Int_t t) { fTRU = t; }
	void SetTRUChannel(Int_t t, Int_t x, Int_t z)
	{
		fTRUCh = t;
		fTRUChX = x;
		fTRUChZ = z;
	}

	Int_t TRU() const { return fTRU; }
	Int_t TRUChannel() const { return fTRUCh; }
	Int_t TRUChannelX() const { return fTRUChX; }
	Int_t TRUChannelZ() const { return fTRUChZ; }
	Int_t IsTrigger() const { return fTrigger; }
	Int_t Module() const {return fModule; }

protected:
	// AliAODCaloCluster * fCluster;
	Int_t fTRU;      // 0-7
	Int_t fTRUCh;    // 0-111
	Int_t fTRUChX;   // 0-55
	Int_t fTRUChZ;   // 0-63
	Int_t fModule;
	Bool_t fTrigger;

private:
	AliPP13AnalysisCluster(const AliPP13AnalysisCluster &);
	AliPP13AnalysisCluster & operator = (const AliPP13AnalysisCluster &);

	ClassDef(AliPP13AnalysisCluster, 1)
};
#endif
