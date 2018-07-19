#include <AliAnalysisTaskCaloCellsPhysQA.h>
#include <AliCaloCellsPhysQA.h>

ClassImp(AliAnalysisTaskCaloCellsPhysQA);

//________________________________________________________________
AliAnalysisTaskCaloCellsPhysQA::AliAnalysisTaskCaloCellsPhysQA(
	const char * name, Int_t nmods, Int_t det, char * outfile
):
	AliAnalysisTaskCaloCellsQA(name, nmods, det, outfile)

{
	// NB: Override the original QA analysis
	//
	if (fCellsQA)
	{
		delete fCellsQA;
		fCellsQA = 0;
	}

	fCellsQA = new AliCaloCellsPhysQA(nmods, AliCaloCellsQA::kPHOS);
}
