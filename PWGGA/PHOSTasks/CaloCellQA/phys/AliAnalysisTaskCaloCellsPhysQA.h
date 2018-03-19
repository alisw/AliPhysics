#ifndef ALIANALYSISTASKCALOCELLSPHYSQA_H
#define ALIANALYSISTASKCALOCELLSPHYSQA_H


#include <AliAnalysisTaskCaloCellsQA.h>

#include <AliCaloCellsQA.h>
#include <AliAnalysisManager.h>
#include <AliInputEventHandler.h>

class AliAnalysisTaskCaloCellsPhysQA: public AliAnalysisTaskCaloCellsQA
{
public:
	AliAnalysisTaskCaloCellsPhysQA(): AliAnalysisTaskCaloCellsQA() {}
	AliAnalysisTaskCaloCellsPhysQA(const char * name, Int_t nmods = 10, Int_t det = AliAnalysisTaskCaloCellsQA::kPHOS, char * outfile = 0);
	
protected:
	ClassDef(AliAnalysisTaskCaloCellsPhysQA, 1);
};
#endif
