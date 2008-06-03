//-------------------------------------------------------------------------
//    Description: 
//    This class include into LRC library for Long-Range Correlation analysis
//    it is the PtPt class
//    calculates PtPt correlations for abs and rel var
//    Origin: Petr Naumenko, SPbSU-CERN, Petr.Naoumenko@cern.ch
//-------------------------------------------------------------------------
#ifndef ALILRCPTPT_H
#define ALILRCPTPT_H
/*  See cxx source for full Copyright notice */


/* $Id$ */

//-------------------------------------------------------------------------
//         LRC library for Long-Range Correlation analysis
//
//    Origin: Petr Naumenko, SPbSU-CERN, Petr.Naoumenko@cern.ch
//-------------------------------------------------------------------------

#include "AliLRCAnalysis.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TProfile.h"
#include "math.h"

class TH2D;

class AliLRCPtPt:public AliLRCAnalysis
{
    public:
	AliLRCPtPt();
	~AliLRCPtPt();
	AliLRCPtPt(char *name, TH2D* sourceHist, double ptd, TH2D* nb);
	AliLRCPtPt(char *fileHistname, char *histname, char *profname, double ptd, char *errhistname);	
	void MakeHistogramm(char *name, TH2D* sourceHist, double ptd, TH2D* nb);
	void MakeHistogramm(char *fileHistname, char *histname, char *profname, double ptd, char *errhistname) ;	
	ClassDef(AliLRCPtPt,0)                 // macro for rootcint

};

#endif

