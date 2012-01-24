//-------------------------------------------------------------------------
//    Description: 
//    This class include into LRC library for Long-Range Correlation analysis
//    it is the PtPt class
//    calculates PtPt correlations for abs and rel var
//    Origin: Petr Naumenko, SPbSU-CERN, Petr.Naoumenko@cern.ch,
//    Andrey Ivanov (SPbSU-CERN), Igor Altsebeev (SPbSU-CERN) 
//-------------------------------------------------------------------------
#ifndef ALILRCPTPT_H
#define ALILRCPTPT_H
/*  See cxx source for full Copyright notice */


/* $Id$ */


#include "AliLRCAnalysis.h"

class TProfile;
class TList;
class TH2D;

class AliLRCPtPt:public AliLRCAnalysis
{
    public:
	AliLRCPtPt();
	~AliLRCPtPt();
	AliLRCPtPt(char *name, TH2D* sourceHist, double ptd, TH2D* nb);
	AliLRCPtPt(char *name, TH2D* sourceHist, double ptd, TProfile* nb);
	AliLRCPtPt(char *fileHistname, char *histname, char *profname, double ptd, char *errhistname);	
	AliLRCPtPt(const TList * const LHist, char *histname, char *profname, char *ptdname, char *errhistname);
	void MakeHistogramm(char *name, TH2D* sourceHist, double ptd, TH2D* nb);
	void MakeHistogramm(char *name, TH2D* sourceHist, double ptd, TProfile* nb);
	void MakeHistogramm(char *fileHistname, char *histname, char *profname, double ptd, char *errhistname) ;
	void MakeHistogramm(const TList * const LHist, char *histname, char *profname, char *ptdname, char *errhistname) ;		
	ClassDef(AliLRCPtPt,0)                 // macro for rootcint

};

#endif

