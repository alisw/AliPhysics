//-------------------------------------------------------------------------
//    Description: 
//    This class include into LRC library for Long-Range Correlation analysis
//    it is the NN class
//    calculates NN correlations for abs and rel var
//    Origin: Petr Naumenko, SPbSU-CERN, Petr.Naoumenko@cern.ch,
//    Andrey Ivanov (SPbSU-CERN), Igor Altsebeev (SPbSU-CERN) 
//-------------------------------------------------------------------------

#ifndef ALILRCNN_H
#define ALILRCNN_H
/*  See cxx source for full Copyright notice */


/* $Id$ */
#include "AliLRCAnalysis.h"
class TList;
class TH2D;
class AliLRCAnalysis;

class AliLRCNN : public AliLRCAnalysis{
    public:
	AliLRCNN();
	virtual ~AliLRCNN();
	AliLRCNN(char *name, TH2D* sourceHist);
	AliLRCNN(char *fileHistname, char *histname, char *profname);
	AliLRCNN(const TList * const LHist, char *histname, char *profname);
	void MakeHistogramm(char *name, TH2D* sourceHist);
	void MakeHistogramm(char *fileHistname, char *histname, char *profname);
	void MakeHistogramm(const TList * const LHist, char *histname, char *profname);
	ClassDef(AliLRCNN,0)                 // macro for rootcint
};



#endif

