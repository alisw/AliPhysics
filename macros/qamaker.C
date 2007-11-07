/*
 *  qamaker.C
 *  
 *
 *  Created by Yves Schutz on 17.10.07.
 *  Copyright 2007 __CERN__. All rights reserved.
 *
 */
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TStopwatch.h"
#include "AliQA.h"
#include "AliQADataMakerSteer.h"
#endif
void qamaker()
{
	AliQADataMakerSteer qas ; 
	TStopwatch timer;
	timer.Start();
	qas.Run(AliQA::kRAWS, "07000008224001.100.root");
	qas.Run(AliQA::kHITS);
	qas.Reset();
	qas.Run(AliQA::kSDIGITS);
	qas.Reset();
	qas.Run(AliQA::kDIGITS);
	qas.Reset();
	qas.Run(AliQA::kRECPOINTS);
	qas.Reset();
	qas.Run(AliQA::kESDS);
	timer.Stop();
	timer.Print();
}



