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
	const char * detectors = "ALL" ; 
	const char * rawFileName = "raw.root" ; 
	AliQADataMakerSteer qas ; 
	TStopwatch timer;
	timer.Start();
	qas.Run(detectors, AliQA::kRAWS, rawFileName);
	qas.Run(detectors, AliQA::kHITS);
	qas.Run(detectors, AliQA::kSDIGITS);
	qas.Run(detectors, AliQA::kDIGITS);
	qas.Run(detectors, AliQA::kRECPOINTS);
	qas.Run(detectors, AliQA::kESDS);
	timer.Stop();
	timer.Print();
}



