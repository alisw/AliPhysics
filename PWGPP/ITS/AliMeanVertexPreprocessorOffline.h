#ifndef ALI_MEANVERTEX_PREPROCESSOROFFLINE_H
#define ALI_MEANVERTEX_PREPROCESSOROFFLINE_H

/* Copyright(c) 1998-2011, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
/* $Id: AliMeanVertexPreprocessor.h $ */


// Mean vertex preprocessor. 
// Davide Caffarri
//
#include "TNamed.h"
class AliCDBStorage;

class AliMeanVertexPreprocessorOffline: public TNamed 
{
  public:
	AliMeanVertexPreprocessorOffline();  
	virtual ~AliMeanVertexPreprocessorOffline();

	void  ProcessOutput(const char *filename, AliCDBStorage *db, Int_t runNb);
	Int_t GetStatus();

	void SetShowPlots(Bool_t showPlots){fShowPlots = showPlots;}

	void ModObject(const char* url, double zv, double zs, const char* commentAdd);

  private:
	AliMeanVertexPreprocessorOffline(const AliMeanVertexPreprocessorOffline & proc); // copy constructor	
	AliMeanVertexPreprocessorOffline& operator=(const AliMeanVertexPreprocessorOffline&); //operator
	
	enum EStatusCode_t {
	  kOk,
	  kInputError, /* open file error, missing histos */
	  kLowStatistics, /* too low statistics */
	  kStoreError, /* problems storing OCDB */
	  kWriteMeanVertexSPD, /*write MeanVertex computed online*/
	  kUseOfflineSPDvtx,  /*write SPD vtx offline*/
	  kLumiRegCovMatrixProblem, /*lumi region or cov matrix computation problems, default values set*/
	  kFitUpdateZFailed, /*problem in fitting the Z for the update in  CPass1*/
	  kNStatusCodes
	};

	Int_t fStatus; /* status code */
	static const Char_t *fgkStatusCodeName[kNStatusCodes];
	Bool_t fShowPlots; /* status code */
	
	ClassDef(AliMeanVertexPreprocessorOffline, 3);
};

#endif
