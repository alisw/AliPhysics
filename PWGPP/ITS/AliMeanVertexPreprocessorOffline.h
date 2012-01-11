#ifndef ALI_MEANVERTEX_PREPROCESSOROFFLINE_H
#define ALI_MEANVERTEX_PREPRECESSOROFFLINE_H

/* Copyright(c) 1998-2011, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
/* $Id: AliMeanVertexPreprocessor.h $ */


// Mean vertex preprocessor. 
// Davide Caffarri
//
#include "TNamed.h"

class AliMeanVertexPreprocessorOffline: public TNamed 
{
  public:
	AliMeanVertexPreprocessorOffline();  
	virtual ~AliMeanVertexPreprocessorOffline();

	void  ProcessOutput(const char *filename, const char *dbString, Int_t runNb);


  private:
	AliMeanVertexPreprocessorOffline(const AliMeanVertexPreprocessorOffline & proc); // copy constructor	
	AliMeanVertexPreprocessorOffline& operator=(const AliMeanVertexPreprocessorOffline&); //operator
	

	ClassDef(AliMeanVertexPreprocessorOffline, 1);
};

#endif
