#ifndef ALISTRUCTFUNCTYPE_H
#define ALISTRUCTFUNCTYPE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

// Helper class to interface pdflib and the TPythia 
// the c++ interface for Pythia
// Author: andreas.morsch@cern.ch

#include <TObject.h>

class AliStructFuncType : public TObject {
    
 public:
    AliStructFuncType(){;}
    virtual ~AliStructFuncType(){;}
    static void PdfSet(char parm[20][20], Double_t value[20]);
    static void StructA(Double_t xx, Double_t qq, Double_t a,
			Double_t& upv, Double_t& dnv, Double_t& usea,
			Double_t& dsea,
			Double_t& str, Double_t& chm, Double_t& bot,
			Double_t& top, Double_t& gl);
    ClassDef(AliStructFuncType,1) // Library for partonic energy loss
};

typedef enum
{
   kCTEQ4L     = 19170,    
   kCTEQ4M     = 19150,    
   kCTEQ5L     = 19070,
   kCTEQ5M     = 19050,
   kGRVLO98    = 80060,    
   kCTEQ6      = 10040,
   kCTEQ61     = 10100,
   kCTEQ6m     = 10050,
   kCTEQ6l     = 10041,
   kCTEQ6ll    = 10042
}
StrucFunc_t;

#endif

