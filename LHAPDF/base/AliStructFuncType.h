#ifndef ALISTRUCTFUNCTYPE_H
#define ALISTRUCTFUNCTYPE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

// Helper class to interface pdflib and the TPythia 
// the c++ interface for Pythia
// Author: andreas.morsch@cern.ch

#include <TObject.h>
#include <TString.h>
typedef enum
{
    kCTEQ4L,  // cteq4l.LHgrid   
    kCTEQ4M,  // cteq4m.LHgrid   
    kCTEQ5L,  // cteq5l.LHgrid
    kCTEQ5M,  // cteq5m.LHgrid
    kGRVLO98, // GRV98lo.LHgrid   
    kCTEQ6,   // cteq6.LHpdf
    kCTEQ61,  // cteq61.LHpdf
    kCTEQ6m,  // cteq6m.LHpdf
    kCTEQ6l,  // cteq6l.LHpdf 
    kCTEQ6ll, // cteq6ll.LHpdf
    kCT10,    // CT10.LHgrid
    kCT10nlo  // CT10nlo.LHgrid
}
StrucFunc_t;


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
    static Int_t   PDFsetIndex(StrucFunc_t pdf);
    static TString PDFsetName(StrucFunc_t pdf);    
    ClassDef(AliStructFuncType,1) // Library for partonic energy loss
};

#endif

