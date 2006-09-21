#ifndef TFLUKASCORINGOPTION
#define TFLUKASCORINGOPTION

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                                                                           //
// Class to store FLUKA specific scoring options                             //
//                                                                           //
//                                                                           //
// Authors: Andreas Morsch <andreas.morsch@cern.ch>                          //
//          Barbara Dalena <Barbara.Dalena@ba.infn.it>                       //
//                                                                           // 
///////////////////////////////////////////////////////////////////////////////


#include <TNamed.h>
class TFlukaMCGeometry;


class TFlukaScoringOption : public TNamed
{
public:
   // Constructors
    TFlukaScoringOption();
    TFlukaScoringOption(const char* name, const char* sdum, Int_t nopfl, char* outfile, Float_t* what);
    TFlukaScoringOption(const char* name, const char* sdum, Int_t nopfl, char* outfile, Float_t* what,
                        const char* det1, const char* det2, const char* det3);
    // Getters
    Float_t     What(Int_t indi) const       {return fWhat[indi - 1];}
    Int_t       Par() const                  {return fNopfl;}
    char*       GetFileName() const          {return fOutFile;}
    Float_t     GetLun()  const              {return fLun;}
    const char* GetRegName(Int_t ndet);
    void        SetPar(Int_t val)            {fNopfl   = val;}
    void        SetFileName(char* outfile)   {fOutFile = outfile;}
    void        SetLun(Float_t lun)          {fLun     = lun;}
 
//
    void        WriteFlukaInputCards();
    void        WriteOpenFlukaFile();
    Int_t       GetRegionByName(const char* detname);

    static void SetStaticInfo(FILE* file, TFlukaMCGeometry* geom)
                {fgFile = file; fgGeom = geom;}

 protected:
    Int_t         fNopfl;        // Number of paramters
    Float_t       fWhat[12];     // WHAT()
    const char*   fName[3];      // Region Name
    char*         fOutFile;      // Output file
    Float_t       fLun;          // Logical Unit Number for Fluka output

    // Static
    static FILE*             fgFile;      // Output file
    static TFlukaMCGeometry* fgGeom;      // Pointer to geometry
 
 private:
    // Copy constructor and operator= declared but not implemented (-Weff++ flag)
    TFlukaScoringOption(const TFlukaScoringOption&);
    TFlukaScoringOption& operator=(const TFlukaScoringOption&);

    ClassDef(TFlukaScoringOption, 1)  // Fluka Scoring Option
};
	
#endif
	
