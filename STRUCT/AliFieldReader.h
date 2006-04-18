#ifndef ALIFIELDREADER_H
#define ALIFIELDREADER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////
//                                                       //
//  Class to generate the particles for the MC           //
//  The base class is empty                              //
//                                                       //
///////////////////////////////////////////////////////////

#include <TObject.h>
class  AliMagFMaps;
class  TNtuple;

class AliFieldReader : public TObject
{

 public:
    AliFieldReader();
    virtual ~AliFieldReader();
    virtual void Init();
    virtual void ReadMap();
    virtual void ReadMapSolenoid();
    virtual void SetCatalogueName(char* name = "goodfiles.list") {fCatalogueName = name;}
    virtual void SetStepSize(Float_t dz = 0.08) {fStepSize = dz;}
    virtual void SetZStart(Float_t zstart = 1383.) {fZStart = zstart;}    
    virtual void SetPolarity(Float_t pol = 1.) {fPolarity = pol;}
 protected:
    void MakeHtmlHeaderMain(FILE* file);
    void MakeHtmlHeaderPict(FILE* file);
    void MakeHtmlPict(FILE* file, char* name);
    void MakeHtmlTableEntry(FILE* htmlmain, char* fileName, char* htmlFile, Float_t x, Float_t y, Int_t i, Float_t bdl, Int_t ifile);
    void MakeHtmlTrailor(FILE* file);
    void ReadRegisterMap();
    void ReadRegisterMapSolenoid();
 protected:
    AliMagFMaps* fField;           // Pointer to calculated map
    TNtuple*     fMap;             // Pointer to measured map
    FILE*        fCatalogue;       // Pointer to file catalogue
    FILE*        fHtmlMain;        // Pointer to the html output file
    Int_t        fRegMap[200][3];  // Mapping between addresses and physical location 
    Float_t      fStepSize;        // Step size in z 
    Float_t      fZStart;          // Starting position in z
    Float_t      fDd;              // Distance between sensors
    Float_t      fDz;              // Distance between sensor planes  
    Float_t      fPolarity;        // Polarity of the field
    char*        fCatalogueName;   // Name of the catalogue
    ClassDef(AliFieldReader,1) 
};

#endif














