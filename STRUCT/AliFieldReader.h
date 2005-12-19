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
 private:
    void MakeHtmlHeaderMain(FILE*);
    void MakeHtmlHeaderPict(FILE*);
    void MakeHtmlPict(FILE*, char*);
    void MakeHtmlTableEntry(FILE* htmlmain, char* fileName, char* htmlFile, Float_t x, Float_t y, Int_t i, Float_t bdl, Int_t ifile);
    void MakeHtmlTrailor(FILE*);
    void ReadRegisterMap();
    void ReadRegisterMapSolenoid();
 protected:
    AliMagFMaps* fField;
    TNtuple*     fMap;
    FILE*        fCatalogue;
    FILE*        fHtmlMain;
    Int_t        fRegMap[200][3];
    Float_t      fStepSize;
    Float_t      fZStart;
    Float_t      fDd;
    Float_t      fDz;
    Float_t      fPolarity;
    char*        fCatalogueName;
    ClassDef(AliFieldReader,1) 
};

#endif














