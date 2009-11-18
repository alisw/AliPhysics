#ifndef ALIITSSURVEYTOALIGN_H
#define ALIITSSURVEYTOALIGN_H
/* Copyright(c) 2008-2010, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//////////////////////////////////////////////////////////////////////////
//   Class to convert survey tables in alignment objects
//   for SSD and SDD
//////////////////////////////////////////////////////////////////////////

#include "AliSurveyToAlignObjs.h"

class TClonesArray;
class TGeoMatrix;
class TSystem;
class TFile;
class AliSurveyObj;
class AliAlignObjParams;
class AliCDBStorage;
class AliCDBEntry;


class AliITSSurveyToAlign : public AliSurveyToAlignObjs
{

public:
    AliITSSurveyToAlign(Int_t run = 0, Int_t repSDD = 845069, Int_t repVerSDD = 1,  Int_t repModSSD = 887877,
	    Int_t repModVerSSD =3, Int_t repLaddSSD = 980521, Int_t repLaddVerSSD = 2);
    AliITSSurveyToAlign(const AliITSSurveyToAlign& align); // copy constructor
    AliITSSurveyToAlign &operator = (const AliITSSurveyToAlign& /* align */); //assignment operator
    virtual ~AliITSSurveyToAlign();

    void Run();
    Bool_t CreateAlignObjs();
    void CreateAlignObjDummySPD();
    void CreateAlignObjSDD();
    void CreateAlignObjDummySDD();
    void CreateAlignObjSSDModules();
    void CreateAlignObjDummySSDModules();
    void CreateAlignObjSSDLadders();
    Bool_t ApplyAlignObjSSDLadders();

private:
    Int_t   fRun;                         // the run number for the OCDB
    Int_t   fSDDrepNumber;
    Int_t   fSDDrepVersion;
    Int_t   fSSDModuleRepNumber;
    Int_t   fSSDModuleRepVersion;
    Int_t   fSSDLadderRepNumber;
    Int_t   fSSDLadderRepVersion;

    Double_t fSDDmeP[6][6];   //measured positions of ref. marks for current module 
    Double_t fSDDidP[6][3];      //ideal positions of ref. marks for current module
    Bool_t     fSDDisMe[6];

    static const Double_t fgkLocR[6][3]; //id. pos. of ref. marks in RS of right oriented modules
    static const Double_t fgkLocL[6][3]; //id. pos. of ref. marks in RS of lefr oriented modules

    void GetIdPosSDD(Int_t uid, Int_t layer, Int_t module, Int_t iPoint);
    void ReadPointNameSDD(const char str[], Int_t &iLayer, Int_t &iLader, Int_t &iModul, Int_t &iPoint) const;
    void ConvertToRSofModulesAndRotSDD(Int_t Layer, Int_t Module);
    void CalcShiftSDD(Double_t &x0,Double_t &y0,Double_t &z0) const;
    void CalcShiftRotSDD(Double_t &tet,Double_t &psi,Double_t &phi,Double_t &x0,Double_t &y0,Double_t &z0);

    ClassDef(AliITSSurveyToAlign,0);
};
#endif

