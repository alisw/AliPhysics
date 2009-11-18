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
#include "TString.h"

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
    AliITSSurveyToAlign(Int_t run = 0, Int_t repModSDD = 845069, Int_t repModVerSDD = 1, 
                        Int_t repLadSDD = 999999, Int_t repLadVerSDD = 1, Int_t repModSSD = 887877,
	    Int_t repModVerSSD =3, Int_t repLaddSSD = 980521, Int_t repLaddVerSSD = 2);
    AliITSSurveyToAlign(const AliITSSurveyToAlign& align); // copy constructor
    AliITSSurveyToAlign &operator = (const AliITSSurveyToAlign& align); //assignment operator
    virtual ~AliITSSurveyToAlign();

    void Run();
    Bool_t CreateAlignObjs();
    void CreateAlignObjDummySPD();

    void CreateAlignObjDummySDDModules();
    void CreateAlignObjDummySDDLadders();
    void CreateAlignObjSDDModules();
    void CreateAlignObjSDDLadders();
    Bool_t ApplyAlignObjSDD();

    void CreateAlignObjSSDModules();
    void CreateAlignObjDummySSDModules();
    void CreateAlignObjSSDLadders();
    Bool_t ApplyAlignObjSSDLadders();

private:
    Int_t   fRun;                         // the run number for the OCDB
    Int_t   fSDDModuleRepNumber;
    Int_t   fSDDModuleRepVersion;
    Int_t   fSDDLadderRepNumber;
    Int_t   fSDDLadderRepVersion;
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

    
    // these are tmp vars. 
    //to be removed later
    Int_t    fuidSDDm[260];      
    TString  fsymnameSDDm[260];
    Double_t fxSDDm[260];
    Double_t fySDDm[260];
    Double_t fzSDDm[260];
    Double_t fpsiSDDm[260];
    Double_t ftetSDDm[260];
    Double_t fphiSDDm[260];
    Int_t    fuidSDDl[36];
    TString  fsymnameSDDl[36];
    Double_t fxSDDl[36];
    Double_t fySDDl[36];
    Double_t fzSDDl[36];
    Double_t fpsiSDDl[36];
    Double_t ftetSDDl[36];
    Double_t fphiSDDl[36];
    
    ClassDef(AliITSSurveyToAlign,0);
};
#endif

