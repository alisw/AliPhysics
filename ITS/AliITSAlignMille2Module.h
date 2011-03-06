#ifndef ALIITSALIGNMILLE2MODULE_H
#define ALIITSALIGNMILLE2MODULE_H 
/* Copyright(c) 2007-2009 , ALICE Experiment at CERN, All rights reserved. * 
 * See cxx source for full Copyright notice                               */  
 
/// \ingroup rec 
/// \class AliITSAlignMille2Module
/// \brief Class for alignment of ITS 
// 
// Authors: Marcello Lunardon 
//
// RS Converted static arrays fSensVolVolumeID and fSensVolIndex
// to TArrays in user transparent way.
//
/* $Id$  */ 
//#include <TString.h> 
//#include <TObject.h> 
#include <TNamed.h> 
#include <TArrayI.h> 
#include <TArrayS.h> 
#include <TObjArray.h> 
class AliITSAlignMille2;

class AliAlignObjParams; 
class TGeoHMatrix; 

class AliITSAlignMille2Module : public TNamed 
{ 
public: 
  enum {kSPD,kSDD,kSSD};
  enum {kMaxParGeom=6,kMaxParTot=9,kSensDefBit=BIT(14),kGlobalGeomBit=BIT(15),kNotInConfBit=BIT(16),kVdSDDSameLRBit=BIT(17)};
  enum {kDOFTX,kDOFTY,kDOFTZ,kDOFPS,kDOFTH,kDOFPH,kDOFT0,kDOFDVL,kDOFDVR};
  //
  AliITSAlignMille2Module(); 
  AliITSAlignMille2Module(UShort_t volid);
  AliITSAlignMille2Module(Int_t index, UShort_t volid, const char* symname, const TGeoHMatrix *m, Int_t nsv=0, const UShort_t *volidsv=NULL);
  AliITSAlignMille2Module(const AliITSAlignMille2Module& rhs); // copy constructor
  AliITSAlignMille2Module& operator=(const AliITSAlignMille2Module& rhs);  
  //
  virtual ~AliITSAlignMille2Module(); 
  //
  // geometry methods  
  Int_t        GetIndex()                             const {return fIndex;} 
  UShort_t     GetVolumeID()                          const {return fVolumeID;}  
  Int_t        GetNSensitiveVolumes()                 const {return fNSensVol;} 
  Int_t        GetSensVolIndex(Int_t at)              const {return fSensVolIndex[at];}
  Short_t      GetSensVolVolumeID(Int_t at)           const {return fSensVolVolumeID[at];}
  TGeoHMatrix *GetMatrix()                            const {return fMatrix;}
  void         GetLocalMatrix(TGeoHMatrix& mat)       const;
  UShort_t    *GetSensitiveVolumeVolumeID()           const {return (UShort_t*)fSensVolVolumeID.GetArray();}
  Float_t      GetSigmaFactor(Int_t i)                const {return fSigmaFactor[i];}
  Float_t      GetSigmaXFactor()                      const {return fSigmaFactor[0];}
  Float_t      GetSigmaYFactor()                      const {return fSigmaFactor[1];}
  Float_t      GetSigmaZFactor()                      const {return fSigmaFactor[2];}
  Int_t        GetNProcessedPoints()                  const {return fNProcPoints;}
  Bool_t       IsFreeDOF(Int_t dof)                   const {return dof<fNParTot && fParCstr[dof]>0;}
  Bool_t       AreSensorsProvided()                   const {return TestBit(kSensDefBit);}
  Bool_t       GeomParamsGlobal()                     const {return TestBit(kGlobalGeomBit);}
  Bool_t       IsNotInConf()                          const {return TestBit(kNotInConfBit);}
  Bool_t       IsVDriftLRSame()                       const {return TestBit(kVdSDDSameLRBit);}
  //
  Bool_t       IsIn(UShort_t volid)                   const;
  Bool_t       IsAlignable()                          const;
  Bool_t       BelongsTo(AliITSAlignMille2Module* parent) const;
  AliITSAlignMille2Module* GetParent()                const {return fParent;}
  AliITSAlignMille2Module* GetChild(Int_t i)          const {return (AliITSAlignMille2Module*)fChildren[i];}
  Int_t        GetNChildren()                         const {return fChildren.GetLast()+1;}
  //
  void         Print(Option_t* opt="")                const; 
  //
  void         EvaluateDOF();
  UShort_t     GetNParTot()                           const {return fNParTot;}
  UShort_t     GetNParFree()                          const {return fNParFree;}
  Float_t     *GetParVals()                           const {return fParVals;}
  Double_t     GetParVal(int par)                     const {return par<fNParTot ? fParVals[par] : 0;}
  Double_t     GetParErr(int par)                     const {return par<fNParTot ? fParErrs[par] : 0;}
  Double_t     GetParConstraint(int par)              const {return par<fNParTot ? fParCstr[par] : 0;}
  Int_t        GetParOffset(Int_t par)                const {return par<fNParTot ? fParOffs[par] : -1;}
  Int_t        GetDetType()                           const {return fDetType;}
  Bool_t       IsParConstrained(Int_t par)            const {return fParCstr[par]>0 && fParCstr[par]<fgkDummyConstraint;}
  Bool_t       IsSPD()                                const {return fDetType == kSPD;}
  Bool_t       IsSDD()                                const {return fDetType == kSDD;}
  Bool_t       IsSSD()                                const {return fDetType == kSSD;}
  Bool_t       IsSensor()                             const {return IsSensor(fVolumeID);}
  void         SetDetType(Int_t tp)                         {fDetType = tp;}
  void         SetParOffset(Int_t par,Int_t offs)           {fParOffs[par] = offs;}
  //
  void         SetParVals(Double_t *vl,Int_t npar);          
  void         SetParVal(Int_t par,Double_t v=0)            {fParVals[par] = v;}
  void         SetParErr(Int_t par,Double_t e=0)            {fParErrs[par] = e;}
  void         SetParConstraint(Int_t par,Double_t s=1e6)   {fParCstr[par] = s>0. ? s:0.0;}
  void         SetSigmaFactor(Int_t i,Float_t v)            {fSigmaFactor[i]=TMath::Max(0.001F,v);}
  void         SetSigmaXFactor(Float_t v)                   {SetSigmaFactor(0,v);}
  void         SetSigmaYFactor(Float_t v)                   {SetSigmaFactor(1,v);}
  void         SetSigmaZFactor(Float_t v)                   {SetSigmaFactor(2,v);}
  void         IncNProcessedPoints(Int_t step=1)            {fNProcPoints += step;}
  void         SetNProcessedPoints(Int_t v)                 {fNProcPoints = v;}
  void         SetParent(AliITSAlignMille2Module* par)      {fParent = par;}
  void         AddChild(AliITSAlignMille2Module* cld)       {fChildren.Add(cld);}
  void         SetFreeDOF(Int_t dof,Double_t cstr);
  void         SetSensorsProvided(Bool_t v=kTRUE)           {SetBit(kSensDefBit,v);}
  void         SetGeomParamsGlobal(Bool_t v=kTRUE)          {SetBit(kGlobalGeomBit,v);}
  void         SetNotInConf(Bool_t v=kTRUE)                 {SetBit(kNotInConfBit,v);}
  void         SetVDriftLRSame(Bool_t v=kTRUE)              {SetBit(kVdSDDSameLRBit,v);}
  Int_t        Set(Int_t index,UShort_t volid, const char* symname, const TGeoHMatrix *m,Int_t nsv=0, const UShort_t *volidsv=0);
  //
  void         AddSensitiveVolume(UShort_t volid);
  void         DelSensitiveVolume(Int_t at);
  void         DelSensitiveVolumes()                        {fNSensVol = 0;}
  //
  void         GetGeomParamsGlo(Double_t *pars);
  void         GetGeomParamsLoc(Double_t *pars);
  //
  TGeoHMatrix *GetSensitiveVolumeMatrix(UShort_t voluid);
  TGeoHMatrix *GetSensitiveVolumeOrigGlobalMatrix(UShort_t voluid);
  TGeoHMatrix *GetSensitiveVolumeModifiedMatrix(UShort_t voluid, const Double_t *delta,Bool_t local=kTRUE); 
  AliAlignObjParams *GetSensitiveVolumeMisalignment(UShort_t voluid, const AliAlignObjParams *a); 
  AliAlignObjParams *GetSensitiveVolumeMisalignment(UShort_t voluid, const Double_t *deltalocal); 
  //
  void         GetGlobalParams(Double_t *t, Double_t *r);
  void         GetGlobalParams(const Double_t *loct, const Double_t *locr,Double_t *t, Double_t *r);
  void         GetLocalParams(const Double_t *loct, const Double_t *locr,Double_t *t, Double_t *r);
  //
  void         GetSensVolGlobalParams(UShort_t volid,Double_t *t, Double_t *r);
  void         GetSensVolLocalParams(UShort_t volid,Double_t *t, Double_t *r);
  void         GetSensVolGlobalParams(UShort_t volid,const Double_t* loct,const Double_t* locr,Double_t *t, Double_t *r);
  void         GetSensVolLocalParams(UShort_t volid,const Double_t* loct,const Double_t* locr,Double_t *t, Double_t *r);
  //
  void         CalcDerivLocGlo(Double_t *deriv);
  void         CalcDerivGloLoc(Int_t idx,Double_t *deriv);
  void         CalcDerivGloLoc(Int_t sensVol,Int_t paridx,Double_t* derivative);
  void         CalcDerivCurLoc(Int_t sensVol,Int_t paridx,Double_t* derivative);
  void         CalcDerivDPosDPar(Int_t sensVol,const Double_t *pl,Double_t *deriv);
  //
  // forse non serve...
  AliAlignObjParams *GetSensitiveVolumeGlobalMisalignment(UShort_t voluid, const Double_t *deltalocal); 
  // mo' proviamo questo
  AliAlignObjParams *GetSensitiveVolumeTotalMisalignment(UShort_t voluid, const Double_t *deltalocal); 
  //
  static Int_t    GetIndexFromVolumeID(UShort_t volid);
  static UShort_t GetVolumeIDFromSymname(const Char_t *symname);
  static UShort_t GetVolumeIDFromIndex(Int_t index);
  static Bool_t   IsSensor(UShort_t vid);
  static Int_t    SensVolMatrix(UShort_t volid, TGeoHMatrix *m); 
  static Int_t    SensVolOrigGlobalMatrix(UShort_t volid, TGeoHMatrix *m); 

  //
protected:
  //
  void         AssignDetType();
  //
protected:
  //
  Int_t          fNSensVol;                       // number of sensor it refers to
  Int_t          fIndex;                          // aliroot index
  Int_t          fDetType;                        // Detector type
  UShort_t       fVolumeID;                       // aliroot volune ID
  UShort_t       fNParTot;                        // total number of parameters
  UShort_t       fNParFree;                       // number of free parameters
  TArrayS        fParOffs;                        // offsets of free params in the fit results
  Int_t          fNProcPoints;                    // number of processed points
  Float_t        fSigmaFactor[3];                 // multiplicative factor for referred sensor X,Y,Z error
  Float_t       *fParVals;                        // values of the fitted params
  Float_t       *fParErrs;                        // errors of the fitted params
  Float_t       *fParCstr;                        // Gaussian type constraint on parameter, 0 means fixed param
  //
  TArrayI        fSensVolIndex;                   // aliroot indices for sensors
  TArrayS        fSensVolVolumeID;                // aliroot indices for sensors volumes
  TGeoHMatrix   *fMatrix;                         // ideal TGeoHMatrix of the supermodule
  TGeoHMatrix   *fSensVolMatrix;                  // sensor's ideal matrices
  TGeoHMatrix   *fSensVolModifMatrix;             // sensor's modified matrices
  //
  AliITSAlignMille2Module* fParent;               // optional parent pointer
  TObjArray      fChildren;                       // array of optional children
  //
  static AliAlignObjParams fgTempAlignObj;        // temp.alignment object used as a buffer               
  static const Float_t fgkDummyConstraint;        // dummy (lose) contraint on parameter
  //
  ClassDef(AliITSAlignMille2Module, 0)
}; 

#endif 
