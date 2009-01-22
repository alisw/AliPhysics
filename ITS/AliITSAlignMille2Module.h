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
class AliITSAlignMille2;

class AliAlignObjParams; 
class TGeoHMatrix; 

class AliITSAlignMille2Module : public TNamed 
{ 
public: 
  AliITSAlignMille2Module(); 
  AliITSAlignMille2Module(UShort_t volid);
  AliITSAlignMille2Module(Int_t index, UShort_t volid, char* symname, TGeoHMatrix *m, Int_t nsv=0, UShort_t *volidsv=NULL);
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
  Bool_t       IsFreeDOF(Int_t dof)                   const {return TestBit(1<<dof);}
  UInt_t       GetFreePattern()                       const {return TestBits(0x3f);}
  Bool_t       AreSensorsProvided()                   const {return TestBit(1<<10);}
  Bool_t       IsIn(UShort_t volid)                   const;
  Bool_t       IsAlignable()                          const;
  Bool_t       BelongsTo(AliITSAlignMille2Module* parent) const;
  AliITSAlignMille2Module* GetParent()                const {return fParent;}
  void         Print(Option_t* opt="")                const; 
  //
  void         SetSigmaFactor(Int_t i,Float_t v)            {fSigmaFactor[i]=v;}
  void         SetSigmaXFactor(Float_t v)                   {fSigmaFactor[0]=v;}
  void         SetSigmaYFactor(Float_t v)                   {fSigmaFactor[1]=v;}
  void         SetSigmaZFactor(Float_t v)                   {fSigmaFactor[2]=v;}
  void         IncNProcessedPoints(Int_t step=1)            {fNProcPoints += step;}
  void         SetNProcessedPoints(Int_t v)                 {fNProcPoints = v;}
  void         SetParent(AliITSAlignMille2Module* par)      {fParent = par;}
  void         SetFreeDOF(Int_t dof,Bool_t free=kTRUE)      {SetBit(1<<dof,free);}
  void         SetSensorsProvided(Bool_t v=kTRUE)            {SetBit(1<<10,v);}
  Int_t        Set(Int_t index,UShort_t volid,char* symname,TGeoHMatrix *m,Int_t nsv=0,UShort_t *volidsv=0);
  //
  void         AddSensitiveVolume(UShort_t volid);
  void         DelSensitiveVolume(Int_t at);
  //
  TGeoHMatrix *GetSensitiveVolumeMatrix(UShort_t voluid);
  TGeoHMatrix *GetSensitiveVolumeOrigGlobalMatrix(UShort_t voluid);
  TGeoHMatrix *GetSensitiveVolumeModifiedMatrix(UShort_t voluid, Double_t *delta,Bool_t local=kTRUE); 
  AliAlignObjParams *GetSensitiveVolumeMisalignment(UShort_t voluid, AliAlignObjParams *a); 
  AliAlignObjParams *GetSensitiveVolumeMisalignment(UShort_t voluid, Double_t *deltalocal); 
  //
  // forse non serve...
  AliAlignObjParams *GetSensitiveVolumeGlobalMisalignment(UShort_t voluid, Double_t *deltalocal); 
  // mo' proviamo questo
  AliAlignObjParams *GetSensitiveVolumeTotalMisalignment(UShort_t voluid, Double_t *deltalocal); 
  //
  static Int_t GetIndexFromVolumeID(UShort_t volid);
  static UShort_t GetVolumeIDFromSymname(const Char_t *symname);
  static UShort_t GetVolumeIDFromIndex(Int_t index);
  //
protected:
  //
  Int_t        SensVolMatrix(UShort_t volid, TGeoHMatrix *m); 
  Int_t        SensVolOrigGlobalMatrix(UShort_t volid, TGeoHMatrix *m); 
  //
protected:
  //
  Int_t          fNSensVol;                       // number of sensor it refers to
  Int_t          fIndex;                          // aliroot index
  UShort_t       fVolumeID;                       // aliroot volune ID
  Int_t          fNProcPoints;                    // number of processed points
  Float_t        fSigmaFactor[3];                 // multiplicative factor for referred sensor X,Y,Z error
  //
  TArrayI        fSensVolIndex;                   // aliroot indices for sensors
  TArrayS        fSensVolVolumeID;                // aliroot indices for sensors volumes
  TGeoHMatrix   *fMatrix;                         // ideal TGeoHMatrix of the supermodule
  TGeoHMatrix   *fSensVolMatrix;                  // sensor's ideal matrices
  TGeoHMatrix   *fSensVolModifMatrix;             // sensor's modified matrices
  //
  AliITSAlignMille2Module* fParent;               // optional parent pointer
  //
  static AliAlignObjParams fgTempAlignObj;         // temp.alignment object used as a buffer               
  //
  ClassDef(AliITSAlignMille2Module, 0)
}; 

#endif 
