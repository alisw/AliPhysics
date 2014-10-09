#ifndef ALIITSALIGNMILLEMODULE_H
#define ALIITSALIGNMILLEMODULE_H 
/* Copyright(c) 2007-2009 , ALICE Experiment at CERN, All rights reserved. * 
 * See cxx source for full Copyright notice                               */  
 
/// \ingroup rec 
/// \class AliITSAlignMilleModule 
/// \brief Class for alignment of ITS 
// 
// Authors: Marcello Lunardon 

/* $Id$  */ 
//#include <TString.h> 
//#include <TObject.h> 
#include <TNamed.h> 


class AliAlignObjParams; 
class TGeoHMatrix; 

class AliITSAlignMilleModule : public TNamed 
{ 
public: 
  AliITSAlignMilleModule(); 
  AliITSAlignMilleModule(UShort_t volid); // basic single volume constructor
  AliITSAlignMilleModule(Int_t index, UShort_t volid, char* symname, const TGeoHMatrix *m, Int_t nsv=0, const UShort_t *volidsv=NULL); // general constructor

  AliITSAlignMilleModule(const AliITSAlignMilleModule& rhs); // copy constructor
  AliITSAlignMilleModule& operator=(const AliITSAlignMilleModule& rhs);  
    
  virtual ~AliITSAlignMilleModule(); 
   
  // geometry methods  
  Int_t     GetIndex() const {return fIndex;} 
  UShort_t  GetVolumeID() const {return fVolumeID;}  
  Int_t     GetNSensitiveVolumes() const {return fNSensVol;} 
  TGeoHMatrix *GetMatrix() const {return fMatrix;} 
  const UShort_t *GetSensitiveVolumeVolumeID() const {return fSensVolVolumeID;}

  Int_t     Set(Int_t index, UShort_t volid, char* symname, const TGeoHMatrix * const m, Int_t nsv=0, const UShort_t *volidsv=NULL); // initialize a super module
  
  // util
  static Int_t GetIndexFromVolumeID(UShort_t volid);
  static UShort_t GetVolumeIDFromSymname(const Char_t *symname);
  static UShort_t GetVolumeIDFromIndex(Int_t index);

  // methods
  Bool_t    IsIn(UShort_t volid) const;
  TGeoHMatrix *GetSensitiveVolumeMatrix(UShort_t voluid);
  TGeoHMatrix *GetSensitiveVolumeOrigGlobalMatrix(UShort_t voluid);
  TGeoHMatrix *GetSensitiveVolumeModifiedMatrix(UShort_t voluid, const Double_t * const deltalocal); 
  AliAlignObjParams *GetSensitiveVolumeMisalignment(UShort_t voluid, const AliAlignObjParams *a); 
  AliAlignObjParams *GetSensitiveVolumeMisalignment(UShort_t voluid, const Double_t * const deltalocal); 
  // forse non serve...
  AliAlignObjParams *GetSensitiveVolumeGlobalMisalignment(UShort_t voluid, const Double_t * const deltalocal); 
  // mo' proviamo questo
  AliAlignObjParams *GetSensitiveVolumeTotalMisalignment(UShort_t voluid, const Double_t * const deltalocal); 
  void      Print(Option_t*) const; 

protected:
  Int_t     SensVolMatrix(UShort_t volid, TGeoHMatrix *m); 
  Int_t     SensVolOrigGlobalMatrix(UShort_t volid, TGeoHMatrix *m); 
  void      AddSensitiveVolume(UShort_t volid);

private:
  static const Int_t fgkSensModules = 2198; // number of sensors
  Int_t          fNSensVol; ///
  Int_t          fIndex; ///
  UShort_t       fVolumeID; ///
  // il symname e' il nome del TNamed...
  Int_t          fSensVolIndex[fgkSensModules]; ///
  UShort_t       fSensVolVolumeID[fgkSensModules]; ///
  TGeoHMatrix   *fMatrix; /// ideal TGeoHMatrix of the supermodule
  TGeoHMatrix   *fSensVolMatrix; ///
  TGeoHMatrix   *fSensVolModifMatrix; ///
  AliAlignObjParams *fTempAlignObj; ///
	
  ClassDef(AliITSAlignMilleModule, 0)

}; 

#endif 
