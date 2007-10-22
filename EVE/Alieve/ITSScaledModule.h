// $Header$

#ifndef ALIEVE_ITSScaledModule_H
#define ALIEVE_ITSScaledModule_H

#include <Reve/Reve.h>
#include <Alieve/ITSModule.h>

#include <TQObject.h>

namespace Alieve {

/**************************************************************************/
// DigitScaleInfo
/**************************************************************************/

class DigitScaleInfo : public TQObject, public Reve::ReferenceBackPtr
{
public:
  enum StatType_e { ST_Occup, ST_Average, ST_Rms };

  // Bool_t           fAutoUpdatePalette;
private:
  DigitScaleInfo(const DigitScaleInfo&);            // Not implemented
  DigitScaleInfo& operator=(const DigitScaleInfo&); // Not implemented

protected:
  Int_t            fScale;    
  Int_t            fStatType;
  
  Bool_t           fSyncPalette;

public:
  DigitScaleInfo();
  virtual ~DigitScaleInfo(){}
    
  Int_t            GetScale() { return fScale; }
  void             ScaleChanged(Int_t s);

  Int_t            GetStatType() { return fStatType; }
  void             StatTypeChanged(Int_t t);

  Bool_t           GetSyncPalette(){return fSyncPalette;}
  void             SetSyncPalette(Bool_t x){fSyncPalette = x;}

  ClassDef(DigitScaleInfo, 1);
};

/**************************************************************************/
// ScaledDigit
/**************************************************************************/

class ScaledDigit : public TObject
{
public:
  Int_t N;
  Float_t sum;
  Float_t sqr_sum;
  Int_t min_i,min_j;
  Int_t max_i,max_j;

  ScaledDigit();
  ScaledDigit(Int_t di, Int_t dj);
  
  virtual void Dump() const;
}; 

/**************************************************************************/
// ITSScaledModule
/**************************************************************************/

class ITSScaledModule : public ITSModule
{
  friend class ITSSDSubEditor;
private:
  map<Int_t, ScaledDigit> fDigitsMap;  
  
  ITSScaledModule(const ITSScaledModule&);            // Not implemented
  ITSScaledModule& operator=(const ITSScaledModule&); // Not implemented

protected:
  Int_t       fNx;  //  per module 
  Int_t       fNz;

  Int_t       fNCx;  // per cell
  Int_t       fNCz;

  DigitScaleInfo* fScaleInfo;

public:
  ITSScaledModule(Int_t gid, ITSDigitsInfo* info, DigitScaleInfo* si );
  virtual ~ITSScaledModule();

  virtual void DigitSelected(Int_t idx);

  virtual void LoadQuads();
  void         SetQuadValues();

  void         SyncPalette();

  void         GetScaleData(Int_t& cnx, Int_t& cnz, Int_t& total);
  DigitScaleInfo*  GetScaleInfo(){ return fScaleInfo; }

  ClassDef(ITSScaledModule, 1);
}; // endclass ITSScaledModule

}

#endif
