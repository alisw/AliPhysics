/// \class AliNanoAODHeader


#ifndef _ALINANOAODHEADER_H_
#define _ALINANOAODHEADER_H_

#include "AliVAODHeader.h"
#include "AliNanoAODStorage.h"
#include <map>


class AliNanoAODHeader : public AliVAODHeader, public AliNanoAODStorage
{
public:
  using AliVHeader::ClassName;
  AliNanoAODHeader();
  AliNanoAODHeader(Int_t size);
  AliNanoAODHeader(Int_t size, Int_t sizeInt);
  virtual ~AliNanoAODHeader(){;}


  // Interface methods
  // AliNanoAODHeader(const AliVHeader& evt); 
  AliNanoAODHeader& operator=(const AliNanoAODHeader& evt);
  
  virtual UShort_t  GetBunchCrossNumber()   const { return GetVarInt(fBunchCrossNumber);}
  virtual UInt_t    GetOrbitNumber()        const { return GetVarInt(fOrbitNumber);}
  virtual UInt_t    GetPeriodNumber()       const { return GetVarInt(fPeriodNumber);}
  virtual UInt_t    GetTimeStamp()          const { NotImplemented();return 0;}
  virtual ULong64_t GetTriggerMask()        const { NotImplemented();return 0;}
  virtual ULong64_t GetTriggerMaskNext50()  const { NotImplemented();return 0;}
  virtual UChar_t   GetTriggerCluster()     const { NotImplemented();return 0;}
  virtual UInt_t    GetEventType()          const { NotImplemented();return 0;}

  virtual Bool_t   InitMagneticField()             const    { NotImplemented(); return 0;};
  virtual void     SetRunNumber(Int_t n)                    { SetVarInt(fRunNumber, n); };
  virtual void     SetMagneticField(Double_t /*magFld*/)        {NotImplemented();};
  virtual void     SetMuonMagFieldScale(Double_t /*magFldScl*/) {NotImplemented();};
  virtual void     SetDiamond(Float_t */*xy[2]*/,Float_t */*cov[3]*/) {NotImplemented();}; 
  virtual void     SetDiamondZ(Float_t /*z*/, Float_t /*sig2z*/)    {NotImplemented();};
  virtual Double_t GetMuonMagFieldScale()          const    {NotImplemented(); return 0;};
  virtual Double_t GetDiamondX()                   const    {NotImplemented(); return 0;};
  virtual Double_t GetDiamondY()                   const    {NotImplemented(); return 0;};
  virtual Double_t GetDiamondZ()                   const    {NotImplemented(); return 0;};
  virtual void     GetDiamondCovXY(Float_t */*cov[3]*/) const    {NotImplemented();};
  virtual Double_t GetSigma2DiamondX()             const    {NotImplemented(); return 0;};
  virtual Double_t GetSigma2DiamondY()             const    {NotImplemented(); return 0;};
  virtual Double_t GetSigma2DiamondZ()             const    {NotImplemented(); return 0;};

  virtual void     SetBunchCrossNumber(UShort_t nBx )       { SetVarInt(fBunchCrossNumber, nBx); };
  virtual void     SetOrbitNumber(UInt_t nOr)               { SetVarInt(fOrbitNumber, nOr); };
  virtual void     SetPeriodNumber(UInt_t nPer)             { SetVarInt(fPeriodNumber, nPer); };
  virtual void     SetTimeStamp(UInt_t /* t */)               {NotImplemented(); };
  virtual void     SetEventType(UInt_t /* evttype */)         {NotImplemented(); };
  virtual void     SetTriggerMask(ULong64_t /* trigMsk */)    {NotImplemented(); };
  virtual void     SetTriggerMaskNext50(ULong64_t /* trigMsk */) {NotImplemented(); };
  virtual void     SetTriggerCluster(UChar_t /* trigClus */)  {NotImplemented(); };
  virtual void     SetFiredTriggerClasses(TString varlist);
  virtual TString  GetFiredTriggerClasses() const;
  virtual Double_t GetZDCN1Energy()         const             {NotImplemented(); return 0;};
  virtual Double_t GetZDCP1Energy()         const             {NotImplemented(); return 0;};
  virtual Double_t GetZDCN2Energy()         const             {NotImplemented(); return 0;};
  virtual Double_t GetZDCP2Energy()         const             {NotImplemented(); return 0;};

  virtual Double_t GetZDCEMEnergy(Int_t /* i */) const            {NotImplemented(); return 0;};
  virtual Int_t    GetNumberOfESDTracks()  const              { return GetVarInt(fNumberOfESDTracks); };
  virtual Int_t    GetNumberOfTPCClusters() const {NotImplemented(); return 0;};
  virtual Int_t    GetNumberOfTPCTracks()   const {NotImplemented(); return 0;};
  
  virtual UInt_t   GetNumberOfITSClusters(Int_t /* ilay */) const {NotImplemented(); return 0;};
  virtual Float_t  GetT0spread(Int_t i)               const { return GetVar(fT0Spread[i]); };

  virtual AliCentrality* GetCentralityP()  const {NotImplemented(); return 0;};
  virtual Double_t GetCentralityV0M() const { return GetVar(fCentr); }
  virtual Double_t GetCentralityTRK() const { return GetVar(fCentrTRK); }
  virtual Double_t GetCentralityCL0() const { return GetVar(fCentrCL0); }
  virtual Double_t GetCentralityCL1() const { return GetVar(fCentrCL1); }
  
  virtual AliEventplane* GetEventplaneP()  const {NotImplemented(); return 0;};
  virtual Double_t       GetEventplane()     const {NotImplemented(); return 0;};
  virtual const Float_t* GetVZEROEqFactors() const {NotImplemented(); return 0;};
  virtual Float_t        GetVZEROEqFactors(Int_t /* i */) const {NotImplemented(); return 0;};
  virtual void           SetVZEROEqFactors(const Float_t* /*factors*/) {NotImplemented(); } 

  virtual UInt_t GetOfflineTrigger()  { return GetVarInt(fOfflineTrigger);};

  virtual void Print(Option_t* /*option = ""*/) const  {Printf("I'm a special header!");}
 
  virtual void Clear(Option_t * opt) ;

  virtual Int_t  GetIRInt2ClosestInteractionMap()                  const {NotImplemented(); return 0;};
  virtual Int_t  GetIRInt1ClosestInteractionMap(Int_t /*gap = 3*/) const {NotImplemented(); return 0;};


  virtual Int_t     GetRefMultiplicity()    const { NotImplemented(); return 0; }

  Double_t  GetMagneticField()      const { return GetVar(fMagField); }
  Double_t  GetCentrality () const; 
  Double_t  GetCentr (const char *x) const; 
  virtual Int_t  GetRunNumber() const; 

  TString GetCentralityMethod() const {return fCentralityMethod;}
  void SetCentralityMethod(const char * method)  { fCentralityMethod = method; } 

  void SetBunchCrossNumberIndex(Int_t var) { fBunchCrossNumber     = var; }
  void SetOrbitNumberIndex(Int_t var) { fOrbitNumber     = var; }
  void SetPeriodNumberIndex(Int_t var) { fPeriodNumber     = var; }
  void SetCentrIndex      (Int_t var) { fCentr     = var; }
  void SetCentrTRKIndex   (Int_t var) { fCentrTRK  = var; }
  void SetCentrCL0Index   (Int_t var) { fCentrCL0  = var; }
  void SetCentrCL1Index   (Int_t var) { fCentrCL1  = var; }
  void SetFiredTriggerClassesIndex   (Int_t var) { fFiredTriggerClasses  = var; }
  void SetMagFieldIndex   (Int_t var) { fMagField  = var; }
  void SetOfflineTriggerIndex   (Int_t var) { fOfflineTrigger  = var; }
  void SetRunNumberIndex  (Int_t var) { fRunNumber = var; }
  void SetT0SpreadIndex  (Int_t i, Int_t var) { fT0Spread[i] = var; }
  void SetNumberOfESDTracksIndex  (Int_t var) { fNumberOfESDTracks = var; }

  
  Int_t GetBunchCrossNumberIndex      () { return fBunchCrossNumber     ; }
  Int_t GetOrbitNumberIndex      () { return fOrbitNumber     ; }
  Int_t GetPeriodNumberIndex      () { return fPeriodNumber     ; }
  Int_t GetCentrIndex      () { return fCentr     ; }
  Int_t GetCentrTRKIndex   () { return fCentrTRK  ; }
  Int_t GetCentrCL0Index   () { return fCentrCL0  ; }
  Int_t GetCentrCL1Index   () { return fCentrCL1  ; }
  Int_t GetFiredTriggerClassesIndex () { return fFiredTriggerClasses  ; }
  Int_t GetMagFieldIndex   () { return fMagField  ; }
  Int_t GetOfflineTriggerIndex () { return fOfflineTrigger  ; }
  Int_t GetRunNumberIndex  () { return fRunNumber ; }
  Int_t GetT0SpreadIndex  (Int_t i) { return fT0Spread[i]; }
  Int_t GetNumberOfESDTracksIndex  () { return fNumberOfESDTracks ; }

  std::map<TString,int> GetMapCstVar () { return fMapCstVar; } 
  void SetMapCstVar (std::map<TString,int> cstmap) { fMapCstVar = cstmap; }
  std::map<TString,int> GetMapFiredTriggerClasses () { return fMapFiredTriggerClasses; } 
  void SetMapFiredTriggerClasses (TString trigClasses);

  Int_t GetVarIndex(TString varName); 
  
  // assume event filtering is done before and that nano AOD don't need the detector staus bits
  virtual Bool_t IsDetectorOn(ULong_t /*detMask*/) const { return kTRUE; }
  
  static void SetFatalMode() { fFatalMode = kTRUE; }

  ClassDef(AliNanoAODHeader, 5)
private:
  void NotImplemented() const;

  TString fCentralityMethod;


  Int_t fBunchCrossNumber;      ///< index of stored variable
  Int_t fOrbitNumber;///< index of stored variable
  Int_t fPeriodNumber;///< index of stored variable
  Int_t fCentr;      ///< index of stored variable
  Int_t fCentrTRK;   ///< index of stored variable
  Int_t fCentrCL0;   ///< index of stored variable
  Int_t fCentrCL1;   ///< index of stored variable
  Int_t fFiredTriggerClasses;   ///< index of stored variable
  Int_t fMagField;   ///< index of stored variable
  Int_t fOfflineTrigger;///< index of stored variable
  Int_t fRunNumber;  ///< index of stored variable
  Int_t fT0Spread[4]; ///< index of stored variable
  Int_t fNumberOfESDTracks;  ///< index of stored variable
  
  static Bool_t fFatalMode; //! throw AliFatal if invalid field is accessed

  std::map<TString,int> fMapCstVar;///< Map of indexes of custom variables: CACHE THIS TO CONST INTs IN YOUR TASK TO AVOID CONTINUOUS STRING COMPARISONS

  std::map<TString,int> fMapFiredTriggerClasses;///< Map of indexes of fired trigger Classes
};

#endif /* _ALINANOAODHEADER_H_ */
