//-*- Mode: C++ -*-
#ifndef ALIESDCentrality_H
#define ALIESDCentrality_H
/* This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

//*****************************************************
//   Class AliCentralitySelectionTask
//   author: Alberica Toia
//*****************************************************

#include "TNamed.h"

/** 
 * Class to hold the centrality for ESDs 
 */
class AliESDCentrality : public TNamed
{
 public:
  /** 
   * constructor
   */
  AliESDCentrality();  
  /** 
   * Destructor
   */
  ~AliESDCentrality();  
  /** 
   * Copy constructor 
   * 
   * @param cnt Object to copy from
   */
  AliESDCentrality(const AliESDCentrality& cnt); 
  /** 
   * Assignment operator
   * 
   * @param cnt Object to assign from
   * 
   * @return Reference to this object
   */
  AliESDCentrality& operator=(const AliESDCentrality& cnt);   

  /** 
   * set centrality result
   * 
   * @param cent Centrality result
   */
  void SetCentralityV0M(Float_t cent) {fCentralityV0M = cent;} 
  /** 
   * set centrality result
   * 
   * @param cent Centrality result
   */
  void SetCentralityFMD(Float_t cent) {fCentralityFMD = cent;}
  /** 
   * set centrality result
   * 
   * @param cent Centrality result
   */
  void SetCentralityTRK(Float_t cent) {fCentralityTRK = cent;}
  /** 
   * set centrality result
   * 
   * @param cent Centrality result
   */
  void SetCentralityTKL(Float_t cent) {fCentralityTKL = cent;}
  /** 
   * set centrality result
   * 
   * @param cent Centrality result
   */
  void SetCentralityCL0(Float_t cent) {fCentralityCL0 = cent;}
  /** 
   * set centrality result
   * 
   * @param cent Centrality result
   */
  void SetCentralityV0MvsFMD(Float_t cent) {fCentralityV0MvsFMD = cent;}
  /** 
   * set centrality result
   * 
   * @param cent Centrality result
   */
  void SetCentralityTKLvsV0M(Float_t cent) {fCentralityTKLvsV0M = cent;}
  /** 
   * set centrality result
   * 
   * @param cent Centrality result
   */
  void SetCentralityZEMvsZDC(Float_t cent) {fCentralityZEMvsZDC = cent;}

  /** 
   * Get the centrality result 
   * 
   * @param method Method to use 
   * 
   * @return Centrality 
   */
  Float_t GetCentralityPercentile(const char *method);
  /** 
   * Get the centrality class for 10% bins 
   * 
   * @param method Method
   * 
   * @return Bin number
   */
  Int_t   GetCentralityClass10(const char *method);
  /** 
   * Get the centrality class for 5% bins 
   * 
   * @param method Method
   * 
   * @return Bin number
   */
  Int_t   GetCentralityClass5(const char *method);
  /** 
   * Check if event is in a given centrality class
   * 
   * @param a       Lower bound
   * @param b       Upper bound
   * @param method  Method
   * 
   * @return true if so. 
   */
  Bool_t  IsEventInCentralityClass(Float_t a, Float_t b, const char *method);

 private:
  Float_t fCentralityV0M;   // Centrality from V0
  Float_t fCentralityFMD;   // Centrality from FMD
  Float_t fCentralityTRK;   // Centrality from tracks
  Float_t fCentralityTKL;   // Centrality from tracklets
  Float_t fCentralityCL0;   // Centrality from Clusters in layer 0
  Float_t fCentralityV0MvsFMD;   // Centrality from V0 vs FMD
  Float_t fCentralityTKLvsV0M;   // Centrality from tracklets vs V0
  Float_t fCentralityZEMvsZDC;   // Centrality from ZEM vs ZDC

  ClassDef(AliESDCentrality, 2)
};

inline
AliESDCentrality::AliESDCentrality() : TNamed("ESDCentrality", "Centrality"),
  fCentralityV0M(0),
  fCentralityFMD(0),
  fCentralityTRK(0),
  fCentralityTKL(0),
  fCentralityCL0(0),
  fCentralityV0MvsFMD(0),
  fCentralityTKLvsV0M(0),
  fCentralityZEMvsZDC(0)
{
  /// constructor
}

inline
AliESDCentrality::AliESDCentrality(const AliESDCentrality& cnt) : 
  TNamed(cnt), 
  fCentralityV0M(cnt.fCentralityV0M),
  fCentralityFMD(cnt.fCentralityFMD),
  fCentralityTRK(cnt.fCentralityTRK),
  fCentralityTKL(cnt.fCentralityTKL),
  fCentralityCL0(cnt.fCentralityCL0),
  fCentralityV0MvsFMD(cnt.fCentralityV0MvsFMD),
  fCentralityTKLvsV0M(cnt.fCentralityTKLvsV0M),
  fCentralityZEMvsZDC(cnt.fCentralityZEMvsZDC)
{
  /// Copy constructor
}

inline
AliESDCentrality& AliESDCentrality::operator=(const AliESDCentrality& c)
{
  /// Assignment operator
  if (this!=&c) {
    TNamed::operator=(c);
    fCentralityV0M = c.fCentralityV0M;
    fCentralityFMD = c.fCentralityFMD;
    fCentralityTRK = c.fCentralityTRK;
    fCentralityTKL = c.fCentralityTKL;
    fCentralityCL0 = c.fCentralityCL0;
    fCentralityV0MvsFMD = c.fCentralityV0MvsFMD;
    fCentralityTKLvsV0M = c.fCentralityTKLvsV0M;
    fCentralityZEMvsZDC = c.fCentralityZEMvsZDC;
  }

  return *this;
}

inline
AliESDCentrality::~AliESDCentrality()
{
  /// destructor
}

inline
Float_t AliESDCentrality::GetCentralityPercentile(const char *x)
{
  TString method = x;
  if(method.CompareTo("V0M")==0)      return fCentralityV0M;
  if(method.CompareTo("FMD")==0)      return fCentralityFMD;
  if(method.CompareTo("TRK")==0)      return fCentralityTRK;
  if(method.CompareTo("TKL")==0)      return fCentralityTKL;
  if(method.CompareTo("CL0")==0)      return fCentralityCL0;
  if(method.CompareTo("V0MvsFMD")==0) return fCentralityV0MvsFMD;
  if(method.CompareTo("TKLvsV0M")==0) return fCentralityTKLvsV0M;
  if(method.CompareTo("ZENvsZDC")==0) return fCentralityZEMvsZDC;
  return -1;
}

inline
Int_t AliESDCentrality::GetCentralityClass10(const char *x)
{
  TString method = x;
  if(method.CompareTo("V0M")==0)      return (Int_t) (fCentralityV0M / 10.0);
  if(method.CompareTo("FMD")==0)      return (Int_t) (fCentralityFMD / 10.0);
  if(method.CompareTo("TRK")==0)      return (Int_t) (fCentralityTRK / 10.0);
  if(method.CompareTo("TKL")==0)      return (Int_t) (fCentralityTKL / 10.0);
  if(method.CompareTo("CL0")==0)      return (Int_t) (fCentralityCL0 / 10.0);
  if(method.CompareTo("V0MvsFMD")==0) return (Int_t) (fCentralityV0MvsFMD / 10.0);
  if(method.CompareTo("TKLvsV0M")==0) return (Int_t) (fCentralityTKLvsV0M / 10.0);
  if(method.CompareTo("ZENvsZDC")==0) return (Int_t) (fCentralityZEMvsZDC / 10.0);
  return -1;
}

Int_t AliESDCentrality::GetCentralityClass5(const char *x)
{
 TString method = x;
  if(method.CompareTo("V0M")==0)      return (Int_t) (fCentralityV0M / 5.0);
  if(method.CompareTo("FMD")==0)      return (Int_t) (fCentralityFMD / 5.0);
  if(method.CompareTo("TRK")==0)      return (Int_t) (fCentralityTRK / 5.0);
  if(method.CompareTo("TKL")==0)      return (Int_t) (fCentralityTKL / 5.0);
  if(method.CompareTo("CL0")==0)      return (Int_t) (fCentralityCL0 / 5.0);
  if(method.CompareTo("V0MvsFMD")==0) return (Int_t) (fCentralityV0MvsFMD / 5.0);
  if(method.CompareTo("TKLvsV0M")==0) return (Int_t) (fCentralityTKLvsV0M / 5.0);
  if(method.CompareTo("ZENvsZDC")==0) return (Int_t) (fCentralityZEMvsZDC / 5.0);
  return -1;
}

inline
Bool_t AliESDCentrality::IsEventInCentralityClass(Float_t a, Float_t b, const char *x)
{
  TString method = x;
  if ((method.CompareTo("V0M")==0) && (fCentralityV0M >=a && fCentralityV0M < b)) return kTRUE;
  if ((method.CompareTo("FMD")==0) && (fCentralityFMD >=a && fCentralityFMD < b)) return kTRUE;
  if ((method.CompareTo("TRK")==0) && (fCentralityTRK >=a && fCentralityTRK < b)) return kTRUE;
  if ((method.CompareTo("TKL")==0) && (fCentralityTKL >=a && fCentralityTKL < b)) return kTRUE;
  if ((method.CompareTo("CL0")==0) && (fCentralityCL0 >=a && fCentralityCL0 < b)) return kTRUE;
  if ((method.CompareTo("V0MvsFMD")==0) && (fCentralityV0MvsFMD >=a && fCentralityV0MvsFMD < b)) return kTRUE;
  if ((method.CompareTo("TKLvsV0M")==0) && (fCentralityTKLvsV0M >=a && fCentralityTKLvsV0M < b)) return kTRUE;
  if ((method.CompareTo("ZEMvsZDC")==0) && (fCentralityZEMvsZDC >=a && fCentralityZEMvsZDC < b)) return kTRUE;
  else return kFALSE;
}




#endif //ALIESDCENTRALITY_H
