#ifndef ALIEMCALTRIGGERPATCHINFO_H
#define ALIEMCALTRIGGERPATCHINFO_H

// $Id$

#include "TObject.h"

#include <TLorentzVector.h>
#include <TMath.h>
#include "AliEMCALTriggerTypes.h"
#include "AliEmcalTriggerSetupInfo.h"

class AliEMCALGeometry;
class TArrayI;

class AliEmcalTriggerPatchInfo: public TObject {
 public:
  AliEmcalTriggerPatchInfo();
  AliEmcalTriggerPatchInfo(const AliEmcalTriggerPatchInfo &p); 
  AliEmcalTriggerPatchInfo &operator=(const AliEmcalTriggerPatchInfo &p);
  virtual ~AliEmcalTriggerPatchInfo();


  Double_t GetPhiGeo() const { return fCenterGeo.Phi(); }
  Double_t GetPhiCM()  const { return fCenterMass.Phi(); }
  Double_t GetPhiMin() const { return fEdge1.Phi(); }
  Double_t GetPhiMax() const { return fEdge2.Phi(); }
  Double_t GetEtaGeo() const { return fCenterGeo.Eta(); }
  Double_t GetEtaCM()  const { return fCenterMass.Eta(); }
  Double_t GetEtaMin() const { return fEdge2.Eta(); }
  Double_t GetEtaMax() const { return fEdge1.Eta(); }
  Double_t GetPatchE() const { return fCenterGeo.E(); }
  Int_t    GetADCAmp() const { return fADCAmp; }
  Double_t GetADCAmpGeVRough() const { return (Double_t)fADCAmp * kEMCL1ADCtoGeV; }
  Int_t    GetTriggerBits() const { return fTriggerBits; }
  Int_t    GetEdgeCellX() const { return fEdgeCell[0]; }
  Int_t    GetEdgeCellY() const { return fEdgeCell[1]; }
  void     GetCellIndices( AliEMCALGeometry *geom, TArrayI *cells );
  
  Bool_t   IsJetLow() const { return (Bool_t)((fTriggerBits >> (fOffSet + kL1JetLow))&(!(fTriggerBits >> 25))&1); }
  Bool_t   IsJetHigh() const { return (Bool_t)((fTriggerBits >> (fOffSet + kL1JetHigh))&(!(fTriggerBits >> 25))&1); }
  Bool_t   IsMainTrigger() const { return (Bool_t)((fTriggerBits >> 24)&(!(fTriggerBits >> 25))&1); }
  Bool_t   IsJetLowSimple() const { return (Bool_t)((fTriggerBits >> (fOffSet + kL1JetLow))&(fTriggerBits >> 25)&1); }
  Bool_t   IsJetHighSimple() const { return (Bool_t)((fTriggerBits >> (fOffSet + kL1JetHigh))&(fTriggerBits >> 25)&1); }
  Bool_t   IsMainTriggerSimple() const { return (Bool_t)((fTriggerBits >> 24)&(fTriggerBits >> 25)&1); }
  Bool_t   IsOfflineSimple() const { return (Bool_t)((fTriggerBits >> 25)&1); }
  
  void SetCenterGeo( TVector3 &v, Double_t e ) { SetLorentzVector( fCenterGeo, v, e ); }
  void SetCenterGeo( TLorentzVector &v ) { fCenterGeo = v; }
  void SetCenterMass( TLorentzVector &v ) { fCenterMass = v; }
  void SetCenterMass( TVector3 &v, Double_t e ) { SetLorentzVector( fCenterMass, v, e ); }
  void SetEdge1( TLorentzVector &v ) { fEdge1 = v; }
  void SetEdge1( TVector3 &v, Double_t e ) { SetLorentzVector( fEdge1, v, e ); }
  void SetEdge2( TLorentzVector &v ) { fEdge2 = v; }
  void SetEdge2( TVector3 &v, Double_t e ) { SetLorentzVector( fEdge2, v, e ); }
  void SetADCAmp( Int_t a ) { fADCAmp = a; }
  void SetEdgeCell( Int_t x, Int_t y ) { fEdgeCell[0] = x; fEdgeCell[1] = y; }

  void SetLorentzVector( TLorentzVector &lv, TVector3 &v, Double_t e );

  void SetTriggerBits( Int_t i ) { fTriggerBits = i; }

  void SetOffSet(Int_t i)        { fOffSet      = i; }


 protected:
  TLorentzVector   &GetLorentzVector(const Double_t *vertex = 0)  const;

  TLorentzVector    fCenterGeo;                     // geometrical center
  TLorentzVector    fCenterMass;                    // CM
  TLorentzVector    fEdge1;                         // max eta/ min phi edge
  TLorentzVector    fEdge2;                         // min eta/ max phi edge
  Int_t             fADCAmp;                        // online ADC amplitude
  Int_t             fTriggerBits;                   //trigger bit mask
  Int_t             fEdgeCell[2];                   // cell "bottom lower" edge (min phi, max eta)
  Int_t             fOffSet;                        // offset of bit (different in data and MC)

  ClassDef(AliEmcalTriggerPatchInfo, 4) // Emcal particle class
};
#endif
