#ifndef ALIFIDUCIALCUT_H
#define ALIFIDUCIALCUT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

//_________________________________________________________________________
/// \class AliFiducialCut
/// \ingroup CaloTrackCorrelationsBase
/// \brief Store the acceptance cuts for clusters and tracks or particle objects
///
/// Class for track/cluster/particle acceptance selection
/// Selection in Central barrel, DCAL and PHOS.
///  
/// Several selection regions possible for the different detectors
///
/// More information can be found in this [twiki](https://twiki.cern.ch/twiki/bin/viewauth/ALICE/PhotonHadronCorrelations).
///
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, LPSC-IN2P3-CNRS
//_________________________________________________________________________

// --- ROOT system ---
#include <TObject.h> 
#include <TArrayF.h> 

class TString ;
//class TLorentzVector ;

class AliFiducialCut : public TObject {
  
public: 
  
  AliFiducialCut() ; // ctor
  virtual  ~AliFiducialCut() ;//virtual dtor
  
  void      InitParameters();

  Bool_t    CheckFiducialRegion(Float_t eta, Float_t phi,
                                const TArrayF* minphi, const TArrayF* maxphi,
                                const TArrayF* mineta, const TArrayF* maxeta) const ;

  Bool_t    IsInFiducialCut    (Float_t eta, Float_t phi, Int_t det) const ;
  
  void      DoCTSFiducialCut  (Bool_t b)     { fCTSFiducialCut   = b    ; }
  void      DoEMCALFiducialCut(Bool_t b)     { fEMCALFiducialCut = b    ; }
  void      DoPHOSFiducialCut (Bool_t b)     { fPHOSFiducialCut  = b    ; }
  void      DoDCALFiducialCut (Bool_t b)     { fDCALFiducialCut  = b    ; }
  
  Bool_t    GetCTSFiducialCutStatus()  const { return fCTSFiducialCut   ; }
  Bool_t    GetEMCALFiducialCut()      const { return fEMCALFiducialCut ; }
  Bool_t    GetPHOSFiducialCutStatus() const { return fPHOSFiducialCut  ; }
  Bool_t    GetDCALFiducialCut()       const { return fDCALFiducialCut  ; }

  void      SetSimpleCTSFiducialCut  (Float_t abseta, Float_t phimin, Float_t phimax) ;
  void      SetSimpleEMCALFiducialCut(Float_t abseta, Float_t phimin, Float_t phimax) ;
  void      SetSimplePHOSFiducialCut (Float_t abseta, Float_t phimin, Float_t phimax) ;
  void      SetSimpleDCALFiducialCut (Float_t abseta, Float_t phimin, Float_t phimax) ;
  void      SetDCALFiducialCut       (Float_t etaminFull , Float_t etamaxFull , Float_t phiminFull , Float_t phimaxFull ,
                                      Float_t etaminThird, Float_t etamaxThird, Float_t phiminThird, Float_t phimaxThird);
  
  void      Print(const Option_t * opt)const;
  
  void      AddCTSFidCutMaxEtaArray(Int_t size, Float_t* array)  
                                              { fCTSFidCutMaxEta->Set(size,array)   ; } 
  TArrayF * GetCTSFidCutMaxEtaArray() const   { return fCTSFidCutMaxEta             ; }
  
  void      AddCTSFidCutMaxPhiArray(Int_t size, Float_t* array)  
                                              { fCTSFidCutMaxPhi->Set(size,array)   ; }
  TArrayF * GetCTSFidCutMaxPhiArray() const   { return fCTSFidCutMaxPhi             ; }
  
  void      AddCTSFidCutMinEtaArray(Int_t size, Float_t* array)  
                                              { fCTSFidCutMinEta->Set(size,array)   ; } 
  TArrayF * GetCTSFidCutMinEtaArray() const   { return fCTSFidCutMinEta             ; }
  
  void      AddCTSFidCutMinPhiArray(Int_t size, Float_t* array)  
                                              { fCTSFidCutMinPhi->Set(size,array)   ; }
  TArrayF * GetCTSFidCutMinPhiArray() const   { return fCTSFidCutMinPhi             ; }
  
  void      AddEMCALFidCutMaxEtaArray(Int_t size, Float_t* array)  
                                              { fEMCALFidCutMaxEta->Set(size,array) ; } 
  TArrayF * GetEMCALFidCutMaxEtaArray() const { return fEMCALFidCutMaxEta           ; }
  
  void      AddEMCALFidCutMaxPhiArray(Int_t size, Float_t* array)  
                                              { fEMCALFidCutMaxPhi->Set(size,array) ; }
  TArrayF * GetEMCALFidCutMaxPhiArray() const { return fEMCALFidCutMaxPhi           ; }
  
  void      AddEMCALFidCutMinEtaArray(Int_t size, Float_t* array)  
                                              { fEMCALFidCutMinEta->Set(size,array) ; } 
  TArrayF * GetEMCALFidCutMinEtaArray() const { return fEMCALFidCutMinEta           ; }
  
  void      AddEMCALFidCutMinPhiArray(Int_t size, Float_t* array)  
                                              { fEMCALFidCutMinPhi->Set(size,array) ; }
  TArrayF * GetEMCALFidCutMinPhiArray() const { return fEMCALFidCutMinPhi           ; }
  
  void      AddPHOSFidCutMaxEtaArray(Int_t size, Float_t* array)  
                                              { fPHOSFidCutMaxEta->Set(size,array)  ; } 
  TArrayF * GetPHOSFidCutMaxEtaArray() const  { return fPHOSFidCutMaxEta            ; }
  
  void      AddPHOSFidCutMaxPhiArray(Int_t size, Float_t* array)  
                                              { fPHOSFidCutMaxPhi->Set(size,array)  ; }
  TArrayF * GetPHOSFidCutMaxPhiArray() const  { return fPHOSFidCutMaxPhi            ; }
  
  void      AddPHOSFidCutMinEtaArray(Int_t size, Float_t* array)  
                                              { fPHOSFidCutMinEta->Set(size,array)  ; } 
  TArrayF * GetPHOSFidCutMinEtaArray() const  { return fPHOSFidCutMinEta            ; }

  void      AddPHOSFidCutMinPhiArray(Int_t size, Float_t* array)  
                                              { fPHOSFidCutMinPhi->Set(size,array)  ; }
  TArrayF * GetPHOSFidCutMinPhiArray() const  { return fPHOSFidCutMinPhi            ; }
  
  void      AddDCALFidCutMaxEtaArray(Int_t size, Float_t* array)
  { fDCALFidCutMaxEta->Set(size,array) ; }
  TArrayF * GetDCALFidCutMaxEtaArray() const { return fDCALFidCutMaxEta           ; }
  
  void      AddDCALFidCutMaxPhiArray(Int_t size, Float_t* array)
  { fDCALFidCutMaxPhi->Set(size,array) ; }
  TArrayF * GetDCALFidCutMaxPhiArray() const { return fDCALFidCutMaxPhi           ; }
  
  void      AddDCALFidCutMinEtaArray(Int_t size, Float_t* array)
  { fDCALFidCutMinEta->Set(size,array) ; }
  TArrayF * GetDCALFidCutMinEtaArray() const { return fDCALFidCutMinEta           ; }
  
  void      AddDCALFidCutMinPhiArray(Int_t size, Float_t* array)
  { fDCALFidCutMinPhi->Set(size,array) ; }
  TArrayF * GetDCALFidCutMinPhiArray() const { return fDCALFidCutMinPhi           ; }

  enum detector { kEMCAL = 0, kPHOS = 1, kCTS = 2, kDCAL = 3, kDCALPHOS = 4 } ;
  
private:
  
  // Detector acceptance cuts
  
  Bool_t    fEMCALFiducialCut  ; ///< Apply fiducial cuts to EMCAL clusters
  Bool_t    fDCALFiducialCut   ; ///< Apply fiducial cuts to DCAL clusters
  Bool_t    fPHOSFiducialCut   ; ///< Apply fiducial cuts to PHOS clusters
  Bool_t    fCTSFiducialCut    ; ///< Apply fiducial cuts to  CTS tracks
  
  TArrayF * fCTSFidCutMinEta   ; ///< Take particles in CTS with eta > fCTSFidCutMinEta
  TArrayF * fCTSFidCutMinPhi   ; ///< Take particles in CTS with phi > fCTSFidCutMinPhi
  TArrayF * fCTSFidCutMaxEta   ; ///< Take particles in CTS with eta < fCTSFidCutMaxEta
  TArrayF * fCTSFidCutMaxPhi   ; ///< Take particles in CTS with phi > fCTSFidCutMaxPhi
  
  TArrayF * fEMCALFidCutMinEta ; ///< Take particles in EMCAL with eta > fEMCALFidCutMinEta
  TArrayF * fEMCALFidCutMinPhi ; ///< Take particles in EMCAL with phi > fEMCALFidCutMinPhi
  TArrayF * fEMCALFidCutMaxEta ; ///< Take particles in EMCAL with eta < fEMCALFidCutMaxEta
  TArrayF * fEMCALFidCutMaxPhi ; ///< Take particles in EMCAL with phi > fEMCALFidCutMaxPhi
  
  TArrayF * fPHOSFidCutMinEta  ; ///< Take particles in PHOS with eta > fPHOSFidCutMinEta
  TArrayF * fPHOSFidCutMinPhi  ; ///< Take particles in PHOS with phi > fPHOSFidCutMinPhi
  TArrayF * fPHOSFidCutMaxEta  ; ///< Take particles in PHOS with eta < fPHOSFidCutMaxEta
  TArrayF * fPHOSFidCutMaxPhi  ; ///< Take particles in PHOS with phi > fPHOSFidCutMaxPhi

  TArrayF * fDCALFidCutMinEta  ; ///< Take particles in DCAL with eta > fDCALFidCutMinEta
  TArrayF * fDCALFidCutMinPhi  ; ///< Take particles in DCAL with phi > fDCALFidCutMinPhi
  TArrayF * fDCALFidCutMaxEta  ; ///< Take particles in DCAL with eta < fDCALFidCutMaxEta
  TArrayF * fDCALFidCutMaxPhi  ; ///< Take particles in DCAL with phi > fDCALFidCutMaxPhi
  
  /// Copy constructor not implemented.
  AliFiducialCut(              const AliFiducialCut & fc) ; 
  
  /// Assignment operator not implemented.
  AliFiducialCut & operator = (const AliFiducialCut & fc) ; 
  
  /// \cond CLASSIMP
  ClassDef(AliFiducialCut,3) ;
  /// \endcond

} ;


#endif //ALIFIDUCIALCUT_H



