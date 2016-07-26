/**
 * \file AliEmcalPythiaInfo.h
 * \brief Declaration of class AliEmcalPythiaInfo
 *
 * In this header file the class AliEmcalPythiaInfo is declared
 *
 * \date Feb 3, 2016
 */

#ifndef ALIPYTHIAINFO_H
#define ALIPYTHIAINFO_H

/* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TNamed.h>
#include "AliTLorentzVector.h"

/**
 * \class AliEmcalPythiaInfo
 * \brief Store some informaion about a Pythia event
 * \ingroup EMCALCOREFW
 *
 * This class is used to store some information
 * about a Pythia event (also for embedding)
 *
 * \author Leticia Cunqueiro <leticia.cunqueiro.mendez@cern.ch>,Muenster University
 * \author Salvatore Aiola <salvatore.aiola@cern.ch>, Yale University
 * \author Davide Caffari <davide.caffarri@cern.ch>,Cern 
*/
class AliEmcalPythiaInfo : public TNamed{

 public:
  AliEmcalPythiaInfo();
  AliEmcalPythiaInfo(const char* name);

  void SetPartonFlag6(Int_t flag6)                                        { fPartonFlag6   = flag6                ; }
  void SetParton6(Float_t pt, Float_t eta, Float_t phi, Float_t mass=0)   { fParton6.SetPtEtaPhiM(pt,eta,phi,mass); }
  void SetPartonFlag7(Int_t flag7)                                        { fPartonFlag7   = flag7                ; }
  void SetParton7(Float_t pt, Float_t eta, Float_t phi, Float_t mass=0)   { fParton7.SetPtEtaPhiM(pt,eta,phi,mass); }
  void SetPythiaEventWeight(Float_t ptWeight)                             { fPythiaEventWeight = ptWeight         ; }
  
  Int_t                    GetPartonFlag6()       const { return fPartonFlag6         ; }
  Float_t                  GetPartonPt6()         const { return fParton6.Pt()        ; }
  Float_t                  GetPartonEta6()        const { return fParton6.Eta()       ; }
  Float_t                  GetPartonPhi6()        const { return fParton6.Phi_0_2pi() ; }

  Int_t                    GetPartonFlag7()       const { return fPartonFlag7         ; }
  Float_t                  GetPartonPt7()         const { return fParton7.Pt()        ; }
  Float_t                  GetPartonEta7()        const { return fParton7.Eta()       ; }
  Float_t                  GetPartonPhi7()        const { return fParton7.Phi_0_2pi() ; }

  const AliTLorentzVector& GetParton6Momentum()   const { return fParton6             ; }
  const AliTLorentzVector& GetParton7Momentum()   const { return fParton7             ; }

  Float_t                  GetPythiaEventWeight() const { return fPythiaEventWeight   ; }

 private: 
  Int_t                fPartonFlag6       ; //!<! Parton 6 flag
  AliTLorentzVector    fParton6           ; //!<! Parton 6 momentum
  Int_t                fPartonFlag7       ; //!<! Parton 7 flag
  AliTLorentzVector    fParton7           ; //!<! Parton 7 momentum
  Float_t              fPythiaEventWeight ; //!<! The Pythia event weight
    
  AliEmcalPythiaInfo(const AliEmcalPythiaInfo&);
  AliEmcalPythiaInfo& operator=(const AliEmcalPythiaInfo&);
  
  /// \cond CLASSIMP
  ClassDef(AliEmcalPythiaInfo, 1);
  /// \endcond
};
#endif
