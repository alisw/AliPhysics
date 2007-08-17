#ifndef ALIGAMMAMCREADER_H
#define ALIGAMMAMCREADER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id$ */

/* History of cvs commits:
 *
 * $Log$
 * Revision 1.1.2.1  2007/07/26 10:32:09  schutz
 * new analysis classes in the the new analysis framework
 *
 *
 */

//_________________________________________________________________________
// Class for reading data (Kinematics) in order to do prompt gamma correlations
//  Class created from old AliPHOSGammaJet
//  (see AliRoot versions previous Release 4-09)

//*-- Author: Gustavo Conesa (INFN-LNF)

// --- ROOT system ---
#include <TParticle.h> 
#include <TClonesArray.h> 
#include "AliStack.h"
#include "AliGammaReader.h" 
 
class TH2F ; 

class AliGammaMCReader : public AliGammaReader {

public: 

  AliGammaMCReader() ; // ctor
  AliGammaMCReader(const AliGammaMCReader & g) ; // cpy ctor
  AliGammaMCReader & operator = (const AliGammaMCReader & g) ;//cpy assignment
  virtual ~AliGammaMCReader() {;} //virtual dtor

  enum decay_t {kNoDecay, kGeantDecay, kDecay};
 
  void InitParameters();

  Bool_t  IsInEMCAL(Double_t phi, Double_t eta) ;
  Bool_t  IsInPHOS(Double_t phi, Double_t eta) ;

  Int_t    GetDecayPi0Flag() const {return fDecayPi0 ; }
  Float_t  GetEMCALIPDistance()  {  return fEMCALIPDistance ; }
  Float_t  GetPHOSIPDistance()  {  return fPHOSIPDistance ; }
  Float_t  GetEMCALMinDistance()  {  return fEMCALMinDistance ; }
  Float_t  GetPHOSMinDistance()  {  return fPHOSMinDistance ; }

  void Print(const Option_t * opt)const;
  
  void SetDecayPi0Flag(Int_t d){ fDecayPi0 = d ; }
  void SetEMCALIPDistance(Float_t  d){  fEMCALIPDistance = d ; }
  void SetPHOSIPDistance(Float_t  d){  fPHOSIPDistance = d ; }
  void SetEMCALMinDistance(Float_t  d){  fEMCALMinDistance = d ; }
  void SetPHOSMinDistance(Float_t  d){  fPHOSMinDistance = d ; }

  void CreateParticleList(TObject * stack, TObject * ,
			  TClonesArray * plCh, TClonesArray * plEMCAL, 
			  TClonesArray * plPHOS, TClonesArray * plParton);
  void MakePi0Decay(TParticle * particle, TClonesArray * plEMCAL, Int_t &indexEMCAL,
		    TClonesArray * plPHOS, Int_t &indexPHOS);
  
  void Pi0Decay(TLorentzVector &p0, TLorentzVector &p1, TLorentzVector &p2, 
		Double_t &angle);

  void SetGeantDecay(TParticle * particle, AliStack * stack,
		     TClonesArray * plEMCAL, Int_t &indexEMCAL,
		     TClonesArray * plPHOS, Int_t &indexPHOS);
 private:

 
  Float_t      fEMCALIPDistance; //Calorimeter IP distance.
  Float_t      fPHOSIPDistance; //Calorimeter IP distance
  Float_t      fEMCALMinDistance; //Gamma decay minimum aperture.
  Float_t      fPHOSMinDistance; //Gamma decay minimum aperture.

  Int_t      fDecayPi0; //Decay Pi0.

  ClassDef(AliGammaMCReader,0)
} ;
 

#endif //ALIGAMMAMCREADER_H



