#ifndef ALIGAMMAMCREADER_H
#define ALIGAMMAMCREADER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id$ */

/* History of cvs commits:
 *
 * $Log$
 *
 */

//_________________________________________________________________________
// Class for reading data (Kinematics) in order to do prompt gamma correlations
//  Class created from old AliPHOSGammaJet
//  (see AliRoot versions previous Release 4-09)

//*-- Author: Gustavo Conesa (INFN-LNF)

// --- ROOT system ---
class TArrayI ;

// --- AliRoot system ---
#include "AliGammaReader.h" 
 
class AliGammaMCReader : public AliGammaReader {

public: 

  AliGammaMCReader() ; // ctor
  AliGammaMCReader(const AliGammaMCReader & g) ; // cpy ctor
  AliGammaMCReader & operator = (const AliGammaMCReader & g) ;//cpy assignment
  virtual ~AliGammaMCReader() ;//virtual dtor

  enum DecayType {kNoDecay, kGeantDecay, kDecay, kDecayGamma};
 
  void InitParameters();

  Bool_t  IsInEMCAL(Double_t phi, Double_t eta) const ;
  Bool_t  IsInPHOS(Double_t phi, Double_t eta) const ;

  Int_t    GetDecayPi0Flag() const {return fDecayPi0 ; }

  void Print(const Option_t * opt) const;
  
  void SetDecayPi0Flag(Int_t d){ fDecayPi0 = d ; }

  void SetCheckOverlapping(Bool_t check){fCheckOverlapping = check ;}
  Bool_t IsCheckOverlappingOn() const {return fCheckOverlapping ;}

  void SetPhotonStatus(TParticle* pphoton, TParticle* mother);

  void AddNeutralParticlesArray(TArrayI & array)  
  { fNeutralParticlesArray   = new TArrayI(array) ; }
  TArrayI * GetNeutralParticlesArray() const   {return  fNeutralParticlesArray;}

  private:
  
  void CaseDecayGamma(Int_t index, TParticle * particle, AliStack * stack,
		      TClonesArray * plEMCAL, Int_t &indexEMCAL,
		      TClonesArray * plPHOS, Int_t &indexPHOS);
  
  void CaseGeantDecay(TParticle * particle, AliStack * stack,
		      TClonesArray * plEMCAL, Int_t &indexEMCAL,
		      TClonesArray * plPHOS, Int_t &indexPHOS);
  
  void CasePi0Decay(TParticle * particle, TClonesArray * plEMCAL, Int_t &indexEMCAL,
		    TClonesArray * plPHOS, Int_t &indexPHOS);
  
  void CreateParticleList(TObject * stack, TObject * ,
			  TClonesArray * plCh, TClonesArray * plEMCAL, 
			  TClonesArray * plPHOS, TClonesArray * plParton,TClonesArray *,TClonesArray *);
  
  void FillListWithDecayGammaOrPi0(TParticle * pPi0, TParticle * pdaug0, TParticle * pdaug1,
				   TClonesArray * plEMCAL, Int_t &indexEMCAL,
				   TClonesArray * plPHOS, Int_t &indexPHOS);  
  void MakePi0Decay(TLorentzVector &p0, TLorentzVector &p1, TLorentzVector &p2);//, Double_t &angle);
  
  Bool_t SkipNeutralParticles(Int_t pdg) const ;
 
  private:
  
  Int_t      fDecayPi0; //Options to study pi0 decay
  Bool_t   fCheckOverlapping; // if True, check if gammas from decay overlapp in calorimeters.
  TArrayI * fNeutralParticlesArray ; //Do not keep neutral particles of the list.
 
  ClassDef(AliGammaMCReader,1)
} ;


#endif //ALIGAMMAMCREADER_H



