#ifndef ALIGAMMAMCDATAREADER_H
#define ALIGAMMAMCDATAREADER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id$ */

/* History of cvs commits:
 *
 * $Log$
 * Revision 1.2  2007/08/17 12:40:04  schutz
 * New analysis classes by Gustavo Conesa
 *
 * Revision 1.1.2.1  2007/07/26 10:32:09  schutz
 * new analysis classes in the the new analysis framework
 *
 *
 */

//_________________________________________________________________________
// Class for reading data (Kinematics and ESDs) in order to do prompt gamma correlations
//  Class created from old AliPHOSGammaJet
//  (see AliRoot versions previous Release 4-09)

//*-- Author: Gustavo Conesa (INFN-LNF)

// --- ROOT system ---
#include <TParticle.h> 
#include <TClonesArray.h> 
#include "AliGammaReader.h" 

class AliESDEvent ;

class AliGammaMCDataReader : public AliGammaReader {

public: 

  AliGammaMCDataReader() ; // ctor
  AliGammaMCDataReader(const AliGammaMCDataReader & g) ; // cpy ctor
  AliGammaMCDataReader & operator = (const AliGammaMCDataReader & g) ;//cpy assignment
  virtual ~AliGammaMCDataReader() {;} //virtual dtor

  void CreateParticleList(TObject * esd, TObject * stack, 
			  TClonesArray * plCh, TClonesArray * plEMCAL, TClonesArray * plPHOS,
			  TClonesArray * plPrimCh, TClonesArray * plPrimEMCAL, TClonesArray * plPrimPHOS);

  TParticle *  GetMotherParticle(Int_t label, AliStack *stack,  TString calo, TLorentzVector momentum) ;

  private:

  ClassDef(AliGammaMCDataReader,1)
} ;
 

#endif //ALIGAMMAMCDATAREADER_H



