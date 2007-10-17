#ifndef ALIGAMMADATAREADER_H
#define ALIGAMMADATAREADER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id$ */

/* History of cvs commits:
 *
 * $Log$
 * Revision 1.2  2007/08/21 08:38:20  schutz
 * Missing heared file added
 *
 * Revision 1.1.2.1  2007/07/26 10:32:09  schutz
 * new analysis classes in the the new analysis framework
 *
 *
 */

//_________________________________________________________________________
// Class for reading data (ESDs) in order to do prompt gamma correlations
//  Class created from old AliPHOSGammaJet
//  (see AliRoot versions previous Release 4-09)

//*-- Author: Gustavo Conesa (INFN-LNF)

// --- ROOT system ---
#include <TParticle.h> 
#include <TClonesArray.h> 
#include "AliGammaReader.h" 

class AliESDEvent ;

class AliGammaDataReader : public AliGammaReader {

public: 

  AliGammaDataReader() ; // ctor
  AliGammaDataReader(const AliGammaDataReader & g) ; // cpy ctor
  AliGammaDataReader & operator = (const AliGammaDataReader & g) ;//cpy assignment
  virtual ~AliGammaDataReader() {;} //virtual dtor

  void CreateParticleList(TObject * esd, TObject *,TClonesArray * plCh, 
			  TClonesArray * plEMCAL, TClonesArray * plPHOS, TClonesArray *,TClonesArray *, TClonesArray *);
  
 private:

  ClassDef(AliGammaDataReader,1)
} ;
 

#endif //ALIGAMMADATAREADER_H



