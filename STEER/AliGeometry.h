#ifndef ALIGEOMETRY_H
#define ALIGEOMETRY_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  AliGeometry Base Class pABC               //
//                                            //
//  Author Yves Schutz     SUBATECH           //
//                                            //  
//                                            //
////////////////////////////////////////////////

// --- ROOT system ---

#include "TNamed.h"
#include "TVector3.h" 

// --- Standard library ---

// --- AliRoot header files ---
 
#include "AliRecPoint.h"

class AliRecPoint;

class AliGeometry : public TNamed {

public:

  AliGeometry() ;          // ctor            
  virtual ~AliGeometry() ; // dtor
 
  virtual void GetGlobal(const AliRecPoint * p, TVector3 & pos, TMatrix & mat) = 0   ; 
  virtual void GetGlobal(const AliRecPoint * p, TVector3 & pos) = 0 ; 

protected:

  AliGeometry(const Text_t* name, const Text_t* title) : TNamed (name,title) {}                                   

public:

  ClassDef(AliGeometry,1)  // description , version 1

};

#endif // ALIGEOMETRY_H
