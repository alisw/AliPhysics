/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// $Id$

///////////////////////////////////////////////////////////////////////////
// Class AliAttribObj
// Generic handling of detector signal (calibration) attributes.
//
// This class is meant to provide an AliAttrib object which is derived
// from TObject such that it can be stored in e.g. TObjArray etc...
// and that it can be written out using the ROOT I/O machinery.
//
// Example :
// ---------
// AliAttrib a;
// a.SetGain(250.7);
// a.SetGain(1340,3);
// a.SetEdgeOn(3);
// a.SetOffset(-22.5,2);
// a.SetDead(1);
// a.Data();
//
// AliAttribObj b(a);
// b.Data();
//
// AliAttribObj c;
// c.Load(a);
// c.Data();
//
//--- Author: Nick van Eijndhoven 18-sep-2003 Utrecht University
//- Modified: NvE $Date$ Utrecht University
///////////////////////////////////////////////////////////////////////////

#include "AliAttribObj.h"
#include "Riostream.h"
 
ClassImp(AliAttribObj) // Class implementation to enable ROOT I/O
 
AliAttribObj::AliAttribObj() : TObject(),AliAttrib()
{
// Creation of an AliAttrib object and initialisation of parameters.
// Several values of the same type (e.g. gain) can be stored in different slots.
// If needed, the storage for values will be expanded automatically
// when entering values.
}
///////////////////////////////////////////////////////////////////////////
AliAttribObj::AliAttribObj(AliAttrib& a) : TObject(),AliAttrib(a)
{
// Creation of an AliAttrib object and initialisation of parameters.
// All attributes are initialised to the values of the input AliAttrib.
}
///////////////////////////////////////////////////////////////////////////
AliAttribObj::~AliAttribObj()
{
// Destructor to delete dynamically allocated memory
}
///////////////////////////////////////////////////////////////////////////
AliAttribObj::AliAttribObj(const AliAttribObj& a) : TObject(a),AliAttrib(a)
{
// Copy constructor
}
///////////////////////////////////////////////////////////////////////////
TObject* AliAttribObj::Clone(const char* name) const
{
// Make a deep copy of the current object and provide the pointer to the copy.
// This memberfunction enables automatic creation of new objects of the
// correct type depending on the object type, a feature which may be very useful
// for containers when adding objects in case the container owns the objects.

 AliAttribObj* att=new AliAttribObj(*this);
 if (name)
 {
  if (strlen(name))
  {
   cout << " *" << ClassName() << "::Clone* No support for SetName." << endl;
  }
 }
 return att;
}
///////////////////////////////////////////////////////////////////////////
