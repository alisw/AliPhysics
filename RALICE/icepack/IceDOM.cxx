/*******************************************************************************
 * Copyright(c) 2003, IceCube Experiment at the South Pole. All rights reserved.
 *
 * Author: The IceCube RALICE-based Offline Project.
 * Contributors are mentioned in the code where appropriate.
 *
 * Permission to use, copy, modify and distribute this software and its
 * documentation strictly for non-commercial purposes is hereby granted
 * without fee, provided that the above copyright notice appears in all
 * copies and that both the copyright notice and this permission notice
 * appear in the supporting documentation.
 * The authors make no claims about the suitability of this software for
 * any purpose. It is provided "as is" without express or implied warranty.
 *******************************************************************************/

// $Id$

///////////////////////////////////////////////////////////////////////////
// Class IceDOM
// Signal/Hit handling of a generic IceCube Digital Optical Module (DOM).
// Basically this class provides an IceCube tailored user interface
// to the functionality of the class AliDevice via the generic IceGOM class.
//
// See IceGOM for some usage examples.
//
//--- Author: Nick van Eijndhoven 23-jun-2004 Utrecht University
//- Modified: NvE $Date$ Utrecht University
///////////////////////////////////////////////////////////////////////////

#include "IceDOM.h"
#include "Riostream.h"
 
ClassImp(IceDOM) // Class implementation to enable ROOT I/O
 
IceDOM::IceDOM() : IceGOM()
{
// Default constructor.
}
///////////////////////////////////////////////////////////////////////////
IceDOM::~IceDOM()
{
// Default destructor.
}
///////////////////////////////////////////////////////////////////////////
IceDOM::IceDOM(const IceDOM& m) : IceGOM(m)
{
// Copy constructor.
}
///////////////////////////////////////////////////////////////////////////
TObject* IceDOM::Clone(const char* name) const
{
// Make a deep copy of the current object and provide the pointer to the copy.
// This memberfunction enables automatic creation of new objects of the
// correct type depending on the object type, a feature which may be very useful
// for containers like AliEvent when adding objects in case the
// container owns the objects. This feature allows e.g. AliEvent
// to store either IceDOM objects or objects derived from IceDOM
// via tha AddDevice memberfunction, provided these derived classes also have
// a proper Clone memberfunction. 

 IceDOM* m=new IceDOM(*this);
 if (name)
 {
  if (strlen(name)) m->SetName(name);
 }
 return m;
}
///////////////////////////////////////////////////////////////////////////
