// $Header$

#include "NLTPointSetGL.h"

using namespace Reve;

//______________________________________________________________________
// NLTPointSetGL
//
// A hack around a bug in fglrx that makes rendering of projected pointsets
// terribly slow with display-lists on and rendering as crosses.

ClassImp(NLTPointSetGL)

NLTPointSetGL::NLTPointSetGL() : TPointSet3DGL()
{
  fDLCache = kFALSE; // Disable display list.
}

NLTPointSetGL::~NLTPointSetGL()
{}
