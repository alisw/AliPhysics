///
/// \file AliFemtoTypes.h
/// \author Frank Laue <laue@bnl.gov>, BNL
///
/// This file simply loades other header files with 'AliFemto' datatypes:
///
///   * AliFemtoString
///   * AliFemtoVector
///   * AliFemtoHelix
///   * AliFemtoEnumeration
///

/*
 * $Log$
 * Revision 1.1.1.1  2007/04/25 15:38:41  panos
 * Importing the HBT code dir
 *
 * Revision 1.1.1.1  2007/03/07 10:14:49  mchojnacki
 * First version on CVS
 *
 * Revision 1.14  2001/06/21 19:15:48  laue
 * Modified files:
 *   CTH.h : new constructor added
 *   AliFemtoEvent, AliFemtoKink, AliFemtoTrack : constructors from the persistent
 *                                   (TTree) classes added
 *   AliFemtoLikeSignAnalysis : minor changes, for debugging
 *   AliFemtoTypes: split into different files
 * Added files: for the new TTree muDst's
 *   StExceptions.cxx StExceptions.h AliFemtoEnumeration.h
 *   AliFemtoHelix.h AliFemtoHisto.h AliFemtoString.h AliFemtoTFile.h
 *   AliFemtoTTreeEvent.cxx AliFemtoTTreeEvent.h AliFemtoTTreeKink.cxx
 *   AliFemtoTTreeKink.h AliFemtoTTreeTrack.cxx AliFemtoTTreeTrack.h
 *   AliFemtoTTreeV0.cxx AliFemtoTTreeV0.h AliFemtoVector.h
 *
 *
 */

//
// I split this up into different files, so that I do not have to
// load/recompile everything all over again.
//

//#include "AliFemtoHisto.h"
#include "AliFemtoString.h"
#include "AliFemtoVector.h"
#include "AliFemtoHelix.h"
#include "AliFemtoEnumeration.h"
