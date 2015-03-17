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

#include "AliMUONAttPainter.h"

#include <Riostream.h>

/// \class AliMUONAttPainter
///
/// Basic attributes shared by all AliMUONVPainter objects
///
/// One basic attribute is "what side" of the detector part we are representing.
/// Take the example of one full chamber. We can represent it as seen from the
/// interaction point, i.e. showing all the "cathode0" of all detection elements.
///
/// Or we might want to see only the bending (or non bending) planes of all the
/// detection elements.
///
/// This is governed by the SetCathode() and SetPlane() methods (and retrieved
/// using IsCathodeDefined() and IsPlaneDefined() methods. In the first case
/// above IsCathodeDefined() would be true and IsPlaneDefined() would be false.
/// The second case would be the contrary.
///
/// \author Laurent Aphecetche, Subatech

using std::cout;
using std::endl;
///\cond CLASSIMP
ClassImp(AliMUONAttPainter)
///\endcond

//_____________________________________________________________________________
AliMUONAttPainter::AliMUONAttPainter()
: TObject(), fName()
{
	/// ctor
	SetSingle(kTRUE);
	SetValid(kTRUE);
	SetCathodeAndPlaneMutuallyExclusive(kFALSE);
	SetCathodeAndPlaneDisabled(kFALSE);
	SetName();
}

//_____________________________________________________________________________
AliMUONAttPainter::~AliMUONAttPainter()
{
	/// dtor
}

//_____________________________________________________________________________
void AliMUONAttPainter::SetName()
{
	/// Build our name

	fName = "Invalid";

	if ( !IsValid() ) return;

	fName = "";

	if ( CathodeName().Length() > 0 ) fName = CathodeName();
	if ( PlaneName().Length() > 0 )
	{
		if ( fName.Length() > 0 ) fName += "-";
		fName += PlaneName();
	}

	//  if ( ViewPointName().Length() > 0 )
	//  {
	//    if ( name.Length() > 0 ) name += "-";
	//    name += ViewPointName();
	//  }
}

//_____________________________________________________________________________
TString
AliMUONAttPainter::CathodeName() const
{
	/// Return cathode name in short form

	if ( IsCathode0() && IsCathode1() ) return "Both";
	else if ( !IsCathode0() && !IsCathode1() ) return "";
	else if ( IsCathode0() ) return "0";
	else if ( IsCathode1() ) return "1";
	return "";
}

//_____________________________________________________________________________
void
AliMUONAttPainter::Invert()
{
	/// Invert our cathode/plane states

	if ( IsCathodeDefined() )
	{
		Bool_t cath0(IsCathode0());
		Bool_t cath1(IsCathode1());
		SetCathode(!cath0,!cath1);
	}

	if ( IsPlaneDefined() )
	{
		Bool_t b(IsBendingPlane());
		Bool_t nb(IsNonBendingPlane());
		SetPlane(!b,!nb);
	}

	SetName();
}

//_____________________________________________________________________________
TString
AliMUONAttPainter::PlaneName() const
{
	/// Return plane name in short form
	if ( IsBendingPlane() && IsNonBendingPlane() ) return "Both";
	else if ( !IsBendingPlane() && !IsNonBendingPlane() ) return "";
	else if ( IsBendingPlane() ) return "Bending";
	else if ( IsNonBendingPlane() ) return "NonBending";
	return "";
}

//_____________________________________________________________________________
TString
AliMUONAttPainter::ViewPointName() const
{
	/// Return name of view point
	if ( IsFrontView() ) return "Front";
	if ( IsBackView() ) return "Back";
	return "";
}

//_____________________________________________________________________________
void
AliMUONAttPainter::Print(Option_t*) const
{
	/// Printout

	if ( !IsValid() ) cout << "INVALID : ";

	if ( IsCathodeDefined() )
	{
		cout << "Cathode-defined " << CathodeName() << ". ";
	}
	if ( IsPlaneDefined() )
	{
		cout << "Plane-defined " << PlaneName() << ". ";
	}
	if ( IsCathodeAndPlaneMutuallyExclusive() )
	{
		cout << "Cathode and Plane mutually exclusive. ";
	}
	cout << ViewPointName() << endl;
}

//_____________________________________________________________________________
void
AliMUONAttPainter::SetCathode(Bool_t cath0, Bool_t cath1)
{
	/// Set our cathode states
	SetBit(kIsCathode0,cath0);
	SetBit(kIsCathode1,cath1);
	SetName();
}

//_____________________________________________________________________________
void
AliMUONAttPainter::SetPlane(Bool_t bending, Bool_t nonBending)
{
	/// Set our plane states
	SetBit(kIsBendingPlane,bending);
	SetBit(kIsNonBendingPlane,nonBending);
	SetName();
}

//_____________________________________________________________________________
void
AliMUONAttPainter::SetSingle(Bool_t value)
{
	/// Set single status
	SetBit(kIsSinglePainter,value);
	SetName();
}

//_____________________________________________________________________________
void
AliMUONAttPainter::SetViewPoint(Bool_t front, Bool_t back)
{
	/// Set view point
	SetBit(kIsFrontView,front);
	SetBit(kIsBackView,back);
	SetName();
}

//_____________________________________________________________________________
void
AliMUONAttPainter::SetCathodeAndPlaneMutuallyExclusive(Bool_t value)
{
	/// Set mutually exclusive flag
	SetBit(kIsCathodeAndPlaneMutuallyExclusive,value);
	SetName();
}

//_____________________________________________________________________________
void
AliMUONAttPainter::SetValid(Bool_t value)
{
	/// Set valid flag
	SetBit(kIsValid,value);
	SetName();
}

//_____________________________________________________________________________
void
AliMUONAttPainter::SetCathodeAndPlaneDisabled(Bool_t value)
{
	/// Set cathode & plane disable flag
	SetBit(kIsCathodeAndPlaneDisabled,value);
	SetName();
}
