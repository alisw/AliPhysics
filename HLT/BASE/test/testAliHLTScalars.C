// $Id: $

/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Artur Szostak <artursz@iafrica.com>                   *
 *                  for The ALICE HLT Project.                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/// @file   testAliHLTScalars.C
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   27 Oct 2010
/// @brief  Test program for the AliHLTScalars class.
///

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliHLTScalars.h"
#include "TObjArray.h"
#include "TString.h"
#include "Riostream.h"
#endif

/**
 * Tests basic functionality of the AliHLTScalars::AliScalar class.
 */
bool CheckScalarItemClass()
{
	AliHLTScalars::AliScalar s1("s1", "scalar one", 1);
	AliHLTScalars::AliScalar s2("s2", "scalar two", 2);
	AliHLTScalars::AliScalar s3("s2", "scalar two", 3);
	if (TString(s1.GetName()) != s1.Name())
	{
		cerr << "ERROR: AliHLTScalars::AliScalar::GetName() returns a different value than AliHLTScalars::AliScalar::Name()." << endl;
		return false;
	}
	if (TString(s1.GetTitle()) != s1.Description())
	{
		cerr << "ERROR: AliHLTScalars::AliScalar::GetTitle() returns a different value than AliHLTScalars::AliScalar::Description()." << endl;
		return false;
	}
	if (s2 == s3)
	{
		cerr << "ERROR: equals operator for AliHLTScalars::AliScalar returns the wrong value." << endl;
		return false;
	}
	if (! s2.IsEqual(&s3))
	{
		cerr << "ERROR: AliHLTScalars::AliScalar::IsEqual returns the wrong value." << endl;
		return false;
	}
	s2.Increment();
	if (! (s2 == s3))
	{
		cerr << "ERROR: equals operator for AliHLTScalars::AliScalar returns the wrong value." << endl;
		return false;
	}
	TObjArray list;
	list.Add(&s2);
	list.Add(&s1);
	list.Sort();
	if (TString(list.At(0)->GetName()) != "s1")
	{
		cerr << "ERROR: Sorting objects of type AliHLTScalars::AliScalar is not working correctly." << endl;
		return false;
	}
	return true;
}

/**
 * Tests functionality of the AliHLTScalars class.
 */
bool CheckScalarsListClass()
{
	AliHLTScalars s;
	s.Add("a", "one", 1);
	s.Add("b", "two", 2);
	if (s.NumberOfScalars() != 2)
	{
		cerr << "ERROR: The number of added scalars is wrong for class AliHLTScalars." << endl;
		return false;
	}
	if (! s.Exists("a"))
	{
		cerr << "ERROR: AliHLTScalars claims scalar 'a' does not exist event though it was added." << endl;
		return false;
	}
	if (! s.Exists("b"))
	{
		cerr << "ERROR: AliHLTScalars claims scalar 'b' does not exist event though it was added." << endl;
		return false;
	}
	s.Remove("a");
	if (s.Exists("a"))
	{
		cerr << "ERROR: AliHLTScalars claims scalar 'a' does not exist event though it was removed." << endl;
		return false;
	}
	s.Add("a", "one", 1);
	const AliHLTScalars& p = s;
	if (p.GetScalar("a").Value() != 1)
	{
		cerr << "ERROR: Constant version of AliHLTScalars::GetScalar(\"a\") returns the wrong scalar object." << endl;
		return false;
	}
	if (TString(p.GetScalar("c").Name()) != "" || TString(p.GetScalar("c").Description()) != "" || p.GetScalar("c").Value() != 0)
	{
		cerr << "ERROR: Constant version of AliHLTScalars::GetScalar(\"c\") does not return a sentinel object." << endl;
		return false;
	}
	if (s.GetScalar("a").Value() != 1)
	{
		cerr << "ERROR: AliHLTScalars::GetScalar(\"a\") returns the wrong scalar object." << endl;
		return false;
	}
	s.GetScalar("c").Value(3);
	if (TString(s.GetScalar("c").Name()) != "c" || TString(s.GetScalar("c").Description()) != "")
	{
		cerr << "ERROR: AliHLTScalars::GetScalar(\"c\") does not create a new object." << endl;
		return false;
	}
	s.Add("c", "three", 33);
	if (s.GetScalar("c").Value() != 33 || TString(s.GetScalar("c").Description()) != "")
	{
		cerr << "ERROR: AliHLTScalars::Add did not update an exisiting scalar correctly." << endl;
		return false;
	}
	if (TString(p.GetScalarN(0).Name()) != "b" || TString(p.GetScalarN(1).Name()) != "a" || TString(p.GetScalarN(2).Name()) != "c")
	{
		cerr << "ERROR: Constant version of AliHLTScalars::GetScalarN(0) returns the wrong scalar object." << endl;
		return false;
	}
	if (TString(s.GetScalarN(0).Name()) != "b" || TString(s.GetScalarN(1).Name()) != "a" || TString(s.GetScalarN(2).Name()) != "c")
	{
		cerr << "ERROR: AliHLTScalars::GetScalarN(0) returns the wrong scalar object." << endl;
		return false;
	}
	if (TString(p.GetScalarN(4).Name()) != "" || TString(p.GetScalarN(4).Description()) != "" || p.GetScalarN(4).Value() != 0)
	{
		cerr << "ERROR: Constant version of AliHLTScalars::GetScalarN(4) returns the wrong scalar object." << endl;
		return false;
	}
	s.GetScalarN(4).Value(5);
	if (TString(s.GetScalarN(4).Name()) != "Scalar4" || TString(s.GetScalarN(4).Description()) != "" || s.GetScalarN(4).Value() != 5)
	{
		cerr << "ERROR: AliHLTScalars::GetScalarN(4) does not create a new scalar object correctly." << endl;
		return false;
	}
	if (TString(p.GetScalarN(3).Name()) != "Scalar3" || TString(p.GetScalarN(3).Description()) != "" || p.GetScalarN(3).Value() != 0)
	{
		cerr << "ERROR: AliHLTScalars::GetScalarN(4) did not initialise the third scalar as expected." << endl;
		return false;
	}
	
	// The following is a special check to check for compilation ambiguity
	// rather than runtime behaviour.
	if (s[4].Value() != s["Scalar4"].Value())
	{
		cerr << "ERROR: AliHLTScalars::operator[](UInt_t) did not return the same value as AliHLTScalars::operator[](const char*)." << endl;
		return false;
	}
	
	// Here we check to see that the AliHLTScalars::GetScalarN class correctly
	// checks and finds an unused name.
	s.Add("Scalar7", "six", 6);
	s.Add("Scalar7_0", "seven", 7);
	s.GetScalarN(7).Value(8);
	if (! s.Exists("Scalar7_1") || s.GetScalarN(7).Value() != 8 || s.GetScalar("Scalar7_1").Value() != 8)
	{
		cerr << "ERROR: AliHLTScalars::GetScalarN is not creating a scalar object with a unique name as expected." << endl;
		return false;
	}
	
	// Check the copying of the object.
	AliHLTScalars* c1 = (AliHLTScalars*) s.Clone();
	AliHLTScalars c2;
	c2 = s;
	AliHLTScalars c3;
	s.Copy(c3);
	if (! (*c1 == s) || *c1 != s)
	{
		cerr << "ERROR: The equals operator of AliHLTScalars is not working as expected." << endl;
		return false;
	}
	if (c2 != s)
	{
		cerr << "ERROR: The assignment operator of AliHLTScalars is not working as expected." << endl;
		return false;
	}
	if (c3 != s)
	{
		cerr << "ERROR: The method AliHLTScalars::Copy is not working as expected." << endl;
		return false;
	}
	delete c1;
	
	// Now check the IsEqual and Reset methods:
	if (! c2.IsEqual(&c3))
	{
		cerr << "ERROR: The AliHLTScalars::IsEqual method is not working as expected." << endl;
		return false;
	}
	
	c3.Reset();
	for (UInt_t i = 0; i < c3.NumberOfScalars(); ++i)
	{
		if (c3[i].Value() != 0)
		{
			cerr << "ERROR: AliHLTScalars::Reset did not reset all scalar values to zero." << endl;
			return false;
		}
		if (TString(c3[i].Name()) != c2[i].Name())
		{
			cerr << "ERROR: AliHLTScalars::Reset modified the name by mistake." << endl;
			return false;
		}
		if (TString(c3[i].Description()) != c2[i].Description())
		{
			cerr << "ERROR: AliHLTScalars::Reset modified the description by mistake." << endl;
			return false;
		}
	}
	if (! c2.IsEqual(&c3))
	{
		cerr << "ERROR: The AliHLTScalars::IsEqual method is not working as expected after call to Reset." << endl;
		return false;
	}
	if (c2 == c3)
	{
		cerr << "ERROR: The equals operator for AliHLTScalars is not working as expected after call to Reset." << endl;
		return false;
	}
	
	c2.Remove("c");
	if (c2.IsEqual(&c3))
	{
		cerr << "ERROR: The AliHLTScalars::IsEqual method is not working as expected after call to Remove." << endl;
		return false;
	}
	if (c2 == c3)
	{
		cerr << "ERROR: The equals operator for AliHLTScalars is not working as expected after call to Remove." << endl;
		return false;
	}
	
	return true;
}

/**
 * Runs the unit test for the AliHLTScalars class.
 * \returns true if the class passed the test and false otherwise.
 */
bool testAliHLTScalars()
{
	if (! CheckScalarItemClass()) return false;
	if (! CheckScalarsListClass()) return false;
	return true;
}

#ifndef __MAKECINT__

int main(int /*argc*/, const char** /*argv*/)
{
	bool resultOk = testAliHLTScalars();
	if (not resultOk) return 1;
	return 0;
}

#endif // __MAKECINT__
