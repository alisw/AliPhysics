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

/// @file   testAliHLTTriggerCounters.C
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   28 Oct 2010
/// @brief  Test program for the AliHLTTriggerCounters class.
///

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliHLTTriggerCounters.h"
#include "TObjArray.h"
#include "TString.h"
#include "Riostream.h"
#endif

/**
 * Tests basic functionality of the AliHLTTriggerCounters::AliCounter class.
 */
bool CheckCounterItemClass()
{
	AliHLTTriggerCounters::AliCounter s1("s1", "counter one", 1, 0.1);
	AliHLTTriggerCounters::AliCounter s2("s2", "counter two", 2, 0.2);
	AliHLTTriggerCounters::AliCounter s3("s2", "counter two", 3, 0.2);
	if (TString(s1.GetName()) != s1.Name())
	{
		cerr << "ERROR: AliHLTTriggerCounters::AliCounter::GetName() returns a different value than AliHLTTriggerCounters::AliCounter::Name()." << endl;
		return false;
	}
	if (TString(s1.GetTitle()) != s1.Description())
	{
		cerr << "ERROR: AliHLTTriggerCounters::AliCounter::GetTitle() returns a different value than AliHLTTriggerCounters::AliCounter::Description()." << endl;
		return false;
	}
	if (s2 == s3)
	{
		cerr << "ERROR: equals operator for AliHLTTriggerCounters::AliCounter returns the wrong value." << endl;
		return false;
	}
	if (! s2.IsEqual(&s3))
	{
		cerr << "ERROR: AliHLTTriggerCounters::AliCounter::IsEqual returns the wrong value." << endl;
		return false;
	}
	s2.Increment();
	if (! (s2 == s3))
	{
		cerr << "ERROR: equals operator for AliHLTTriggerCounters::AliCounter returns the wrong value." << endl;
		return false;
	}
	TObjArray list;
	list.Add(&s2);
	list.Add(&s1);
	list.Sort();
	if (TString(list.At(0)->GetName()) != "s1")
	{
		cerr << "ERROR: Sorting objects of type AliHLTTriggerCounters::AliCounter is not working correctly." << endl;
		return false;
	}
	return true;
}

/**
 * Tests functionality of the AliHLTTriggerCounters class.
 */
bool CheckCountersListClass()
{
	AliHLTTriggerCounters s;
	s.Add("a", "one", 1);
	s.Add("b", "two", 2, 5);
	if (s.NumberOfScalars() != 2)
	{
		cerr << "ERROR: The number of added counters is wrong for class AliHLTTriggerCounters." << endl;
		return false;
	}
	if (! s.Exists("a"))
	{
		cerr << "ERROR: AliHLTTriggerCounters claims counter 'a' does not exist event though it was added." << endl;
		return false;
	}
	if (! s.Exists("b"))
	{
		cerr << "ERROR: AliHLTTriggerCounters claims counter 'b' does not exist event though it was added." << endl;
		return false;
	}
	s.Remove("a");
	if (s.Exists("a"))
	{
		cerr << "ERROR: AliHLTTriggerCounters claims counter 'a' does not exist event though it was removed." << endl;
		return false;
	}
	s.Add("a", "one", 1);
	const AliHLTTriggerCounters& p = s;
	if (p.GetCounter("a").Rate() != 1)
	{
		cerr << "ERROR: Constant version of AliHLTTriggerCounters::GetCounter(\"a\") returns the wrong counter object." << endl;
		return false;
	}
	if (TString(p.GetCounter("c").Name()) != "" || TString(p.GetCounter("c").Description()) != "" || p.GetCounter("c").Rate() != 0 || p.GetCounter("c").Counter() != 0)
	{
		cerr << "ERROR: Constant version of AliHLTTriggerCounters::GetCounter(\"c\") does not return a sentinel object." << endl;
		return false;
	}
	if (s.GetCounter("a").Rate() != 1)
	{
		cerr << "ERROR: AliHLTTriggerCounters::GetCounter(\"a\") returns the wrong counter object." << endl;
		return false;
	}
	s.GetCounter("c").Value(3);
	if (TString(s.GetCounter("c").Name()) != "c" || TString(s.GetCounter("c").Description()) != "")
	{
		cerr << "ERROR: AliHLTTriggerCounters::GetCounter(\"c\") does not create a new object." << endl;
		return false;
	}
	s.Add("c", "three", 33, 7);
	if (s.GetCounter("c").Rate() != 33 || s.GetCounter("c").Counter() != 7 || TString(s.GetCounter("c").Description()) != "")
	{
		cerr << "ERROR: AliHLTTriggerCounters::Add did not update an exisiting counter correctly." << endl;
		return false;
	}
	if (TString(p.GetCounterN(0).Name()) != "b" || TString(p.GetCounterN(1).Name()) != "a" || TString(p.GetCounterN(2).Name()) != "c")
	{
		cerr << "ERROR: Constant version of AliHLTTriggerCounters::GetCounterN(0) returns the wrong counter object." << endl;
		return false;
	}
	if (TString(s.GetCounterN(0).Name()) != "b" || TString(s.GetCounterN(1).Name()) != "a" || TString(s.GetCounterN(2).Name()) != "c")
	{
		cerr << "ERROR: AliHLTTriggerCounters::GetCounterN(0) returns the wrong counter object." << endl;
		return false;
	}
	if (TString(p.GetCounterN(4).Name()) != "" || TString(p.GetCounterN(4).Description()) != "" || p.GetCounterN(4).Rate() != 0)
	{
		cerr << "ERROR: Constant version of AliHLTTriggerCounters::GetCounterN(4) returns the wrong counter object." << endl;
		return false;
	}
	s.GetCounterN(4).Value(5);
	if (TString(s.GetCounterN(4).Name()) != "Scalar4" || TString(s.GetCounterN(4).Description()) != "" || s.GetCounterN(4).Rate() != 5)
	{
		cerr << "ERROR: AliHLTTriggerCounters::GetCounterN(4) does not create a new counter object correctly." << endl;
		return false;
	}
	if (TString(p.GetCounterN(3).Name()) != "Scalar3" || TString(p.GetCounterN(3).Description()) != "" || p.GetCounterN(3).Rate() != 0)
	{
		cerr << "ERROR: AliHLTTriggerCounters::GetCounterN(4) did not initialise the third counter as expected." << endl;
		return false;
	}
	
	// The following is a special check to check for compilation ambiguity
	// rather than runtime behaviour.
	if (s[4].Rate() != s["Scalar4"].Rate())
	{
		cerr << "ERROR: AliHLTTriggerCounters::operator[](UInt_t) did not return the same value as AliHLTTriggerCounters::operator[](const char*)." << endl;
		return false;
	}
	
	// Here we check to see that the AliHLTTriggerCounters::GetCounterN class correctly
	// checks and finds an unused name.
	s.Add("Scalar7", "six", 6, 2);
	s.Add("Scalar7_0", "seven", 7, 3);
	s.GetCounterN(7).Value(8);
	if (! s.Exists("Scalar7_1") || s.GetCounterN(7).Rate() != 8 || s.GetCounter("Scalar7_1").Rate() != 8)
	{
		cerr << "ERROR: AliHLTTriggerCounters::GetCounterN is not creating a counter object with a unique name as expected." << endl;
		return false;
	}
	
	// Check the copying of the object.
	AliHLTTriggerCounters* c1 = (AliHLTTriggerCounters*) s.Clone();
	AliHLTTriggerCounters c2;
	c2 = s;
	AliHLTTriggerCounters c3;
	s.Copy(c3);
	if (! (*c1 == s) || *c1 != s)
	{
		cerr << "ERROR: The equals operator of AliHLTTriggerCounters is not working as expected." << endl;
		return false;
	}
	if (c2 != s)
	{
		cerr << "ERROR: The assignment operator of AliHLTTriggerCounters is not working as expected." << endl;
		return false;
	}
	if (c3 != s)
	{
		cerr << "ERROR: The method AliHLTTriggerCounters::Copy is not working as expected." << endl;
		return false;
	}
	c1->UpdateTimeStamp();
	if (*c1 == s)
	{
		cerr << "ERROR: Modification of the time stamp for AliHLTTriggerCounters did not work as expected or comparison operator is not working." << endl;
		return false;
	}
	delete c1;
	
	// Now check the IsEqual and Reset methods:
	if (! c2.IsEqual(&c3))
	{
		cerr << "ERROR: The AliHLTTriggerCounters::IsEqual method is not working as expected." << endl;
		return false;
	}
	
	c3.Reset();
	for (UInt_t i = 0; i < c3.NumberOfScalars(); ++i)
	{
		if (c3[i].Rate() != 0 || c3[i].Counter() != 0)
		{
			cerr << "ERROR: AliHLTTriggerCounters::Reset did not reset all counter values to zero." << endl;
			return false;
		}
		if (TString(c3[i].Name()) != c2[i].Name())
		{
			cerr << "ERROR: AliHLTTriggerCounters::Reset modified the name by mistake." << endl;
			return false;
		}
		if (TString(c3[i].Description()) != c2[i].Description())
		{
			cerr << "ERROR: AliHLTTriggerCounters::Reset modified the description by mistake." << endl;
			return false;
		}
	}
	if (! c2.IsEqual(&c3))
	{
		cerr << "ERROR: The AliHLTTriggerCounters::IsEqual method is not working as expected after call to Reset." << endl;
		return false;
	}
	if (c2 == c3)
	{
		cerr << "ERROR: The equals operator for AliHLTTriggerCounters is not working as expected after call to Reset." << endl;
		return false;
	}
	
	c2.Remove("c");
	if (c2.IsEqual(&c3))
	{
		cerr << "ERROR: The AliHLTTriggerCounters::IsEqual method is not working as expected after call to Remove." << endl;
		return false;
	}
	if (c2 == c3)
	{
		cerr << "ERROR: The equals operator for AliHLTTriggerCounters is not working as expected after call to Remove." << endl;
		return false;
	}
	
	return true;
}

/**
 * Runs the unit test for the AliHLTTriggerCounters class.
 * \returns true if the class passed the test and false otherwise.
 */
bool testAliHLTTriggerCounters()
{
	if (! CheckCounterItemClass()) return false;
	if (! CheckCountersListClass()) return false;
	return true;
}

#ifndef __MAKECINT__

int main(int /*argc*/, const char** /*argv*/)
{
	bool resultOk = testAliHLTTriggerCounters();
	if (not resultOk) return 1;
	return 0;
}

#endif // __MAKECINT__
