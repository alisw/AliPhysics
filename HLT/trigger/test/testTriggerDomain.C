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

/**
 * This macro is used to test the behaviour of the AliHLTTriggerDomain class.
 * We specifically check that the AliHLTTriggerDomain operators behave correctly
 * as set operations.
 */

#if defined(__CINT__) && (! defined(__MAKECINT__))
#error This macro must be compiled. Try running as testTriggerDomain.C++
#endif

#include "Riostream.h"
#include "TSystem.h"
#include "TClassTable.h"
#include "TRandom3.h"
#include "TObjArray.h"
#include "AliHLTTriggerDomain.h"
#include "AliHLTDomainEntry.h"

const int N = 6;
const char* datatype[N] = {"*******", "AAAAAAAA", "BBBBBBBB", "CCCCCCCC", "DDDDDDDD", "EEEEEEEE"};
const char* origin[N] = {"***", "xxxx", "yyyy", "zzzz", "wwww", "uuuu"};
int spec[N] = {0, 0x11111111, 0x22222222, 0x33333333, 0x44444444, 0x55555555};

/**
 * Randomly builds a trigger domain from the various data types, origins and specifications.
 * \param domain <i>[out]</i> The new constructed domain is filled into this output variable.
 * \param opcount <i>[in]</i> The number of Add / Remove operations to apply to the domain.
 */
void BuildTriggerDomain(AliHLTTriggerDomain& domain, int opcount = 100)
{
	domain.Clear();
	for (int i = 0; i < opcount; i++)
	{
		// Make sure not to use the last data type, origin or specification in
		// the respective lists because we want to use them later to test the
		// matching against wild card values.
		int itype = gRandom->Integer(N-1);
		int iorigin = gRandom->Integer(N-1);
		int ispec = gRandom->Integer(N-1);
		if (gRandom->Integer(2) == 0)
		{
			if (spec[ispec] == 0)
				domain.Add(datatype[itype], origin[iorigin]);
			else
				domain.Add(datatype[itype], origin[iorigin], spec[ispec]);
		}
		else
		{
			if (spec[ispec] == 0)
				domain.Remove(datatype[itype], origin[iorigin]);
			else
				domain.Remove(datatype[itype], origin[iorigin], spec[ispec]);
		}
	}
}

/**
 * Returns the string representation of the data type.
 */
const char* TypeToString(const AliHLTComponentDataType& type)
{
	static char str[kAliHLTComponentDataTypefIDsize+1];
	for (int i = 0; i < kAliHLTComponentDataTypefIDsize; i++)
	{
		str[i] = type.fID[i];
	}
	return str;
}

/**
 * Returns the string representation of the data origin.
 */
const char* OriginToString(const AliHLTComponentDataType& type)
{
	static char str[kAliHLTComponentDataTypefOriginSize+1];
	for (int i = 0; i < kAliHLTComponentDataTypefOriginSize; i++)
	{
		str[i] = type.fOrigin[i];
	}
	return str;
}

/**
 * Checks to see if the domain is correctly constructed compared to the list of entries
 * that were applied to the domain.
 * \param entryList <i>[out]</i> The list of entries that were used to construct the domain.
 * \param domain <i>[in]</i> The trigger domain to check.
 */
bool DomainOK(const TObjArray& entryList, const AliHLTTriggerDomain& domain)
{
	for (int i = 1; i < N; i++)
	for (int j = 1; j < N; j++)
	for (int k = 1; k < N; k++)
	{
		AliHLTDomainEntry entry(datatype[i], origin[j], spec[k]);
		AliHLTComponentBlockData block;
		block.fDataType = entry;
		block.fSpecification = spec[k];
		
		bool result = false;
		bool result2 = false;
		for (int n = 0; n < entryList.GetEntriesFast(); n++)
		{
			const AliHLTDomainEntry* rule = static_cast<const AliHLTDomainEntry*>( entryList[n] );
			if (*rule == entry) result = rule->Inclusive();
			if (*rule == &block) result2 = rule->Inclusive();
		}
		
		bool containsResult = domain.Contains(entry);
		bool includeResult = domain.IncludeInReadout(&block);
		if (containsResult != result or includeResult != result2)
		{
			cout << "FAILED DomainOK test!" << endl;
			cout << "==============================================================================" << endl;
			cout << "Dump of how AliHLTTriggerDomain was built:" << endl;
			cout << "  AliHLTTriggerDomain domain;" << endl;
			for (int n = 0; n < entryList.GetEntriesFast(); n++)
			{
				const AliHLTDomainEntry* rule = static_cast<const AliHLTDomainEntry*>( entryList[n] );
				const char* opstr = (rule->Exclusive() ? "Remove" : "Add");
				if (rule->IsValidSpecification())
				{
					cout << "  domain." << opstr << "(\"" << TypeToString(rule->DataType())
						<< "\", \"" << OriginToString(rule->DataType())
						<< "\", 0x" << hex << rule->Specification() << dec
						<< ");" << endl;
				}
				else
				{
					cout << "  domain." << opstr << "(\"" << TypeToString(rule->DataType())
						<< "\", \"" << OriginToString(rule->DataType())
						<< "\");" << endl;
				}
			}
			cout << "==============================================================================" << endl;
			cout << "Dump of the trigger domain contents:" << endl;
			domain.Print();
			cout << "==============================================================================" << endl;
		}
		if (containsResult != result)
		{
			cout << "Failed for entry = ";
			entry.Print();
			cout << "Result of domain.Contains(entry) == " << containsResult << endl;
			cout << " Expected domain.Contains(entry) == " << result << endl;
			return false;
		}
		if (includeResult != result2)
		{
			cout << "Failed for block: type = " << TypeToString(block.fDataType)
				<< ", origin = " << OriginToString(block.fDataType)
				<< ", specification = " << hex << block.fSpecification << dec;
			cout << "Result of domain.IncludeInReadout(&block) == " << includeResult << endl;
			cout << " Expected domain.IncludeInReadout(&block) == " << result2 << endl;
			return false;
		}
	}
	
	return true;
}

/**
 * Randomly builds a trigger domain and tests the Add / Remove operations of the
 * AliHLTTriggerDomain class.
 * \param opcount <i>[in]</i> The number of Add / Remove operations to use to
 *     build the domain.
 */
bool AddRemoveOK(int opcount = 100)
{
	TObjArray entryList;
	entryList.SetOwner(kTRUE);
	AliHLTTriggerDomain domain;
	
	for (int i = 0; i < opcount; i++)
	{
		// Make sure not to use the last data type, origin or specification in
		// the respective lists because we want to use them later to test the
		// matching against wild card values.
		int itype = gRandom->Integer(N-1);
		int iorigin = gRandom->Integer(N-1);
		int ispec = gRandom->Integer(N-1);
		if (gRandom->Integer(2) == 0)
		{
			if (spec[ispec] == 0)
			{
				entryList.Add(new AliHLTDomainEntry(kFALSE, datatype[itype], origin[iorigin]));
				domain.Add(datatype[itype], origin[iorigin]);
			}
			else
			{
				entryList.Add(new AliHLTDomainEntry(kFALSE, datatype[itype], origin[iorigin], spec[ispec]));
				domain.Add(datatype[itype], origin[iorigin], spec[ispec]);
			}
		}
		else
		{
			if (spec[ispec] == 0)
			{
				entryList.Add(new AliHLTDomainEntry(kTRUE, datatype[itype], origin[iorigin]));
				domain.Remove(datatype[itype], origin[iorigin]);
			}
			else
			{
				entryList.Add(new AliHLTDomainEntry(kTRUE, datatype[itype], origin[iorigin], spec[ispec]));
				domain.Remove(datatype[itype], origin[iorigin], spec[ispec]);
			}
		}
	}
	
	if (not DomainOK(entryList, domain)) return false;
	
	return true;
}


#define DefOperation(name, expr, logicExpr) \
	class name \
	{ \
	public: \
		static const char* Name() \
		{ \
			return #name; \
		} \
		static const char* Expression() \
		{ \
			return #expr; \
		} \
		static AliHLTTriggerDomain Apply(const AliHLTTriggerDomain& a, const AliHLTTriggerDomain& b) \
		{ \
			return expr; \
		} \
		static bool ExpectedResult(bool a, const bool b) \
		{ \
			return logicExpr; \
		} \
	};


DefOperation( OpComplement, ~a,    !a       );
DefOperation( OpUnion,      a | b, a | b    );
DefOperation( OpIntersect,  a & b, a & b    );
DefOperation( OpXor,        a ^ b, a ^ b    );
DefOperation( OpPlus,       a + b, a | b    );
DefOperation( OpMinus,      a - b, a & (!b) );

/**
 * Randomly builds two trigger domains and tests a overloaded operator of the
 * AliHLTTriggerDomain class.
 * \param opcount <i>[in]</i> The number of Add / Remove operations to use to
 *     build the domains.
 */
template <class Op>
bool OperatorOK(int opcount = 100)
{
	AliHLTTriggerDomain d1;
	BuildTriggerDomain(d1, opcount);
	AliHLTTriggerDomain d2;
	BuildTriggerDomain(d2, opcount);
	AliHLTTriggerDomain d3 = Op::Apply(d1, d2);
	
	for (int i = 1; i < N; i++)
	for (int j = 1; j < N; j++)
	for (int k = 1; k < N; k++)
	{
		AliHLTDomainEntry entry(datatype[i], origin[j], spec[k]);
		bool d1Result = d1.Contains(entry);
		bool d2Result = d2.Contains(entry);
		bool d3Result = d3.Contains(entry);
		bool d3ExpectedResult = Op::ExpectedResult(d1Result, d2Result);
		if (d3Result != d3ExpectedResult)
		{
			cout << "FAILED OperatorOK<" << Op::Name() << "> test!" << endl;
			cout << "==============================================================================" << endl;
			cout << "Dump of the trigger domain a contents:" << endl;
			d1.Print();
			cout << "==============================================================================" << endl;
			cout << "Dump of the trigger domain b contents:" << endl;
			d2.Print();
			cout << "==============================================================================" << endl;
			cout << "Dump of the trigger domain c = " << Op::Expression() << " contents:" << endl;
			d3.Print();
			cout << "==============================================================================" << endl;
			cout << "Failed for entry = ";
			entry.Print();
			cout << "Result of c.Contains(entry) == " << d3Result << endl;
			cout << " Expected c.Contains(entry) == " << d3ExpectedResult << endl;
			return false;
		}
	}
	
	return true;
}


#define DefEqualExprCheck(name, expr1, expr2) \
	class name \
	{ \
	public: \
		static const char* Name() \
		{ \
			return #name; \
		} \
		static const char* Expr1() \
		{ \
			return #expr1; \
		} \
		static const char* Expr2() \
		{ \
			return #expr2; \
		} \
		static AliHLTTriggerDomain Apply1(const AliHLTTriggerDomain& a, const AliHLTTriggerDomain& b) \
		{ \
			return expr1; \
		} \
		static AliHLTTriggerDomain Apply2(const AliHLTTriggerDomain& a, const AliHLTTriggerDomain& b) \
		{ \
			return expr2; \
		} \
	};


DefEqualExprCheck( XorExprCheck1,  a ^ b, (a | b) - (a & b)             );
DefEqualExprCheck( XorExprCheck2,  a ^ b, (a - (a & b)) | (b - (a & b)) );
DefEqualExprCheck( MinusExprCheck1, a - b, a & (a ^ b)                  );
DefEqualExprCheck( MinusExprCheck2, a - b, a & ~(a & b)                 );

/**
 * Randomly builds two trigger domains and tests two expressions applied to the
 * domains that should be equivalent.
 * \param opcount <i>[in]</i> The number of Add / Remove operations to use to
 *     build the domains.
 */
template <class Expr>
bool EquivalentExpressionsOK(int opcount = 100)
{
	AliHLTTriggerDomain d1;
	BuildTriggerDomain(d1, opcount);
	AliHLTTriggerDomain d2;
	BuildTriggerDomain(d2, opcount);
	AliHLTTriggerDomain d3 = Expr::Apply1(d1, d2);
	AliHLTTriggerDomain d4 = Expr::Apply2(d1, d2);
	
	for (int i = 1; i < N; i++)
	for (int j = 1; j < N; j++)
	for (int k = 1; k < N; k++)
	{
		AliHLTDomainEntry entry(datatype[i], origin[j], spec[k]);
		bool d3Result = d3.Contains(entry);
		bool d4Result = d4.Contains(entry);
		if (d3Result != d4Result)
		{
			cout << "FAILED EquivalentExpressionsOK<" << Expr::Name() << "> test!" << endl;
			cout << "==============================================================================" << endl;
			cout << "Dump of the trigger domain a contents:" << endl;
			d1.Print();
			cout << "==============================================================================" << endl;
			cout << "Dump of the trigger domain b contents:" << endl;
			d2.Print();
			cout << "==============================================================================" << endl;
			cout << "Dump of the trigger domain c = " << Expr::Expr1() << " contents:" << endl;
			d3.Print();
			cout << "==============================================================================" << endl;
			cout << "Dump of the trigger domain d = " << Expr::Expr2() << " contents:" << endl;
			d3.Print();
			cout << "==============================================================================" << endl;
			cout << "Failed for entry = ";
			entry.Print();
			cout << "Result of c.Contains(entry) == " << d3Result << endl;
			cout << "Result of d.Contains(entry) == " << d4Result << endl;
			cout << "But the results should be the same." << endl;
			return false;
		}
	}
	
	return true;
}


#define DefTest(name, routine, description) \
	class name \
	{ \
	public: \
		static bool Run(int opcount) \
		{ \
			return routine(opcount); \
		} \
		static const char* Description() \
		{ \
			return description; \
		} \
	}

// Declarations of the different tests.
DefTest(AddRemoveTest,           AddRemoveOK,              "AliHLTTriggerDomain::Add and AliHLTTriggerDomain::Remove");
DefTest(ComplementOperationTest, OperatorOK<OpComplement>, "operator ~a");
DefTest(UnionOperationTest,      OperatorOK<OpUnion>,      "operator a | b");
DefTest(IntersectOperationTest,  OperatorOK<OpIntersect>,  "operator a & b");
DefTest(XorOperationTest,        OperatorOK<OpXor>,        "operator a ^ b");
DefTest(PlusOperationTest,       OperatorOK<OpPlus>,       "operator a + b");
DefTest(MinusOperationTest,      OperatorOK<OpMinus>,      "operator a - b");
DefTest(XorExpressionTest1,   EquivalentExpressionsOK<XorExprCheck1>,   "expression a ^ b == (a | b) - (a & b)");
DefTest(XorExpressionTest2,   EquivalentExpressionsOK<XorExprCheck2>,   "expression a ^ b == (a - (a & b)) | (b - (a & b))");
DefTest(MinusExpressionTest1, EquivalentExpressionsOK<MinusExprCheck1>, "expression a - b == a & (a ^ b)");
DefTest(MinusExpressionTest2, EquivalentExpressionsOK<MinusExprCheck2>, "expression a - b == a & ~(a & b)");

/**
 * Routine to run an individual test.
 * \param print  Should information be printed showing the progress of testing.
 * \param numOfTests  The number of test iterations to run.
 */
template <class Test>
bool Run(bool print = true, int numOfTests = 100)
{
	if (print) cout << "Testing " << Test::Description() << " ..." << endl;
	for (int i = 0; i < numOfTests; i++)
	{
		if (not Test::Run(10)) return false;
		if (not Test::Run(100)) return false;
		if (not Test::Run(1000)) return false;
		if (print and ((i+1) % (numOfTests / 10) == 0))
		{
			cout << "  Completed " << ((i+1) * 100 / numOfTests) << "%" << endl;
		}
	}
	return true;
}

/**
 * This is the top level testing method which calls individual tests.
 * \param print  Should information be printed showing the progress of testing.
 * \param numOfTests  The number of test iterations to run.
 * \param seed  The random number seed to use.
 */
bool testTriggerDomain(bool print = true, int numOfTests = 100, int seed = 0)
{
	if (gClassTable->GetID("AliHLTDomainEntry") < 0)
	{
		gSystem->Load("libAliHLTTrigger.so");
	}
	
	gRandom->SetSeed(seed);
	
	if (not Run<AddRemoveTest>(print, numOfTests)) return false;
	if (not Run<ComplementOperationTest>(print, numOfTests)) return false;
	if (not Run<UnionOperationTest>(print, numOfTests)) return false;
	if (not Run<IntersectOperationTest>(print, numOfTests)) return false;
	if (not Run<XorOperationTest>(print, numOfTests)) return false;
	if (not Run<PlusOperationTest>(print, numOfTests)) return false;
	if (not Run<MinusOperationTest>(print, numOfTests)) return false;
	if (not Run<XorExpressionTest1>(print, numOfTests)) return false;
	if (not Run<XorExpressionTest2>(print, numOfTests)) return false;
	if (not Run<MinusExpressionTest1>(print, numOfTests)) return false;
	if (not Run<MinusExpressionTest2>(print, numOfTests)) return false;
	
	return true;
}

#ifndef __MAKECINT__

int main(int /*argc*/, const char** /*argv*/)
{
	bool resultOk = testTriggerDomain();
	if (not resultOk) return 1;
	return 0;
}

#endif // __MAKECINT__

