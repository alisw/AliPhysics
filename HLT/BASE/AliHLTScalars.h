//-*- Mode: C++ -*-
// $Id$
#ifndef AliHLTSCALARS_H
#define AliHLTSCALARS_H
/* This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

///  @file   AliHLTScalars.h
///  @author Artur Szostak <artursz@iafrica.com>
///  @date   28 Sep 2010
///  @brief  Declares the a base class for named scalar values.

#include "TObject.h"
#include "TNamed.h"
#include "TTimeStamp.h"
#include "TClonesArray.h"
#include "THashTable.h"

/**
 * @class AliHLTScalars
 * @brief Container for named scalar values.
 *
 * This class contains a list of named scalars for an event as summary information.
 * These can be used by the trigger components to perform event selection or used
 * for monitoring purposes.
 *
 * \ingroup alihlt_base
 */
class AliHLTScalars : public TObject
{
public:
	/**
	 * This class stores a single scalar value and name.
	 */
	class AliScalar : public TNamed
	{
	public:
		/// Default constructor
		AliScalar() : TNamed(), fValue(0) {}
		
		/// Constructor to set the initial value.
		AliScalar(const char* name, const char* description, Double_t value) :
			TNamed(name, description), fValue(value)
		{}
		
		/// Default destructor
		virtual ~AliScalar() {}
		
		/// Inherited from TObject. Compares two scalar names.
		virtual Int_t Compare(const TObject *obj) const
		{
			return fName.CompareTo(obj->GetName());
		}
		
		/// Inherited from TObject. Returns true.
		virtual Bool_t IsSortable() const { return kTRUE; }
		
		/**
		 * Inherited from TObject.
		 * Returns true if the names of the scalars are the same.
		 */
		virtual Bool_t IsEqual(const TObject *obj) const
		{
			return fName == obj->GetName();
		}
		
		/// Resets the scalar value to zero.
		virtual void Clear(Option_t* /*option*/ = "") { fValue = 0; }
	
		/// Inherited from TObject. Performs a deep copy.
		virtual void Copy(TObject& object) const;
		
		/// Returns the value of the scalar.
		Double_t Value() const { return fValue; }

		/// Sets a new value for the scalar.
		void Value(Double_t value) { fValue = value; }

		/**
		 * Increments the scalar by a value of 'count'.
		 * \param count  The number to increment the scalar by. The default is 1.
		 */
		void Increment(UInt_t count = 1) { fValue += count; }
		
		/// Returns the name of the scalar.
		const char* Name() const { return fName.Data(); }
		
		/// Returns the description string for the scalar.
		const char* Description() const { return fTitle.Data(); }
		
		/// Checks if two scalar objects are identical.
		bool operator == (const AliScalar& x) const
		{
			return fValue == x.fValue and fName == x.fName and fTitle == x.fTitle;
		}
		
		/// Checks if two scalar objects are not identical.
		bool operator != (const AliScalar& x) const
		{
			return not (this->operator == (x));
		}
		
		/// Typecast operator for returning the value directly.
		operator Double_t () { return fValue; }
		
	private:
		Double_t fValue; // The scalar's value.
		
		ClassDef(AliScalar, 1);  // HLT scalar value.
	};
	
	/// Default constructor.
	AliHLTScalars();
	
	/// The copy constructor performs a deep copy.
	AliHLTScalars(const AliHLTScalars& obj);
	
	/// Default destructor.
	virtual ~AliHLTScalars();
	
	/**
	 * Adds a new scalar to the end of the scalars list.
	 * If the scalar already exists then its values are updated instead.
	 * \param name  The name of the scalar.
	 * \param description  A short description of the scalar.
	 * \param value  The value of the new scalar.
	 * \returns true if the scalar already exists and false otherwise.
	 */
	virtual bool Add(const char* name, const char* description = NULL, Double_t value = 0);

	/**
	 * Removes a named scalar from the scalars list.
	 * \param name  The name of the scalar to remove.
	 * \returns true if the scalar existed and false otherwise.
	 * \note The scalars list is compressed so this method will be slow.
	 *    In addition, scalar positions will change if not removing from the end.
	 */
	virtual bool Remove(const char* name);
	
	/// Checks to see if the named scalar exists.
	bool Exists(const char* name) const { return fMap.FindObject(name) != NULL; }

	/**
	 * Fetches the specified scalar object.
	 * \param name  The name of the scalar object.
	 * \returns the found scalar object, otherwise an empty sentinel object with
	 *    zeros. One can tell it is a sentinel because the name will be empty.
	 */
	const AliScalar& GetScalar(const char* name) const;

	/**
	 * Fetches the specified scalar object for editing.
	 * \param name  The name of the scalar object.
	 * \returns the found scalar object. If the scalar does not already
	 *     exist then a new one is created and returned.
	 */
	AliScalar& GetScalar(const char* name);

	/// Returns the number of scalar values.
	UInt_t NumberOfScalars() const { return UInt_t(fScalars.GetEntriesFast()); }

	// Note: the following GetScalarN methods do not use the same name as
	// GetScalar above because the parameter type would unfortunately be
	// ambiguous to an ISO c++ compiler.
	
	/**
	 * Fetches the n'th scalar object.
	 * \param n  The number of the scalar object.
	 * \returns the found scalar object, otherwise an empty sentinel object with
	 *    zeros. One can tell it is a sentinel because the name will be empty.
	 */
	const AliScalar& GetScalarN(UInt_t n) const;

	/**
	 * Fetches the n'th scalar object for editing.
	 * \param n  The number of the scalar object.
	 * \returns the found scalar object. If the scalar does not already
	 *     exist then a new one is created and returned.
	 */
	AliScalar& GetScalarN(UInt_t n);
	
	/// Resets all scalar values to zero.
	virtual void Reset();
	
	/**
	 * Removes all the scalars from the internal array.
	 * \param option  This is passed onto the internal Delete method.
	 */
	virtual void Clear(Option_t* option = "");
	
	/// Inherited form TObject. Performs a deep copy.
	virtual void Copy(TObject& object) const;
	
	/// Finds the scalar object by name.
	virtual TObject* FindObject(const char* name) const
	{
		return fMap.FindObject(name);
	}
	
	/// Finds the scalar object with the same name as obj->GetName().
	virtual TObject* FindObject(const TObject* obj) const
	{
		return fMap.FindObject(obj->GetName());
	}
	
	/**
	 * Inherited from TObject, this prints the contents of all the scalars.
	 * \param option  Can be "compact", which will just print all the values on one line.
	 */
	virtual void Print(Option_t* option = "") const;
	
	/**
	 * The assignment operator performs a deep copy.
	 */
	AliHLTScalars& operator = (const AliHLTScalars& obj);
	
	/// Returns the n'th scalar or a zero sentinel if n is out of range.
	const AliScalar& operator [] (UInt_t n) const { return GetScalarN(n); }

	/// Returns the n'th scalar for editing. A new scalar is created if n is out of range.
	AliScalar& operator [] (UInt_t n) { return GetScalarN(n); }

	/// Returns the named scalar or a zero sentinel if no such scalar is found.
	const AliScalar& operator [] (const TString& name) const { return GetScalar(name.Data()); }

	/// Returns the named scalar for editing. A new scalar is created if n is out of range.
	AliScalar& operator [] (const TString& name) { return GetScalar(name.Data()); }

	/**
	 * Inherited from TObject.
	 * Returns true if the names of the two sets of scalars are the same.
	 * \note The actual values are not checked. Use the comparison operator for that.
	 */
	virtual Bool_t IsEqual(const TObject *obj) const;
	
	/**
	 * Comparison operator to check if two sets of scalars have the same values.
	 * \note The description strings are not checked so they could be different
	 *   and the order of the scalars does not matter either.
	 */
	bool operator == (const AliHLTScalars& obj) const;
	
	/**
	 * Comparison operator to check if two sets of scalars are different.
	 * \note The description strings are not checked, only the values are.
	 *   In addition, the order of the scalars does not matter.
	 */
	bool operator != (const AliHLTScalars& obj) const
	{
		return not (this->operator == (obj));
	}

protected:
	
	/**
	 * Constructor that can be used by deriving classes to overload the class stored
	 * in the fScalars TClonesArray.
	 * \param cl  The class to use in the fScalars as passed to the TClonesArray constructor.
	 * \param initSize  The initial approximate number of elements in fScalars. (Default = 128).
	 * \note The class used in <i>cl</i> must derive from AliHLTScalars::AliScalar.
	 */
	AliHLTScalars(const TClass* cl, Int_t initSize = 128);
	
	/**
	 * This method creates a new scalar object in the fScalars TClonesArray.
	 * \param i  Location of the new object to construct in the TClonesArray.
	 * \param name  The name of the new scalar.
	 * \param description  The description of the new scalar.
	 * \param value  The value of the new scalar.
	 * \returns the pointer to the new object created.
	 * \note This method must be overridden by classes inheriting from this class if
	 *    the protected AliHLTScalars(const TClass*, Int_t) constructor is used to
	 *    change the class stored in the fScalars TClonesArray.
	 *    One should use the method ScalarForConstructor to get the location where
	 *    the new scalar object will be constructed.
	 */
	virtual AliScalar* NewScalar(UInt_t i, const char* name, const char* description, Double_t value);
	
	/**
	 * Returns a pointer to the memory where a new scalar object should be constructed.
	 * \param i  The position of the new object.
	 */
	TObject*& ScalarForConstructor(UInt_t i) { return fScalars[Int_t(i)]; }
	
	/**
	 * This method should return an empty sentinel object to mark that a scalar was
	 * not found in the list.
	 * \note This method must be overridden by classes inheriting from this class if
	 *    the protected AliHLTScalars(const TClass*, Int_t) constructor is used to
	 *    change the class stored in the fScalars TClonesArray.
	 */
	virtual const AliScalar& Sentinel() const;
	
	/**
	 * This is an internal Add method which can be faster to use than the public Add method
	 * directly for classes derived from AliHLTScalars.
	 * [out] \param scalar This gets filled with the pointer of the new scalar created or
	 *    the existing one found.
	 * [in] \param name The name of the scalar.
	 * [in] \param description  A short description of the scalar.
	 * [in] \param value  The value of the new scalar.
	 * \returns true if the scalar already exists and false otherwise.
	 */
	bool Add(AliScalar*& scalar, const char* name, const char* description, Double_t value);
	
	/**
	 * Utility method for classes deriving from AliHLTScalars to fetch the i'th scalar
	 * from the TClonesArray without checking that the index is valid.
	 * \param i  The index number of the scalar to fetch.
	 */
	AliScalar* ScalarUncheckedAt(UInt_t i) const { return static_cast<AliScalar*>(fScalars.UncheckedAt(Int_t(i))); }
	
private:
	
	TClonesArray fScalars;  // List of scalar objects.
	THashTable fMap;        //! Hash table of pointers to the scalars for fast lookup.
	
	ClassDef(AliHLTScalars, 1);  // Set of HLT scalars.
};

#endif // AliHLTSCALARS_H
