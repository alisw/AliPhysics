#ifndef EMCALHLTGUI_SHARED_PTR
#define EMCALHLTGUI_SHARED_PTR
/****************************************************************************************
 *  shared pointer implementation for ALICE                                             *
 *  Copyright (C) 2016 The ALICE collaboration                                          *
 *                                                                                      *
 *  This program is free software: you can redistribute it and/or modify                *
 *  it under the terms of the GNU General Public License as published by                *
 *  the Free Software Foundation, either version 3 of the License, or                   *
 *  (at your option) any later version.                                                 *
 *                                                                                      *
 *  This program is distributed in the hope that it will be useful,                     *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of                      *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                       *
 *  GNU General Public License for more details.                                        *
 *                                                                                      *
 *  You should have received a copy of the GNU General Public License                   *
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.               *
 ****************************************************************************************/
// Ported from EMCALHLTGui - a simple monitoring program for EMCAL data from the HLT
// see https://gitlab.cern.ch/mfasel/EMCALHLTgui for details

#include <cstdlib>

/**
 * @class shared_ptr_base
 * @brief Data handler class for the shared ptr
 * @author Markus Fasel <markus.fasel@cern.ch>, Oak Ridge National Laboratory
 * @since Oct 11, 2016
 *
 * This class is the backbone of the shared ptr as it handles the
 * data, destruction and the reference count. All shared ptrs contain
 * a pointer to the same shared_ptr_base instance. Once the copy contructor
 * or assingment operator is invoked, the shared_pointer Connects to the
 * shared_ptr_base, meaning its reference count is increased. Once the
 * shared_ptr is destroyed, the shared_ptr_base is released meaning the
 * reference count is increased. Once the reference count reaches 0, meaning
 * non shared_ptr is pointing to the data, the object is destroyed. The destructor
 * can never be called from outside, but destruction should happen via the
 * release function.
 */
template<typename T>
class AliSharedPtrBase{
	T			*fData;				///< Data content of the shared ptr
	int			fRefCount;			///< Number of pointer to the data package

	/**
 	 * Destructor private, public methods are Connect and Release
 	 */
	~AliSharedPtrBase() { if(fData) delete fData; }

	// Copy constructor and assignment operator forbidden (would not make
	// sense as it would create new independent ref counts)
	AliSharedPtrBase(const AliSharedPtrBase<T> &ref);
	AliSharedPtrBase &operator=(const AliSharedPtrBase<T> &ref);
public:
	/**
 	 * Constructor - must always be called with the
 	 * data. The constructor
 	 */
	AliSharedPtrBase(T* data):
		fData(data),
		fRefCount(1)
	{
	}

	/**
 	 * Connect new shared pointer to the data.
 	 * Connecting a shared pointer means increasing the ref count
 	 */
	void Connect(){
		fRefCount++;
	}

	/**
 	 * Releasing a pointer from the data. In this case the
 	 * ref count is removed. A ref count of 0 means no pointer
 	 * is pointing to the data, consequently the data package
 	 * can be deleted.
 	 */
	void Release(){
		fRefCount--;
		if(!fRefCount) delete this;
	}

	/**
 	 * Return pointer to the underlying data of the shared_ptr
 	 * @return Pointer to the underlying data
 	 */
	T* GetData() const {
		return fData;
	}
};

/**
 * @class AliSharedPtr
 * @brief shared pointer, a smart pointer containing a reference count
 * @author Markus Fasel <markus.fasel@cern.ch>, Oak Ridge National Laboratory
 * @since Oct 11, 2016
 *
 * shared_pointers are integral part of automatic memory management via
 * reference counting. Shared pointers can be shared among several places.
 * This is implemented by copying the shared pointer, which automatically
 * increases its reference count. Once the shared pointer gets out-of-scope,
 * the reference count is decreased. As soon as the reference count reaches
 * 0 the data handled by the shared pointer is removed.
 */
template<typename T>
class AliSharedPtr {
	AliSharedPtrBase<T>			*fSharedData;				///< underlying shared_ptr data and ref count
public:
	/**
 	 * Default constructor. This is needed as it is called by
 	 * array definitions. However this doesn't fully initialize
 	 * the shared pointer for which it is not yet fully operational.
 	 */
	AliSharedPtr():
		fSharedData(new AliSharedPtrBase<T>(NULL))
	{
	}

	/**
 	 * Main constuctor, called with the data
 	 * for the first time. At this step an underlying
 	 * data handler is created, which sets the reference
 	 * count to 1 (first pointer is pointing to the data)
 	 *
 	 * @param[in] data Input data to which the shared pointers point to.
 	 */
	AliSharedPtr(T * data):
		fSharedData(new AliSharedPtrBase<T>(data))
	{
	}

	/**
 	 * Copy constructor. At this step the shared pointer
 	 * really becomes "shared", as it is copied from another
 	 * existing shared pointer. Connecting to the shared data.
 	 *
 	 * @param[in] ref Pointer from which this pointer is initialized.
 	 */
	AliSharedPtr(const AliSharedPtr<T> &ref):
		fSharedData(ref.fSharedData)
	{
		fSharedData->Connect();
	}

	/**
 	 * Assigment operator. At this step the shared pointer
 	 * really becomes "shared", as it is copied from another
 	 * existing shared pointer. Connecting to the shared data.
 	 *
 	 * The asignment operator has to handle a special case
 	 * which is that the existing shared pointe could be connected
 	 * to a different set of data. In this case it has to release
 	 * the connection first. Afterward it can connect to the new
 	 * data. As self-assingment is forbidden it is guaranteed that
 	 * the assignment operator never produces a 0 ref count.
 	 *
 	 * @param[in] ref Pointer from which this pointer is initialized.
 	 */
 	AliSharedPtr<T> &operator=(const AliSharedPtr<T> &ref){
		if(this != &ref){
			if(fSharedData) fSharedData->Release();
			fSharedData = ref.fSharedData;
			fSharedData->Connect();
		}
		return *this;
	}

	/**
 	 * Destructor. The destructor does not delete the
 	 * underlying object but only reduces the reference
 	 * count as there might be pointers still pointing
 	 * to the data. Once the reference count becomes 0
 	 * the data is deleted.
 	 */
	~AliSharedPtr(){ fSharedData->Release(); }

	/**
 	 * Access to data in form of a pointer.
 	 * @return Pointer to underlying data
 	 */
	T * get() const { return fSharedData->GetData(); }

	/**
 	 * Redefine pointer. Initialise with different address
 	 * In this case the memory it was pointing to previously
 	 * is released.
 	 *
 	 * @param[in] data new data to be used
 	 */
	void set(T *data){
		fSharedData->Release();
		fSharedData = new AliSharedPtrBase<T>(data);
	}

	/**
 	 * Dereferencing the shared pointer. Attention: Data
 	 * can be modified.
 	 * @return Reference to underlying data.
 	 */
	T& operator*() const { return *(fSharedData->GetData()); }

	/**
 	 * pointer operator. Mimics actual pointer
 	 * functionality
 	 *
 	 * @return Pointer to the underlying data
 	 */
	T* operator->() const { return fSharedData->GetData(); }

	/**
 	 * operator bool, checks if the shared pointer
 	 * is defined (meaning the memory it points to
 	 * is not NULL.
 	 *
 	 * @return True if not null, false otherwise
 	 */
	operator bool() const {
		return get() != NULL;
	}
};

/***************************
 * operators
 ***************************/

template<typename T>
bool operator==(const AliSharedPtr<T> &lhs, const AliSharedPtr<T> &rhs){
	return lhs.get() == rhs.get();
}

template<typename T>
bool operator!=(const AliSharedPtr<T> &lhs, const AliSharedPtr<T> &rhs){
	return !(lhs == rhs);
}

template<typename T>
bool operator<(const AliSharedPtr<T> &lhs, const AliSharedPtr<T> &rhs){
	return lhs.get() < rhs.get();
}

template<typename T>
bool operator<=(const AliSharedPtr<T> &lhs, const AliSharedPtr<T> &rhs){
	return lhs.get() <= rhs.get();
}

template<typename T>
bool operator>(const AliSharedPtr<T> &lhs, const AliSharedPtr<T> &rhs){
	return lhs.get() > rhs.get();
}

template<typename T>
bool operator>=(const AliSharedPtr<T> &lhs, const AliSharedPtr<T> &rhs){
	return lhs.get() >= rhs.get();
}

#endif
