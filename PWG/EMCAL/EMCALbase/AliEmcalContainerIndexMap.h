#if !(defined(__CINT__) || defined(__MAKECINT__))
#ifndef ALIEMCALCONTAINERINDEXMAP_H
#define ALIEMCALCONTAINERINDEXMAP_H

#include <map>
#include <typeinfo>

#include <TObject.h>
#include <TClonesArray.h>
#include <TObjArray.h>
#include <TClass.h>
#include <AliLog.h>

#include "AliEmcalContainer.h"

/**
 * @class AliEmcalContainerIndexMap
 * @ingroup EMCALCOREFW
 * @brief Container or array to index mapping for AliEmcalContainer derived classes
 *
 * AliEmcalContainerIndexMap contains index mapping for AliEmcalContainer derived classes. The mapping is from
 * a single global linear array to an index in a collection of a container. For example, if there are two
 * containers, then container a would be at offset 0 and container b at offset 1000. Then, the third particle
 * in container b would be indexed as 1002, while the third particle in container a would be indexed as 2.
 *
 * Below, input object (or inputObject) should be read as a TClonesArray or AliEmcalContainer derived class.
 * An object of interest will usually be what we are trying to retrieve, so usually a AliVCluster or AliVParticle
 * derived class, but could be any type derived from TObject.
 *
 * Generally, U corresponds to the input object and V corresponds to the object of interest.
 *
 * Actually using the index mapping requires a few steps. The mapping between the TClonesArrays used by cluster
 * or particle containers and indices is stored in each container class. Thus, the first step is to copy that
 * mapping, adapting it to map between the containers in your task and the indices. To do so, create an 
 * instance of AliEmcalContainerIndexMap for your class and then copy the mapping with ``CopyMappingFrom()``.
 * For example:
 *
 * ~~~{cxx}
 * // This example will focus on particle containers.
 * // The particle containers are assumed to be stored in fParticleCollArray, which is
 * // a TObjArray defined in AliAnalysisTaskEmcal (assuming your task inherits from it).
 * 
 * // Define the object to container the container to index mapping.
 * AliEmcalContainerIndexMap<AliParticleContainer, AliVParticle> particleMap;
 * // It is recommended to call CopyMappingFrom() in ExecOnce().
 * // At the very least, it must be after the arrays are set in the containers!
 * particleMap.CopyMappingFrom(AliParticleContainer::GetEmcalContainerIndexMap(), particleCollArray);
 * ~~~
 *
 * Now, we can access the objects as desired. Some example might include:
 *
 * ~~~{cxx}
 * // We must already have the global index from somewhere, such as the id of a
 * // track matched to a cluster.
 *
 * // Get the particle
 * AliVParticle * particleOne = particleMap.GetObjectFromGlobalIndex(globalIndex);
 * // Or get the index and the associated container
 * auto result = particleMap.LocalIndexFromGlobalIndex(globalIndex);
 * Int_t localIndex = result.first;
 * AliParticleContainer * partCont = result.second;
 * // Can use the result to check whether a particle is accepted.
 * AliVParticle * particleTwo = partCont.GetAcceptParticle(localIndex);
 * ~~~
 *
 * For some additional information, see the framework details in the
 * [embedding documentation](\ref READMEemcEmbedding).
 *
 * @author Raymond Ehlers <raymond.ehlers@cern.ch>, Yale University
 * @date Nov 21, 2016
 */

template <class U, class V>
class AliEmcalContainerIndexMap {
 public:
  AliEmcalContainerIndexMap();
  // Implemented using copy-and-swap mechanism
  AliEmcalContainerIndexMap(const AliEmcalContainerIndexMap & map);
  AliEmcalContainerIndexMap &operator=(AliEmcalContainerIndexMap map);
  // Change the template variables to avoid shadowing U and V.
  // See: https://stackoverflow.com/a/4039892
  template<class X, class Y>
  friend void swap(AliEmcalContainerIndexMap<X, Y> & first, AliEmcalContainerIndexMap<X, Y> & second);
  AliEmcalContainerIndexMap(AliEmcalContainerIndexMap<U, V> && other);

  // Don't delete fClass, as it was returned by TClass and deleting it will cause segfaults in ROOT
  virtual ~AliEmcalContainerIndexMap() {}

  // Setup index map for AliEmcalContainer or TClonesArray
  int RegisterArray(U * inputObject);

  // Copy the mapping from the source "map"; the mapping is between the global index found in "map" for object "cont" and the object "cont" itself
  template<class U2>
  void CopyMappingFrom(const AliEmcalContainerIndexMap<U2, V>& map, U* cont);

  // Copy the mapping from the source "map"; the mapping is between the global index found in "map" for each object in "containers" and the objects in "containers"
  template<class U2>
  void CopyMappingFrom(const AliEmcalContainerIndexMap<U2, V>& map, TCollection & containers);

  // Index operations
  // Returns the offset of an input object 
  template<class U2>
  int GetOffset(const U2 * inputObject) const;
  // Returns the global index of an input object and the local index of the object of interest
  int GlobalIndexFromLocalIndex(const U * inputObject, const int localIndex) const;
  // Returns a pair with the local index and the array for a given global index
  std::pair<int, U *> LocalIndexFromGlobalIndex(const int globalIndex) const;
  // Get object directly
  V * GetObjectFromGlobalIndex(const int globalIndex) const;

 protected:
  //template<class Z = U>
  const TClonesArray * GetObject(const AliEmcalContainer * inputObject) const;
  const TClonesArray * GetObject(const TClonesArray * inputObject) const;
  bool IsUnderlyingInputObjectTypeCompatible(const U * inputObject) const;

  std::map <int, U *> fGlobalIndexMap; //!<! Map between index and input object

  int fOffset;                         ///< Offset between each TClonesArray
  TClass* fClass;                      ///< Used to compare the type of V against the underlying input object type
};

/**
 * Standard constructor
 */
template<class U, class V>
AliEmcalContainerIndexMap<U, V>::AliEmcalContainerIndexMap():
  fGlobalIndexMap(),
  fOffset(100000),
  fClass(TClass::GetClass(typeid(V)))
{
}

/**
 * Copy Constructor
 */
template<class U, class V>
AliEmcalContainerIndexMap<U, V>::AliEmcalContainerIndexMap(const AliEmcalContainerIndexMap<U, V> & map):
  fGlobalIndexMap(map.fGlobalIndexMap),
  fOffset(map.fOffset),
  fClass(map.fClass) // Should be okay to be the same pointer, as it just points to the class type
{
}

/**
 * Move constructor
 */
template<class U, class V>
AliEmcalContainerIndexMap<U, V>::AliEmcalContainerIndexMap(AliEmcalContainerIndexMap<U, V> && other):
  AliEmcalContainerIndexMap()
{
  swap(*this, other);
}

/**
 * Assignment operator. Note that we pass by _value_, so a copy is created and it is
 * fine to swap the values with the created object!
 */
template<class U, class V>
AliEmcalContainerIndexMap<U, V> & AliEmcalContainerIndexMap<U, V>::operator=(AliEmcalContainerIndexMap<U, V> map)
{
  swap(*this, map);

  return *this;
}

/**
 * Swap function. Created using guide described here: https://stackoverflow.com/a/3279550
 */
template<class X, class Y>
void swap(AliEmcalContainerIndexMap<X, Y> & first, AliEmcalContainerIndexMap<X, Y> & second)
{
  using std::swap;

  swap(first.offset, second.offset);
  swap(first.fGlobalIndexMap, second.fGlobalIndexMap);
  swap(first.fClass, second.fClass);
}

/**
 * Get the proper object to use in comparison. This function handles AliEmcalContainer derived classes.
 *
 * Note: If we use enable_if<> in the future, we should check that it is an AliEmcalContainer derived class with std::is_conertible<>.
 *
 * @param inputObject AliEmcalContainer derived class.
 *
 * @return TClonesArray to be mapped.
 */
template<class U, class V>
inline const TClonesArray * AliEmcalContainerIndexMap<U, V>::GetObject(const AliEmcalContainer * inputObject) const
{
  return inputObject->GetArray();
}

/**
 * Get the proper object to use in comparison. This function handles a TClonesArray.
 * It is really just a dummy function, but is needed so that AliEmcalContainer derived classes work correctly.
 *
 * @param inputObject The TClonesArray to be mapped.
 *
 * @return TClonesArray to be mapped.
 */
template<class U, class V>
inline const TClonesArray * AliEmcalContainerIndexMap<U, V>::GetObject(const TClonesArray * inputObject) const
{
  return inputObject;
}

/**
 * Check if the type to cast the object of interest (ie V) is compatible with the underlying object
 * type of the input object. Note that the input object must implement GetClass() (which for AliEmcalContainer
 * derived classes, it passes through the GetClass() of the underlying TClonesArray) in order to check the type.
 *
 * If the types are incompatible (ie returns false), then an AliError is thrown!
 *
 * @param inputObject Input object to be checked.
 *
 * @return True if the underlying type is **compatible**.
 */
template<class U, class V>
inline bool AliEmcalContainerIndexMap<U, V>::IsUnderlyingInputObjectTypeCompatible(const U * inputObject) const
{
  if (!(inputObject->GetClass()->InheritsFrom(fClass))) {
    AliErrorGeneral("AliEmcalContainerIndexMap", Form("Cannot register array %s. This map can only accept arrays of type %s.", inputObject->GetName(), fClass->GetName()));
    return false;
  }

  return true;
}

/**
 * Add arrays to the map.
 * We are just comparing the address that the TClonesArray points to. Given that we should be in control of
 * everyting, this should always be a fine way to make the comparison, and much faster than comparing the
 * entire object.
 *
 * @param inputObject Input object to be registered in the map. Can be either AliEmcalContainer or TClonesArray
 *
 * @return Index where the input object is registered.
 */
template<class U, class V>
int AliEmcalContainerIndexMap<U, V>::RegisterArray(U * inputObject)
{
  if (IsUnderlyingInputObjectTypeCompatible(inputObject) == false) {
    return -1;
  }

  int index = 0;
  bool addToMap = true;
  for (auto val : fGlobalIndexMap)
  {
    // Check if the array is already added
    // We want to compare the address of the TClonesArrays to simplify the comparison
    if (GetObject(val.second) == GetObject(inputObject)) {
      addToMap = false;
      index = val.first;
      break;
    }

    // Find the max index so that we can add at the proper place
    if (val.first >= index) {
      index = val.first + fOffset;
    }
  }

  // Add to map
  if (addToMap) {
    // We want to compare the address of the TClonesArrays to simplify the comparison
    fGlobalIndexMap[index] = inputObject;
  }

  return index;
}

/**
 * Get the global index from the local index of an object of interest in the input object.
 * The input object is needed to determine the offset.
 *
 * @param inputObject Input object where the object of interest resides.
 * @param localIndex Local index of the object of interest in the input object.
 *
 * @return globalIndex The global index of the object of interest.
 */
template<class U, class V>
int AliEmcalContainerIndexMap<U, V>::GlobalIndexFromLocalIndex(const U * inputObject, const int localIndex) const
{
  int globalIndex = GetOffset(inputObject);

  // Only add the local index if we found a match!
  if (globalIndex >= 0) {
    globalIndex += localIndex;
  }
  else {
    AliWarningGeneral("AliEmcalContainerIndexMap", TString::Format("Unable to retrieve global index for input object %s. Was the input object registered?", inputObject->GetName()));
  }

  return globalIndex;
}

/**
 * Get the offset assigned to a particular input object.
 *
 * @param inputObject Input object which we want to get the offset of in the map.
 *
 * @return Global index corresponding to the offset of the input object in the map.
 */
template<class U, class V>
template<class U2>
int AliEmcalContainerIndexMap<U, V>::GetOffset(const U2 * inputObject) const
{
  int globalIndex = -1;
  for (auto val : fGlobalIndexMap)
  {
    // We want to compare the address of the TClonesArrays to simplify the comparison
    if (GetObject(val.second) == GetObject(inputObject)) {
      globalIndex = val.first;
      break;
    }
  }

  return globalIndex;
}

/**
 * Get the local index of an object of interest from the global index.
 *
 * @param globalIndex Global index of the object of interest
 *
 * @return std::pair of the local index and the associated input object (could be TClonesArray * or AliEmcalContainer *). The input object can be used to retrieve the object of interest.
 */
template<class U, class V>
std::pair<int, U *> AliEmcalContainerIndexMap<U, V>::LocalIndexFromGlobalIndex(const int globalIndex) const
{
  // Round to nearest offset
  int index = (globalIndex + fOffset/2)/fOffset * fOffset;

  // Used at() to allow the function to be const (operator[] is not const for std::map)
  U * array = fGlobalIndexMap.at(index);

  return std::pair<int, U*> (globalIndex - index, array);
}

/**
 * Gets the object of interest at the global index. Leverages LocalIndexFromGlobalIndex() to get the array
 * location of the object, as well as the input object.
 *
 * @param globalIndex
 *
 * @return The object of interest, cast to the proper type.
 */
template<class U, class V>
V * AliEmcalContainerIndexMap<U, V>::GetObjectFromGlobalIndex(const int globalIndex) const
{
  auto res = LocalIndexFromGlobalIndex(globalIndex);

  // We don't need to call GetObject() if we access using operator[]!
  return static_cast<V *>((*res.second)[res.first]);
}

/**
 * Copy a given index mapping from one AliEmcalContainerIndexMap object to the current one. The map values which
 * are copied are those corresponding to the input objects collection passed to this function. For a single
 * input object, see the other signature for this function.
 *
 * @param map AliEmcalContainerIndexMap object containing the index map that should be copied.
 * @param containers The collection of input objects whose index values should be copied from the index map.
 *
 * @return
 */
template<class U, class V>
template<class U2>
void AliEmcalContainerIndexMap<U, V>::CopyMappingFrom(const AliEmcalContainerIndexMap<U2, V>& map, TCollection & containers)
{
  TIter next(&containers);
  TObject* obj = 0;
  while((obj = next())) {
    U* cont = dynamic_cast<U*>(obj);
    if (!cont) continue;
    if (IsUnderlyingInputObjectTypeCompatible(cont) == false) {
      continue;
    }

    int gindex = map.GetOffset(cont);
    if (gindex >= 0) fGlobalIndexMap[gindex] = cont;
  }
}

/**
 * Copy a given index mapping from one AliEmcalContainerIndexMap object to the current one. The map values which
 * are copied are those corresponding to the input object passed to this function. For a collection of input
 * objects, see the other signature for this function.
 *
 * @param map AliEmcalContainerIndexMap object containing the index map that should be copied.
 * @param cont The input object whose index value should be copied from the index map.
 *
 * @return
 */
template<class U, class V>
template<class U2>
void AliEmcalContainerIndexMap<U, V>::CopyMappingFrom(const AliEmcalContainerIndexMap<U2, V>& map, U* cont)
{
  if (IsUnderlyingInputObjectTypeCompatible(cont) == false) {
    return;
  }

  int gindex = map.GetOffset(cont);
  if (gindex >= 0) fGlobalIndexMap[gindex] = cont;
}

#endif /* AliEmcalContainerIndexMap.h */
#endif /* Hiding from CINT */
