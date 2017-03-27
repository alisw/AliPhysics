#ifndef ALIEMCALCONTAINERUTILS_H
#define ALIEMCALCONTAINERUTILS_H

#include <string>

#include <TObjArray.h>

class AliVEvent;

/**
 * @namespace AliEmcalContainerUtils
 * @ingroup EMCALCOREFW
 * @brief Helper functions related to AliEmcalContainer derived objects
 *
 * This namespace includes a variety of helper functions related to AliEmcalContainer derived objects,
 * as well as configuring input objects.
 *
 * Functions include:
 *  - Automatically determine the proper input object branch name based on the input object and file
 *  (AOD or ESD) types. This is determined by the "usedefault" pattern (see the code for precise
 *  implementation). It can also return the proper type input object to be stored in a collection
 *  such as a TClonesArray.
 *  - Centralizing the code to add and get containers from the collections of containers often
 *  stored in analysis tasks. 
 *
 * Usage of the add, and get containers code should be something like (for cluster containers):
 * ~~~{.cxx}
 * AliClusterContainer * AddClusterContainer(const char * name)    { AddContainer<AliClusterContainer>(name, fClusterCollArray); }
 * AliClusterContainer * GetClusterContainer(Int_t i)              { GetContainer<AliClusterContainer>(i, fClusterCollArray); }
 * AliClusterContainer * GetClusterContainer(const char * name)    { GetContainer<AliClusterContainer>(name, fClusterCollArray); }
 * ~~~
 * The code is similar for other containers. See AliEmcalCorrectionTask for an example.
 *
 * @author Raymond Ehlers <raymond.ehlers@cern.ch>, Yale University
 * @date Dec 6, 2016
 */

class AliEmcalContainerUtils {
 public:
  /** 
   * @enum InputObject_t
   * @brief %Type of input object to be created
   */
  enum InputObject_t {
    kNoDefinedInputObject = -1,    //!<! Not initialized type
    kCaloCells = 0,                //!<! Calo cells
    kCluster = 1,                  //!<! Cluster container
    kTrack = 2,                    //!<! Track container
  };

  // Utility functions
  static std::string DetermineUseDefaultName(InputObject_t objType, bool esdMode, bool returnObjectType = false);
  static AliVEvent * GetEvent(AliVEvent * inputEvent, bool isEmbedding = false);

#if !(defined(__CINT__) || defined(__MAKECINT__))
  template <class T>
  static inline T * AddContainer(const char *n, TObjArray & collection);
  template <class T>
  static T * GetContainer(Int_t i, const TObjArray & collection);
  template <class T>
  static T * GetContainer(const char *name, const TObjArray & collection);

  // Then to use these functions, add below (for example) to each class (or you can use it directly if desired)
  /*AliMCParticleContainer * AddMCParticleContainer(const char * name) { AddContainer<AliMCPartileContainer>(name, fParticleCollArray); }
  AliTrackContainer * AddTrackContainer(const char * name)        { AddContainer<AliTrackContainer>(name, fParticleCollArray); }
  AliParticleContainer * AddTrackContainer(const char * name)     { AddContainer<AliParticleContainer>(name, fParticleCollArray); }
  */
#endif /* Hide from CINT */
};

#if !(defined(__CINT__) || defined(__MAKECINT__))
/**
 * Create new container and attach it to the passed collection. The name
 * provided to this function must match the name of the array attached
 * to the new container inside the input event.
 *
 * @param[in] n Name of the container and the array the container points to
 * @param[in] collection Collection into which the new container will be added
 * @return Pointer to the new container
 */
template <class T>
inline T * AliEmcalContainerUtils::AddContainer(const char *n, TObjArray & collection)
{
  if (TString(n).IsNull()) return 0;

  T * cont = new T(n);

  collection.Add(cont);

  return cont;
}

/**
 * Get \f$ i^{th} \f$ container contained in the collection
 * @param[in] i Index of the container
 * @param[in] collection Collection which stores available containers
 * @return Container found for the given index (NULL if no container exists for that index)
 */
template <class T>
inline T * AliEmcalContainerUtils::GetContainer(Int_t i, const TObjArray & collection)
{
  if (i<0 || i>collection.GetEntriesFast()) return 0;
  T *cont = static_cast<T *>(collection.At(i));
  return cont;
}

/**
 * Find container container in the collection according to its name
 * @param[in] name Name of the container
 * @param[in] collection Collection which stores available containers
 * @return Container found under the given name
 */
template <class T>
inline T * AliEmcalContainerUtils::GetContainer(const char *name, const TObjArray & collection)
{
  T *cont = static_cast<T *>(collection.FindObject(name));
  return cont;
}
#endif /* Hide from CINT */

#endif /* AliEmcalContainerUtils.h */
