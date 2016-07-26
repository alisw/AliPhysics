#ifndef ALIEMCALITERABLECONTAINER_H
#define ALIEMCALITERABLECONTAINER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <iterator>
#include <TArrayI.h>
#include <TObject.h>

class AliEmcalContainer;

/**
 * @class AliEmcalIterableContainer
 * @brief Container implementing iterable functionality of the EMCAL containers
 * @ingroup EMCALCOREFW
 * @author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 * @author Salvatore Aiola, Yale University
 * @date March 23rd, 2016
 *
 * Providing an interface to iterator functionality for the AliEmcalContainer and
 * inheriting objects, iterating over either all or only accepted objects inside
 * the container. The content is specified in the constructor.
 *
 * EMCAL iterable containers should not be created by hand. Instead, the EMCAL
 * container provides the functionality to create the interface for both cases:
 *
 * ~~~{.cxx}
 * AliEmcalContainer *cont;
 * AliEmcalIterableContainer *accepted = cont->accepted(), // iterative container over accepted entries
 *                           *all = cont->all();           // iterative container over all entries
 * ~~~
 *
 * Once created, EMCAL iterable containers implement the functions begin(), end(),
 * rbegin() and rend() creating stl iterators (type AliEmcalIterableContainer::iterator).
 * These can be used as normal stl iterators
 *
 * ~~~{.cxx}
 * for(AliEmcalIterableContainer::iterator iter = all.begin(); iter != all.end(); ++iter){
 *   // Do something with the object
 * }
 * ~~~
 *
 * In case c++11 is used this code simplifies to
 *
 * ~~~{.cxx}
 * for(auto en : all){
 *   // Do something with the object
 * }
 * ~~~
 */
class AliEmcalIterableContainer : public TObject {
public:
  /**
   * @class iterator
   * @brief bidirectional stl iterator over the EMCAL iterable container
   * @author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
   * @date March 23rd, 2016
   *
   * stl iterator corresponding to the EMCAL iterable container. The iterator
   * iterates over all object in the EMCAL iterable container as specified in
   * its constructor (all or accepted). It can be both forward or backward iterator.
   *
   * As stl iterator it implements the operators required for an iterator
   * - operator!= to determine the end of the iteration
   * - prefix and postfix operator++ and operator-- forwarding the position
   * - operator* to access the content of the iteration
   *
   * In case of c++11 the iterator also allows range-based iteration.
   */
  class iterator : public std::iterator<std::bidirectional_iterator_tag,
                                        TObject,std::ptrdiff_t,
                                        TObject **, TObject *>{
  public:
    iterator(const AliEmcalIterableContainer *cont, int currentpos, bool forward = true);
    iterator(const iterator &ref);
    iterator &operator=(const iterator &ref);
    virtual ~iterator(){}

    bool operator!=(const iterator &ref) const;

    iterator &operator++();
    iterator operator++(int);
    iterator &operator--();
    iterator operator--(int);

    TObject *operator*() const;

  private:
    iterator();

    const AliEmcalIterableContainer         *fkData;    ///< container with data
    int                                      fCurrent;  ///< current index in the container
    bool                                     fForward;  ///< use forward or backward direction
  };

  AliEmcalIterableContainer();
  AliEmcalIterableContainer(const AliEmcalContainer *cont, bool useAccept);
  AliEmcalIterableContainer(const AliEmcalIterableContainer &ref);
  AliEmcalIterableContainer &operator=(const AliEmcalIterableContainer &cont);

  /**
   * Destructor
   */
  virtual ~AliEmcalIterableContainer() {}

  TObject *operator[](int index) const;

  /**
   * Integer conversion operator: Returning the size if the container (number of entries)
   * @return Number of entries in the container
   */
  operator int() const { return GetEntries(); }

  /**
   * Access to underlying EMCAL container
   * @return Underlying EMCAL container
   */
  const AliEmcalContainer *GetContainer() const { return fkContainer; }

  int GetEntries() const;

  /**
   * Creating forward iterator at the beginning of the container
   * (first entry).
   * @return Iterator at the beginning of the container.
   */
  iterator begin() const { return iterator(this, 0, true); }

  /**
   * Creating forward iterator behind the last entry of the
   * container.
   * @return Iterator behind the container.
   */
  iterator end() const { return iterator(this, GetEntries(), true); }

  /**
   * Creating backward iterator at the end of the container
   * (last entry).
   * @return Iterator at the end of the container
   */
  iterator rbegin() const { return iterator(this, GetEntries()-1, false); }

  /**
   * Creating backward iterator before the beginning of the
   * container.
   * @return Iterator before the container.
   */
  iterator rend() const { return iterator(this, -1, false); }

protected:
  void BuildAcceptIndices();

private:
  const AliEmcalContainer     *fkContainer;         ///< Container to be iterated over
  TArrayI                     fAcceptIndices;       ///< Array of accepted indices
  Bool_t                      fUseAccepted;         ///< Switch between accepted and all objects

  /// \cond CLASSIMP
  ClassDef(AliEmcalIterableContainer, 1);
  /// \endcond
};

#endif /* ALIEMCALITERABLECONTAINER_H */
