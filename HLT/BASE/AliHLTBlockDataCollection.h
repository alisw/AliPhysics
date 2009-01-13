// -*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTBLOCKDATACOLLECTION_H
#define ALIHLTBLOCKDATACOLLECTION_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTBlockDataCollection.h
    @author Matthias Richter
    @date   
    @brief  A collection of AliHLTComponentBlockData descriptors providing
            argument parsing and basic selection.
*/

#include "AliHLTLogging.h"
#include "vector"
#include "TObject.h"

/**
 * @class AliHLTBlockDataCollection
 * Class handles a list of AliHLTComponentBlockData entries and parsing of
 * argument list to fill it. Originally taken from AliHLTBlickFilterComponent,
 * but decided to be commonly of benefit.
 *
 * See ScanArgument() function for description of available arguments
 * <pre>
 * -datatype ID ORIGIN
 * -typeid ID
 * -origin ORIGIN
 * -dataspec SPEC
 * </pre>
 *
 * @ingroup alihlt_base
 */
class AliHLTBlockDataCollection : public TObject, public AliHLTLogging
{
 public:
  /** standard constructor */
  AliHLTBlockDataCollection();
  /** destructor */
  virtual ~AliHLTBlockDataCollection();

  /**
   * Add data block descriptor.
   */
  int Add(const AliHLTComponentBlockData& block);

  /**
   * Check if the data block is selected by the filter rules.
   * @return 1 if selected
   */
  int IsSelected(const AliHLTComponentBlockData& block);

  /**
   * Scan argument and read block descriptor data.
   * The function is invoked by components in the course of argument
   * scan.
   *
   * Scan the list for known arguments, terminates at the first unknown argument.
   * Recognized arguments:
   * \li -datatype     <i> id origin      </i>                            <br>
   *      e.g. <tt> -datatype 'ESD_TREE' 'TPC ' </tt>                     <br>
   *      \b Note: due to the 4-character data origin it might be necessary to
   *      append a blank to the detectorname, e.g. <tt>TPC -> 'TPC '</tt>
   *
   * \li -origin  <i> origin  </i>                                        <br>
   *      e.g. -origin 'TPC ', \b Note:  the filter rule has type id 'ANY'
   *
   * \li -typeid  <i> id      </i>                                        <br>
   *      e.g. -typeid ESD_TREE, \b Note: the filter rule has origin 'ANY'
   *
   * \li -dataspec     <i> specification </i>                             <br>
   *      data specification treated as decimal number or hex number if
   *      prepended by '0x'
   * 
   * @return number of arguments which have been treated.
   */
  int ScanArgument(int argc, const char** argv );

  /**
   * Check collection for content.
   * @return 1 if empty, 0 if content available
   */
  int IsEmpty();
 protected:

 private:
  /** copy constructor prohibited */
  AliHLTBlockDataCollection(const AliHLTBlockDataCollection&);
  /** assignment operator prohibited */
  AliHLTBlockDataCollection& operator=(const AliHLTBlockDataCollection&);

  /** filtering rules, only the data type and specification members are use */
  vector<AliHLTComponentBlockData> fFilterRules;                       //! transient

  ClassDef(AliHLTBlockDataCollection, 0)
};
#endif
