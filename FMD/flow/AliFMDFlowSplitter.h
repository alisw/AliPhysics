// -*- mode: C++ -*-
/* Copyright (C) 2007 Christian Holm Christensen <cholm@nbi.dk>
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1 of
 * the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 * USA
 */
/** @file 
    @brief Declaration of an Axis in a Flow "histogram" */
#ifndef ALIFMDFLOWSPLITTER_H
#define ALIFMDFLOWSPLITTER_H
#include "flow/AliFMDFlowAxis.h"
#include <TArrayI.h>
#ifndef ROOT_TObject
# include <TObject.h>
#endif
class AliFMDFlowAxis;
class TProfile;
//
// Object used by AliFMDFlowBinned1D to split an event into sub-events.
// Default is to split randomly.  
// User defined derived classes can do other stuff.
//
//______________________________________________________
/** @class AliFMDFlowSplitter flow/AliFMDFlowSplitter.h <flow/AliFMDFlowSplitter.h>

    @brief Object used by AliFMDFlowBinned1D to split an event into
    sub-events. Default is to split randomly.  User defined derived
    classes can do other stuff.  

    @ingroup c_binned 
*/
class AliFMDFlowSplitter : public TObject
{
public:
  /** Constructor 
      @param axis Axis object */
  AliFMDFlowSplitter() {}
  /** Copy constructor 
      @param o Other splitter to copy from */
  AliFMDFlowSplitter(const AliFMDFlowSplitter& o) {}
  /** Destructor */
  virtual ~AliFMDFlowSplitter() {}
  /** Assignment operator  
      @param o Other splitter to assign from. 
      @return Reference to this object */
  AliFMDFlowSplitter& operator=(const AliFMDFlowSplitter& o) { return *this; }
  /** Called at the beginning of an event */
  virtual void Begin() {}
  /** Prepare for an event 
      @param phis List of phis. 
      @param xs   List of bin variable 
      @param n    Number of entries in @a phis and @a n 
      @return true on success, false otherwise */ 
  virtual void Event(Double_t* phis, Double_t* xs, ULong_t n) {}
  /** Decide whether entry should go in A or B sub-event. 
      @param entry The entry number 
      @return true if this should go in sub-event A */
  virtual Bool_t Select(ULong_t entry) const;
  /** Called at the end of an event */
  virtual void End() {}
protected:
  /** Reference to axis object */
  ClassDef(AliFMDFlowSplitter,1) // Split events
};

//____________________________________________________________________
class AliFMDFlowShuffle : public AliFMDFlowSplitter
{
public:
  /** Constuctor */
  AliFMDFlowShuffle() {}
  /** Prepare for an event 
      @param phis List of phis. 
      @param xs   List of bin variable 
      @param n    Number of entries in @a phis and @a n 
      @return true on success, false otherwise */ 
  virtual void Event(Double_t* phis, Double_t* xs, ULong_t n);
  /** Decide whether entry should go in A or B sub-event. 
      @param entry The entry number 
      @return true if this should go in sub-event A */
  virtual Bool_t Select(ULong_t entry) const;
protected:
  /** Function to shuffle */
  void Shuffle();
  /** Array of indicies */ 
  TArrayI fIdx;
  /** Curren number of entries */ 
  ULong_t fN;
  ClassDef(AliFMDFlowShuffle,1) // Split events
};


#endif
//
// EOF
//
