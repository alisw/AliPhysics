//STARTHEADER
// $Id: NestedDefsPlugin.hh 1394 2009-01-17 05:08:31Z soyez $
//
// Copyright (c) 2007-2008, Matteo Cacciari, Gavin Salam and Gregory Soyez
//
//----------------------------------------------------------------------
// This file is part of FastJet.
//
//  FastJet is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  The algorithms that underlie FastJet have required considerable
//  development and are described in hep-ph/0512210. If you use
//  FastJet as part of work towards a scientific publication, please
//  include a citation to the FastJet paper.
//
//  FastJet is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with FastJet; if not, write to the Free Software
//  Foundation, Inc.:
//      59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//----------------------------------------------------------------------
//ENDHEADER

#ifndef __NESTEDALGSPLUGIN_HH__
#define __NESTEDALGSPLUGIN_HH__

#include "fastjet/JetDefinition.hh"
#include <list>
#include <memory>
#include <cmath>

// questionable whether this should be in fastjet namespace or not...
namespace fastjet {      // defined in fastjet/internal/base.hh

// another forward declaration to reduce includes
class PseudoJet;

//----------------------------------------------------------------------
//
/// NestedAglsPlugin is a plugin for fastjet (v2.4 upwards) that, given
/// a list of jet definitions, performs the clustering by feeding the 
/// particles to the first algorithm and then, successively feeding the 
/// output to the next algorithm in the list.
//
class NestedDefsPlugin : public JetDefinition::Plugin {
public:
  /// Main constructor for the NestedDefs Plugin class.  
  ///
  /// The argument is an initialised list of jet algorithms
  NestedDefsPlugin (std::list<JetDefinition> &defs) :
    _defs(defs){}

  /// copy constructor
  NestedDefsPlugin (const NestedDefsPlugin & plugin) {
    *this = plugin;
  }

  // the things that are required by base class
  virtual std::string description () const;
  virtual void run_clustering(ClusterSequence &) const;
  /// the plugin mechanism's standard way of accessing the jet radius
  /// here we return the R of the last alg in the list
  virtual double R() const {return _defs.rbegin()->R();}

private:
  std::list<JetDefinition> _defs;
};

} // fastjet namespace 

#endif // __SISCONEPLUGIN_HH__

