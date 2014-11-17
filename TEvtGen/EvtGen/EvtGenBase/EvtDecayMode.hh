//-----------------------------------------------------------------------
// File and Version Information: 
//      $Id: EvtDecayMode.hh,v 1.2 2009-03-16 16:44:33 robbep Exp $
// 
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations. If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information:
//      Copyright (C) 1998 Caltech, UCSB
//
// Module creator:
//      Alexei Dvoretskii, Caltech, 2001-2002.
//-----------------------------------------------------------------------

#ifndef EVT_DECAY_MODE_HH
#define EVT_DECAY_MODE_HH

#include <string>
#include <vector>
#include "EvtGenBase/EvtCyclic3.hh"

#include <iosfwd>

class EvtDecayMode {

public:

  EvtDecayMode(const char* decay);
  EvtDecayMode(const EvtDecayMode& other);
  EvtDecayMode(std::string mother,std::vector<std::string> dau);
  ~EvtDecayMode();

  const char* mother() const;
  int nD() const;
  const char* dau(int i) const; 

  std::ostream& print(std::ostream&) const;


  // Frequent name combinations

  std::string m(EvtCyclic3::Pair i) const;
  std::string q(EvtCyclic3::Pair i) const;
  std::string dal(EvtCyclic3::Pair i, EvtCyclic3::Pair j) const;
  std::string mode() const;

private:

  std::string _mother;
  std::vector<std::string> _dau;

};


std::ostream& operator<<(std::ostream&,const EvtDecayMode&);

#endif
