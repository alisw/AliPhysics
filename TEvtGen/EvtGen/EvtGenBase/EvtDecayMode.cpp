//-----------------------------------------------------------------------
// File and Version Information: 
//      $Id: EvtDecayMode.cpp,v 1.3 2009-03-16 15:53:27 robbep Exp $
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
#include "EvtGenBase/EvtPatches.hh"

// Parses a decay string to identify the name
// of the mother and of the daughters. The string should 
// be in standard format e.g. "B+ -> pi+ pi+ pi-"

#include "EvtGenBase/EvtPatches.hh"
#include <assert.h>
#include <iostream>
#include "EvtGenBase/EvtDecayMode.hh"
#include "EvtGenBase/EvtReport.hh"
using std::endl;
using std::ostream;

using std::string;
using std::vector;


EvtDecayMode::EvtDecayMode(std::string mother,vector<string> dau)
  : _mother(mother),
    _dau(dau)
{
}


EvtDecayMode::EvtDecayMode(const EvtDecayMode& other)
  : _mother(other._mother),
    _dau(other._dau)
{
}


EvtDecayMode::EvtDecayMode(const char* decay)
{
  // Parse the decay string, it should be in a standard 
  // format, e.g. "B+ -> pi+ pi+ pi-" with all spaces
  
  string s(decay);

  // mother

  string::size_type i = s.find_first_not_of(" ");
  string::size_type j = s.find_first_of(" ",i);

  if(i == string::npos) {

    report(Severity::Info,"EvtGen") << "No non-space character found" << endl;
    assert(0);
  }

  if(j == string::npos) {
    
    report(Severity::Info,"EvtGen") << "No space before -> found" << endl;
    assert(0);
  }

  _mother = string(s,i,j-i);

  i = s.find_first_not_of(" ",j);
  j = s.find_first_of("->",j);
  if(i != j) {

    report(Severity::Info,"EvtGen") << "Multiple mothers?" << i << "," << j << endl;
    assert(0);
  }
  j += 2;

  while(1) {

    i = s.find_first_not_of(" ",j);
    j = s.find_first_of(" ",i);

    if(i == string::npos) break;
    if(j == string::npos) {
      _dau.push_back(string(s,i,s.size()-i+1));
      break;
    } else {
      _dau.push_back(string(s,i,j-i));
    }
  }
}



EvtDecayMode::~EvtDecayMode()
{}


const char* EvtDecayMode::mother() const
{
  return _mother.c_str();
}


int EvtDecayMode::nD() const
{
  return _dau.size();
}


const char* EvtDecayMode::dau(int i) const
{
  assert(0<=i && i< (int) _dau.size());
  return _dau[i].c_str();
}

std::string EvtDecayMode::mode() const
{
  string ret = _mother + string(" -> ");

  for(size_t i=0;i<_dau.size()-1;i++) {
    ret += string(_dau[i]) + string(" ");
  }
  ret += _dau[_dau.size()-1];
  return ret;
}


ostream& EvtDecayMode::print(ostream& os) const
{
  os << _mother.c_str() << " ->";
  for(size_t i=0;i<_dau.size();i++) {
    os << " " << _dau[i].c_str();
  }
  return os;
}


std::string EvtDecayMode::m(EvtCyclic3::Pair i) const
{
  string s("m(");
  s.append(dau(EvtCyclic3::first(i)));
  s.append(",");
  s.append(dau(EvtCyclic3::second(i)));
  s.append(")");
  return s;
}


std::string EvtDecayMode::q(EvtCyclic3::Pair i) const
{
  string s("q(");
  s.append(dau(EvtCyclic3::first(i)));
  s.append(",");
  s.append(dau(EvtCyclic3::second(i)));
  s.append(")");
  return s;
}


std::string EvtDecayMode::dal(EvtCyclic3::Pair i, EvtCyclic3::Pair j) const
{
  string s(q(i));
  s.append(":");
  s.append(q(j));
  return s;
}


ostream& operator<<(ostream& os, const EvtDecayMode& mode)
{
  mode.print(os);
  return os;
}
