//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//      Copyright (C) 1998      Caltech, UCSB
//
// Module: EvtGen/EvtId.hh
//
// Description:Class for particle Id used in EvtGen.
//
// Modification history:
//
//    DJL/RYD     May 26, 1998         Module created
//
//------------------------------------------------------------------------

#ifndef EVTID_HH
#define EVTID_HH

#include <iostream>
#include <string>

class EvtId {

public:

  //need a default constructor
  EvtId():_id(-1),_alias(-1){}

  EvtId(int id,int alias):_id(id),_alias(alias){} 

  friend std::ostream& operator<<(std::ostream& s, const EvtId& v);  

  int operator==(const EvtId& id) const { return _id==id._id; }
  int operator!=(const EvtId& id) const { return _id!=id._id; }
  int operator<(const EvtId& id) const { return _id<id._id; }

  int isConjugate(const EvtId & id) const;

  int getId() const { return _id;}

  int getAlias() const { return _alias;}

  int isAlias() const { return _id!=_alias;}

  std::string getName() const;

private:

  //particle number 0..n. The order of particles are determined 
  //by the order in pdt.table
  int _id;
  //if the particle is an alias to another particle alias!=id
  //The only place where the alias should be used is for looking
  //up decays in the decay table.
  int _alias;
  
}; 

#endif
