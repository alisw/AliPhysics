///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2010
//
//    This file is part of starlight.
//
//    starlight is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    starlight is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with starlight. If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////
//
// File and Version Information:
// $Rev:: 263                         $: revision of last commit
// $Author:: butter                   $: author of last commit
// $Date:: 2016-06-05 00:03:58 +0200 #$: date of last commit
//
// Description:
//
//
//
///////////////////////////////////////////////////////////////////////////


#ifndef STARLIGHTPARTICLE_H
#define STARLIGHTPARTICLE_H


#include "lorentzvector.h"


class starlightParticle : public lorentzVector
{
   public:
      
      starlightParticle();
      starlightParticle ( double px, double py, double pz, double e, double mass, int pdgCode, short charge, 
			  double vx = 0., double vy = 0, double vz = 0, double vt = 0,
			  int firstParent = 0, int lastParent = 0, int firstDaughter = 0, int lastDaughter = 0, int status = 0);
      virtual ~starlightParticle();
   
      void setPdgCode(int pdgCode) { _pdgCode = pdgCode; }
      int getPdgCode() const { return _pdgCode; }
      
      void setCharge(short charge) { _charge = charge; }
      short getCharge() const { return _charge; }

      void setMass(short mass) { _mass = mass; }
      short getMass() const { return _mass; }
      
      void setFirstParent(int parent) { _firstParent = parent; }
      void setLastParent(int parent) { _lastParent = parent; }
      int getFirstParent() const { return _firstParent; }
      int getLastParent() const { return _lastParent; }
      
      void setFirstDaughter(int first) { _firstDaughter = first; }
      int getFirstDaughter() const { return _firstDaughter; }
      
      void setLastDaughter(int first) { _lastDaughter = first; }
      int getLastDaughter() const { return _lastDaughter; }
      
      void setStatus(int status) { _status = status; }
      int getStatus() const { return _status; }
      
      void setVertex(lorentzVector v) { _vertex = v; }
      lorentzVector getVertex() const { return _vertex; }
      
   private:
     
    lorentzVector _vertex;
    
    int _pdgCode;
    short _charge;
    double _mass;
    
    int _firstParent;
    int _lastParent;
    int _firstDaughter;
    int _lastDaughter;

    int _status;
  
};


#endif  // STARLIGHTPARTICLE_H
