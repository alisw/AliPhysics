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
// $Rev:: 133                         $: revision of last commit
// $Author:: odjuvsla                 $: author of last commit
// $Date:: 2013-09-05 22:08:42 +0200 #$: date of last commit
//
// Description:
//
//
//
///////////////////////////////////////////////////////////////////////////


#include "starlightparticle.h"


starlightParticle::starlightParticle() :
   lorentzVector()
   ,_vertex(0,0,0,0)
   ,_pdgCode(0)
   ,_charge(-999)
   ,_mass(-1)
   ,_firstParent(0)
   ,_lastParent(0)
   ,_firstDaughter(0)
   ,_lastDaughter(0)
   ,_status(0)
{ }


starlightParticle::starlightParticle(double px, double py, double pz, double e, double mass, int pdgCode, short charge,
				     double vx, double vy, double vz, double vt,
				     int firstParent, int lastParent, int firstDaughter, int lastDaughter, int status) :
lorentzVector(px, py, pz, e)
,_vertex(vx,vy,vz,vt)
,_pdgCode(pdgCode)
,_charge(charge)
,_mass(mass)
,_firstParent(firstParent)
,_lastParent(lastParent)
,_firstDaughter(firstDaughter)
,_lastDaughter(lastDaughter)
,_status(status)
{ }


starlightParticle::~starlightParticle()
{ }
