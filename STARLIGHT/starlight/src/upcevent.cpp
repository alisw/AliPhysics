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


#include "upcevent.h"


upcEvent::upcEvent() :
        _particles(0)
        ,_vertices(0)
	,_isGammaavm(false)
	,_bSlope(0)
	,_t(0)
	,_Egam(0)
	,_Epom(0)
{
  for (int i=0; i<2; ++i)
    _ptGam[i] = _ptPom[i] = 0.0;
}

upcEvent::upcEvent(starlightConstants::event &ev) :
        _particles(0)
        ,_vertices(0)
	,_isGammaavm(false)
	,_bSlope(0)
	,_t(0)
	,_Egam(0)
	,_Epom(0)
{
  for (int i=0; i<2; ++i)
    _ptGam[i] = _ptPom[i] = 0.0;

  for(int i = 0; i < ev._numberOfTracks; i++)
    {
      starlightParticle p(
			  ev.px[i], 
			  ev.py[i], 
			  ev.pz[i], 
			  starlightConstants::UNKNOWN, 
			  starlightConstants::UNKNOWN, 
			  ev._fsParticle[i],
			  ev._charge[i]
			  );
      addParticle(p);
    }
}

upcEvent::~upcEvent()
{ }


upcEvent& upcEvent::operator=(const upcEvent& rhs)
{
  if(this != &rhs)
  {
    _particles     = rhs._particles;
    _vertices      = rhs._vertices;
    _gammaEnergies = rhs._gammaEnergies;
    _isGammaavm    = rhs._isGammaavm;
    _bSlope        = rhs._bSlope;
    _t             = rhs._t;
    _Egam          = rhs._Egam;
    _Epom          = rhs._Epom;
    for (int i=0; i<2; ++i) {
      _ptGam[i] = rhs._ptGam[i];
      _ptPom[i] = rhs._ptPom[i];
    }
  }
  return *this;
}

upcEvent& upcEvent::operator+(const upcEvent& ev)
{
  for(unsigned int n = 0; n < ev._particles.size(); n++)
  {
    this->_particles.push_back(ev._particles.at(n));
  }
  for(unsigned int n = 0; n < ev._vertices.size(); n++)
  {
    this->_vertices.push_back(ev._vertices.at(n));
  }
 for(unsigned int n = 0; n < ev._gammaEnergies.size(); n++)
  {
    this->_gammaEnergies.push_back(ev._gammaEnergies.at(n));
  }
  // what to do with _bSlope, _t2, ...?
  return *this;
}

void upcEvent::boost(double rapidity)
{
    vector3 boostVector(0, 0, tanh(rapidity));
    std::vector<starlightParticle>::iterator part = _particles.begin();
      
    for (part = _particles.begin(); part != _particles.end(); part++)
    {
      (*part).Boost(boostVector);
    }
}
