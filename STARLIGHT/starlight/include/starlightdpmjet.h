
/*
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) <year>  <name of author>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

*/

#ifndef STARLIGHTDPMJET_H
#define STARLIGHTDPMJET_H

#include <eventchannel.h>
#include <spectrum.h>

class starlightDpmJet : public eventChannel
{

public:
   
    starlightDpmJet(const inputParameters& input,beamBeamSystem& beamsystem);

    int init();
    
    virtual upcEvent produceEvent();
    
    virtual upcEvent produceSingleEvent(int zdirection, float egamma);
    
    virtual upcEvent produceDoubleEvent();
    
    virtual starlightConstants::event produceEvent(int& /*ievent*/) { return starlightConstants::event(); }
    
    void setSingleMode() { _doDoubleEvent = false; }
    
    void setDoubleMode() { _doDoubleEvent = true; }
    
    void setMinGammaEnergy(double energy) { _minGammaEnergy = energy; }
    
    void setMaxGammaEnergy(double energy) { _maxGammaEnergy = energy; }
    
    void setProtonMode(bool v = true) { _protonMode = v; }
    
    
   private:
     
      /** Contains the photon spectrum */
      spectrum *_spectrum;

      /** Should we produce a double event? */
      bool _doDoubleEvent;
      
      /** Min gamma energy */
      double _minGammaEnergy;
      
      /** Max gamma energy */
      double _maxGammaEnergy;
      

      /** Proton-nucleus mode */
      double _protonMode;
      
      /** Default constructor not implemented */
    starlightDpmJet();
      
};

#endif // STARLIGHTDPMJET_H
