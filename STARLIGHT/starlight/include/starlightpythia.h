/*
    <one line to give the library's name and an idea of what it does.>
    Copyright (C) 2011  Oystein Djuvsland <oystein.djuvsland@gmail.com>

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with this library; if not, write to the Free Software
    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/


#ifndef STARLIGHTPYTHIA_H
#define STARLIGHTPYTHIA_H

#include "upcevent.h"
#include "inputParameters.h"
#include "beambeamsystem.h"
#include "eventchannel.h"

class spectrum;

class starlightPythia : public eventChannel
{

public:

    starlightPythia(const inputParameters& input,randomGenerator* randy,beamBeamSystem& bbsystem);
    virtual ~starlightPythia();

    int init(std::string pythiaParams, bool fullEventRecord = false);

    virtual upcEvent produceEvent();

    virtual upcEvent produceSingleEvent(int /*zdirection*/, float /*egamma*/){return upcEvent();}

    virtual upcEvent produceDoubleEvent(){return upcEvent();}

    virtual starlightConstants::event produceEvent(int& /*ievent*/){ return starlightConstants::event();}

    void setSingleMode() {
        _doDoubleEvent = false;
    }

    void setDoubleMode() {
        _doDoubleEvent = true;
    }

    void setMinGammaEnergy(double energy) {
        _minGammaEnergy = energy;
    }

    void setMaxGammaEnergy(double energy) {
        _maxGammaEnergy = energy;
    }

    void setFullEventRecord(bool fer = true) { _fullEventRecord = fer; }
    
private:

    /** Contains the photon spectrum */
    spectrum *_spectrum;

    /** Should we produce a double event? */
    bool _doDoubleEvent;

    /** Min gamma energy */
    double _minGammaEnergy;

    /** Max gamma energy */
    double _maxGammaEnergy;

    /** Full event record or not */
    bool _fullEventRecord;
    
    /** Prohibited */
    starlightPythia();
    starlightPythia(const starlightPythia& other);
    starlightPythia& operator=(const starlightPythia& other);
    bool operator==(const starlightPythia& other) const;

};

#endif // STARLIGHTPYTHIA_H
