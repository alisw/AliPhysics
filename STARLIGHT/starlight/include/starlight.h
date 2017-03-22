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


#ifndef STARLIGHT_H
#define STARLIGHT_H


#include <string>

#include "upcevent.h"
#include "eventchannel.h"


class inputParameters;
class beam;
class beamBeamSystem;


class starlight {

public:
      
	starlight();
	~starlight();
      
	bool init();

	upcEvent produceEvent();
      
        std::string   baseFileName  () const { return _baseFileName;		  }
	unsigned long nmbAttempts   () const { return _eventChannel->nmbAttempts(); }
	unsigned long nmbAccepted   () const { return _eventChannel->nmbAccepted(); } 
	double getTotalCrossSection () const { return _eventChannel->getTotalChannelCrossSection(); } 
	void setInputParameters(inputParameters* inputParams) { _inputParameters = inputParams; }   

private:
      
	bool luminosityTableIsValid() const;

	bool createEventChannel();
      
	beamBeamSystem*    _beamSystem;
	eventChannel*      _eventChannel;
	unsigned int       _nmbEventsPerFile;
	std::string        _baseFileName;
	std::string        _lumLookUpTableFileName;
	bool               _isInitialised;
	inputParameters*   _inputParameters;
};


#endif  // STARLIGHT_H
