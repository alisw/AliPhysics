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
// $Rev:: 293                         $: revision of last commit
// $Author:: butter                   $: author of last commit
// $Date:: 2017-11-11 15:46:05 +0100 #$: date of last commit
//
// Description:
//
//
//
///////////////////////////////////////////////////////////////////////////


#ifndef STARLIGHTSTANDALONE_H
#define STARLIGHTSTANDALONE_H


#include <string>

class upcEvent;

class starlight;
class inputParameters;


class starlightStandalone {

public:

	starlightStandalone();
	~starlightStandalone();

	bool init();  ///< reads configuration file and initializes startlight
	bool run ();  ///< creates output file and runs starlight

        const std::string& baseFileName     () const { return _baseFileName;        }
	const std::string& configFileName   () const { return _configFileName;      }  ///< returns path to config file
	const std::string& eventDataFileName() const { return _eventDataFileName;   }  ///< returns path to output file

        void setBaseFileName     (const std::string& baseFileName     ) { _baseFileName      = baseFileName;    }  ///< sets path to base file
	void setConfigFileName   (const std::string& configFileName   ) { _configFileName    = configFileName;    }  ///< sets path to config file
	void setEventDataFileName(const std::string& eventDataFileName) { _eventDataFileName = eventDataFileName; }  ///< sets path to output file
	
	void boostEvent(upcEvent &e); ///< Boost event from beam CMS to lab system

private:

        std::string      _baseFileName;       ///< path to base filename
	std::string      _configFileName;     ///< path to configuration file
	std::string      _eventDataFileName;  ///< path to output file

	starlight*       _starlight;         ///< pointer to starlight instance
        inputParameters* _inputParameters;   ///< pointer to parameter instance
	randomGenerator* _randomGenerator;   ///< pointer to random instance
	unsigned int     _nmbEventsTot;      ///< total number of events to generate (taken from configuration file)
	unsigned int     _nmbEventsPerFile;  ///< maximum number of events written to a single file (not yet implemented)

};


#endif  // STARLIGHTSTANDALONE_H
