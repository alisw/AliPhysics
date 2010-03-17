//-*- Mode: C++ -*-
// $Id$

#ifndef ALI_HLT_PENDOLINO_LOGGER_H
#define ALI_HLT_PENDOLINO_LOGGER_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTPendolinoLogger.h
    @author Sebastian Bablok
    @date   
    @brief  
*/


#include <TObject.h>




/**
 * Class that defines the interface, that a Logger for the Pendolino has to 
 * implement. 
 *
 * @author Sebastian Bablok
 *
 * @date 2007-10-24
 */
class AliHLTPendolinoLogger : public TObject {
    public:

		/**
		 * Constructor for AliHLTPendolinoLogger
		 */
		AliHLTPendolinoLogger();

		/**
		 * Destructor for AliHLTPendolinoLogger
		 */
		virtual ~AliHLTPendolinoLogger();

		/**
		 * logging function for the Pendolino. Pure virtual, is implemented in 
		 * inherited class
		 *
		 * @param detector the detector, which is making the log entry
		 * @param msg the log message
		 */
		virtual void log(const char* detector, const char* msg) = 0;

	protected:

	private:

		ClassDef(AliHLTPendolinoLogger, 0);	
};

#endif

