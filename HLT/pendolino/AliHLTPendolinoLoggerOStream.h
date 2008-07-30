//-*- Mode: C++ -*-
// $Id$

#ifndef ALI_HLT_PENDOLINO_LOGGER_OSTREAM_H
#define ALI_HLT_PENDOLINO_LOGGER_OSTREAM_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTPendolinoLoggerOStream.h
    @author Sebastian Bablok
    @date   
    @brief  
*/


#include "AliHLTPendolinoLogger.h"



/**
 * Class that implements the interface for a Pendolino Logger and simply 
 * streams the log message to the command line. 
 *
 * @author Sebastian Bablok
 *
 * @date 2007-10-24
 */
class AliHLTPendolinoLoggerOStream : public AliHLTPendolinoLogger {
    public:

		/**
		 * Constructor for AliHLTPendolinoLoggerOStream
		 */
		AliHLTPendolinoLoggerOStream();

		/**
		 * Destructor for AliHLTPendolinoLoggerOStream
		 */
		virtual ~AliHLTPendolinoLoggerOStream();

		/**
		 * Implementated logging interface
		 *
		 * @param detector the detector for which the log entry shall be made.
		 * @param msg the log message
		 */
		virtual void log(const char* detector, const char* msg);

	protected:

	private:

		ClassDef(AliHLTPendolinoLoggerOStream, 0);
};

#endif

