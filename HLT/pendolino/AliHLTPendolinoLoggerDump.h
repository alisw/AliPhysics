//-*- Mode: C++ -*-
// $Id$

#ifndef ALI_HLT_PENDOLINO_LOGGER_DUMP_H
#define ALI_HLT_PENDOLINO_LOGGER_DUMP_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTPendolinoLoggerDump.h
    @author Sebastian Bablok
    @date   
    @brief  
*/


#include "AliHLTPendolinoLogger.h"



/**
 * Class that implements the interface for a Pendolino Logger and dumps 
 * the log messages. No output is made! 
 *
 * @author Sebastian Bablok
 *
 * @date 2007-10-29
 */
class AliHLTPendolinoLoggerDump : public AliHLTPendolinoLogger {
    public:

		/**
		 * Constructor for AliHLTPendolinoLoggerDump
		 */
		AliHLTPendolinoLoggerDump();

		/**
		 * Destructor for AliHLTPendolinoLoggerDump
		 */
		virtual ~AliHLTPendolinoLoggerDump();

		/**
		 * Implementated logging interface
		 *
		 * @param detector the detector for which the log entry shall be made.
		 * @param msg the log message
		 */
		virtual void log(const char* detector, const char* msg);

	protected:

	private:

		ClassDef(AliHLTPendolinoLoggerDump, 0);
};

#endif

