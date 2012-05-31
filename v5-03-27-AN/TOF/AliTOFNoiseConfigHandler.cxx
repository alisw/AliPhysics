/*************************************************************************
 * * Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * *                                                                        *
 * * Author: The ALICE Off-line Project.                                    *
 * * Contributors are mentioned in the code where appropriate.              *
 * *                                                                        *
 * * Permission to use, copy, modify and distribute this software and its   *
 * * documentation strictly for non-commercial purposes is hereby granted   *
 * * without fee, provided that the above copyright notice appears in all   *
 * * copies and that both the copyright notice and this permission notice   *
 * * appear in the supporting documentation. The authors make no claims     *
 * * about the suitability of this software for any purpose. It is          *
 * * provided "as is" without express or implied warranty.                  *
 * **************************************************************************/

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  The SAX XML file handler used in the TOFnoiseda                       //
//                                                                        //
//  Author:                                                               //
//    Chiara Zampolli (Chiara.Zampolli@cern.ch)                           //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include <cstdlib>
#include <Riostream.h>

#include <TList.h>
#include <TObject.h>
#include <TXMLAttr.h>
#include <TSAXParser.h>

#include "AliLog.h"
#include "AliTOFNoiseConfigHandler.h"

ClassImp(AliTOFNoiseConfigHandler)

  
//_____________________________________________________________________________
AliTOFNoiseConfigHandler::AliTOFNoiseConfigHandler()
	:TObject(),
	 fDebugFlag(0)
{
	//
	// AliTOFNoiseConfigHandler default constructor
	//
}

//_____________________________________________________________________________
AliTOFNoiseConfigHandler::AliTOFNoiseConfigHandler(const AliTOFNoiseConfigHandler &sh)
	:TObject(sh),
	 fDebugFlag(sh.fDebugFlag)
{
	//
	// AliTOFNoiseConfigHandler copy constructor
	//
}

//_____________________________________________________________________________
AliTOFNoiseConfigHandler &AliTOFNoiseConfigHandler::operator=(const AliTOFNoiseConfigHandler &sh)
{
	//
	// Assignment operator
	//
	if (&sh == this) return *this;
	
	new (this) AliTOFNoiseConfigHandler(sh);
	return *this;
}

//_____________________________________________________________________________
AliTOFNoiseConfigHandler::~AliTOFNoiseConfigHandler()
{
	//
	// AliTOFNoiseConfigHandler destructor
	//	
}

//_____________________________________________________________________________
void AliTOFNoiseConfigHandler::OnStartDocument()
{
	// if something should happen right at the beginning of the
	// XML document, this must happen here
	AliInfo("Reading XML file for TOFnoiseda Config");
}

//_____________________________________________________________________________
void AliTOFNoiseConfigHandler::OnEndDocument()
{
	// if something should happen at the end of the XML document
	// this must be done here
}

//_____________________________________________________________________________
void AliTOFNoiseConfigHandler::OnStartElement(const char *name, const TList *attributes)
{
	// when a new XML element is found, it is processed here

	// set the current system if necessary
	TString strName(name);
	AliDebug(2,Form("name = %s",strName.Data()));
	TXMLAttr* attr;
	TIter next(attributes);
	while ((attr = (TXMLAttr*) next())) {
		TString attrName = attr->GetName();
		AliDebug(2,Form("Name = %s",attrName.Data())); 
		if (attrName == "DebugFlag"){
			TString debugFlag = (TString)(attr->GetValue());
			if (debugFlag == "ON" || debugFlag == "On" || debugFlag == "on"){	
				fDebugFlag = 1;
			}
			else if (debugFlag == "OFF" || debugFlag == "Off"|| debugFlag == "off"){
				fDebugFlag = 0;
			}
			else {
				AliWarning("Invalid Debug Flag. Keeping debug off");
				fDebugFlag = 0;
			}
		}
	}	
	AliDebug(2,Form("Debug Flag = %i",fDebugFlag)); 
	return;
}
//_____________________________________________________________________________
void AliTOFNoiseConfigHandler::OnEndElement(const char *name)
{
	// do everything that needs to be done when an end tag of an element is found
	TString strName(name);
	AliDebug(2,Form("name = %s",strName.Data()));
}

//_____________________________________________________________________________
void AliTOFNoiseConfigHandler::OnCharacters(const char *characters)
{
	// copy the text content of an XML element
	//fContent = characters;
	TString strCharacters(characters);
	AliDebug(2,Form("characters = %s",strCharacters.Data()));
}

//_____________________________________________________________________________
void AliTOFNoiseConfigHandler::OnComment(const char* /*text*/)
{
	// comments within the XML file are ignored
}

//_____________________________________________________________________________
void AliTOFNoiseConfigHandler::OnWarning(const char *text)
{
	// process warnings here
	AliInfo(Form("Warning: %s",text));
}

//_____________________________________________________________________________
void AliTOFNoiseConfigHandler::OnError(const char *text)
{
	// process errors here
	AliError(Form("Error: %s",text));
}

//_____________________________________________________________________________
void AliTOFNoiseConfigHandler::OnFatalError(const char *text)
{
	// process fatal errors here
	AliFatal(Form("Fatal error: %s",text));
}

//_____________________________________________________________________________
void AliTOFNoiseConfigHandler::OnCdataBlock(const char* /*text*/, Int_t /*len*/)
{
	// process character data blocks here
	// not implemented and should not be used here
}

