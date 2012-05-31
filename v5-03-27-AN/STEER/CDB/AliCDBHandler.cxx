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
//  The SAX XML file handler used in the CDBManager                       //
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
#include "AliCDBHandler.h"

ClassImp(AliCDBHandler)

  
//_____________________________________________________________________________
AliCDBHandler::AliCDBHandler()
	:TObject(),
	 fRun(-1),
	 fStartRunRange(-1),
	 fEndRunRange(-1),
	 fOCDBFolder("")
{
	//
	// AliCDBHandler default constructor
	//
}

//_____________________________________________________________________________
AliCDBHandler::AliCDBHandler(Int_t run)
	:TObject(),
	 fRun(run),
	 fStartRunRange(-1),
	 fEndRunRange(-1),
	 fOCDBFolder("")
{
	//
	// AliCDBHandler constructor with requested run
	//
}

//_____________________________________________________________________________
AliCDBHandler::AliCDBHandler(const AliCDBHandler &sh)
	:TObject(sh),
	 fRun(sh.fRun),
	 fStartRunRange(sh.fStartRunRange),
	 fEndRunRange(sh.fEndRunRange),
	 fOCDBFolder(sh.fOCDBFolder)
{
	//
	// AliCDBHandler copy constructor
	//
}

//_____________________________________________________________________________
AliCDBHandler &AliCDBHandler::operator=(const AliCDBHandler &sh)
{
	//
	// Assignment operator
	//
	if (&sh == this) return *this;
	
	new (this) AliCDBHandler(sh);
	return *this;
}

//_____________________________________________________________________________
AliCDBHandler::~AliCDBHandler()
{
	//
	// AliCDBHandler destructor
	//	
}

//_____________________________________________________________________________
void AliCDBHandler::OnStartDocument()
{
	// if something should happen right at the beginning of the
	// XML document, this must happen here
	AliInfo("Reading XML file for LHCPerdiod <-> Run Range correspondance");
}

//_____________________________________________________________________________
void AliCDBHandler::OnEndDocument()
{
	// if something should happen at the end of the XML document
	// this must be done here
}

//_____________________________________________________________________________
void AliCDBHandler::OnStartElement(const char *name, const TList *attributes)
{
	// when a new XML element is found, it is processed here

	// set the current system if necessary
	TString strName(name);
	AliDebug(2,Form("name = %s",strName.Data()));
	Int_t startRun=-1;
	Int_t endRun=-1;
	TXMLAttr* attr;
	TIter next(attributes);
	while ((attr = (TXMLAttr*) next())) {
		TString attrName = attr->GetName();
		AliDebug(2,Form("Name = %s",attrName.Data())); 
		if (attrName == "StartRunRange"){
			startRun = (Int_t)(((TString)(attr->GetValue())).Atoi());
			AliDebug(2,Form("startRun = %d",startRun));
		}
		if (attrName == "EndRunRange"){
			endRun = (Int_t)(((TString)(attr->GetValue())).Atoi());
			AliDebug(2,Form("endRun = %d",endRun));
		}				
		if (attrName == "OCDBFolder"){
			if (fRun>=startRun && fRun<=endRun && startRun!=-1 && endRun!=-1){
				fOCDBFolder = (TString)(attr->GetValue());
				AliDebug(2,Form("OCDBFolder = %s",fOCDBFolder.Data()));
				fStartRunRange = startRun;
				fEndRunRange = endRun;
			}
		}
	}	
	return;
}
//_____________________________________________________________________________
void AliCDBHandler::OnEndElement(const char *name)
{
	// do everything that needs to be done when an end tag of an element is found
	TString strName(name);
	AliDebug(2,Form("name = %s",strName.Data()));
}

//_____________________________________________________________________________
void AliCDBHandler::OnCharacters(const char *characters)
{
	// copy the text content of an XML element
	//fContent = characters;
	TString strCharacters(characters);
	AliDebug(2,Form("characters = %s",strCharacters.Data()));
}

//_____________________________________________________________________________
void AliCDBHandler::OnComment(const char* /*text*/)
{
	// comments within the XML file are ignored
}

//_____________________________________________________________________________
void AliCDBHandler::OnWarning(const char *text)
{
	// process warnings here
	AliInfo(Form("Warning: %s",text));
}

//_____________________________________________________________________________
void AliCDBHandler::OnError(const char *text)
{
	// process errors here
	AliError(Form("Error: %s",text));
}

//_____________________________________________________________________________
void AliCDBHandler::OnFatalError(const char *text)
{
	// process fatal errors here
	AliFatal(Form("Fatal error: %s",text));
}

//_____________________________________________________________________________
void AliCDBHandler::OnCdataBlock(const char* /*text*/, Int_t /*len*/)
{
	// process character data blocks here
	// not implemented and should not be used here
}

