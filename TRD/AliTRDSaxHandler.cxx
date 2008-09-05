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

/* $Id: AliTRDSaxHandler.cxx 26327 2008-06-02 15:36:18Z cblume $ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  The SAX XML file handler used in the preprocessor                     //
//                                                                        //
//  Author:                                                               //
//    Frederick Kramer (kramer@ikf.uni-frankfurt.de)                      //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include <cstdlib>
#include <Riostream.h>
#include <TList.h>
#include <TXMLAttr.h>
#include <TSAXParser.h>
#include <TObjArray.h>

#include "AliLog.h"

#include "AliTRDSaxHandler.h"
#include "AliTRDgeometry.h"

#include "Cal/AliTRDCalDCS.h"
#include "Cal/AliTRDCalDCSFEE.h"
#include "Cal/AliTRDCalDCSPTR.h"
#include "Cal/AliTRDCalDCSGTU.h"

ClassImp(AliTRDSaxHandler)

//_____________________________________________________________________________
AliTRDSaxHandler::AliTRDSaxHandler()
  :TObject()
  ,fHandlerStatus(0)
  ,fNDCSPTR(0)
  ,fNDCSGTU(0)
  ,fFEEArr(new TObjArray(540))
  ,fPTRArr(new TObjArray(6))
  ,fGTUArr(new TObjArray(19))
  ,fSystem(0)
  ,fCurrentSM(0)
  ,fCurrentStack(0)
  ,fContent(0)
  ,fDCSFEEObj(0)
  ,fDCSPTRObj(0)
  ,fDCSGTUObj(0)
  ,fCalDCSObj(new AliTRDCalDCS())
{
  //
  // AliTRDSaxHandler default constructor
  //
  fFEEArr->SetOwner();
  fPTRArr->SetOwner();
  fGTUArr->SetOwner();
}

//_____________________________________________________________________________
AliTRDSaxHandler::AliTRDSaxHandler(const AliTRDSaxHandler &sh)
  :TObject(sh)
  ,fHandlerStatus(0)
  ,fNDCSPTR(0)
  ,fNDCSGTU(0)
  ,fFEEArr(0)
  ,fPTRArr(0)
  ,fGTUArr(0)
  ,fSystem(0)
  ,fCurrentSM(0)
  ,fCurrentStack(0)
  ,fContent(0)
  ,fDCSFEEObj(0)
  ,fDCSPTRObj(0)
  ,fDCSGTUObj(0)
  ,fCalDCSObj(0)
{
  //
  // AliTRDSaxHandler copy constructor
  //
}

//_____________________________________________________________________________
AliTRDSaxHandler &AliTRDSaxHandler::operator=(const AliTRDSaxHandler &sh)
{
  //
  // Assignment operator
  //
  if (&sh == this) return *this;

  new (this) AliTRDSaxHandler(sh);
  return *this;
}

//_____________________________________________________________________________
AliTRDSaxHandler::~AliTRDSaxHandler()
{
  //
  // AliTRDSaxHandler destructor
  //

  if (fFEEArr) {
    delete fFEEArr;
    fFEEArr    = 0x0;
  }
  if (fPTRArr) {
    delete fPTRArr;
    fPTRArr    = 0x0;
  }
  if (fGTUArr) {
    delete fGTUArr;
    fGTUArr    = 0x0;
  }
  if (fCalDCSObj) {
    delete fCalDCSObj;
    fCalDCSObj = 0x0;
  }

}

//_____________________________________________________________________________
AliTRDCalDCS* AliTRDSaxHandler::GetCalDCSObj()
{
  // put the arrays in the global calibration object and return this
  fCalDCSObj->SetNumberOfTimeBins(0); //test test test
  fCalDCSObj->SetFEEArr(fFEEArr);
  fCalDCSObj->SetPTRArr(fPTRArr);
  fCalDCSObj->SetGTUArr(fGTUArr);
  return fCalDCSObj;
}

//_____________________________________________________________________________
void AliTRDSaxHandler::OnStartDocument()
{
  // if something should happen right at the beginning of the
  // XML document, this must happen here
}

//_____________________________________________________________________________
void AliTRDSaxHandler::OnEndDocument()
{
  // if something should happen at the end of the XML document
  // this must be done here
}

//_____________________________________________________________________________
void AliTRDSaxHandler::OnStartElement(const char *name, const TList *attributes)
{
  // when a new XML element is found, it is processed here
  fContent    = "";
  Int_t dcsId = 0;
  TString strName  = name;
  TString dcsTitle = "";

  // set the current system if necessary
  if (strName.Contains("FEE")) fSystem = kInsideFEE;
  if (strName.Contains("PTR")) fSystem = kInsidePTR;
  if (strName.Contains("GTU")) fSystem = kInsideGTU;

  // get the attributes of the element
  TXMLAttr *attr;
  TIter next(attributes);
  while ((attr = (TXMLAttr*) next())) {
    TString attribName = attr->GetName();
    if (attribName.Contains("id") && strName.Contains("DCS")) {
      dcsTitle = name;
      dcsId = atoi(attr->GetValue());
    }
    if (attribName.Contains("roc") && strName.Contains("ack")) {
      if (atoi(attr->GetValue()) != fDCSFEEObj->GetDCSid())
	fDCSFEEObj->SetStatusBit(3); // consistence check
    }
    if (attribName.Contains("sm") && strName.Contains("DCS")) {
      fCurrentSM = atoi(attr->GetValue()); // only for GTU/PTR
    }
    if (attribName.Contains("id") && strName.Contains("STACK")) {
      fCurrentStack = atoi(attr->GetValue()); // only for GTU/PTR
    }
  }

  // if there is a new DCS element put it in the correct array
  if (strName.Contains("DCS")) {
    if (fSystem == kInsideFEE) {
      fDCSFEEObj = new AliTRDCalDCSFEE(name,dcsTitle);
      fDCSFEEObj->SetDCSid(dcsId);
    }
    if (fSystem == kInsidePTR) {
      fDCSPTRObj = new AliTRDCalDCSPTR(name,dcsTitle);
      fDCSPTRObj->SetDCSid(dcsId);
    }
    if (fSystem == kInsideGTU) {
      fDCSGTUObj = new AliTRDCalDCSGTU(name,dcsTitle);
      fDCSGTUObj->SetDCSid(dcsId);
    }
  }
}

//_____________________________________________________________________________
void AliTRDSaxHandler::OnEndElement(const char *name)
{
  // do everything that needs to be done when an end tag of an element is found
  TString strName = name;
  
  // if done with this DCS board, put it in the correct array
  // no check for </ack> necessary since this check is done during XML validation
  if (strName.Contains("DCS")) {
    if (fSystem == kInsideFEE) {
      Int_t detID = 0;
      if (fDCSFEEObj->GetStatusBit() == 0) {
	// if there were no errors (StatusBit==0) the following should match
        detID = fDCSFEEObj->GetDCSid();
        AliTRDgeometry aliGeo;
        Int_t calDetID = aliGeo.GetDetector(fDCSFEEObj->GetLayer(),
	  				    fDCSFEEObj->GetStack(),
					    fDCSFEEObj->GetSM());
	if (detID != calDetID) fDCSFEEObj->SetStatusBit(4);
      } else {
	// if the dcs board didn't properly respond, don't compare
        detID = fDCSFEEObj->GetDCSid();
      }
      fFEEArr->AddAt(fDCSFEEObj,detID);
    }
    if (fSystem == kInsidePTR) {
      fPTRArr->AddAt(fDCSPTRObj,fNDCSPTR);
      fNDCSPTR++;
    }
    if (fSystem == kInsideGTU) {
      fGTUArr->AddAt(fDCSGTUObj,fNDCSGTU);
      fNDCSGTU++;
    }
    fCurrentSM = 99; // 99 for no SM set
    fDCSFEEObj = 0;  // just to be sure
    return;
  }

  // done with this stack?
  if (strName.Contains("STACK")) {
    fCurrentStack = 99; // 99 for no stack set
  }

  // store informations of the FEE DCS-Board
  if (fSystem == kInsideFEE) {
    if (strName.Contains("DNR")) fDCSFEEObj->SetStatusBit(fContent.Atoi());
    if (strName.Contains("CFGNME")) fDCSFEEObj->SetConfigName(fContent);
    if (strName.Contains("CFGTAG")) fDCSFEEObj->SetConfigTag(fContent.Atoi());
    if (strName.Contains("CFGVRSN")) fDCSFEEObj->SetConfigVersion(fContent);
    if (strName.Contains("NTB")) fDCSFEEObj->SetNumberOfTimeBins(fContent.Atoi());
    if (strName.Contains("SM-ID")) fDCSFEEObj->SetSM(fContent.Atoi());
    if (strName.Contains("STACK-ID")) fDCSFEEObj->SetStack(fContent.Atoi());
    if (strName.Contains("LAYER-ID")) fDCSFEEObj->SetLayer(fContent.Atoi());
    if (strName.Contains("SINGLEHITTHRES")) fDCSFEEObj->SetSingleHitThres(fContent.Atoi());
    if (strName.Contains("THRPADCLSTHRS")) fDCSFEEObj->SetThreePadClustThres(fContent.Atoi());
    if (strName.Contains("SELNOZS")) fDCSFEEObj->SetSelectiveNoZS(fContent.Atoi());
    if (strName.Contains("FASTSTATNOISE")) fDCSFEEObj->SetFastStatNoise(fContent.Atoi());
    if (strName.Contains("FILTWEIGHT")) fDCSFEEObj->SetTCFilterWeight(fContent.Atoi());
    if (strName.Contains("FILTSHRTDCYPRM")) fDCSFEEObj->SetTCFilterShortDecPar(fContent.Atoi());
    if (strName.Contains("FILTLNGDCYPRM")) fDCSFEEObj->SetTCFilterLongDecPar(fContent.Atoi());
    if (strName.Contains("FLTR")) fDCSFEEObj->SetFilterType(fContent);
    if (strName.Contains("READOUTPAR")) fDCSFEEObj->SetReadoutParam(fContent);
    if (strName.Contains("TESTPATTERN")) fDCSFEEObj->SetTestPattern(fContent);
    if (strName.Contains("TRCKLTMODE")) fDCSFEEObj->SetTrackletMode(fContent);
    if (strName.Contains("TRCKLTDEF")) fDCSFEEObj->SetTrackletDef(fContent);
    if (strName.Contains("TRIGGERSETUP")) fDCSFEEObj->SetTriggerSetup(fContent);
    if (strName.Contains("ADDOPTIONS")) fDCSFEEObj->SetAddOptions(fContent);
  }

  
  // store pretrigger informations
  if (fSystem == kInsidePTR) {
    // no informations available yet
  }
  // store GTU informations
  if (fSystem == kInsideGTU) {
    if (strName.Contains("SMMASK"))
      fHandlerStatus = fDCSGTUObj->SetSMMask(fContent);
    if (strName.Contains("LINKMASK")) 
      fHandlerStatus = fDCSGTUObj->SetLinkMask(fCurrentSM, fCurrentStack, fContent);
    if (strName.Contains("STMASK"))
      fDCSGTUObj->SetStackMaskBit(fCurrentSM, fCurrentStack, fContent.Atoi());
  }
}

//_____________________________________________________________________________
void AliTRDSaxHandler::OnCharacters(const char *characters)
{
  // copy the the text content of an XML element
  fContent = characters;
}

//_____________________________________________________________________________
void AliTRDSaxHandler::OnComment(const char* /*text*/)
{
  // comments within the XML file are ignored
}

//_____________________________________________________________________________
void AliTRDSaxHandler::OnWarning(const char *text)
{
  // process warnings here
  AliInfo(Form("Warning: %s",text));
}

//_____________________________________________________________________________
void AliTRDSaxHandler::OnError(const char *text)
{
  // process errors here
  AliError(Form("Error: %s",text));
}

//_____________________________________________________________________________
void AliTRDSaxHandler::OnFatalError(const char *text)
{
  // process fatal errors here
  AliError(Form("Fatal error: %s",text)); // use AliFatal?
}

//_____________________________________________________________________________
void AliTRDSaxHandler::OnCdataBlock(const char* /*text*/, Int_t /*len*/)
{
  // process character data blocks here
  // not implemented and should not be used here
}

