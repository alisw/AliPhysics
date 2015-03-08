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
//  The SAX XML file handler used in the TRD preprocessor                 //
//                                                                        //
//  Authors:                                                              //
//    Frederick Kramer (kramer@ikf.uni-frankfurt.de)                      //
//    Thomas Bird      (thomas@thomasbird.com)                            //
//                                                                        //
////////////////////////////////////////////////////////////////////////////


#include <TList.h>
#include <TMath.h>
#include "AliLog.h"
#include "AliTRDCalDCSGTUTgu.h"
#include "AliTRDCalDCSPTR.h"

#include <TXMLAttr.h>
#include <TObjArray.h>
#include "AliTRDSaxHandler.h"
#include "AliTRDgeometry.h"
#include "AliTRDcalibDB.h"
#include "AliTRDCalDCSv2.h"
#include "AliTRDCalDCSFEEv2.h"
#include "AliTRDCalDCSGTU.h"

ClassImp(AliTRDSaxHandler)
  
//_____________________________________________________________________________
AliTRDSaxHandler::AliTRDSaxHandler()
  :TObject()
  ,fHandlerStatus(0)
  ,fNDCSPTR(0)
  ,fNDCSGTU(0)
  ,fFEEArr(new TObjArray(540))
  ,fPTRArr(new TObjArray(6))
  ,fSystem(0)
  ,fInsideRstate(0)
  ,fCurrentSM(0)
  ,fCurrentStack(0)
  ,fCurrentROB(-1)
  ,fCurrentMCM(-1)
  ,fCurrentADC(-1)
  ,fContent(0)
  ,fDCSFEEObj(0)
  ,fDCSPTRObj(0)
  ,fDCSGTUObj(0)
  ,fCalDCSObj(new AliTRDCalDCSv2())
  ,fLevel1Tag(-2)
  ,fLevel2Tag(-2)
  ,fInsideBoardInfo(false)
  ,fTmu(0)
  ,fCtpOpc(0)
  ,fSegment(0)
  ,fBoardInfo(0)
{
  // AliTRDSaxHandler default constructor
  fFEEArr->SetOwner();
  fPTRArr->SetOwner();
}

//_____________________________________________________________________________
AliTRDSaxHandler::AliTRDSaxHandler(const AliTRDSaxHandler &sh)
  :TObject(sh)
  ,fHandlerStatus(0)
  ,fNDCSPTR(0)
  ,fNDCSGTU(0)
  ,fFEEArr(0)
  ,fPTRArr(0)
  ,fSystem(0)
  ,fInsideRstate(0)
  ,fCurrentSM(0)
  ,fCurrentStack(0)
  ,fCurrentROB(-1)
  ,fCurrentMCM(-1)
  ,fCurrentADC(-1)
  ,fContent(0)
  ,fDCSFEEObj(0)
  ,fDCSPTRObj(0)
  ,fDCSGTUObj(0)
  ,fCalDCSObj(0)
  ,fLevel1Tag(-2)
  ,fLevel2Tag(-2)
  ,fInsideBoardInfo(false)
  ,fTmu(0)
  ,fCtpOpc(0)
  ,fSegment(0)
  ,fBoardInfo(0)
{
  // AliTRDSaxHandler copy constructor
}

//_____________________________________________________________________________
AliTRDSaxHandler &AliTRDSaxHandler::operator=(const AliTRDSaxHandler &sh)
{
  // Assignment operator
  if (&sh == this) return *this;
  new (this) AliTRDSaxHandler(sh);
  return *this;
}

//_____________________________________________________________________________
AliTRDSaxHandler::~AliTRDSaxHandler()
{
  // AliTRDSaxHandler destructor
  if (fFEEArr) {
    delete fFEEArr;
    fFEEArr    = 0x0;
  }
  if (fPTRArr) {
    delete fPTRArr;
    fPTRArr    = 0x0;
  }
  if (fDCSGTUObj) {
    delete fDCSGTUObj;
    fDCSGTUObj    = 0x0;
  }
  if (fCalDCSObj) {
    delete fCalDCSObj;
    fCalDCSObj = 0x0;
  }
}

//_____________________________________________________________________________
AliTRDCalDCSv2* AliTRDSaxHandler::GetCalDCSObj()
{
  // put the arrays in the global calibration object and return this

  fCalDCSObj->SetFEEArr(fFEEArr);
  fCalDCSObj->SetPTRArr(fPTRArr);
  fCalDCSObj->SetGTUObj(fDCSGTUObj);
  return fCalDCSObj;
}

//_____________________________________________________________________________
void AliTRDSaxHandler::ParseConfigName(TString cfgname) const
{
  // Evaluate the config name and set the individual parameters

  TString cfg = "", par = "", pars = "";
  Int_t nPar = AliTRDcalibDB::GetNumberOfParsDCS(cfgname);
  if (nPar == 0) return;

  for (Int_t i=1; i<=nPar; i++) {
    // Get the configuration parameter
    AliTRDcalibDB::GetDCSConfigParOption(cfgname, i, 0, cfg);

    // Set Parameters accordingly
    if (i == AliTRDcalibDB::kFltrSet) fDCSFEEObj->SetFilterType(cfg);
    if (i == AliTRDcalibDB::kTrigSet) fDCSFEEObj->SetTriggerSetup(cfg);
    if (i == AliTRDcalibDB::kAddOpti) fDCSFEEObj->SetAddOptions(cfg);
    if (i == AliTRDcalibDB::kTimebin) fDCSFEEObj->SetNumberOfTimeBins(AliTRDcalibDB::ExtractTimeBinsFromString(cfg));
    if (i == AliTRDcalibDB::kReadout) fDCSFEEObj->SetReadoutParam(cfg);
    if (i == AliTRDcalibDB::kTrkMode) fDCSFEEObj->SetTrackletMode(cfg);

    // Set options of parameters accordingly
    Int_t nOpt = AliTRDcalibDB::GetNumberOfOptsDCS(cfgname, i);
    if (nOpt == 0) continue;

    for (Int_t j=1; j<=nOpt; j++) {
      // Get the parameter option
      AliTRDcalibDB::GetDCSConfigParOption(cfgname, i, j, par);

      if (i == AliTRDcalibDB::kReadout) {
	if (par.EqualTo("stat")) fDCSFEEObj->SetFastStatNoise(1);
      }
      if (i == AliTRDcalibDB::kTrkMode) {
	if ((j > 1) && (par.Length() != 0)) pars += "-";
	pars += par;
      }
      // SetTCFilterWeight, SetTCFilterShortDecPar, SetTCFilterLongDecPar might be filled here, too
      // SetSingleHitThres, SetThreePadClustThres, SetSelectiveNoZS, SetTestPattern might be filled here, too
    }

    fDCSFEEObj->SetTrackletDef(pars);
  }
}

//_____________________________________________________________________________
void AliTRDSaxHandler::OnStartDocument() const
{
  // if something should happen right at the beginning of the
  // XML document, this must happen here
}

//_____________________________________________________________________________
void AliTRDSaxHandler::OnEndDocument() const
{
  // if something should happen at the end of the XML document
  // this must be done here
}

//_____________________________________________________________________________
bool AliTRDSaxHandler::CompareString(TString str, const char *str2)
{
  // compre strings, ignoring case
  // returns true if they are the same, else false
  return !(bool)str.CompareTo(str2,str.kIgnoreCase);
}


//_____________________________________________________________________________
void AliTRDSaxHandler::OnStartElement(const char *name, const TList *attributes)
{
  // when a new XML element is found, it is processed here
  fContent    = ""; // Technically <p> This <em>is</em> ok but would be a problem here</p>
  Int_t dcsId = 0;
  TString tagName  = name;
  TString dcsTitle = "";

  // set the current system if necessary
  if (CompareString(tagName, "FEE")) fSystem = kInsideFEE;
  if (CompareString(tagName, "PTR")) fSystem = kInsidePTR;
  if (CompareString(tagName, "GTU")) {
    fSystem = kInsideGTU;
    fDCSGTUObj = new AliTRDCalDCSGTU(tagName,tagName);
  }

  if (fSystem == kInsideGTU) {
    if (CompareString(tagName, "tgu")) fLevel1Tag = kInsideTgu;
    if (CompareString(tagName, "board_info")) {
      fInsideBoardInfo = true;
      fBoardInfo = new AliTRDCalDCSGTUBoardInfo(tagName,tagName);
    }
    if (CompareString(tagName(0,tagName.Length()-3), "segment")) { 
      fSegment = new AliTRDCalDCSGTUSegment(tagName,tagName);
      fSegment->SetId(TString(tagName(8,2)).Atoi());
      fLevel1Tag = kInsideSegment;
    }
    if (fLevel1Tag == kInsideTgu) {
      if (CompareString(tagName, "ctp_opc"))   fCtpOpc = new AliTRDCalDCSGTUCtpOpc(tagName,tagName);
    } else if (fLevel1Tag == kInsideSegment) {
      if (CompareString(tagName, "smu")) {
	fLevel2Tag = kInsideSmu;
      }
      if (CompareString(tagName(0,3), "tmu")) {
	fTmu = new AliTRDCalDCSGTUTmu(tagName,tagName);
	fTmu->SetId(TString(tagName(4,2)).Atoi());
	fLevel2Tag = kInsideTmu;
      }
    }
  } else if (fSystem == kInsideFEE) {
    if (CompareString(tagName, "gaintbl")) fLevel1Tag = kInsideGainTable;
  }

  // set if we are inside rstate 
  // (in principle not necessary - just to be more safe against stupid tags)
  if (CompareString(tagName, "rstate")) fInsideRstate = 1;

  // get the attributes of the element
  TXMLAttr *attr;
  TIter next(attributes);
  while ((attr = (TXMLAttr*) next())) {
    TString attribName = attr->GetName();
    TString attribValue = attr->GetValue();
    if (fSystem == kInsideFEE && fLevel1Tag == kInsideNone) {
      if (CompareString(attribName, "id") && CompareString(tagName, "DCS")) {
	dcsTitle = name;
	dcsId = atoi(attr->GetValue());
      }
      if (CompareString(attribName, "roc") && CompareString(tagName, "ack")) {
	if (attribValue.Atoi() != fDCSFEEObj->GetDCSid())
	  fDCSFEEObj->SetStatusBit(4); // consistency check
      }
      if (CompareString(attribName, "rob") && CompareString(tagName, "ro-board") && (fInsideRstate == 1)) {
	fCurrentROB = attribValue.Atoi();
      }
      if (CompareString(attribName, "mcm") && CompareString(tagName, "m") && (fInsideRstate == 1)) {
	fCurrentMCM = attribValue.Atoi();
      }
      if (CompareString(attribName, "sm") && CompareString(tagName, "DCS")) {
	fCurrentSM = attribValue.Atoi(); // only for GTU/PTR
      }
      if (CompareString(attribName, "id") && CompareString(tagName, "STACK")) {// hmmmm not exist?
	fCurrentStack = attribValue.Atoi(); // only for GTU/PTR
      }
    } else if (fSystem == kInsideFEE && fLevel1Tag == kInsideGainTable) {
      if (CompareString(tagName, "roc") && CompareString(attribName, "type"))    fDCSFEEObj->SetGainTableRocType(attribValue);
      if (CompareString(tagName, "roc") && CompareString(attribName, "serial"))  fDCSFEEObj->SetGainTableRocSerial(attribValue.Atoi());
      if (CompareString(tagName, "mcm") && CompareString(attribName, "rob"))     fCurrentROB = attribValue.Atoi();
      if (CompareString(tagName, "mcm") && CompareString(attribName, "pos"))     fCurrentMCM = attribValue.Atoi();
      if (CompareString(tagName, "adc") && CompareString(attribName, "id"))      fCurrentADC = attribValue.Atoi();
      
    } else if (fSystem == kInsideGTU && fLevel1Tag == kInsideNone) {
      if (CompareString(tagName, "publisher")) {
	if (CompareString(attribName, "at"))         fDCSGTUObj->SetSORFlag(attribValue.Atoi());
	if (CompareString(attribName, "serial"))     fDCSGTUObj->SetSerial(attribValue.Atoi());
	if (CompareString(attribName, "runnr"))      fDCSGTUObj->SetRunNumber(attribValue.Atoi());
      }
    } else if (fSystem == kInsideGTU && fLevel1Tag == kInsideTgu) {
      if (CompareString(tagName, "from")) {
	if (CompareString(attribName, "at"))         fDCSGTUObj->GetTgu()->SetFromSORFlag(attribValue.Atoi());
	if (CompareString(attribName, "runnr"))      fDCSGTUObj->GetTgu()->SetFromRunNumber(attribValue.Atoi());
	if (CompareString(attribName, "child"))      fDCSGTUObj->GetTgu()->SetFromChild(attribValue.Atoi());
      }
      if (CompareString(tagName, "segmentmask") && CompareString(attribName, "value"))  fDCSGTUObj->GetTgu()->SetSegmentMask(attribValue);
      if (CompareString(tagName, "busymask") && CompareString(attribName, "value"))     fDCSGTUObj->GetTgu()->SetBusyMask(attribValue);
      if (CompareString(tagName, "contribmask") && CompareString(attribName, "value"))  fDCSGTUObj->GetTgu()->SetContribMask(attribValue);
      
      if (CompareString(tagName, "ctp_opc") && CompareString(attribName, "id"))         fCtpOpc->SetId(attribValue.Atoi());
      if (CompareString(tagName, "ctp_opc") && CompareString(attribName, "opcode"))     fCtpOpc->SetOpcode(attribValue.Atoi());
      if (CompareString(tagName, "ctp_opc") && CompareString(attribName, "direction"))  fCtpOpc->SetDirection(attribValue.Atoi());
      if (CompareString(tagName, "ctp_opc") && CompareString(attribName, "inverted"))   fCtpOpc->SetInverted(attribValue.Atoi());
      if (CompareString(tagName, "ctp_opc") && CompareString(attribName, "delay"))      fCtpOpc->SetDelay(attribValue.Atoi());
      if (CompareString(tagName, "ctp_opc") && CompareString(attribName, "connected"))  fCtpOpc->SetConnected(attribValue.Atoi());
      
    } else if (fSystem == kInsideGTU && fLevel1Tag == kInsideSegment) {
      if (CompareString(tagName, "from")) {
	if (CompareString(attribName, "at"))         fSegment->SetFromSORFlag(attribValue.Atoi());
	if (CompareString(attribName, "runnr"))      fSegment->SetFromRunNumber(attribValue.Atoi());
	if (CompareString(attribName, "child"))      fSegment->SetFromChild(attribValue.Atoi());
      }
      if (fLevel2Tag == kInsideSmu) {
	if (CompareString(tagName, "stackmask") && CompareString(attribName, "value"))     fSegment->SetSmuStackMask(attribValue);
	if (CompareString(tagName, "tracklets") && CompareString(attribName, "send"))      fSegment->SetSmuTracklets(attribValue.Atoi());
	if (CompareString(tagName, "tracks") && CompareString(attribName, "send"))         fSegment->SetSmuTracks(attribValue.Atoi());
	if (CompareString(tagName, "idelay") && CompareString(attribName, "value"))        fSegment->SetSmuIdelay(attribValue.Atoi());
	if (CompareString(tagName, "ttc_emulator") && CompareString(attribName, "enable")) fSegment->SetSmuTtcEmulatorEnable(attribValue.Atoi());
	
	if (CompareString(tagName, "trigger_window") && CompareString(attribName, "l1_low"))  
	  fSegment->SetSmuTriggerWindowL1Low(attribValue.Atoi());
	if (CompareString(tagName, "trigger_window") && CompareString(attribName, "l1_high"))  
	  fSegment->SetSmuTriggerWindowL1High(attribValue.Atoi());
	if (CompareString(tagName, "trigger_window") && CompareString(attribName, "l2_low"))  
	  fSegment->SetSmuTriggerWindowL2Low(attribValue.Atoi());
	if (CompareString(tagName, "trigger_window") && CompareString(attribName, "l2_high"))  
	  fSegment->SetSmuTriggerWindowL2High(attribValue.Atoi());
	
      } else if (fLevel2Tag == kInsideTmu) {
	if (CompareString(tagName, "linkmask") && CompareString(attribName, "value"))      fTmu->SetLinkMask(attribValue);
	if (CompareString(tagName, "pattern_generator") && CompareString(attribName, "enable")) 
	  fTmu->SetPatternGeneratorEnable(attribValue.Atoi());
	if (CompareString(tagName, "pattern_generator") && CompareString(attribName, "datawords")) 
	  fTmu->SetPatternGeneratorDataWords(attribValue.Atoi());
	if (CompareString(tagName, "pattern_generator") && CompareString(attribName, "trackletwords")) 
	  fTmu->SetPatternGeneratorTrackletWords(attribValue.Atoi());
      }
    }
    
    if (fInsideBoardInfo) {
      if (CompareString(tagName, "board_info") && CompareString(attribName, "board_id"))    fBoardInfo->SetId(attribValue);
      if (CompareString(tagName, "board_info") && CompareString(attribName, "design_type")) fBoardInfo->SetType(attribValue.Atoi());
      if (CompareString(tagName, "board_info") && CompareString(attribName, "pci_ga"))      fBoardInfo->SetPciGa(attribValue.Atoi());
      if (CompareString(tagName, "hardware") && CompareString(attribName, "date"))          fBoardInfo->SetHwDate(attribValue);
      if (CompareString(tagName, "hardware") && CompareString(attribName, "rev"))           fBoardInfo->SetHwRev(attribValue.Atoi());
      if (CompareString(tagName, "hardware") && CompareString(attribName, "clean"))         fBoardInfo->SetHwClean(attribValue.Atoi());
      if (CompareString(tagName, "software") && CompareString(attribName, "date"))          fBoardInfo->SetSwDate(attribValue);
      if (CompareString(tagName, "software") && CompareString(attribName, "rev"))           fBoardInfo->SetSwRev(attribValue.Atoi());
      if (CompareString(tagName, "software") && CompareString(attribName, "clean"))         fBoardInfo->SetSwClean(attribValue.Atoi());
    }
  }

  // if there is a new DCS element put it in the correct array
  if (CompareString(tagName, "DCS")) {
    if (fSystem == kInsideFEE) {
      fDCSFEEObj = new AliTRDCalDCSFEEv2();
      fDCSFEEObj->SetDCSid(dcsId);
    }
    if (fSystem == kInsidePTR) {
//       fDCSPTRObj = new AliTRDCalDCSPTR(name,dcsTitle);
//       fDCSPTRObj->SetDCSid(dcsId);
    }
    if (fSystem == kInsideGTU) {
//       fDCSGTUObj = new AliTRDCalDCSGTU(name,dcsTitle);
//       fDCSGTUObj->SetDCSid(dcsId);
    }
  }
}

//_____________________________________________________________________________
void AliTRDSaxHandler::OnEndElement(const char *name)
{
  // do everything that needs to be done when an end tag of an element is found
  TString tagName = name;
  
  // if done with this DCS board, put it in the correct array
  // no check for </ack> necessary since this check is done during XML validation
  if (CompareString(tagName, "DCS")) {
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
//     if (fSystem == kInsideGTU) {
//       fGTUArr->AddAt(fDCSGTUObj,fNDCSGTU);
//       fNDCSGTU++;
//     }
    fCurrentSM = 99; // 99 for no SM set
    fDCSFEEObj = 0;  // just to be sure
    return;
  }

  // done with this stack? 
  if (CompareString(tagName, "STACK")) {// TODO: errrrm ???? always 99?
    fCurrentStack = 99; // 99 for no stack set
  }

  // outside of rstate again?
  if (CompareString(tagName, "rstate")) {
    fInsideRstate = 0;
    fCurrentROB   = -1;
    fCurrentMCM   = -1;
  }
  if (CompareString(tagName, "ro-board")) fCurrentROB = -1;
  
  // store informations of the FEE DCS-Board
  if (fSystem == kInsideFEE) {
    if (CompareString(tagName, "DNR"))            fDCSFEEObj->SetStatusBit(fContent.Atoi());
    if (CompareString(tagName, "CFGTAG"))         fDCSFEEObj->SetConfigTag(fContent.Atoi());
    if (CompareString(tagName, "CFGVRSN"))        fDCSFEEObj->SetConfigVersion(fContent);
    if (CompareString(tagName, "SM-ID"))          fDCSFEEObj->SetSM(fContent.Atoi());
    if (CompareString(tagName, "STACK-ID"))       fDCSFEEObj->SetStack(fContent.Atoi());
    if (CompareString(tagName, "LAYER-ID"))       fDCSFEEObj->SetLayer(fContent.Atoi());
    if (CompareString(tagName, "CFGNME")) {
      fDCSFEEObj->SetConfigName(fContent);
      ParseConfigName(fContent);
    }
    if (CompareString(tagName, "gaintbl")) {
      fLevel1Tag = kInsideNone;
      fCurrentROB = -1;
      fCurrentMCM = -1;
      fCurrentADC = -1;
    }
    if (fLevel1Tag == kInsideGainTable) {
      if (CompareString(tagName, "name"))   fDCSFEEObj->SetGainTableName(fContent);
      if (CompareString(tagName, "desc"))   fDCSFEEObj->SetGainTableDesc(fContent);
      if (fCurrentROB>=0 && fCurrentMCM>=0) {
	if (CompareString(tagName, "adcdac")) fDCSFEEObj->SetGainTableAdcdac(fCurrentROB, fCurrentMCM, fContent.Atoi());
	if (fCurrentADC>=0) {
	  if (CompareString(tagName, "fgfn"))   fDCSFEEObj->SetGainTableFgfn(fCurrentROB, fCurrentMCM, fCurrentADC, fContent.Atoi());
	  if (CompareString(tagName, "fgan"))   fDCSFEEObj->SetGainTableFgan(fCurrentROB, fCurrentMCM, fCurrentADC, fContent.Atoi());
	}
      }
    }
    if (fInsideRstate == 1) {
      if (fCurrentROB>=0 && fCurrentMCM>=0) {
	if (CompareString(tagName, "gsm")) fDCSFEEObj->SetMCMGlobalState(fCurrentROB, fCurrentMCM, fContent.Atoi());
        if (CompareString(tagName, "ni")) fDCSFEEObj->SetMCMStateNI(fCurrentROB, fCurrentMCM, fContent.Atoi());
	if (CompareString(tagName, "ev")) fDCSFEEObj->SetMCMEventCnt(fCurrentROB, fCurrentMCM, fContent.Atoi());
	if (CompareString(tagName, "ptrg")) fDCSFEEObj->SetMCMPtCnt(fCurrentROB, fCurrentMCM, fContent.Atoi());
      }
    }
  }

  if (fSystem == kInsideGTU) {
//     if (CompareString(tagName, "run")) { 
//       fDCSGTUObj->SetSORFlag(TString(fContent(fContent.Length()-1,1)).Atoi());
//       fDCSGTUObj->SetRunNumber(TString(fContent(0,fContent.Length()-2)).Atoi());
//     }
//     if (CompareString(tagName, "serial"))         fDCSGTUObj->SetSerial(fContent.Atoi());
    if (CompareString(tagName, "board_info")) {
      fInsideBoardInfo = false;
      if (fLevel1Tag == kInsideTgu)                                  fDCSGTUObj->GetTgu()->SetBoardInfo(fBoardInfo);
      if (fLevel1Tag == kInsideSegment && fLevel2Tag == kInsideSmu)  fSegment->SetSmuBoardInfo(fBoardInfo);
      if (fLevel1Tag == kInsideSegment && fLevel2Tag == kInsideTmu)  fTmu->SetBoardInfo(fBoardInfo);
    }
    if (CompareString(tagName, "dnr"))            fDCSGTUObj->SetDNR(fContent.Atoi());
    if (CompareString(tagName, "tgu"))            fLevel1Tag = kInsideNone;
    if (CompareString(tagName(0,tagName.Length()-3), "segment")) { 
      fDCSGTUObj->GetSegmentArray()->Add(fSegment);
      fLevel1Tag = kInsideNone;
    }
    if (fLevel1Tag == kInsideTgu) {
      if (CompareString(tagName, "ctp_opc"))        fDCSGTUObj->GetTgu()->GetCtpOpcArray()->Add(fCtpOpc);
    } else if (fLevel1Tag == kInsideSegment) {
      if (CompareString(tagName, "smu"))          fLevel2Tag = kInsideNone;
      if (CompareString(tagName(0,3), "tmu")) {
	fSegment->GetTmuArray()->Add(fTmu);
	fLevel2Tag = kInsideNone;
      }
    }
  }

  
  // store pretrigger informations
  if (fSystem == kInsidePTR) {
    // no informations available yet
  }
//   // store GTU informations
//   if (fSystem == kInsideGTU) {
//     if (CompareString(tagName, "SMMASK"))
//       fHandlerStatus = fDCSGTUObj->SetSMMask(fContent);
//     if (CompareString(tagName, "LINKMASK")) 
//       fHandlerStatus = fDCSGTUObj->SetLinkMask(fCurrentSM, fCurrentStack, fContent);
//     if (CompareString(tagName, "STMASK"))
//       fDCSGTUObj->SetStackMaskBit(fCurrentSM, fCurrentStack, fContent.Atoi());
//   }
}

//_____________________________________________________________________________
void AliTRDSaxHandler::OnCharacters(const char *characters)
{
  // copy the the text content of an XML element
  fContent = characters;
}

//_____________________________________________________________________________
void AliTRDSaxHandler::OnComment(const char* /*text*/) const
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
void AliTRDSaxHandler::OnCdataBlock(const char* /*text*/, Int_t /*len*/) const
{
  // process character data blocks here
  // not implemented and should not be used here
}

