//-*- Mode: C++ -*-
// $Id: AliHLTMultiplicityCorrelationsComponent.cxx $
/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Jochen Thaeder <jochen@thaeder.de>                    *
 *                  for The ALICE HLT Project.                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/** @file    AliHLTMultiplicityCorrelationsComponent.cxx
    @author  Jochen Thaeder <jochen@thaeder.de>
    @brief   Component for Multiplicty Correlations
*/

#if __GNUC__>= 3
using namespace std;
#endif

#include "TMap.h"
#include "TSystem.h"
#include "TTimeStamp.h"
#include "TObjString.h"
#include "AliESDVZERO.h"
#include "AliESDtrackCuts.h"
#include "AliHLTMultiplicityCorrelations.h"
#include <AliHLTDAQ.h>

#include "AliHLTErrorGuard.h"
#include "AliHLTDataTypes.h"
#include "AliHLTMultiplicityCorrelationsComponent.h"
#include "AliHLTITSClusterDataFormat.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTMultiplicityCorrelationsComponent)

/*
 * ---------------------------------------------------------------------------------
 *                            Constructor / Destructor
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
AliHLTMultiplicityCorrelationsComponent::AliHLTMultiplicityCorrelationsComponent() :
  AliHLTProcessor(),
  fESDTrackCuts(NULL),  
  fCorrObj(NULL),
  fUID(0) {
  // an example component which implements the ALICE HLT processor
  // interface and does some analysis on the input raw data
  //
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
  //
  // NOTE: all helper classes should be instantiated in DoInit()
}

// #################################################################################
AliHLTMultiplicityCorrelationsComponent::~AliHLTMultiplicityCorrelationsComponent() {
  // see header file for class documentation
}

/*
 * ---------------------------------------------------------------------------------
 * Public functions to implement AliHLTComponent's interface.
 * These functions are required for the registration process
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
const Char_t* AliHLTMultiplicityCorrelationsComponent::GetComponentID() { 
  // see header file for class documentation
  return "MultiplicityCorrelations";
}

// #################################################################################
void AliHLTMultiplicityCorrelationsComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list) {
  // see header file for class documentation
  list.push_back(kAliHLTDataTypeESDObject|kAliHLTDataOriginAny);
  list.push_back(kAliHLTDataTypeClusters|kAliHLTDataOriginITSSPD);
  list.push_back(kAliHLTDataTypeESDContent|kAliHLTDataOriginVZERO);
}

// #################################################################################
AliHLTComponentDataType AliHLTMultiplicityCorrelationsComponent::GetOutputDataType() {
  // see header file for class documentation
  return kAliHLTDataTypeTObject|kAliHLTDataOriginHLT;
}

// #################################################################################
void AliHLTMultiplicityCorrelationsComponent::GetOutputDataSize( ULong_t& constBase, Double_t& inputMultiplier ) {
  // see header file for class documentation
  constBase = 1000;
  inputMultiplier = 0.5;
}

// #################################################################################
void AliHLTMultiplicityCorrelationsComponent::GetOCDBObjectDescription( TMap* const targetMap) {
  // see header file for class documentation

  if (!targetMap) return;
  targetMap->Add(new TObjString("HLT/ConfigGlobal/MultiplicityCorrelations"),
		 new TObjString("configuration object"));

  return;
}

// #################################################################################
AliHLTComponent* AliHLTMultiplicityCorrelationsComponent::Spawn() {
  // see header file for class documentation
  return new AliHLTMultiplicityCorrelationsComponent;
}

/*
 * ---------------------------------------------------------------------------------
 * Protected functions to implement AliHLTComponent's interface.
 * These functions provide initialization as well as the actual processing
 * capabilities of the component. 
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
Int_t AliHLTMultiplicityCorrelationsComponent::DoInit( Int_t argc, const Char_t** argv ) {
  // see header file for class documentation

  Int_t iResult=0;

  // -- Initialize members
  // -----------------------
  do {

    fCorrObj = new AliHLTMultiplicityCorrelations;
    if (!fCorrObj) {
      iResult=-ENOMEM;
      break;
    }

    fESDTrackCuts = new AliESDtrackCuts("AliESDtrackCuts","HLT");
    if (!fESDTrackCuts) {
      iResult=-ENOMEM;
      break;
    }

    fUID = 0;
    
    // implement further initialization
  } while (0);

  if (iResult<0) {
    // implement cleanup

    if (fCorrObj) 
      delete fCorrObj;
    fCorrObj = NULL;

    if (fESDTrackCuts) 
      delete fESDTrackCuts;
    fESDTrackCuts = NULL;
  }

  if (iResult>=0) {
    SetDefaultConfiguration();

    // -- Read configuration object : HLT/ConfigGlobal/MultiplicityCorrelations
    TString cdbPath="HLT/ConfigGlobal/";
    cdbPath+=GetComponentID();
    iResult=ConfigureFromCDBTObjString(cdbPath);
    
    // -- Read the component arguments
    if (iResult>=0) {
      iResult=ConfigureFromArgumentString(argc, argv);
    }
  }

  if (iResult>=0) {
    HLTInfo("ESD track cuts : %s",fESDTrackCuts->GetTitle() );

    fCorrObj->SetESDTrackCuts(fESDTrackCuts);
    fCorrObj->Initialize();
  }

  return iResult;
}

// #################################################################################
void AliHLTMultiplicityCorrelationsComponent::SetDefaultConfiguration() {
  // see header file for class documentation

  if (fESDTrackCuts) {
    fESDTrackCuts->SetEtaRange(-0.9,0.9);
    fESDTrackCuts->SetPtRange(0.2,200);
    fESDTrackCuts->SetMinNClustersTPC(80);
    
    fESDTrackCuts->SetDCAToVertex2D(kFALSE);
    fESDTrackCuts->SetRequireSigmaToVertex(kFALSE);
    
    fESDTrackCuts->SetMaxDCAToVertexXY(3.0);
    fESDTrackCuts->SetMaxDCAToVertexZ(3.0);
    
    fESDTrackCuts->SetRequireTPCRefit(kFALSE);
    fESDTrackCuts->SetRequireITSRefit(kFALSE);
  }
 
  return;
}

// #################################################################################
Int_t AliHLTMultiplicityCorrelationsComponent::ScanConfigurationArgument(Int_t argc, const Char_t** argv) {
  // Scan configuration arguments
  // Return the number of processed arguments
  //        -EPROTO if argument format error (e.g. number expected but not found)
  //
  // The AliHLTComponent base class implements a parsing loop for argument strings and
  // arrays of strings which is invoked by ConfigureFromArgumentString/ConfigureFromCDBTObjString
  // The component needs to implement ScanConfigurationArgument in order to decode the arguments.

  if (argc<=0) return 0;
  Int_t ii =0;
  TString argument=argv[ii];
  
  if (argument.IsNull()) return 0;

  if( !fESDTrackCuts){
    HLTError("No ESD track cuts availible");
    return -ENOMEM;
  }

  // ---------------------

 // -maxpt
  if (argument.CompareTo("-maxpt")==0) {
    if (++ii>=argc) return -EPROTO;
    argument=argv[ii];

    Float_t minPt, maxPt;
    fESDTrackCuts->GetPtRange(minPt,maxPt);
    maxPt = argument.Atof(); 
    fESDTrackCuts->SetPtRange(minPt,maxPt);

    TString title = fESDTrackCuts->GetTitle();
    if (!title.CompareTo("No track cuts")) title = "";
    else title += " && ";
    title += Form("p_t < %f", maxPt);
    fESDTrackCuts->SetTitle(title);
    return 2;
  }    

  // -minpt
  if (argument.CompareTo("-minpt")==0) {
    if (++ii>=argc) return -EPROTO;
    argument=argv[ii];

    Float_t minPt, maxPt;
    fESDTrackCuts->GetPtRange(minPt,maxPt);
    minPt = argument.Atof(); 
    fESDTrackCuts->SetPtRange(minPt,maxPt);

    TString title = fESDTrackCuts->GetTitle();
    if (!title.CompareTo("No track cuts")) title = "";
    else title += " && ";
    title += Form("p_t > %f", minPt);
    fESDTrackCuts->SetTitle(title);
    return 2;
  }    

  // -min-ldca
  // minimum longitudinal dca to vertex
  if (argument.CompareTo("-min-ldca")==0) {
    if (++ii>=argc) return -EPROTO;
    argument=argv[ii];

    fESDTrackCuts->SetMinDCAToVertexZ(argument.Atof());
    TString title = fESDTrackCuts->GetTitle();
    if (!title.CompareTo("No track cuts")) title = "";
    else title += " && ";
    title += Form("DCAz > %f", argument.Atof());
    fESDTrackCuts->SetTitle(title);
    return 2;
  }
  
  // -max-ldca
  // maximum longitudinal dca to vertex
  if (argument.CompareTo("-max-ldca")==0) {
    if (++ii>=argc) return -EPROTO;
    argument=argv[ii];

    fESDTrackCuts->SetMaxDCAToVertexZ(argument.Atof());
    TString title = fESDTrackCuts->GetTitle();
    if (!title.CompareTo("No track cuts")) title = "";
    else title += " && ";
    title += Form("DCAz < %f", argument.Atof());
    fESDTrackCuts->SetTitle(title);
    return 2;
  }

  // -min-tdca
  // minimum transverse dca to vertex
  if (argument.CompareTo("-min-tdca")==0) {
    if (++ii>=argc) return -EPROTO;
    argument=argv[ii];

    fESDTrackCuts->SetMinDCAToVertexXY(argument.Atof());
    TString title = fESDTrackCuts->GetTitle();
    if (!title.CompareTo("No track cuts")) title = "";
    else title += " && ";
    title += Form("DCAr > %f", argument.Atof());
    fESDTrackCuts->SetTitle(title);
    return 2;
  }
  
  // -max-tdca
  // maximum transverse dca to vertex
  if (argument.CompareTo("-max-tdca")==0) {
    if (++ii>=argc) return -EPROTO;
    argument=argv[ii];

    fESDTrackCuts->SetMaxDCAToVertexXY(argument.Atof());
    TString title = fESDTrackCuts->GetTitle();
    if (!title.CompareTo("No track cuts")) title = "";
    else title += " && ";
    title += Form("DCAr < %f", argument.Atof());
    fESDTrackCuts->SetTitle(title);
    return 2;
  }

  // -etarange
  // +/- eta 
  if (argument.CompareTo("-etarange")==0) {
    if (++ii>=argc) return -EPROTO;
    argument=argv[ii];
    Float_t eta = argument.Atof();

    fESDTrackCuts->SetEtaRange(-eta,eta);     
    TString title = fESDTrackCuts->GetTitle();
    if (!title.CompareTo("No track cuts")) title = "";
    else title += " && ";
    title += Form("Eta[%f,%f]", argument.Atof(),argument.Atof());
    fESDTrackCuts->SetTitle(title);
    return 2;
  }

  // -- BINNING --------------

  // binningVzero
  if (argument.CompareTo("-binningVzero")==0) {
    if (++ii>=argc) return -EPROTO;
    argument=argv[ii];
    Int_t binning = argument.Atoi();
    if (++ii>=argc) return -EPROTO;
    argument=argv[ii];
    Float_t min = argument.Atof();
    if (++ii>=argc) return -EPROTO;
    argument=argv[ii];
    Float_t max = argument.Atof();

    fCorrObj->SetBinningVzero(binning, min, max);
    return 4;
  }

  // binningTpc
  if (argument.CompareTo("-binningTpc")==0) {
    if (++ii>=argc) return -EPROTO;
    argument=argv[ii];
    Int_t binning = argument.Atoi();
    if (++ii>=argc) return -EPROTO;
    argument=argv[ii];
    Float_t min = argument.Atof();
    if (++ii>=argc) return -EPROTO;
    argument=argv[ii];
    Float_t max = argument.Atof();

    fCorrObj->SetBinningTpc(binning, min, max);
    return 4;
  }

  // binningSpd
  if (argument.CompareTo("-binningSpd")==0) {
    if (++ii>=argc) return -EPROTO;
    argument=argv[ii];
    Int_t binning = argument.Atoi();
    if (++ii>=argc) return -EPROTO;
    argument=argv[ii];
    Float_t min = argument.Atof();
    if (++ii>=argc) return -EPROTO;
    argument=argv[ii];
    Float_t max = argument.Atof();

    fCorrObj->SetBinningSpd(binning, min, max);
    return 4;
  }

  // binningZdc
  if (argument.CompareTo("-binningZdc")==0) {
    if (++ii>=argc) return -EPROTO;
    argument=argv[ii];
    Int_t binning = argument.Atoi();
    if (++ii>=argc) return -EPROTO;
    argument=argv[ii];
    Float_t min = argument.Atof();
    if (++ii>=argc) return -EPROTO;
    argument=argv[ii];
    Float_t max = argument.Atof();

    fCorrObj->SetBinningZdc(binning, min, max);
    return 4;
  }

  // binningZnp
  if (argument.CompareTo("-binningZnp")==0) {
    if (++ii>=argc) return -EPROTO;
    argument=argv[ii];
    Int_t binning = argument.Atoi();
    if (++ii>=argc) return -EPROTO;
    argument=argv[ii];
    Float_t min = argument.Atof();
    if (++ii>=argc) return -EPROTO;
    argument=argv[ii];
    Float_t max = argument.Atof();

    fCorrObj->SetBinningZnp(binning, min, max);
    return 4;
  }

  // binningZem
  if (argument.CompareTo("-binningZem")==0) {
    if (++ii>=argc) return -EPROTO;
    argument=argv[ii];
    Int_t binning = argument.Atoi();
    if (++ii>=argc) return -EPROTO;
    argument=argv[ii];
    Float_t min = argument.Atof();
    if (++ii>=argc) return -EPROTO;
    argument=argv[ii];
    Float_t max = argument.Atof();

    fCorrObj->SetBinningZem(binning, min, max);
    return 4;
  }
  // binningZem
  if (argument.CompareTo("-binningCalo")==0) {
    if (++ii>=argc) return -EPROTO;
    argument=argv[ii];
    Int_t binning = argument.Atoi();
    if (++ii>=argc) return -EPROTO;
    argument=argv[ii];
    Float_t min = argument.Atof();
    if (++ii>=argc) return -EPROTO;
    argument=argv[ii];
    Float_t max = argument.Atof();

    fCorrObj->SetBinningCalo(binning, min, max);
    return 4;
  }
  if (argument.CompareTo("-enablePhos")==0) {
    if (++ii>=argc) return -EPROTO;
    argument=argv[ii];
    Int_t enabled = argument.Atoi();
    fCorrObj->SetProcessPhos(enabled);
    return 2;
  }
  if (argument.CompareTo("-enableEmcal")==0) {
    if (++ii>=argc) return -EPROTO;
    argument=argv[ii];
    Int_t enabled = argument.Atoi();
    fCorrObj->SetProcessEmcal(enabled);
    return 2;
  }

  // -- enable
  if (argument.CompareTo("-enableCALO")==0) {
    fCorrObj->SetProcessCALO(kTRUE);
    return 1;
  }
  if (argument.CompareTo("-enableVZERO")==0) {
    fCorrObj->SetProcessVZERO(kTRUE);
    return 1;
  }
  if (argument.CompareTo("-enableZDC")==0) {
    fCorrObj->SetProcessZDC(kTRUE);
    return 1;
  }
  if (argument.CompareTo("-enableTPC")==0) {
    fCorrObj->SetProcessTPC(kTRUE);
    return 1;
  }
  if (argument.CompareTo("-enableSPD")==0) {
    fCorrObj->SetProcessSPD(kTRUE);
    return 1;
  }

  // -- disable
  if (argument.CompareTo("-disableCALO")==0) {
    fCorrObj->SetProcessCALO(kFALSE);
    return 1;
  }
  if (argument.CompareTo("-disableVZERO")==0) {
    fCorrObj->SetProcessVZERO(kFALSE);
    return 1;
  }
  if (argument.CompareTo("-disableZDC")==0) {
    fCorrObj->SetProcessZDC(kFALSE);
    return 1;
  }
  if (argument.CompareTo("-disableTPC")==0) {
    fCorrObj->SetProcessTPC(kFALSE);
    return 1;
  }
  if (argument.CompareTo("-disableSPD")==0) {
    fCorrObj->SetProcessSPD(kFALSE);
    return 1;
  }

  // unknown argument
  return -EINVAL;
}

// #################################################################################
Int_t AliHLTMultiplicityCorrelationsComponent::DoDeinit() {
  // see header file for class documentation

  if (fCorrObj) 
    delete fCorrObj;
  fCorrObj = NULL;
  
  if (fESDTrackCuts) 
    delete fESDTrackCuts;
  fESDTrackCuts = NULL;
  
  fUID = 0;

  return 0;
}

// #################################################################################
Int_t AliHLTMultiplicityCorrelationsComponent::DoEvent(const AliHLTComponentEventData& evtData,
					AliHLTComponentTriggerData& /*trigData*/) {
  // see header file for class documentation

  Int_t iResult=0;

  // -- Only use data event
  if (!IsDataEvent()) 
    return 0;
  
  if( fUID == 0 ){
    TTimeStamp t;
    fUID = ( gSystem->GetPid() + t.GetNanoSec())*10 + evtData.fEventID;
  }

  // -- Get ESD object 
  AliESDEvent *esdEvent = NULL;
  for ( const TObject *iter = GetFirstInputObject(kAliHLTDataTypeESDObject); iter != NULL; iter = GetNextInputObject() ) {
    esdEvent = dynamic_cast<AliESDEvent*>(const_cast<TObject*>( iter ) );
    if( !esdEvent ){ 
      HLTWarning("Wrong ESDEvent object received");
      iResult = -1;
      continue;
    }
    esdEvent->GetStdContent();
  }

  // -- Get VZEROESD object 
  AliESDVZERO *esdVZERO = NULL;
  for ( const TObject *iter = GetFirstInputObject(kAliHLTDataTypeESDContent|kAliHLTDataOriginVZERO); 
	iter != NULL; iter = GetNextInputObject() ) {
    esdVZERO = dynamic_cast<AliESDVZERO*>(const_cast<TObject*>( iter ) );
    if( !esdVZERO ){ 
      HLTWarning("Wrong VZERO ESDEvent object received");
      iResult = -1;
      continue;
    }
  }

  // -- Get SPD clusters
  // ---------------------
  const AliHLTComponentBlockData* iter = NULL;
  Int_t totalClusters = 0;
  Int_t innerClusters = 0;
  Int_t outerClusters = 0;

  for ( iter = GetFirstInputBlock(kAliHLTDataTypeClusters|kAliHLTDataOriginITSSPD); 
	iter != NULL; iter = GetNextInputBlock() ) {

    AliHLTITSClusterData *clusterData = reinterpret_cast<AliHLTITSClusterData*>(iter->fPtr);
    Int_t nClusters = clusterData->fSpacePointCnt;
    for( int icl=0; icl<nClusters; icl++ ) {
      AliHLTITSSpacePointData &d = clusterData->fSpacePoints[icl];
      if ( d.fLayer == 0 ) ++innerClusters;
      else if ( d.fLayer == 1 ) ++outerClusters;
    }  
    
    totalClusters += nClusters;
  }

  fCorrObj->SetSPDClusters(innerClusters,outerClusters);

  // -- Process Event
  // ------------------
  if (esdEvent)
    iResult = fCorrObj->ProcessEvent(esdEvent,esdVZERO,totalClusters);

  if (iResult) {
    HLTError("Error while processing event inside multiplicity correlation object");
    return iResult;
  }

  // -- Send histlist
  PushBack(dynamic_cast<TObject*>(fCorrObj->GetHistList()),
	   kAliHLTDataTypeTObject|kAliHLTDataOriginHLT,fUID);
 
  return iResult;
}

// #################################################################################
Int_t AliHLTMultiplicityCorrelationsComponent::Reconfigure(const Char_t* cdbEntry, const Char_t* chainId) {
  // see header file for class documentation

  Int_t iResult=0;
  TString cdbPath;
  if (cdbEntry) {
    cdbPath=cdbEntry;
  } else {
    cdbPath="HLT/ConfigGlobal/";
    cdbPath+=GetComponentID();
  }

  AliInfoClass(Form("reconfigure '%s' from entry %s%s", chainId, cdbPath.Data(), cdbEntry?"":" (default)"));
  iResult=ConfigureFromCDBTObjString(cdbPath);

  return iResult;
}

// #################################################################################
Int_t AliHLTMultiplicityCorrelationsComponent::ReadPreprocessorValues(const Char_t* /*modules*/) {
  // see header file for class documentation
  ALIHLTERRORGUARD(5, "ReadPreProcessorValues not implemented for this component");
  return 0;
}
