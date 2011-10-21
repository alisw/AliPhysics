// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Federico Ronchetti                                    *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/** @file   AliHLTEMCALAgent.cxx
 @author Federico Ronchetti
 @date   
 @brief  Agent of the libAliHLTEMCAL library
 */

#include "AliHLTEMCALAgent.h"
#include "AliHLTConfiguration.h"
#include "AliHLTEMCALDefinitions.h"
#include "AliHLTOUT.h"
#include "AliHLTOUTHandlerChain.h"
#include "AliRunLoader.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"


// #include "AliHLTEMCALConstant.h"
#include "AliHLTEMCALConstants.h"
#include "AliHLTEMCALMapper.h"

/** global instance for agent registration */
AliHLTEMCALAgent gAliHLTEMCALAgent;

// component headers
//#include "AliHLTEMCALCalibrationComponent.h"
#include "AliHLTCaloClusterAnalyser.h"
#include "AliHLTEMCALClusterizerComponent.h"
#include "AliHLTEMCALDigitMakerComponent.h"
//#include "AliHLTEMCALESDEntriesMakerComponent.h"
//#include "AliHLTEMCALHistogramProducerComponent.h"
//#include "AliHLTEMCALModuleCalibrationProcessorComponent.h"
//#include "AliHLTEMCALMonitorTriggerComponent.h"
#include "AliHLTEMCALRawAnalyzerComponent.h"
#include "AliHLTEMCALRawAnalyzerCrudeComponent.h"
#include "AliHLTEMCALRawAnalyzerPeakFinderComponent.h"
//#include "AliHLTEMCALRcuCalibrationProcessorComponent.h"
//#include "AliHLTEMCALRcuDAComponent.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTEMCALAgent)

AliHLTEMCALAgent::AliHLTEMCALAgent() : AliHLTModuleAgent("EMCAL"),
fRawDataHandler(NULL)
{
    // see header file for class documentation
    // or
    // refer to README to build package
    // or
    // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTEMCALAgent::~AliHLTEMCALAgent()
{
    // see header file for class documentation
}

int AliHLTEMCALAgent::CreateConfigurations(AliHLTConfigurationHandler* handler,
                                           AliRawReader* /*rawReader*/,
                                           AliRunLoader* /*runloader*/) const
{
    // see header file for class documentation
    if (handler) 
    {
        //      const char* cdbEntry="EMCAL/Calib/Parameters";
        //AliCDBManager* pMan=AliCDBManager::Instance();
        //AliEMCALParam* pEMCALParam=NULL;
        
        //       if (pMan) 
        // 	{
        // 	  AliCDBEntry *pEntry = pMan->Get(cdbEntry);
        // 	  if (pEntry && 
        // 	      pEntry->GetObject() &&
        // 	      (pEMCALParam=dynamic_cast<AliEMCALParam*>(pEntry->GetObject()))) 
        // 	    {
        // 	    } else 
        // 	    {
        // 	      HLTWarning("can not load AliEMCALParam from CDB entry %s", cdbEntry);
        // 	    }
        // 	}
        
        int moduleStart = 0;
        int moduleEnd = 9;
        
        int rcuStart = 0;
        int rcuEnd = 1;
        
        Int_t rcusPerModule = 2;
        Int_t ddlOffset = 4608; 
        
        TString mergerInput;
        TString sinkClusterInput;
        TString emInput;
        
        for (int module = moduleStart; module <= moduleEnd; module++) 
        {
            TString clInput;
            
            for(int rcu = rcuStart; rcu < rcuEnd; rcu++) 
            {
                TString arg, publisher, ra, dm;
                // raw data publisher components
                publisher.Form("EMCAL-RP_%02d_%d", module, rcu);
                arg.Form("-verbose -minid %d -datatype 'DDL_RAW ' 'EMCA'  -dataspec 0x%x ", ddlOffset + module*(rcusPerModule) + rcu, 0x1 << (module*rcusPerModule + rcu));
                
                handler->CreateConfiguration(publisher.Data(), "AliRawReaderPublisher", NULL , arg.Data());
                
                // Raw analyzer
                arg = "";
                ra.Form("EMCAL-RA_%02d_%d", module, rcu);
                handler->CreateConfiguration(ra.Data(), "EmcalRawCrude", publisher.Data(), arg.Data());
                
                // digit maker components
                dm.Form("EMCAL-DM_%02d_%d", module, rcu);
                arg="";
                arg.Form("-sethighgainfactor 0.0153 -setlowgainfactor 0.2448 -setdigitthresholds 0.005 0.002");
                handler->CreateConfiguration(dm.Data(), "EmcalDigitMaker", ra.Data(), arg.Data());
                
                if(clInput.Length() > 0) clInput += " ";
                clInput+=dm;
            }
            
            TString arg, cl, ca;
            
            cl.Form("EMCAL-CF_%02d", module);
            arg = "";
            arg.Form("-digitthreshold 0.005 -recpointthreshold 0.1 -modulemode");
            handler->CreateConfiguration(cl.Data(), "EmcalClusterizer", clInput.Data(), arg.Data());
            
            //ca.Form("EMCAL-CA_%02d", module);
            //arg = " ";
            //handler->CreateConfiguration(ca.Data(), "CaloClusterAnalyser", cl.Data(), arg.Data());
            
            if(emInput.Length() > 0) emInput += " ";
            emInput += ca;
        }
        
        
        TString arg, em;
        
        // tracker finder components
        
        // later
        //      em.Form("EMCAL-EM");
        //arg = " ";
        //handler->CreateConfiguration(em.Data(), "EmcalEsdEntriesMaker", emInput.Data(), " ");
        
    }
    
    return 0;
}

const char* AliHLTEMCALAgent::GetReconstructionChains(AliRawReader* /*rawReader*/,
                                                      AliRunLoader* runloader) const
{
    // see header file for class documentation
    if (runloader) {
        // reconstruction chains for AliRoot simulation
        // Note: run loader is only available while running embedded into
        // AliRoot simulation
        
        // if (runloader->GetLoader("EMCALLoader") != NULL)
        //     return "EMCAL-EM";
    }
    return NULL;
}

const char* AliHLTEMCALAgent::GetRequiredComponentLibraries() const
{
    // see header file for class documentation
    return NULL;
}

int AliHLTEMCALAgent::RegisterComponents(AliHLTComponentHandler* pHandler) const
{
    // see header file for class documentation
    if (!pHandler) return -EINVAL;
    
    pHandler->AddComponent(new AliHLTEMCALRawAnalyzerCrudeComponent);
    pHandler->AddComponent(new AliHLTEMCALDigitMakerComponent);
    pHandler->AddComponent(new AliHLTEMCALClusterizerComponent);
    //pHandler->AddComponent(new AliHLTCaloClusterAnalyserComponent);			 
    //pHandler->AddComponent(new AliHLTEMCALESDEntriesMakerComponent);
    
    return 0;
}


int AliHLTEMCALAgent::GetHandlerDescription(AliHLTComponentDataType dt,
                                            AliHLTUInt32_t spec,
                                            AliHLTOUTHandlerDesc& desc) const
{
    // see header file for class documentation

    // FIXME: there is memory allocated in the mapper, this happens for every event
    // and is time consuming, think about initializing the mappers only once, or make
    // the mapper class more flexible to handle more than one specification
    AliHLTEMCALMapper mapper(spec);
    
    // raw data blocks to be fed into offline reconstruction
    if (dt==(kAliHLTDataTypeDDLRaw|kAliHLTDataOriginEMCAL)) 
    {
        if(mapper.GetDDLFromSpec(spec) >= 0)
        {
            desc=AliHLTOUTHandlerDesc(kRawReader, dt, GetModuleId());
            return 1;
        } 
        else 
        {
            HLTWarning("Handler can not process data inconsistent with a single EMCAL DDL from specification % d", spec);
            return 0;
        }
    }
    return 0;
}

AliHLTOUTHandler* AliHLTEMCALAgent::GetOutputHandler(AliHLTComponentDataType dt,
                                                     AliHLTUInt32_t /*spec*/)
{
    // see header file for class documentation
    
    // raw data blocks to be fed into offline reconstruction
    if (dt==(kAliHLTDataTypeDDLRaw|kAliHLTDataOriginEMCAL)) 
    {
        if (!fRawDataHandler) 
        {
            fRawDataHandler = new AliHLTEMCALAgent::AliHLTEMCALRawDataHandler;
        }
        return fRawDataHandler;
    }
    
    return NULL;
}

int AliHLTEMCALAgent::DeleteOutputHandler(AliHLTOUTHandler* pInstance)
{
    // see header file for class documentation
    if (pInstance==NULL) return -EINVAL;
    
    if (pInstance==fRawDataHandler) {
        delete fRawDataHandler;
        fRawDataHandler=NULL;
    }
    return 0;
}

AliHLTEMCALAgent::AliHLTEMCALRawDataHandler::AliHLTEMCALRawDataHandler()
{
    // see header file for class documentation
}

AliHLTEMCALAgent::AliHLTEMCALRawDataHandler::~AliHLTEMCALRawDataHandler()
{
    // see header file for class documentation
}

int AliHLTEMCALAgent::AliHLTEMCALRawDataHandler::ProcessData(AliHLTOUT* pData)
{
    // see header file for class documentation
    if (!pData) return -EINVAL;
    
    AliHLTComponentDataType dt = kAliHLTVoidDataType;
    AliHLTUInt32_t spec=kAliHLTVoidDataSpec;
    
    AliHLTEMCALMapper mapper(spec);
    mapper.InitDDLSpecificationMapping();
    
    int iResult = pData->GetDataBlockDescription(dt, spec);
    if (iResult>=0) 
    {
        int ddl = -1;
        if((ddl = mapper.GetDDLFromSpec(spec)) >=0)
        {
            iResult = ddl;
        }
    } 
    else 
    {
        HLTError("Handler can not process data inconsistent with a single EMCAL DDL from specification % d", spec);
        iResult=-EBADMSG;
    }
    return iResult;
}
