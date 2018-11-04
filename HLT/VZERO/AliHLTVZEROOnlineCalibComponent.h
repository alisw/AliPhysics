//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTVZEROONLINECALIBCOMPONENT_H
#define ALIHLTVZEROONLINECALIBCOMPONENT_H

/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file    AliHLTVZEROOnlineCalibComponent.h
 @author  David Rohr <drohr@cern.ch>
 @brief   VZERO online calibration component
 */

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt


#include "AliHLTProcessor.h"
#include "AliOptionParser.h"
#include "TH1F.h"

class TTree;

class AliRunInfo;
class AliESDVZERO;
class AliVZERORecoParam;

class AliHLTVZEROOnlineCalibComponent : public AliHLTProcessor, public AliOptionParser {
public:
    
    /*
     * ---------------------------------------------------------------------------------
     *                            Constructor / Destructor
     * ---------------------------------------------------------------------------------
     */
    
    /** constructor */
    AliHLTVZEROOnlineCalibComponent();
    
    /** destructor */
    virtual ~AliHLTVZEROOnlineCalibComponent();
    
    /*
     * ---------------------------------------------------------------------------------
     * Public functions to implement AliHLTComponent's interface.
     * These functions are required for the registration process
     * ---------------------------------------------------------------------------------
     */
    
    /** interface function, see @ref AliHLTComponent for description */
    const Char_t* GetComponentID();
    
    /** interface function, see @ref AliHLTComponent for description */
    AliHLTComponentDataType GetOutputDataType();
    void GetInputDataTypes(AliHLTComponentDataTypeList& tgtList);
    
    /** interface function, see @ref AliHLTComponent for description */
    void GetOutputDataSize( ULong_t& constBase, Double_t& inputMultiplier );
    
    /** interface function, see @ref AliHLTComponent for description */
    void GetOCDBObjectDescription( TMap* const targetMap);
    
    /** interface function, see @ref AliHLTComponent for description */
    AliHLTComponent* Spawn();
    
    //overload from AliOptionParser
    int ProcessOption(TString option, TString value);
    
protected:
    
    /*
     * ---------------------------------------------------------------------------------
     * Protected functions to implement AliHLTComponent's interface.
     * These functions provide initialization as well as the actual processing
     * capabilities of the component.
     * ---------------------------------------------------------------------------------
     */
    
    // AliHLTComponent interface functions
    
    /** interface function, see @ref AliHLTComponent for description */
    Int_t DoInit( Int_t argc, const Char_t** argv );
    
    /** interface function, see @ref AliHLTComponent for description */
    Int_t DoDeinit();
    
    /** interface function, see @ref AliHLTComponent for description */
    Int_t DoEvent( const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData);
    
    using AliHLTProcessor::DoEvent;
    
    /** interface function, see @ref AliHLTComponent for description */
    Int_t ScanConfigurationArgument(Int_t argc, const Char_t** argv);
    
    /** interface function, see @ref AliHLTComponent for description */
    Int_t Reconfigure(const Char_t* cdbEntry, const Char_t* chainId);
    
    ///////////////////////////////////////////////////////////////////////////////////
    
private:
    
    /*
     * ---------------------------------------------------------------------------------
     * Private functions to implement AliHLTComponent's interface.
     * These functions provide initialization as well as the actual processing
     * capabilities of the component.
     * ---------------------------------------------------------------------------------
     */
    
    /** copy constructor prohibited */
    AliHLTVZEROOnlineCalibComponent(const AliHLTVZEROOnlineCalibComponent&);
    
    /** assignment operator prohibited */
    AliHLTVZEROOnlineCalibComponent& operator=(const AliHLTVZEROOnlineCalibComponent&);
    
    /*
     * ---------------------------------------------------------------------------------
     *                             Members - private
     * ---------------------------------------------------------------------------------
     */
    
    /** runInfo Object */
    AliRunInfo            *fRunInfo;            // see above
    
    /** VZERO reco param instance */
    AliVZERORecoParam     *fVZERORecoParam;     //! transient
    
    //reference trigger
    TString fRefTrigger;
    TString fTrigger1;
    TString fTrigger2;
    TString fTrigger3;
    
    TH1F fHistMult[24]; //all histograms we'll need
    
    ClassDef(AliHLTVZEROOnlineCalibComponent, 0)
};
#endif
