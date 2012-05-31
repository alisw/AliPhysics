// -*- Mode: C++ -*-

// $Id: AliHLTMCGeneratorComponent.h  $

#ifndef ALIHLTMCGENERATORCOMPONENT_H
#define ALIHLTMCGENERATORCOMPONENT_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTMCGeneratorComponent.h
    @author Jochen Thaeder <thaeder@kip.uni-heidelberg.de>
    @date   
    @brief  Component for generating MC events
    @note   The class is used in Offline (AliRoot) context
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliGenerator.h"
#include "AliRunLoader.h"
#include "AliMCEvent.h"

#include "AliHLTDataSource.h"
#include "AliHLTMCEvent.h"

/**
 * @class AliHLTMCGeneratorComponent
 * An HLT data source component which simulates PYTHIA events
 * right now set up for Jet Production, this can be extended by
 * a person in need.
 *
 * The events are pulished in form of AliHLTMCEvent.
 *
 * <h2>General properties:</h2>
 *
 * Component ID: \b MCGenerator <br>
 * Library:      \b libAliHLTUtil.so<br>
 * Input Data Types: none <br>
 * Output Data Types: kAliHLTDataTypeMCObject|kAliHLTDataOriginHLT<br>
 *                    <i>AliHLTMCEvent</i><br>
 *
 * <h2>Mandatory arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 * \li -seed         <i> Seed for generation </i> <br>
 *
 * \li -runtype      <i> Run type as a string </i> <br>
 *
 * \li -nevents      <i> Number of events to generate </i> <br>
 *
 * \li -coneRadius   <i> Cone radius for PyCell </i> <br>
 *
 * \li -jetCutMinEt  <i> Final state min Et cut for PyCell </i> <br>
 *
 * \li -applyParticleCuts  <i> Apply particle cuts before filling in AliHLTMCEvent </i> <br>
 *
 * 
 * <h2>Optional arguments:</h2>
 *
 * <h2>Configuration:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 * Configuration by component arguments.
 *
 * <h2>Default CDB entries:</h2>
 * The component loads no CDB entries.
 *
 * <h2>Performance:</h2>
 * The component does not process any event data.
 *
 * <h2>Memory consumption:</h2>
 * The component does not process any event data.
 *
 * <h2>Output size:</h2>
 * According to the available data. The component is an AliHLTDataSource
 * and inteded to be used in the AliHLTSystem framework only. The component
 * implements the standard AliHLTSystem adaptive buffer allocation. 
 *
 * @ingroup alihlt_util_components
 */

class AliHLTMCGeneratorComponent : public AliHLTDataSource  {
  
 public:

  /*
   * ---------------------------------------------------------------------------------
   *                            Constructor / Destructor
   * ---------------------------------------------------------------------------------
   */
  
  /** standard constructor */
  AliHLTMCGeneratorComponent();

  /** destructor */
  virtual ~AliHLTMCGeneratorComponent();
  
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

  /** interface function, see @ref AliHLTComponent for description */
  virtual void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier );

  /** interface function, see @ref AliHLTComponent for description */
  AliHLTComponent* Spawn();

 protected:
  /*
   * ---------------------------------------------------------------------------------
   * Protected functions to implement AliHLTComponent's interface.
   * These functions provide initialization as well as the actual processing
   * capabilities of the component. 
   * ---------------------------------------------------------------------------------
   */

  /**
   * Init method. 
   * Overwrites the AliHLTFilePublisher::DoInit() method.  
   * @param argc           size of the argument array
   * @param argv           agument array for component initialization
   * @return number of processed members of the argv <br>
   *         -EINVAL unknown argument <br>
   *         -EPROTO parameter for argument missing <br>
   */
  Int_t DoInit( int argc, const char** argv );

  /**
   * Deinit method. Calls also the one of AliHLTFilePublisher.
   */
  Int_t DoDeinit();

  /**
   * Data processing method for the component.
   * The component uses the @ref alihltcomponent-high-level-interface
   * to put serialized Root object into the output stream. Despite of that it
   * implements the low-level DumpEvent method in order to allow child classes
   * to use the low-level method.
   * @param evtData       event data structure
   * @param trigData	  trigger data structure
   * @param outputPtr     not used
   * @param size          not used
   * @param outputBlocks  not used
   * @return
   */
  Int_t GetEvent( const AliHLTComponentEventData& evtData,
		AliHLTComponentTriggerData& trigData,
		AliHLTUInt8_t* outputPtr, 
		AliHLTUInt32_t& size,
		vector<AliHLTComponentBlockData>& outputBlocks);

  using AliHLTDataSource::GetEvent;

 private:

  /*
   * ---------------------------------------------------------------------------------
   * Private functions to implement AliHLTComponent's interface.
   * These functions provide initialization as well as the actual processing
   * capabilities of the component. 
   * ---------------------------------------------------------------------------------
   */

  /** copy constructor prohibited */
  AliHLTMCGeneratorComponent(const AliHLTMCGeneratorComponent&);

  /** assignment operator prohibited */
  AliHLTMCGeneratorComponent& operator=(const AliHLTMCGeneratorComponent&);

  /** Create Generator Factory
   *  @return          Ptr to generator
   */
  AliGenerator* GeneratorFactory();

  /*
   * ---------------------------------------------------------------------------------
   *                             Members - private
   * ---------------------------------------------------------------------------------
   */

  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // ++ Enum / static const
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  /** Run Type */
  enum PprRun_t {
    kPythia6Jets20_24,   kPythia6Jets24_29,   kPythia6Jets29_35,
    kPythia6Jets35_42,   kPythia6Jets42_50,   kPythia6Jets50_60,
    kPythia6Jets60_72,   kPythia6Jets72_86,   kPythia6Jets86_104,
    kPythia6Jets104_125, kPythia6Jets125_150, kPythia6Jets150_180,
    kPyJetJet, kPyGammaJetPHOS, kRunMax
  };
  
  /** Run type name */
  static const Char_t *fgkPprRunName[];            //! transient

  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // ++ Data Members 
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  /* Ptr to current AliMCEvent, to be shipped out*/
  AliMCEvent          *fpMC;                       //! transient
  
  /* Ptr to current AliHLTMCEvent, to be shipped out*/
  AliHLTMCEvent       *fpHLTMC;                    //! transient
  
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // ++ Event Handling
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  /** Number of events */
  Int_t                fNEvents;                   // see above

  /** Current event */
  Int_t                fCurrentEvent;              // see above

  /** Event number */
  Int_t                fEventNumber;               // see above

  /** Run number */
  Int_t                fRunNumber;                 // see above

  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // ++ MC Generation
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  /** RunLoader */
  AliRunLoader        *fRunLoader;                 //! transient

  /** Event venerator base class */
  AliGenerator        *fGenerator;                 //! transient

  /** Comment line */
  TString              fComment;                   // see above

  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // ++ MC Generation
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  /** Seed for random generation - default 0 */
  UInt_t               fSeed;                      // see above
  
  /** Run type - default kPythia6Jets104*/
  PprRun_t             fRunType;                   // see above

  /** Center of mass energy - default 14 TeV*/
  Float_t              fEcms;                      // see above
  
  /** Max Eta - Final state kinematic cut for PyCell - default 0.2 */
  Float_t              fJetEtaMax;                 // see above

  /** Min Eta - Final state kinematic cut for PyCell - default -0.2 */
  Float_t              fJetEtaMin;                 // see above

  /** Max Et - Final state kinematic cut for PyCell - default 1000 GeV*/ 
  Float_t              fJetEtMax;                  // see above

  /** Min Et - Final state kinematic cut for PyCell - default 10 GeV */  
  Float_t              fJetEtMin;                  // see above

  /** Cone radius for PyCell - default 0.4 */
  Float_t              fJetConeRadius;             // see above
  
  /** Min Pt transfer of the hard scattering - default 10 GeV */
  Double_t             fPtHardMin;                 // see above

  /** Min Pt transfer of the hard scattering - default infinity */
  Double_t             fPtHardMax;                 // see above
  
  /** Quench mode 0 = no, 1 = AM, 2 = IL - default 0*/
  Int_t                fQuenching ;                // see above

  /** Q Hat - default 20 */
  Float_t              fQhat;                      // see above
 
  /** Apply particle cuts, before filling in AliHLTMCEvent */
  Bool_t               fApplyParticleCuts;         // see above

  ClassDef(AliHLTMCGeneratorComponent, 0)
};
#endif
