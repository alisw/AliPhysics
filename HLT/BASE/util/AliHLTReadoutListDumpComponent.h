// -*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTREADOUTLISTDUMPCOMPONENT_H
#define ALIHLTREADOUTLISTDUMPCOMPONENT_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTReadoutListDumpComponent.h
    @author Matthias Richter
    @date   
    @brief  Collect ESDs of multiple events and write toi file
*/

#include "AliHLTFileWriter.h"

class TH1I;
class TH2I;
class AliHLTReadoutList;

/**
 * @class AliHLTReadoutListDumpComponent
 * The ReadoutListDump component fetches the DAQ readout list object
 * and can store the information in different ways, like e.g. in a histogram
 * or a tree.
 *
 * <h2>General properties:</h2>
 *
 * Component ID: \b ReadoutListDump                                      <br>
 * Library: \b libAliHLTUtil.so					         <br>
 * Input Data Types: {HLTRDLST:HLT },                                    <br>
 * Output Data Types: none					         <br>
 *
 * <h2>Mandatory arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *
 * <h2>Optional arguments:</h2>
 * \li -binary                                                           <br>
 *      fetch the binary readout list block (default)
 * \li -decision                                                         <br>
 *      fetch the readout list from the HLT decision object
 * The only AliHLTFileWriter argument of relevance is the \em -directory
 * argument. See AliHLTFileWriter for full list of arguments.
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
 * No data published (AliHLTDataSink).
 *
 * @ingroup alihlt_util_components
 */
class AliHLTReadoutListDumpComponent : public AliHLTFileWriter
{
 public:
  /** standard constructor */
  AliHLTReadoutListDumpComponent();
  /** destructor */
  virtual ~AliHLTReadoutListDumpComponent();

  /**
   * The id of the component.
   * @return component id (string)
   */
  const char* GetComponentID() {return "ReadoutListDump";};

  /**
   * Spawn function.
   * @return new class instance
   */
  AliHLTComponent* Spawn() {return new AliHLTReadoutListDumpComponent;}

  enum {
    kModeBinaryList = 1, // fetch the readout list block
    kModeHLTDecision = 2 // fetch the readout list from the HLT decision object
  };

 protected:
  // interface functions

  /// inherited form AliHLTFileWriter
  int InitWriter();

  /// inherited form AliHLTFileWriter
  int CloseWriter();

  /// inherited form AliHLTDataSink
  int DumpEvent( const AliHLTComponentEventData& evtData,
		 const AliHLTComponentBlockData* blocks, 
		 AliHLTComponentTriggerData& trigData );
  
  using AliHLTFileWriter::DumpEvent;

  /// inherited form AliHLTFileWriter
  int ScanArgument(int argc, const char** argv);

  /**
   * Create the histogram for monitoring of the readout list.
   * Each bin corresponds to a bit in the bitfield of the readout
   * list.
   * The object has to be deleted by the caller. 
   */
  static TH1I* CreateReadoutListHistogram();

  /**
   * Create the histogram for monitoring of the readout list.
   * Plot readout list bits vs. CTP trigger bit
   * The object has to be deleted by the caller. 
   */
  static TH2I* CreateReadoutListVsCTPHistogram();

  /**
   * Fill histogram from the readout list.
   */
  static int FillReadoutListHistogram(TH1I* histo, const AliHLTReadoutList* list);

  /**
   * Fill histogram from the readout list.
   */
  static int FillReadoutListVsCTP(TH2I* histo, const AliHLTReadoutList* list, const AliHLTComponentTriggerData* trigData);

private:
  /** copy constructor prohibited */
  AliHLTReadoutListDumpComponent(const AliHLTReadoutListDumpComponent&);
  /** assignment operator prohibited */
  AliHLTReadoutListDumpComponent& operator=(const AliHLTReadoutListDumpComponent&);

  unsigned fMode; /// how to get the readout list (kModeBinaryList, kModeHLTDecision)
  TH1I*    fBitsHisto; /// the histogram to be filled
  TH2I*    fBitsVsCTP; /// histogram of bits vs. ctp triggers

  ClassDef(AliHLTReadoutListDumpComponent, 0)
};
#endif
