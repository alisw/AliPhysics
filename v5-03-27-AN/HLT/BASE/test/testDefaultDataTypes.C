// $Id$

/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
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

/** @file   testDefaultDataTypes.C
    @author Matthias Richter
    @date   
    @brief  Test program for default data types
 */

#include "AliHLTDataTypes.h"
#include "AliHLTComponent.h"

int main(int /*argc*/, const char** /*argv*/)
{
  /** multiple output data types */
  const AliHLTComponentDataType kAliHLTMultipleDataTypeTest = {
    sizeof(AliHLTComponentDataType),
    {'M','U','L','T','I','P','L','E'},
    {'P','R','I','V'}
  };

  /** data to file exchange subscriber */
  const AliHLTComponentDataType kAliHLTDataTypeFXSCalibTest = {
    sizeof(AliHLTComponentDataType),
    {'F','X','S','_','C','A','L',' '},
    {'H','L','T',' '}
  };

  /** DDL list data type */
  const AliHLTComponentDataType kAliHLTDataTypeDDLTest = {
    sizeof(AliHLTComponentDataType),
    {'D','D','L','L','I','S','T',' '},
    {'H','L','T',' '}
  };

  /** SOR data type */
  const AliHLTComponentDataType kAliHLTDataTypeSORTest = {
    sizeof(AliHLTComponentDataType),
    {'S','T','A','R','T','O','F','R'},
    {'P','R','I','V'}
  };

  /** EOR data type */
  const AliHLTComponentDataType kAliHLTDataTypeEORTest = {
    sizeof(AliHLTComponentDataType),
    {'E','N','D','O','F','R','U','N'},
    {'P','R','I','V'}
  };

  /** run type data block */
  const AliHLTComponentDataType kAliHLTDataTypeRunTypeTest = {
    sizeof(AliHLTComponentDataType),
    {'R','U','N','T','Y','P','E',' '},
    {'P','R','I','V'}
  };

  /** Event type specification */
  const AliHLTComponentDataType kAliHLTDataTypeEventTest = {
    sizeof(AliHLTComponentDataType),
    {'E','V','E','N','T','T','Y','P'},
    {'P','R','I','V'}
  };

  /** Configuration event data type */
  const AliHLTComponentDataType kAliHLTDataTypeComConfTest = {
    sizeof(AliHLTComponentDataType),
    {'C','O','M','_','C','O','N','F'},
    {'P','R','I','V'}
  };

  /** DCS value update event */
  const AliHLTComponentDataType kAliHLTDataTypeUpdtDCSTest = {
    sizeof(AliHLTComponentDataType),
    {'U','P','D','T','_','D','C','S'},
    {'P','R','I','V'}
  };

  /** RAW DDL data specification, data publisher will set type id and origin correctly */
  const AliHLTComponentDataType kAliHLTDataTypeDDLRawTest = {
    sizeof(AliHLTComponentDataType),
    {'D','D','L','_','R','A','W',' '},
    {'*','*','*','\0'}
  };

  /** ESD data specification */
  const AliHLTComponentDataType kAliHLTDataTypeESDObjectTest = {
    sizeof(AliHLTComponentDataType),
    {'A','L','I','E','S','D','V','0'},
    {'*','*','*','\0'}
  };

  /** ESD tree data specification */
  const AliHLTComponentDataType kAliHLTDataTypeESDTreeTest = {
    sizeof(AliHLTComponentDataType),
    {'E','S','D','_','T','R','E','E'},
    {'*','*','*','\0'}
  };

  /** AliRoot TreeD data specification */
  const AliHLTComponentDataType kAliHLTDataTypeAliTreeDTest = {
    sizeof(AliHLTComponentDataType),
    {'A','L','I','T','R','E','E','D'},
    {'*','*','*','\0'}
  };

  /** AliRoot TreeR data specification */
  const AliHLTComponentDataType kAliHLTDataTypeAliTreeRTest = {
    sizeof(AliHLTComponentDataType),
    {'A','L','I','T','R','E','E','R'},
    {'*','*','*','\0'}
  };

  /** 16 bit Hardware address selection data specification, origin is 'any' */
  const AliHLTComponentDataType kAliHLTDataTypeHwAddr16Test = {
    sizeof(AliHLTComponentDataType),
    {'H','W','A','D','D','R','1','6'},
    {'*','*','*','\0'}
  };

  /** Event statistics */
  const AliHLTComponentDataType kAliHLTDataTypeEventStatisticsTest = {
    sizeof(AliHLTComponentDataType),
    {'E','V','_','S','T','A','T','I'},
    {'*','*','*','\0'}
  };

  /** Event summary */
  const AliHLTComponentDataType kAliHLTDataTypeEventSummaryTest = {
    sizeof(AliHLTComponentDataType),
    {'E','V','_','S','U','M','M','A'},
    {'H','L','T',' '}
  };

  /** Run statistics */
  const AliHLTComponentDataType kAliHLTDataTypeRunStatisticsTest = {
    sizeof(AliHLTComponentDataType),
    {'R','U','N','S','T','A','T','I'},
    {'*','*','*','\0'}
  };

  /** Run summary */
  const AliHLTComponentDataType kAliHLTDataTypeRunSummaryTest = {
    sizeof(AliHLTComponentDataType),
    {'R','U','N','S','U','M','M','A'},
    {'H','L','T',' '}
  };

  /** Component statistics */
  const AliHLTComponentDataType kAliHLTDataTypeComponentStatisticsTest = {
    sizeof(AliHLTComponentDataType),
    {'C','O','M','P','S','T','A','T'},
    {'P','R','I','V'}
  };

  /** Component table */
  const AliHLTComponentDataType kAliHLTDataTypeComponentTableTest = {
    sizeof(AliHLTComponentDataType),
    {'C','O','M','P','T','A','B','L'},
    {'P','R','I','V'}
  };

  /** general ROOT TObject */
  const AliHLTComponentDataType kAliHLTDataTypeTObjectTest = {
    sizeof(AliHLTComponentDataType),
    {'R','O','O','T','T','O','B','J'},
    {'*','*','*','\0'}
  };

  /** ROOT TObjArray */
  const AliHLTComponentDataType kAliHLTDataTypeTObjArrayTest = {
    sizeof(AliHLTComponentDataType),
    {'R','O','O','T','O','B','A','R'},
    {'*','*','*','\0'}
  };

  /** ROOT TTree */
  const AliHLTComponentDataType kAliHLTDataTypeTTreeTest = {
    sizeof(AliHLTComponentDataType),
    {'R','O','O','T','T','R','E','E'},
    {'*','*','*','\0'}
  };

  /** ROOT TH1 (can be used for all histograms, they derive from TH1) */
  const AliHLTComponentDataType kAliHLTDataTypeHistogramTest = {
    sizeof(AliHLTComponentDataType),
    {'R','O','O','T','H','I','S','T'},
    {'*','*','*','\0'}
  };

  /** ROOT TNtuple */
  const AliHLTComponentDataType kAliHLTDataTypeTNtupleTest = {
    sizeof(AliHLTComponentDataType),
    {'R','O','O','T','T','U','P','L'},
    {'*','*','*','\0'}
  };

  ////////////////////////


  if (kAliHLTMultipleDataTypeTest!=kAliHLTMultipleDataType) {
    cerr << "missmatch comparing kAliHLTMultipleDataType: ";
    cerr << AliHLTComponent::DataType2Text(kAliHLTMultipleDataType).c_str() << " vs. ";
    cerr << AliHLTComponent::DataType2Text(kAliHLTMultipleDataTypeTest).c_str() << endl;
    return -1;
  }

  if (kAliHLTDataTypeFXSCalibTest!=kAliHLTDataTypeFXSCalib) {
    cerr << "missmatch comparing kAliHLTDataTypeFXSCalib: ";
    cerr << AliHLTComponent::DataType2Text(kAliHLTDataTypeFXSCalib).c_str() << " vs. ";
    cerr << AliHLTComponent::DataType2Text(kAliHLTDataTypeFXSCalibTest).c_str() << endl;
    return -1;
  }

  if (kAliHLTDataTypeDDLTest!=kAliHLTDataTypeDDL) {
    cerr << "missmatch comparing kAliHLTDataTypeDDL: ";
    cerr << AliHLTComponent::DataType2Text(kAliHLTDataTypeDDL).c_str() << " vs. ";
    cerr << AliHLTComponent::DataType2Text(kAliHLTDataTypeDDLTest).c_str() << endl;
    return -1;
  }

  if (kAliHLTDataTypeSORTest!=kAliHLTDataTypeSOR) {
    cerr << "missmatch comparing kAliHLTDataTypeSOR: ";
    cerr << AliHLTComponent::DataType2Text(kAliHLTDataTypeSOR).c_str() << " vs. ";
    cerr << AliHLTComponent::DataType2Text(kAliHLTDataTypeSORTest).c_str() << endl;
    return -1;
  }

  if (kAliHLTDataTypeEORTest!=kAliHLTDataTypeEOR) {
    cerr << "missmatch comparing kAliHLTDataTypeEOR: ";
    cerr << AliHLTComponent::DataType2Text(kAliHLTDataTypeEOR).c_str() << " vs. ";
    cerr << AliHLTComponent::DataType2Text(kAliHLTDataTypeEORTest).c_str() << endl;
    return -1;
  }

  if (kAliHLTDataTypeRunTypeTest!=kAliHLTDataTypeRunType) {
    cerr << "missmatch comparing kAliHLTDataTypeRunType: ";
    cerr << AliHLTComponent::DataType2Text(kAliHLTDataTypeRunType).c_str() << " vs. ";
    cerr << AliHLTComponent::DataType2Text(kAliHLTDataTypeRunTypeTest).c_str() << endl;
    return -1;
  }

  if (kAliHLTDataTypeEventTest!=kAliHLTDataTypeEvent) {
    cerr << "missmatch comparing kAliHLTDataTypeEvent: ";
    cerr << AliHLTComponent::DataType2Text(kAliHLTDataTypeEvent).c_str() << " vs. ";
    cerr << AliHLTComponent::DataType2Text(kAliHLTDataTypeEventTest).c_str() << endl;
    return -1;
  }

  if (kAliHLTDataTypeComConfTest!=kAliHLTDataTypeComConf) {
    cerr << "missmatch comparing kAliHLTDataTypeComConf: ";
    cerr << AliHLTComponent::DataType2Text(kAliHLTDataTypeComConf).c_str() << " vs. ";
    cerr << AliHLTComponent::DataType2Text(kAliHLTDataTypeComConfTest).c_str() << endl;
    return -1;
  }

  if (kAliHLTDataTypeUpdtDCSTest!=kAliHLTDataTypeUpdtDCS) {
    cerr << "missmatch comparing kAliHLTDataTypeUpdtDCS: ";
    cerr << AliHLTComponent::DataType2Text(kAliHLTDataTypeUpdtDCS).c_str() << " vs. ";
    cerr << AliHLTComponent::DataType2Text(kAliHLTDataTypeUpdtDCSTest).c_str() << endl;
    return -1;
  }

  /*, data publisher will set type id and origin correctly */
  if (kAliHLTDataTypeDDLRawTest!=kAliHLTDataTypeDDLRaw) {
    cerr << "missmatch comparing kAliHLTDataTypeDDLRaw: ";
    cerr << AliHLTComponent::DataType2Text(kAliHLTDataTypeDDLRaw).c_str() << " vs. ";
    cerr << AliHLTComponent::DataType2Text(kAliHLTDataTypeDDLRawTest).c_str() << endl;
    return -1;
  }

  if (kAliHLTDataTypeESDObjectTest!=kAliHLTDataTypeESDObject) {
    cerr << "missmatch comparing kAliHLTDataTypeESDObject: ";
    cerr << AliHLTComponent::DataType2Text(kAliHLTDataTypeESDObject).c_str() << " vs. ";
    cerr << AliHLTComponent::DataType2Text(kAliHLTDataTypeESDObjectTest).c_str() << endl;
    return -1;
  }

  if (kAliHLTDataTypeESDTreeTest!=kAliHLTDataTypeESDTree) {
    cerr << "missmatch comparing kAliHLTDataTypeESDTree: ";
    cerr << AliHLTComponent::DataType2Text(kAliHLTDataTypeESDTree).c_str() << " vs. ";
    cerr << AliHLTComponent::DataType2Text(kAliHLTDataTypeESDTreeTest).c_str() << endl;
    return -1;
  }

  if (kAliHLTDataTypeAliTreeDTest!=kAliHLTDataTypeAliTreeD) {
    cerr << "missmatch comparing kAliHLTDataTypeAliTreeD: ";
    cerr << AliHLTComponent::DataType2Text(kAliHLTDataTypeAliTreeD).c_str() << " vs. ";
    cerr << AliHLTComponent::DataType2Text(kAliHLTDataTypeAliTreeDTest).c_str() << endl;
    return -1;
  }

  if (kAliHLTDataTypeAliTreeRTest!=kAliHLTDataTypeAliTreeR) {
    cerr << "missmatch comparing kAliHLTDataTypeAliTreeR: ";
    cerr << AliHLTComponent::DataType2Text(kAliHLTDataTypeAliTreeR).c_str() << " vs. ";
    cerr << AliHLTComponent::DataType2Text(kAliHLTDataTypeAliTreeRTest).c_str() << endl;
    return -1;
  }

  /*, origin is 'any' */
  if (kAliHLTDataTypeHwAddr16Test!=kAliHLTDataTypeHwAddr16) {
    cerr << "missmatch comparing kAliHLTDataTypeHwAddr16: ";
    cerr << AliHLTComponent::DataType2Text(kAliHLTDataTypeHwAddr16).c_str() << " vs. ";
    cerr << AliHLTComponent::DataType2Text(kAliHLTDataTypeHwAddr16Test).c_str() << endl;
    return -1;
  }

  if (kAliHLTDataTypeEventStatisticsTest!=kAliHLTDataTypeEventStatistics) {
    cerr << "missmatch comparing kAliHLTDataTypeEventStatistics: ";
    cerr << AliHLTComponent::DataType2Text(kAliHLTDataTypeEventStatistics).c_str() << " vs. ";
    cerr << AliHLTComponent::DataType2Text(kAliHLTDataTypeEventStatisticsTest).c_str() << endl;
    return -1;
  }

  if (kAliHLTDataTypeEventSummaryTest!=kAliHLTDataTypeEventSummary) {
    cerr << "missmatch comparing kAliHLTDataTypeEventSummary: ";
    cerr << AliHLTComponent::DataType2Text(kAliHLTDataTypeEventSummary).c_str() << " vs. ";
    cerr << AliHLTComponent::DataType2Text(kAliHLTDataTypeEventSummaryTest).c_str() << endl;
    return -1;
  }

  if (kAliHLTDataTypeRunStatisticsTest!=kAliHLTDataTypeRunStatistics) {
    cerr << "missmatch comparing kAliHLTDataTypeRunStatistics: ";
    cerr << AliHLTComponent::DataType2Text(kAliHLTDataTypeRunStatistics).c_str() << " vs. ";
    cerr << AliHLTComponent::DataType2Text(kAliHLTDataTypeRunStatisticsTest).c_str() << endl;
    return -1;
  }

  if (kAliHLTDataTypeRunSummaryTest!=kAliHLTDataTypeRunSummary) {
    cerr << "missmatch comparing kAliHLTDataTypeRunSummary: ";
    cerr << AliHLTComponent::DataType2Text(kAliHLTDataTypeRunSummary).c_str() << " vs. ";
    cerr << AliHLTComponent::DataType2Text(kAliHLTDataTypeRunSummaryTest).c_str() << endl;
    return -1;
  }

  if (kAliHLTDataTypeComponentStatisticsTest!=kAliHLTDataTypeComponentStatistics) {
    cerr << "missmatch comparing kAliHLTDataTypeComponentStatistics: ";
    cerr << AliHLTComponent::DataType2Text(kAliHLTDataTypeComponentStatistics).c_str() << " vs. ";
    cerr << AliHLTComponent::DataType2Text(kAliHLTDataTypeComponentStatisticsTest).c_str() << endl;
    return -1;
  }

  if (kAliHLTDataTypeComponentTableTest!=kAliHLTDataTypeComponentTable) {
    cerr << "missmatch comparing kAliHLTDataTypeComponentTable: ";
    cerr << AliHLTComponent::DataType2Text(kAliHLTDataTypeComponentTable).c_str() << " vs. ";
    cerr << AliHLTComponent::DataType2Text(kAliHLTDataTypeComponentTableTest).c_str() << endl;
    return -1;
  }

  if (kAliHLTDataTypeTObjectTest!=kAliHLTDataTypeTObject) {
    cerr << "missmatch comparing kAliHLTDataTypeTObject: ";
    cerr << AliHLTComponent::DataType2Text(kAliHLTDataTypeTObject).c_str() << " vs. ";
    cerr << AliHLTComponent::DataType2Text(kAliHLTDataTypeTObjectTest).c_str() << endl;
    return -1;
  }

  if (kAliHLTDataTypeTObjArrayTest!=kAliHLTDataTypeTObjArray) {
    cerr << "missmatch comparing kAliHLTDataTypeTObjArray: ";
    cerr << AliHLTComponent::DataType2Text(kAliHLTDataTypeTObjArray).c_str() << " vs. ";
    cerr << AliHLTComponent::DataType2Text(kAliHLTDataTypeTObjArrayTest).c_str() << endl;
    return -1;
  }

  if (kAliHLTDataTypeTTreeTest!=kAliHLTDataTypeTTree) {
    cerr << "missmatch comparing kAliHLTDataTypeTTree: ";
    cerr << AliHLTComponent::DataType2Text(kAliHLTDataTypeTTree).c_str() << " vs. ";
    cerr << AliHLTComponent::DataType2Text(kAliHLTDataTypeTTreeTest).c_str() << endl;
    return -1;
  }

  if (kAliHLTDataTypeHistogramTest!=kAliHLTDataTypeHistogram) {
    cerr << "missmatch comparing kAliHLTDataTypeHistogram: ";
    cerr << AliHLTComponent::DataType2Text(kAliHLTDataTypeHistogram).c_str() << " vs. ";
    cerr << AliHLTComponent::DataType2Text(kAliHLTDataTypeHistogramTest).c_str() << endl;
    return -1;
  }

  if (kAliHLTDataTypeTNtupleTest!=kAliHLTDataTypeTNtuple) {
    cerr << "missmatch comparing kAliHLTDataTypeTNtuple: ";
    cerr << AliHLTComponent::DataType2Text(kAliHLTDataTypeTNtuple).c_str() << " vs. ";
    cerr << AliHLTComponent::DataType2Text(kAliHLTDataTypeTNtupleTest).c_str() << endl;
    return -1;
  }

  return 0;
}
