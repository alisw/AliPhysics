#include "AliHLTPHOSProcessor.h"


const AliHLTComponentDataType AliHLTPHOSProcessor::fgkInputDataTypes[]={kAliHLTVoidDataType,{0,"",""}}; //'zero' terminated array


AliHLTPHOSProcessor::AliHLTPHOSProcessor():AliHLTProcessor(), fModuleID(0), fPrintInfoFrequncy(1000)
{

}


AliHLTPHOSProcessor::~AliHLTPHOSProcessor()
{

}


