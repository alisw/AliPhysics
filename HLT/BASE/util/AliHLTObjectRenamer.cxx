#include "AliHLTObjectRenamer.h"
#include "AliHLTOUT.h"
#include "AliHLTMisc.h"
#include "AliHLTSystem.h"
#include "AliHLTComponentHandler.h"
#include "AliHLTPluginBase.h"
#include "AliHLTConfiguration.h"
#include "AliHLTHOMERLibManager.h"
#include "AliHLTHOMERBlockDesc.h"
#include "AliHLTDataSource.h"
#include "AliHLTProcessor.h"
#include "AliHLTDataSink.h"
#include "AliHLTMessage.h"
#include "AliHLTGlobalTriggerDecision.h"
#include "AliHLTCTPData.h"
#include "AliCDBManager.h"
#include "AliDAQ.h"
#include "TServerSocket.h"
#include "TObjArray.h"
#include "TSystem.h"
#include "TTimeStamp.h"
#include <vector>
#include <iostream>
#include <iomanip>
using namespace std;

#include <sys/ioctl.h>
#include <arpa/inet.h>

ClassImp(AliHLTObjectRenamer)

int AliHLTObjectRenamer::DoInit(int argc, const char** argv)
{
	fSuffix = "";
	bool suffixSet = false;
	for (int i = 0; i < argc; i++)
	{
		if (strcmp( argv[i], "-suffix" ) == 0)
		{
			if (suffixSet)
			{
				HLTWarning("The suffix has already been specified."
					" Will replace previous value given by -suffix."
				);
			}
			if ( argc <= i+1 )
			{
				HLTError("The suffix value was not specified.");
				return -EINVAL;
			}
			fSuffix = argv[i+1];
			i++;
			continue;
		}
		
		HLTError("Unknown option '%s'.", argv[i]);
		return -EINVAL;
	}
	
	return 0;
}

int AliHLTObjectRenamer::DoEvent(const AliHLTComponentEventData& /*evtData*/, AliHLTComponentTriggerData& /*trigData*/)
{
	const TObject* obj = NULL;
	for (obj = GetFirstInputObject(); obj != NULL; obj = GetNextInputObject())
	{
		TObject* newObj = NULL; CloneList<TObjArray>(obj, fSuffix);
		if (newObj == NULL) newObj = CloneList<TList>(obj, fSuffix);
		if (newObj == NULL)
		{
			TString name = obj->GetName();
			name += fSuffix;
			newObj = obj->Clone(name);
		}
		PushBack(newObj, GetDataType(), GetSpecification());
		delete newObj;
	}
	return 0;
}

