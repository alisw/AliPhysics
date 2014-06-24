#include <TSystem.h>
#include <TString.h>

#include "AliEventServerUtil.h"

const char* AliEventServerUtil::GetPathToServerConf()
{
	return Form("%s/MONITOR/%s", gSystem->Getenv("ALICE_ROOT"), ALIEVENTSERVER_CONF);
}
