#ifndef __ALIONLINERECONSTRUCTION_UTIL_H__
#define __ALIONLINERECONSTRUCTION_UTIL_H__

// Default configuration file
#define ALIEVENTSERVER_CONF "onlinereco.conf"

//______________DEFAULT SETTINGS________________
#define DEFAULT_SERVER_HOST "tcp://*"
#define DEFAULT_SERVER_PORT 5024
#define DEFAULT_SERVER_SAVE_RECO_DIR "/local/reco"
#define DEFAULT_CDB_STORAGE "local:///local/cdb"
#define DEFAULT_CDB_SPEC_STORAGE_PATH1 "GRP/GRP/Data"
#define DEFAULT_CDB_SPEC_STORAGE_VALUE1 "/local/reco/"
#define DEFAULT_CDB_SPEC_STORAGE_PATH2 "GRP/CTP/Config"
#define DEFAULT_CDB_SPEC_STORAGE_VALUE2 "/local/reco/"
#define DEFAULT_CDB_SPEC_STORAGE_PATH3 ""
#define DEFAULT_CDB_SPEC_STORAGE_VALUE3 ""
#define DEFAULT_QA_RUN ":"
#define DEFAULT_QAREF_STORAGE "local://$ALICE_ROOT/QAref"
#define DEFAULT_QA_RUN_GLOBAL 1
#define DEFAULT_RECO_RUN_PLANE_EFF 1
#define DEFAULT_RECO_WRITE_ESDF 0
#define DEFAULT_RECO_WRITE_ALIGN 1
#define DEFAULT_RECO_CLEAN_ESD 0
#define DEFAULT_RECO_DETECTORS "ALL -PHOS -EMCAL"
#define DEFAULT_LOGBOOK_HOST "host"
#define DEFAULT_LOGBOOK_PORT 3306
#define DEFAULT_LOGBOOK_DB "database"
#define DEFAULT_LOGBOOK_USER "user"
#define DEFAULT_LOGBOOK_PASS "pass123"
#define DEFAULT_DATA_SOURCE "local"
#define DEFAULT_DATA_ONLINE_SOURCE "mem@*:"

#include <TSystem.h>
#include <TString.h>

namespace AliOnlineReconstructionUtil
{
// return full path to the server configuration file
	inline const char* GetPathToServerConf()
	{
		return Form("%s/%s",
			    gSystem->Getenv("HOME"),
			    ALIEVENTSERVER_CONF);
	}
}

#endif
