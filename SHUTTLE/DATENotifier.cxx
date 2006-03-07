#include "DATENotifier.h"

#include "AliLog.h"

#include "AliShuttleTrigger.h"

void DATENotifier::infoHandler() {
	AliInfoGeneral("DATENotifier::infoHandler()",
			"DATE notification received ...");
	fTrigger->Notify();
}

void DATENotifier::errorHandler(int severity, int code, char *msg) {

	AliInfoGeneral("DATENotifier::errorHandler()",
		Form("DIM Error: severity<%d>, code<%d> , message<%s>",
		severity, code, msg));
}
