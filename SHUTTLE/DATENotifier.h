#ifndef DATE_NOTIFIER_H
#define DATE_NOTIFIER_H

#include <dic.hxx>

class AliShuttleTrigger;

class DATENotifier: public DimInfo
{
	AliShuttleTrigger* fTrigger;
public:
        DATENotifier(AliShuttleTrigger* trigger, const char* service): 
		DimInfo(service, -1), fTrigger(trigger) {}

	void infoHandler();

	void errorHandler(int severity, int code, char *msg); 
};

#endif
