#include "AliStorageAdministratorPanel.h"
#include <TApplication.h>


int main(int argc,char **argv)
{
	TApplication app("AliStorageAdmin", &argc, argv);
	AliStorageAdministratorPanel *mainWindow = new AliStorageAdministratorPanel();
	app.Run(kTRUE);
	if(mainWindow){delete mainWindow;}

	
	return 0;
}
