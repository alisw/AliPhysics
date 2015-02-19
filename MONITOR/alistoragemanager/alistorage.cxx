#include "AliStorageClientThread.h"
#include "AliStorageServerThread.h"

#include <TThread.h>
#include <iostream>
#include <string.h>
#include <sstream>
#include <cstdlib>

using namespace std;

void *ClientThreadHandle(void*)
{
	cout<<"ALICE Storage Manager -- Starting client thread"<<endl;
	AliStorageClientThread *client = new AliStorageClientThread();
	if(client){delete client;}
	return 0;
}

void *ServerThreadHandle(void*)
{
	cout<<"\nALICE Storage Manager -- Starting server thread"<<endl;
	AliStorageServerThread *server = new AliStorageServerThread();
	if(server){delete server;}
	return 0;
}

bool isStorageRunning()
{
 // check if there is events server already running

  // first mathod with pidof
  const char *pid = gSystem->GetFromPipe("pidof alistorage").Data();
  int pidSize = gSystem->GetFromPipe("pidof alistorage").Sizeof();
  std::string pidOfAll(pid,pidSize);
  std::stringstream pidStream(pidOfAll);
  int word_count=0; 
  std::string word;
  while( pidStream >> word ) ++word_count;
  if(word_count != 1){return true;}
  
  // second method with ps
  cout<<"checking if alistorage is running with ps"<<endl;
  TString psName = gSystem->GetFromPipe("ps -e -o comm | grep alistorage");
  if(strcmp(psName.Data(),"alistorage")>0){return true;}

  return false;
}

int main()
{
  if(isStorageRunning())
    {
      std::cout<<"There is other Storage Manager running on this machine.\n Cannot start multiple managers on the same machine. Quitting..."<<std::endl;
      return 0;
    }


	TThread *clientThread = new TThread("clientThread", ClientThreadHandle,NULL);
	TThread *serverThread = new TThread("serverThread", ServerThreadHandle,NULL);
    
	clientThread->Run();
	serverThread->Run();
	
	clientThread->Join();
	serverThread->Kill();//if client thread if finished, server thread is killed
	return 0;
}
