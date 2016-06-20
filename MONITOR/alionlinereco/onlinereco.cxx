// Author: Mihai Niculescu 2013

/**************************************************************************
 * Copyright(c) 1998-2013, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliOnlineReconstruction.h"

#include <TSystem.h>

#include <iostream>
#include <string.h>
#include <sstream>
#include <cstdlib>

bool isReconstructionRunning()
{
    // check if there is events server already running
    const char *pid = gSystem->GetFromPipe("pidof alionlinereco").Data();
    int pidSize = gSystem->GetFromPipe("pidof alionlinereco").Sizeof();
    std::string pidOfAll(pid,pidSize);
    std::stringstream pidStream(pidOfAll);
    int word_count=0;
    std::string word;
    while( pidStream >> word ) ++word_count;
    if(word_count != 1){return true;}
    return false;
}

int main(int argc, char **argv)
{
    if(argc!=2)
    {
        std::cout<<"Start Online Reconstruction with run number as a parameter"<<std::endl;
        return 0;
    }
    if(isReconstructionRunning())
    {
        std::cout<<"There are other servers. Cannot start multiple servers on the same machine. Quitting..."<<std::endl;
        return 0;
    }
    if(atoi(argv[1])<=0)
    {
        std::cout<<"Incorrect run number"<<std::endl;
        return 0;
    }
    
    gSystem->Setenv("DATE_RUN_NUMBER",Form("%d",atoi(argv[1])));
    
    std::cout<<"Starting Online Reconstruction for run:"<<atoi(argv[1])<<std::endl;
    AliOnlineReconstruction *onlineReconstruction = new AliOnlineReconstruction(atoi(argv[1]));
    if(onlineReconstruction){delete onlineReconstruction;}
    
    return 0;
}
