/**************************************************************************
 * This file is property of and copyright by the Experimental Nuclear     *
 * Physics Group, Dep. of Physics                                         *
 * University of Oslo, Norway, 2007                                       *
 *                                                                        * 
 * Author: Per Thomas Hille perthi@fys.uio.no for the ALICE HLT Project.  *
 * Contributors are mentioned in the code where appropriate.              *
 * Please report bugs to perthi@fys.uio.no                                * 
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/


#include "AliHLTPHOSOnlineDisplay.h"
#ifndef __CINT__
# include <stdexcept>
# include <TSystem.h>
# include <TApplication.h>
# include "Rtypes.h"
#include <TString.h>
#endif




int
main(int argc, char** argv) 
{
  Bool_t hostnameIsSet = kFALSE;
  Bool_t portIsSet = kFALSE;

  int iResult=0;
  TString argument="";
  //  TString    hostname;  
  char *hostname;
  //  char *port;
  int port;

  int bMissingParam=0;
  printf("Main: the number of argumnets is %d \n", argc);
  for (int i=0; i<argc && iResult>=0; i++) 
    {
      argument=argv[i];
      if (argument.IsNull()) 
	{
	  continue;
	}

      if (argument.CompareTo("-hostname")==0) 
	{
	  if ((bMissingParam=(++i>=argc))) 
	    {
	      break;
	    }
	  hostname = argv[i];
	  hostnameIsSet = kTRUE;
	  printf("\nsetting host to: %s \n", hostname);
	} 
	  
      if (argument.CompareTo("-port")==0) 
	{
	  if ((bMissingParam=(++i>=argc))) 
	    {
	      break;
	    }
	  port = atoi(argv[i]);
	  portIsSet = kTRUE;
	  printf("\nsetting port to: %d \n", port);
	}
   }

  if(hostnameIsSet != kTRUE ||  portIsSet != kTRUE)
    {
      if(hostnameIsSet == kFALSE)
	{
	  printf("\nERROR: no hostname is specified\n");
	}
      
      if( portIsSet == kFALSE)
	{
	  printf("ERROR: no port spcified\n");
	}
      printf("\nYou must specify hostname & port as command line arguments\n\n");
      printf("*****************************************************************\n");
      printf("\nUsage: ./onlinedisplay  -hostname  <hostname>   -port  <port>\n\n");
      printf("*****************************************************************\n\n\n");
    }
  else
    {


      try {
	TApplication app("app", 0, 0);
	AliHLTPHOSOnlineDisplay* phosDisplayPtr = AliHLTPHOSOnlineDisplay::Instance(hostname, port);
	app.Run();
      }
  
      catch (std::exception& e) {
	//   std::cerr << e.what() << std::endl;
	return 1;
      }
      return 0;
      
    }

}


