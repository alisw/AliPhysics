/**************************************************************************
 * Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
/////////////////////////////////////////////////////////////////////////
// ALICE HLT DISPLAY CLASS                                             //
// Author: Mayeul   ROUSSELET                                          //
// e-mail: Mayeul.Rousselet@cern.ch                                    //
// Last update:26/08/2003                                              //
/////////////////////////////////////////////////////////////////////////

#include <TPolyMarker3D.h>

#include "AliDisplayHLT.h"
#include "AliDisplay2.h"
#include "AliModuleInfo.h"
#ifdef ALI_HLT
#include <stdlib.h>
#include "AliL3MemHandler.h"
#include "AliL3SpacePointData.h"
#include "AliL3Transform.h"
#include "AliL3Logging.h"
#endif

ClassImp(AliDisplayHLT)

//_____________________________________________________________
AliDisplayHLT::AliDisplayHLT()
{
  // Default constructor
  fPoints = new TPolyMarker3D[gAliDisplay2->GetNbModules()];
  fName = new char*[gAliDisplay2->GetNbModules()];
  fNb=0;
  for(Int_t i=0;i<gAliDisplay2->GetNbModules();i++){
    fPoints[i].SetMarkerSize(0.2); 
    fPoints[i].SetMarkerColor(2); 
    fPoints[i].SetMarkerStyle(1);
  }
}

//_____________________________________________________________
AliDisplayHLT::~AliDisplayHLT()
{
  delete [] fPoints;
}

//_____________________________________________________________
void AliDisplayHLT::LoadHLT(const char *name,Int_t nevent)
{
  // Loads HLT clusters of TPC
  if(strstr(name,"TPC")) LoadHLTTPC(nevent);
}

//_____________________________________________________________
void AliDisplayHLT::LoadHLTTPC(Int_t nevent)
{
  //load TPC Clusters from the raw data
  //raw data must be in the directorie $ALICE_ROOT/raw
  //First we read the data from files
  //do_mc MUST BE DEFINED AND USED FOR RAW DATA GENERATION
#ifdef ALI_HLT

  fName[fNb]=new char[strlen("TPC")];
  strcpy(fName[fNb],"TPC");
  Char_t fname[256];
  AliL3MemHandler *clusterfile[36][6];
  AliL3SpacePointData *fClusters[36][6];
  UInt_t fNcl[36][6];
  memset(fClusters,0,36*6*sizeof(AliL3SpacePointData*));
  //  strcpy(path,gSystem->Getenv("ALICE_ROOT"));
  //strcat(path,"/raw");
  //printf("\nRaw data path %s",path);
  char path[128];
  strcpy(path,gAliDisplay2->GetRawDataPath());
  for(Int_t s=0; s<36; s++)
    {
      for(Int_t p=0; p<AliL3Transform::GetNPatches(); p++)
	{
	  Int_t patch;
	  patch=-1;
	  clusterfile[s][p] = new AliL3MemHandler();
	  if(nevent<0)
	    sprintf(fname,"%s/points_%d_%d.raw",path,s,patch);
	  else
	    sprintf(fname,"%s/points_%d_%d_%d.raw",path,nevent,s,patch);
	  if(!clusterfile[s][p]->SetBinaryInput(fname))
	    {
	      LOG(AliL3Log::kError,"AliL3Evaluation::Setup","File Open")
		<<"Inputfile "<<fname<<" does not exist"<<ENDLOG; 
	      delete clusterfile[s][p];
              clusterfile[s][p] = 0; 
	      continue;
	    }
	  fClusters[s][p] = (AliL3SpacePointData*)clusterfile[s][p]->Allocate();
	  clusterfile[s][p]->Binary2Memory(fNcl[s][p],fClusters[s][p]);
	  clusterfile[s][p]->CloseBinaryInput();
	  break;
	}
    }

  //Second step: we assign the clusters to the fPoints array
  Int_t nbc=0;
  for(Int_t s=0; s<36; s++){
      for(Int_t p=0;p<6;p++){
	  AliL3SpacePointData *points = fClusters[s][p];
	  if(!points) continue;
	  Float_t xyz[3];
	  for(UInt_t i=0; i<fNcl[s][p]; i++){
	      xyz[0] = points[i].fX;
	      xyz[1] = points[i].fY;
	      xyz[2] = points[i].fZ;
	      fPoints[fNb].SetPoint(i+nbc,xyz[0],xyz[1],xyz[2]);
 	    }
	    nbc += fNcl[s][p];
	}
    }
  fNb++;
#else
  printf("This is event %d\n",nevent);
#endif
}

//_____________________________________________________________
void AliDisplayHLT::Draw()
{
  // Draws HLT clusters
  for(Int_t i=0;i<fNb;i++){
   if(gAliDisplay2->GetModuleInfo()->IsEnabled(fName[i])) fPoints[i].Draw();
  }
}

