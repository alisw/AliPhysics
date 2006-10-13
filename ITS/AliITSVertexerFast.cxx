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
#include <Riostream.h>
#include <TArrayF.h>
#include <TRandom.h>
#include "AliESDVertex.h"
#include <AliITSVertexerFast.h>
#include "AliHeader.h"
#include "AliGenEventHeader.h"
#include "AliRun.h"
#include "AliITSLoader.h"
#include "AliRunLoader.h"

/////////////////////////////////////////////////////////////////////////
//                                                                     //
// Fast vertexer - True (i.e. generated) vertex coordinates            //
// are smeared with gaussians of given width                           //
// Origin: masera@to.infn.it     25/09/2003                            //
//                                                                     //
/////////////////////////////////////////////////////////////////////////
ClassImp(AliITSVertexerFast)



//______________________________________________________________________
AliITSVertexerFast::AliITSVertexerFast():AliITSVertexer(),
fSmear(0) 
{
  // Default Constructor
  fSmear = 0;
}

//______________________________________________________________________
AliITSVertexerFast::AliITSVertexerFast(Double_t *smear):AliITSVertexer(),
fSmear(0)
{
  // Standard constructor
  fSmear = new Double_t[3];
  for(Int_t i=0;i<3;i++)fSmear[i]=smear[i];
  Info("AliITSVertexerFast","Gaussian smaring of the generated vertex. Parameters %f12.5 , %f12.5 , %f12.5 \n",fSmear[0],fSmear[1],fSmear[2]);
}

//______________________________________________________________________
AliITSVertexerFast::~AliITSVertexerFast(){
  // Destructor
  if(fSmear)delete [] fSmear;
  fSmear = 0;
}

//______________________________________________________________________
AliESDVertex* AliITSVertexerFast::FindVertexForCurrentEvent(Int_t evnumb){
  // Defines the AliITSVertex for the current event
  fCurrentVertex = 0;
  AliRunLoader *rl =AliRunLoader::GetRunLoader();
  rl->GetEvent(evnumb);
  TArrayF primaryVertex(3);  // true vertex
  AliHeader* header = gAlice->GetHeader();
  AliGenEventHeader* genEventHeader = header->GenEventHeader();   
  genEventHeader->PrimaryVertex(primaryVertex); 

  // Smearing
  Double_t vrttrue[3],vrtx[3];
  for(Int_t k=0; k<3;k++){
    vrttrue[k] = static_cast<Double_t>(primaryVertex[k]);
    vrtx[k] = gRandom->Gaus(vrttrue[k],fSmear[k]);
  }
  char name[30];
  sprintf(name,"Vertex_%d",evnumb);
  fCurrentVertex = new AliESDVertex(vrtx,fSmear,name);
  fCurrentVertex->SetTruePos(vrttrue);
  return fCurrentVertex;
}

//______________________________________________________________________
void AliITSVertexerFast::FindVertices(){
  // computes the vertices of the events in the range FirstEvent - LastEvent

  AliRunLoader *rl = AliRunLoader::GetRunLoader();
  AliITSLoader* iTSloader =  (AliITSLoader*) rl->GetLoader("ITSLoader");
  iTSloader->ReloadRecPoints();
  for(Int_t i=fFirstEvent;i<=fLastEvent;i++){
    rl->GetEvent(i);
    FindVertexForCurrentEvent(i);   
    if(fCurrentVertex) WriteCurrentVertex();
    else {
      cout<<"Vertex not found for event "<<i<<endl;

    }

  }

}

//________________________________________________________
void AliITSVertexerFast::PrintStatus() const {
  // Print current status
  cout <<"=======================================================\n";

  cout<<"First event to be processed "<<fFirstEvent;
  cout<<"\n Last event to be processed "<<fLastEvent<<endl;
  cout<<"RMS for gaussian smearing: ";
  for(Int_t k=0;k<3;k++)cout<<" "<<fSmear[k];
  cout<<endl;
}

//______________________________________________________________________
AliITSVertexerFast::AliITSVertexerFast(const AliITSVertexerFast &vtxr) : 
                    AliITSVertexer(vtxr) {
  // Copy constructor
  // Copies are not allowed. The method is protected to avoid misuse.
  Error("AliITSVertexerFast","Copy constructor not allowed\n");
}

//______________________________________________________________________
AliITSVertexerFast& AliITSVertexerFast::operator=(const 
                    AliITSVertexerFast& /* vtxr */){
  // Assignment operator
  // Assignment is not allowed. The method is protected to avoid misuse.
  Error("= operator","Assignment operator not allowed\n");
  return *this;
}
