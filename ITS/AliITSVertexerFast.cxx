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
  AliRunLoader *rl =AliRunLoader::Instance();
  TTree *trK=(TTree*)rl->TreeK();
  if(!trK)AliFatal("This class should be used only with simulated events!!");
  rl->LoadHeader(); 
}

//______________________________________________________________________
AliITSVertexerFast::AliITSVertexerFast(Double_t *smear):AliITSVertexer(),
fSmear(0)
{
  // Standard constructor
  fSmear = new Double_t[3];
  for(Int_t i=0;i<3;i++)fSmear[i]=smear[i];
  Info("AliITSVertexerFast","Gaussian smaring of the generated vertex. Parameters %f12.5 , %f12.5 , %f12.5 \n",fSmear[0],fSmear[1],fSmear[2]);
  AliRunLoader *rl =AliRunLoader::Instance();
  TTree *trK=(TTree*)rl->TreeK();
  if(!trK)AliFatal("This class should be used only with simulated events!!");
  rl->LoadHeader(); 

}


//______________________________________________________________________
AliITSVertexerFast::~AliITSVertexerFast(){
  // Destructor
  if(fSmear)delete [] fSmear;
  fSmear = 0;
}

//______________________________________________________________________
AliESDVertex* AliITSVertexerFast::FindVertexForCurrentEvent(TTree *itsClusterTree){
  // Defines the AliITSVertex for the current event
  AliWarning(Form("This class should be used only with simulated events!! Input cluster tree (%p) will not be used!!",itsClusterTree));

  fCurrentVertex = 0;
  AliRunLoader *rl =AliRunLoader::Instance();
  TArrayF primaryVertex(3);  // true vertex
  AliHeader* header = rl->GetHeader();
  AliGenEventHeader* genEventHeader = header->GenEventHeader();   
  genEventHeader->PrimaryVertex(primaryVertex); 
  
  // Smearing
  Double_t vrttrue[3],vrtx[3];
  for(Int_t k=0; k<3;k++){
    vrttrue[k] = static_cast<Double_t>(primaryVertex[k]);
    vrtx[k] = gRandom->Gaus(vrttrue[k],fSmear[k]);
  }
  fCurrentVertex = new AliESDVertex(vrtx,fSmear,"Smeared Generated Vertex");
  return fCurrentVertex;
  
}

//________________________________________________________
void AliITSVertexerFast::PrintStatus() const {
  // Print current status
  cout <<"=======================================================\n";

  cout<<"RMS for gaussian smearing: ";
  for(Int_t k=0;k<3;k++)cout<<" "<<fSmear[k];
  cout<<endl;
}

