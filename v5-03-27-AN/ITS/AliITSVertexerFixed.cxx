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
#include "AliESDVertex.h"
#include "AliLog.h"
#include <TString.h>
#include "AliITSVertexerFixed.h"

/////////////////////////////////////////////////////////////////////////
//                                                                     //
// Fixed vertexer - creates a vertex in a defined postion (x,y,z)      //
// the standard contructor takes a sting to specify the case           //
// Useful for reconstruction of injection tests with beam on TDI       //
//                                                                     //
/////////////////////////////////////////////////////////////////////////

ClassImp(AliITSVertexerFixed)

/* $Id$ */

//______________________________________________________________________
AliITSVertexerFixed::AliITSVertexerFixed():AliITSVertexer()
{
  // Default Constructor
  AliWarning("This contructor sets the vertex in (0,0,0)");
  for(Int_t k=0; k<3;k++){ 
    fVtxPos[k]=0.;
    fVtxErr[k]=0.5;  
  }
}

//______________________________________________________________________
AliITSVertexerFixed::AliITSVertexerFixed(TString option):AliITSVertexer()
{
  // Standard constructor
  if(option.Contains("TDI")){
    fVtxPos[0]=0.;
    fVtxPos[1]=0.;
    fVtxPos[2]=8000.;  // TDI at z=80 m
    fVtxErr[0]=1.;
    fVtxErr[1]=1.;
    fVtxErr[2]=100.;   
   }
  else if(option.Contains("TED")){
    fVtxPos[0]=0.;
    fVtxPos[1]=0.;
    fVtxPos[2]=34000.;  // TED at z=+340 m
    fVtxErr[0]=1.;
    fVtxErr[1]=1.;
    fVtxErr[2]=100.;   
  }else{
    AliError(Form("%s is invalid, sets the vertex in (0,0,0)",option.Data()));
    for(Int_t k=0; k<3;k++){ 
      fVtxPos[k]=0.;
      fVtxErr[k]=0.5;  
    }
  }
}


//______________________________________________________________________
AliESDVertex* AliITSVertexerFixed::FindVertexForCurrentEvent(TTree * /*itsClusterTree */){
  // Defines the AliITSVertex for the current event
  
  fCurrentVertex = new AliESDVertex(fVtxPos,fVtxErr,"Fixed Vertex");
  return fCurrentVertex;
  
}

//________________________________________________________
void AliITSVertexerFixed::PrintStatus() const {
  // Print current status
  cout <<"=======================================================\n";

  cout<<"Fixed positions: ";
  for(Int_t k=0;k<3;k++)cout<<" "<<fVtxPos[k]<<"+-"<<fVtxErr[k];
  cout<<endl;
}

