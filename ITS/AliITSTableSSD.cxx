/**************************************************************************
 * Copyright(c) 2002-2003, ALICE Experiment at CERN, All rights reserved. *
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
#include "AliITSTableSSD.h"

ClassImp(AliITSTableSSD)
////////////////////////////////////////////////////////////////////////
// Version: 0                                                         //
// Origin: Massimo Masera                                             //
// March 2002                                                         //
//                                                                    //
// AliITSTableSSD is used by AliITSsimulationSSD class to fill the AliITSpList
// object starting from the map with energy depositions
////////////////////////////////////////////////////////////////////////
//----------------------------------------------------------------------
AliITSTableSSD::AliITSTableSSD() : TObject(),
fDim(0),
fArray(0){
  // Default Constructor
  for(Int_t i=0;i<2;i++){
    fCurrUse[i]=0;
    fCurrFil[i]=0;
  }
}
//----------------------------------------------------------------------
AliITSTableSSD::AliITSTableSSD(const AliITSTableSSD & source):TObject(source),
fDim(source.fDim),
fArray(source.fArray){
  // Copy constructor

    if(this == &source) return;
    fArray = new Int_t [fDim];
    fCurrUse[0]=(source.fCurrUse)[0];
    fCurrUse[1]=(source.fCurrUse)[1];
    fCurrFil[0]=(source.fCurrFil)[0];
    fCurrFil[1]=(source.fCurrFil)[1];
    for(Int_t i=0;i<fDim;i++)fArray[i]=(source.fArray)[i];
}
//----------------------------------------------------------------------
AliITSTableSSD& AliITSTableSSD::operator=(const AliITSTableSSD & source){
  // = opporator constructor

    if(this == &source) return *this;
    fDim=source.fDim;
    fArray = new Int_t [fDim];
    fCurrUse[0]=(source.fCurrUse)[0];
    fCurrUse[1]=(source.fCurrUse)[1];
    fCurrFil[0]=(source.fCurrFil)[0];
    fCurrFil[1]=(source.fCurrFil)[1];
    for(Int_t i=0;i<fDim;i++)fArray[i]=(source.fArray)[i];
    return *this;
}
//----------------------------------------------------------------------
AliITSTableSSD::AliITSTableSSD(Int_t noelem) : TObject(),
fDim(0),
fArray(0){
  // Standard constructor
  fDim=noelem*2;
  fArray = new Int_t [fDim];
  Clear();
}
//----------------------------------------------------------------------
AliITSTableSSD::~AliITSTableSSD(){
  // Destructor
  delete [] fArray;
  fArray=0;
}
//----------------------------------------------------------------------
void AliITSTableSSD::Add(Int_t side,Int_t strip){
  // Add an element to the table
  if(strip>=fDim/2){
    cerr<<" Error in AliITSTableSSD::Add. "<<strip<<" is out of range\n";
    return;
  }
  if((side!=0) && (side!=1)){
    cerr<<" Error in AliITSTableSSD::Add. side="<<side<<endl;
    cerr<<" side must be 0 or 1\n";
    return;
  }
  if(fCurrFil[side]>(fDim/2-1)){
    cerr<<" Error in AliITSTableSSD::Add. Trying to fill an element out of range\n";
    cerr<<" Element="<<fCurrFil[side]<<" while max limit is "<<(fDim/2-1)<<endl;
    fCurrFil[side]++;
    return;
  }
  // P side = 0 and N side =1
  Int_t index=(fDim/2)*side+fCurrFil[side];
  Int_t * ptr= &fArray[(fDim/2)*side];
  if(SearchValue(ptr,strip,fCurrFil[side])>=0)return;
  fArray[index]=strip;
  fCurrFil[side]++;
}
//----------------------------------------------------------------------
void AliITSTableSSD::Clear(){
  //clear arrays
  fCurrUse[0]= 0;
  fCurrUse[1] = 0;
  fCurrFil[0]= 0;
  fCurrFil[1] = 0;
  for(Int_t i=0;i<fDim;i++)fArray[i]=-1;
}
//----------------------------------------------------------------------
void AliITSTableSSD::DumpTable(){
  // Dumps the contents of the table
  cout<<"==============================================================\n";
  cout<<" AliITSTableSSD::DumpTable \n";
  cout<<" Dimension of the table "<<fDim<<" ("<<fDim/2<<" per side)\n";
  cout<<" Current element to be filled for P side "<<fCurrFil[0]<<endl;
  cout<<" Current element to be filled for N side "<<fCurrFil[1]<<endl;
  cout<<" Current element in use for P side "<<fCurrUse[0]<<endl;
  cout<<" Current element in use for N side "<<fCurrUse[1]<<endl;
  cout<<"\n Elements for P side: \n";
  for(Int_t i=0; i<fCurrFil[0];i++){
    printf("%6d ",fArray[i]);
    if(i%6==0 && i>0)printf("\n");
  }
  printf("\n");
  cout<<"\n Elements for N side: \n";
  for(Int_t i=0; i<fCurrFil[1];i++){
    printf("%6d ",fArray[fDim/2+i]);
    if(i%6==0 && i>0)printf("\n");
  }
  printf("\n");
}

//----------------------------------------------------------------------
Int_t AliITSTableSSD::Use(Int_t side){
  // uses the current element. This means that the current element is returned
  // and its content is replaced by -1. Hence each element can be used only 
  // once.
  Int_t elem=-1;
  if((side!=0) && (side!=1)){
    cerr<<" Error in AliITSTableSSD::Use. side="<<side<<endl;
    cerr<<" side must be 0 or 1\n";
    return elem;
  }
  if(fCurrUse[side]>(fDim/2-1)){
    cerr<<" Error in AliITSTableSSD::Use. Trying to use an element out of range\n";
    cerr<<" Element="<<fCurrUse[side]<<" while max limit is "<<(fDim/2-1)<<endl;
    fCurrUse[side]++;
    return elem;
  }
  Int_t index=(fDim/2)*side+fCurrUse[side];
  elem=fArray[index];
  fArray[index]=-1;
  fCurrUse[side]++;
  return elem;
}
