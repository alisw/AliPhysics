// $Id$

// Author: Uli Frankenfeld <uli.frankenfeld@gsi.de>
//*-- Copyright &copy ALICE HLT Group


#include "speedtest.h"

int main(int arg,char **arc){
  int n = 0;
  if(arg!=2) {cerr<<"usage: speedtest #loops \n";return -1;}  
  n = atoi(arc[1]);
//  cerr<<"allocate: "<<n*sizeof(AliHLTConfMapPoint)<<" Bytes"<<endl;
//  AliHLTConfMapPoint *array = new AliHLTConfMapPoint[n];
  cerr<<"allocate: "<<n*sizeof(double)<<" Bytes"<<endl;
  double *array = new double[n];
//  cerr<<"allocate: "<<n*sizeof(int)<<" Bytes"<<endl;
//  int *array = new int[n];
//  cerr<<"allocate: "<<n*sizeof(AliHLTSpacePointData)<<" Bytes"<<endl;
//  AliHLTSpacePointData *array = new AliHLTSpacePointData[n];

  AliHLTSpacePointData hit;
  hit.fX=103.55;
  hit.fY=22.33;
  hit.fZ=95.312;
  hit.fID=(4<<25)|14042;
  hit.fPadRow=77;
  hit.fXYErr=1;
  hit.fZErr=2.5;

  cerr<<"start loop"<<endl;
  double initCpuTime,cpuTime;
  initCpuTime = CpuTime();

/*  
  for(int i=0;i<n;i++){
    array[i].SetHitNumber(hit.fID);
    array[i].SetPadRow(hit.fPadRow);
    Int_t slice = (hit.fID>>25) & 0x7f;
    array[i].SetSector(slice);
    array[i].SetX(hit.fX);
    array[i].SetY(hit.fY);
    array[i].SetZ(hit.fZ);
    array[i].SetXerr(sqrt((double)hit.fXYErr));
//    array[i].SetYerr(sqrt(hit.fXYErr));
//    array[i].SetZerr(sqrt(hit.fZErr));
  }
*/
//  for(int i=0;i<n;i++) array[i].ReadHits(&hit);
  for(int i=0;i<n;i++) array[i]=sqrt(i); 
//  for(int i=0;i<n;i++) array[i]=i; 
//  for(int i=0;i<n;i++) array[i].fX=hit.fX; 

  cpuTime = CpuTime() - initCpuTime;
//  cerr<<"cpuTime: "<<cpuTime<<endl;
  printf("cpuTime: %d ms\n",int(cpuTime*1000));
  return 0;
}

