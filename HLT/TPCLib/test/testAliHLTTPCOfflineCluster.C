// $Id$

/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Kalliopi Kanaki <Kalliopi.Kanaki@ift.uib.no>        *
 *                  for The ALICE HLT Project.                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/** @file   testAliHLTTPCOfflineCluster.C
    @author Kalliopi Kanaki
    @date   
    @brief  Test macro/program for the AliHLTTPCOfflineCluster class
 */

#ifndef __CINT__
#include "TSystem.h"

#include "AliHLTTPCTransform.h"
#include "AliHLTTPCClusterData.h"
#include "AliHLTTPCSpacePointData.h"
#include "AliHLTTPCOfflineCluster.h"
#include "TObjArray.h"

#include <ostream>
#include <istream>
#endif //__CINT__

using namespace std;

/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
//
// configuration of the test program
//
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
 


#include "/scratch/prog/HLT/TPCLib/AliHLTTPCSpacePointData.h";

int testAliHLTTPCOfflineCluster(){
//#ifdef __CINT__
  //gSystem->Load("libAliHLTUtil.so");
  gSystem->Load("libAliHLTTPC.so");
//#endif

  int iResult=0;
  
  



 AliHLTTPCSpacePointData spacePoint1 = { 5.,3.,2.,7,8,11.2,5.4.,1,0,kFALSE,0 };

//   AliHLTTPCSpacePointData *spacePoint = new AliHLTTPCSpacePointData(); //[1];
//   spacePoint->fX = 3.;
//   spacePoint->fY = 2.;
//   spacePoint->fZ = 5.;
   
  cout << "hlt X: "      <<        spacePoint1.fX      << endl;
  cout << "hlt padrow: " << (Int_t)spacePoint1.fPadRow << endl;
  cout << "hlt Z: "      <<        spacePoint1.fZ      << endl;
  
  
  //Float_t xyz[3] = {spacePoint1.fX , spacePoint1.fPadRow, spacePoint1.fZ };
  //AliHLTTPCTransform::Local2Global(xyz,36);  
  //cout << xyz[0] << " " << xyz[1] << " " << xyz[2] << endl;

  AliTPCclusterMI *offCl1 = new AliTPCclusterMI();  
  AliHLTTPCOfflineCluster *p1 = new AliHLTTPCOfflineCluster();
   
  offCl1 = p1->ConvertHLTToOffline(spacePoint1);
  //offCl1 = new AliHLTTPCOfflineCluster(spacePoint1); // ??? to be tested

  cout << "off pad: "      << offCl1->GetPad() << endl;
  cout << "off row: "      << offCl1->GetRow() << endl;
  cout << "off time bin: " << offCl1->GetTimeBin() << endl;
  
  cout << "================================" << endl;  
  
  AliTPCclusterMI *offCl2 = new AliTPCclusterMI();
  offCl2->SetPad(2.5);
  offCl2->SetRow(4);

  cout << "off pad: " << offCl2->GetPad() << endl;
  cout << "off row: " << offCl2->GetRow() << endl;
  
  AliHLTTPCOfflineCluster *p2 = new AliHLTTPCOfflineCluster();
  AliHLTTPCSpacePointData spacePoint2 = p2->ConvertOfflineToHLT(offCl2);
  
  cout << "hlt X: "      << spacePoint2.fX             << endl;
  cout << "hlt padrow: " << (Int_t)spacePoint2.fPadRow << endl;

  return iResult;
  //delete [] spacePoint;
 
}

int main(int /*argc*/, const char** /*argv*/){

  return testAliHLTTPCOfflineCluster();
}
