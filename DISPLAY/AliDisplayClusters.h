#ifndef ALIDISPLAYCLUSTERS_H
#define ALIDISPLAYCLUSTERS_H
/* Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
/////////////////////////////////////////////////////////////////////////
// ALICE DISPLAY CLUSTERS CLASS                                        //
// Author: Mayeul   ROUSSELET                                          //
// e-mail: Mayeul.Rousselet@cern.ch                                    //
// Last update:26/08/2003                                              //
/////////////////////////////////////////////////////////////////////////

#include <Rtypes.h>
#include <RQ_OBJECT.h>

class TPolyMarker3D;

class AliDisplayClusters{
  //This class is an interface to the clusters data

public:
 AliDisplayClusters();
 virtual ~AliDisplayClusters();

 void          LoadClusters(const char * name,Int_t nevent);
 void          LoadITSClusters(Int_t nevent);
 void          LoadTPCClusters(Int_t nevent);
 void          Draw();
 Int_t         GetNbClusters();

private: 
 TPolyMarker3D *fPoints; //fPoints[i]=set of cluster coordinates in detector i;
 Int_t         fNb;      //Nimber of clusters
 char          **fName; //fName[i]=name of the detector i 

 RQ_OBJECT("AliDisplayClusters")

 ClassDef(AliDisplayClusters,0);
};
#endif
