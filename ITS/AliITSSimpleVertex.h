#ifndef ALIITSSIMPLEVERTEX_H
#define ALIITSSIMPLEVERTEX_H
/* Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


//-------------------------------------------------------
//                    Secondary Vertex Class
//          as calculated by AliITSVertexerTracks
//   Origin: F. Prino, Torino, prino@to.infn.it
//-------------------------------------------------------

/*****************************************************************************
 *                                                                           *
 * This class deals with secondary vertices.                                 *
 * AliITSSimpleVertex objects are created by the class AliITSVertexerTracks  *
 *                                                                           *
 *****************************************************************************/

//---- Root headers -----
#include <TNamed.h>

class AliITSSimpleVertex : public TNamed {
 
 public:
 
  AliITSSimpleVertex();
  AliITSSimpleVertex(Double_t position[3],Double_t dispersion,
		Int_t nContributors);
  virtual ~AliITSSimpleVertex();


  void   SetXYZ(Double_t pos[3]) {for(Int_t j=0; j<3; j++) fPosition[j]=pos[j];}
  void   SetXv(Double_t xVert) {fPosition[0]=xVert; }
  void   SetYv(Double_t yVert) {fPosition[1]=yVert; }
  void   SetZv(Double_t zVert) {fPosition[2]=zVert; }
  void   SetDispersion(Double_t disp) { fSigma=disp; }
  void   SetNContributors(Int_t nContr) {fNContributors=nContr; }

  void     GetXYZ(Double_t position[3]) const;
  Double_t GetXv() const { return fPosition[0]; }
  Double_t GetYv() const { return fPosition[1]; }
  Double_t GetZv() const { return fPosition[2]; }
  Double_t GetDispersion() const { return fSigma; }
  Int_t    GetNContributors() const { return fNContributors; }

  void     Print(Option_t* option = "") const;


 protected:

  Double_t fPosition[3];  // vertex position
  Double_t fSigma;  // track dispersion around found vertex
  Int_t    fNContributors;  // # of tracklets/tracks used for the estimate 


  ClassDef(AliITSSimpleVertex,1)  // Class for Primary Vertex
};

#endif






    
