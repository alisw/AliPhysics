#ifndef ALIESDVERTEX_H
#define ALIESDVERTEX_H
/* Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------
//                    Primary Vertex Class
//          for the Event Data Summary Class
//   Origin: A.Dainese, Padova, andrea.dainese@pd.infn.it
//-------------------------------------------------------

/*****************************************************************************
 *                                                                           *
 * This class deals with primary vertex.                                     *
 * AliESDVertex objects are created by the class AliITSVertexer and its      *
 * derived classes.                                                          *
 * Different constructors are provided:                                      *
 * - for primary vertex determined with pixels in pp (only Z)                *
 * - for primary vertex determined with pixels in ion-ion (X,Y,Z)            *
 * - for primary vertex determined with ITS tracks in pp (X,Y,Z)             *
 * This class replaces the previous version of AliESDVertex, designed        *
 * originally only for A-A collisions. The functionalities of the old class  *
 * are maintained in the AliITSVertexerIons class                            *
 *                                                                           *
 *****************************************************************************/

#include <TMath.h>

#include "AliVertex.h"

class AliESDVertex : public AliVertex {
 
 public:
 
  AliESDVertex();
  AliESDVertex(Double_t positionZ,Double_t sigmaZ,Int_t nContributors,
	       const Char_t *vtxName="Vertex");
  AliESDVertex(Double_t position[3],Double_t covmatrix[6],
	       Double_t chi2,Int_t nContributors,
	       const Char_t *vtxName="Vertex");
  AliESDVertex(Double_t position[3],Double_t sigma[3],
	       const Char_t *vtxName="Vertex");
  AliESDVertex(Double_t position[3],Double_t sigma[3],Double_t snr[3],
	       const Char_t *vtxName="Vertex");

  virtual ~AliESDVertex(){;}


  void     GetSigmaXYZ(Double_t sigma[3]) const;
  void     GetCovMatrix(Double_t covmatrix[6]) const;
  void     GetSNR(Double_t snr[3]) const;

  Double_t GetXRes() const {return TMath::Sqrt(fCovXX);}
  Double_t GetYRes() const {return TMath::Sqrt(fCovYY);}
  Double_t GetZRes() const {return TMath::Sqrt(fCovZZ);}
  Double_t GetXSNR() const { return fSNR[0]; }
  Double_t GetYSNR() const { return fSNR[1]; }
  Double_t GetZSNR() const { return fSNR[2]; }

  Double_t GetChi2() const { return fChi2; }
  Double_t GetChi2toNDF() const 
    { return fChi2/(2.*(Double_t)fNContributors-3.); }

  void     Print(Option_t* option = "") const;
  void     PrintStatus() const {Print();}
  void SetTruePos(Double_t *tp){for(Int_t i=0;i<3;i++)fTruePos[i]=tp[i];}
  void GetTruePos(Double_t *tp) const 
       {for(Int_t i=0;i<3;i++)tp[i]=fTruePos[i];}

  void Reset() { SetToZero(); SetName("Vertex"); }

 protected:

  Double_t fCovXX,fCovXY,fCovYY,fCovXZ,fCovYZ,fCovZZ;  // vertex covariance matrix
  Double_t fSNR[3];  // S/N ratio
  Double_t fChi2;  // chi2 of vertex fit
  Double_t fTruePos[3];   //true vertex position (for comparison purposes)

 private:

  void SetToZero();

  ClassDef(AliESDVertex,5)  // Class for Primary Vertex
};

#endif






    
