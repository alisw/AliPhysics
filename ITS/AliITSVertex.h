#ifndef ALIITSVERTEX_H
#define ALIITSVERTEX_H
/* Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


//-------------------------------------------------------
//                    Primary Vertex Class
//
//   Origin: A.Dainese, Padova, andrea.dainese@pd.infn.it
//-------------------------------------------------------

/*****************************************************************************
 *                                                                           *
 * This class deals with primary vertex.                                     *
 * AliITSVertex objects are created by the class AliITSVertexer and its      *
 * derived classes.                                                          *
 * Different constructors are provided:                                      *
 * - for primary vertex determined with pixels in pp (only Z)                *
 * - for primary vertex determined with pixels in ion-ion (X,Y,Z)            *
 * - for primary vertex determined with ITS tracks in pp (X,Y,Z)             *
 * This class replaces the previous version of AliITSVertex, designed        *
 * originally only for A-A collisions. The functionalities of the old class  *
 * are maintained in the AliITSVertexerIons class                            *
 *                                                                           *
 *****************************************************************************/

//---- Root headers -----
#include <TNamed.h>

class AliITSVertex : public TNamed {
 
 public:
 
  AliITSVertex();
  AliITSVertex(Double_t positionZ,Double_t sigmaZ,Int_t nContributors,
	       const char *vtxName="Vertex");
  AliITSVertex(Double_t phi,Double_t position[3],Double_t covmatrix[6],
	       Double_t chi2,Int_t nContributors,
	       const char *vtxName="Vertex");
  AliITSVertex(Double_t position[3],Double_t sigma[3],
	       const char *vtxName="Vertex");
  AliITSVertex(Double_t position[3],Double_t sigma[3],Double_t snr[3],
	       const char *vtxName="Vertex");

  virtual ~AliITSVertex();


  void     GetXYZ(Double_t position[3]) const;
  void     GetSigmaXYZ(Double_t sigma[3]) const;
  void     GetCovMatrix(Double_t covmatrix[6]) const;
  void     GetSNR(Double_t snr[3]) const;

  void     GetXYZ_ThrustFrame(Double_t position[3]) const;
  void     GetSigmaXYZ_ThrustFrame(Double_t sigma[3]) const;
  void     GetCovMatrix_ThrustFrame(Double_t covmatrix[6]) const;

  Double_t GetXv() const;
  Double_t GetYv() const;
  Double_t GetZv() const;
  Double_t GetXRes() const;
  Double_t GetYRes() const;
  Double_t GetZRes() const;
  Double_t GetXSNR() const { return fSNR[0]; }
  Double_t GetYSNR() const { return fSNR[1]; }
  Double_t GetZSNR() const { return fSNR[2]; }

  Double_t GetPhi() const { return fPhi; }
  Double_t GetChi2() const { return fChi2; }
  Double_t GetChi2toNDF() const 
    { return fChi2/(2.*(Double_t)fNContributors-3.); }
  Int_t    GetNContributors() const { return fNContributors; }

  void     PrintStatus() const;
  void     SetDebug(Int_t dbg = 0) { fDebug = dbg; return; }
  void SetTruePos(Double_t *tp){for(Int_t i=0;i<3;i++)fTruePos[i]=tp[i];}
  void GetTruePos(Double_t &x,Double_t &y,Double_t &z) const 
       {x=fTruePos[0]; y=fTruePos[1]; z=fTruePos[2];}
        
 private:

  void SetToZero();

 protected:

  Double_t fPosition[3];  // vertex position
  Double_t fCovXX,fCovXY,fCovYY,fCovXZ,fCovYZ,fCovZZ;  // vertex covariance matrix
  Double_t fSNR[3];  // S/N ratio
  Double_t fPhi;     // position is given in a frame rotated of an angle fPhi around z axis: but in most cases fPhi = 0 ! 
  Double_t fChi2;  // chi2 of vertex fit
  Int_t    fNContributors;  // # of tracklets/tracks used for the estimate 
  Int_t    fDebug;  //! for debugging
  Double_t fTruePos[3];   //true vertex position (for comparison purposes)
  ClassDef(AliITSVertex,3)  // Class for Primary Vertex
    };

#endif






    
