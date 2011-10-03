/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id$ */

//-------------------------------------------------------------------------
//     AOD track base class
//     Base class for Analysis Object Data
//     Generic version
//     Author: Markus Oldenburg, CERN
//     Inheritance from AliVVertex: A. Dainese
//-------------------------------------------------------------------------

#include "AliAODVertex.h"
#include "AliAODTrack.h"

ClassImp(AliAODVertex)

//______________________________________________________________________________
AliAODVertex::AliAODVertex() : 
  AliVVertex(),
  fChi2perNDF(-999.),
  fID(-1),
  fBCID(AliVTrack::kTOFBCNA),
  fType(kUndef),
  fNprong(0),
  fIprong(0),
  fNContributors(0),
  fCovMatrix(NULL),
  fParent(),
  fDaughters(),
  fProngs(NULL)
  {
  // default constructor

  fPosition[0] = fPosition[1] = fPosition[2] = -999.;
}

//______________________________________________________________________________
AliAODVertex::AliAODVertex(const Double_t position[3], 
			   const Double_t covMatrix[6],
			   Double_t  chi2perNDF,
			   TObject  *parent,
			   Short_t id,
			   Char_t vtype, 
			   Int_t  nprong) :
  AliVVertex(),
  fChi2perNDF(chi2perNDF),
  fID(id),
  fBCID(AliVTrack::kTOFBCNA),
  fType(vtype),
  fNprong(nprong),
  fIprong(0),
  fNContributors(0),
  fCovMatrix(NULL),
  fParent(parent),
  fDaughters(),
  fProngs(0)
{
  // constructor

  SetPosition(position);
  if (covMatrix) SetCovMatrix(covMatrix);
  MakeProngs();
}

//______________________________________________________________________________
AliAODVertex::AliAODVertex(const Float_t position[3], 
			   const Float_t  covMatrix[6],
			   Double_t  chi2perNDF,
			   TObject  *parent,
			   Short_t id,
			   Char_t vtype,
			   Int_t nprong) :

  AliVVertex(),
  fChi2perNDF(chi2perNDF),
  fID(id),
  fBCID(AliVTrack::kTOFBCNA),
  fType(vtype),
  fNprong(nprong),
  fIprong(0),
  fNContributors(0),
  fCovMatrix(NULL),
  fParent(parent),
  fDaughters(),
  fProngs(0)
{
  // constructor

  SetPosition(position);
  if (covMatrix) SetCovMatrix(covMatrix);
  MakeProngs();
}

//______________________________________________________________________________
AliAODVertex::AliAODVertex(const Double_t position[3], 
			   Double_t  chi2perNDF,
			   Char_t vtype, 
			   Int_t nprong) :
  AliVVertex(),
  fChi2perNDF(chi2perNDF),
  fID(-1),
  fBCID(AliVTrack::kTOFBCNA),
  fType(vtype),
  fNprong(nprong),
  fIprong(0),
  fNContributors(0),
  fCovMatrix(NULL),
  fParent(),
  fDaughters(),
  fProngs(0)
{
  // constructor without covariance matrix

  SetPosition(position);  
  MakeProngs();
}

//______________________________________________________________________________
AliAODVertex::AliAODVertex(const Float_t position[3], 
			   Double_t  chi2perNDF,
			   Char_t vtype, Int_t nprong) :
  AliVVertex(),
  fChi2perNDF(chi2perNDF),
  fID(-1),
  fBCID(AliVTrack::kTOFBCNA),
  fType(vtype),
  fNprong(nprong),
  fIprong(0),
  fNContributors(0),
  fCovMatrix(NULL),
  fParent(),
  fDaughters(),
  fProngs(0)
{
  // constructor without covariance matrix

  SetPosition(position);  
  MakeProngs();
}

//______________________________________________________________________________
AliAODVertex::~AliAODVertex() 
{
  // Destructor

  delete fCovMatrix;
  if (fNprong > 0) delete[] fProngs;
}

//______________________________________________________________________________
AliAODVertex::AliAODVertex(const AliAODVertex& vtx) :
  AliVVertex(vtx),
  fChi2perNDF(vtx.fChi2perNDF),
  fID(vtx.fID),
  fBCID(vtx.fBCID),
  fType(vtx.fType),
  fNprong(vtx.fNprong),
  fIprong(vtx.fIprong),
  fNContributors(vtx.fNContributors),
  fCovMatrix(NULL),
  fParent(vtx.fParent),
  fDaughters(vtx.fDaughters),
  fProngs(0)
{
  // Copy constructor.
  
  for (int i = 0; i < 3; i++) 
    fPosition[i] = vtx.fPosition[i];

  if (vtx.fCovMatrix) fCovMatrix=new AliAODRedCov<3>(*vtx.fCovMatrix);
  MakeProngs();
  for (int i = 0; i < fNprong; i++) {
      fProngs[i] = vtx.fProngs[i];
  }
}

//______________________________________________________________________________
AliAODVertex* AliAODVertex::CloneWithoutRefs() const
{
  // Special method to copy all but the refs 
  
  Double_t cov[6] = { 0.0 };
      
  if (fCovMatrix) fCovMatrix->GetCovMatrix(cov);
  
  AliAODVertex* v = new AliAODVertex(fPosition,
                                     cov,
                                     fChi2perNDF,
                                     0x0,
                                     fID,
                                     fType,
                                     0);
  
  v->SetNContributors(fNContributors);  
  
  return v;
}

//______________________________________________________________________________
AliAODVertex& AliAODVertex::operator=(const AliAODVertex& vtx) 
{
  // Assignment operator
  if (this != &vtx) {

    // name and type
    AliVVertex::operator=(vtx);

    //momentum
    for (int i = 0; i < 3; i++) 
      fPosition[i] = vtx.fPosition[i];
    
    fChi2perNDF = vtx.fChi2perNDF;
    fID = vtx.fID;
    fType = vtx.fType;

    //covariance matrix
    delete fCovMatrix;
    fCovMatrix = NULL;   
    if (vtx.fCovMatrix) fCovMatrix=new AliAODRedCov<3>(*vtx.fCovMatrix);
    
    //other stuff
    fParent = vtx.fParent;
    fDaughters = vtx.fDaughters;
    fNprong    = vtx.fNprong;
    fIprong    = vtx.fIprong;  

    MakeProngs();
    for (int i = 0; i < fNprong; i++) {
	fProngs[i] = vtx.fProngs[i];
    }
  }
  
  return *this;
}

//______________________________________________________________________________
void AliAODVertex::AddDaughter(TObject *daughter)
{
  // Add reference to daughter track
    if (!fProngs) {
	if (fDaughters.GetEntries()==0) {
	    TRefArray* arr = &fDaughters;
	    new(arr)TRefArray(TProcessID::GetProcessWithUID(daughter));  	
	}
	fDaughters.Add(daughter);	
    } else {
	if (fIprong < fNprong) {
	    fProngs[fIprong++] = daughter;
	} else {
	    AliWarning("Number of daughters out of range !\n");
	}
    }
  return;
}


//______________________________________________________________________________
template <class T> void AliAODVertex::GetSigmaXYZ(T sigma[3]) const
{
  // Return errors on vertex position in thrust frame
  
  if(fCovMatrix) {
    sigma[0]=fCovMatrix[3]; //GetCovXZ
    sigma[1]=fCovMatrix[4]; //GetCovYZ
    sigma[2]=fCovMatrix[5]; //GetCovZZ
  } else 
    sigma[0]=sigma[1]=sigma[2]=-999.;

  /*
  for (int i = 0, j = 6; i < 3; i++) {
    j -= i+1;
    sigma[2-i] = fCovMatrix ? TMath::Sqrt(fCovMatrix[j]) : -999.;
  }
  */
}

//______________________________________________________________________________
Int_t AliAODVertex::GetNContributors() const 
{
  // Returns the number of tracks used to fit this vertex.
  Int_t cont  = 0;

  TString vtitle = GetTitle();
  if (!vtitle.Contains("VertexerTracks")) {
    cont = fNContributors;
  } else {
    for (Int_t iDaug = 0; iDaug < GetNDaughters(); iDaug++) {
	AliAODTrack* aodT = dynamic_cast<AliAODTrack*>(fDaughters.At(iDaug));
	if (!aodT) continue;
	if (aodT->GetUsedForPrimVtxFit()) cont++;
    } 
    // the constraint adds another DOF
    if(vtitle.Contains("VertexerTracksWithConstraint"))cont++;
  }
  return cont;
}

//______________________________________________________________________________
Bool_t AliAODVertex::HasDaughter(TObject *daughter) const 
{
  // Checks if the given daughter (particle) is part of this vertex.
    if (!fProngs) {
	TRefArrayIter iter(&fDaughters);
	while (TObject *daugh = iter.Next()) {
	    if (daugh == daughter) return kTRUE;
	}
	return kFALSE;
    } else {
	Bool_t has = kFALSE;
	for (int i = 0; i < fNprong; i++) {
	    if (fProngs[i].GetObject() == daughter) has = kTRUE;
	}
	return has;
    }
}

//______________________________________________________________________________
Double_t AliAODVertex::RotatedCovMatrixXX(Double_t phi, Double_t theta) const
{
  // XX term of covariance matrix after rotation by phi around z-axis
  // and, then, by theta around new y-axis

  if (!fCovMatrix) {
    //AliFatal("Covariance matrix not set");
    return -999.;
  }

  Double_t covMatrix[6];

  GetCovMatrix(covMatrix);

  Double_t cp = TMath::Cos(phi);
  Double_t sp = TMath::Sin(phi);
  Double_t ct = TMath::Cos(theta);
  Double_t st = TMath::Sin(theta);
  return
     covMatrix[0]*cp*cp*ct*ct  // GetCovXX
    +covMatrix[1]*2.*cp*sp*ct*ct  // GetCovXY
    +covMatrix[3]*2.*cp*ct*st  // GetCovXZ
    +covMatrix[2]*sp*sp*ct*ct  // GetCovYY
    +covMatrix[4]*2.*sp*ct*st  // GetCovYZ
    +covMatrix[5]*st*st;  // GetCovZZ
}

//______________________________________________________________________________
Double_t AliAODVertex::RotatedCovMatrixXY(Double_t phi, Double_t theta) const
{
  // XY term of covariance matrix after rotation by phi around z-axis
  // and, then, by theta around new y-axis

  if (!fCovMatrix) {
    //AliFatal("Covariance matrix not set");
    return -999.;
  }

  Double_t covMatrix[6];

  GetCovMatrix(covMatrix);

  Double_t cp = TMath::Cos(phi);
  Double_t sp = TMath::Sin(phi);
  Double_t ct = TMath::Cos(theta);
  Double_t st = TMath::Sin(theta);
  return 
    -covMatrix[0]*cp*sp*ct  // GetCovXX
    +covMatrix[1]*ct*(cp*cp-sp*sp)  // GetCovXY
    -covMatrix[3]*sp*st  // GetCovXZ
    +covMatrix[2]*cp*sp*ct  // GetCovYY
    +covMatrix[4]*cp*st;  // GetCovYZ
}

//______________________________________________________________________________
Double_t AliAODVertex::RotatedCovMatrixXZ(Double_t phi, Double_t theta) const
{
  // XZ term of covariance matrix after rotation by phi around z-axis
  // and, then, by theta around new y-axis

  if (!fCovMatrix) {
    //AliFatal("Covariance matrix not set");
    return -999.;
  }

  Double_t covMatrix[6];

  GetCovMatrix(covMatrix);

  Double_t cp = TMath::Cos(phi);
  Double_t sp = TMath::Sin(phi);
  Double_t ct = TMath::Cos(theta);
  Double_t st = TMath::Sin(theta);
  return 
    -covMatrix[0]*cp*cp*ct*st  // GetCovXX
    -covMatrix[1]*2.*cp*sp*ct*st  // GetCovXY
    +covMatrix[3]*cp*(ct*ct-st*st)  // GetCovXZ
    -covMatrix[2]*sp*sp*ct*st  // GetCovYY
    +covMatrix[4]*sp*(ct*ct-st*st)  // GetCovYZ
    +covMatrix[5]*ct*st;  // GetCovZZ
}

//______________________________________________________________________________
Double_t AliAODVertex::RotatedCovMatrixYY(Double_t phi) const
{
  // YY term of covariance matrix after rotation by phi around z-axis
  // and, then, by theta around new y-axis

  if (!fCovMatrix) {
    //AliFatal("Covariance matrix not set");
    return -999.;
  }

  Double_t covMatrix[6];

  GetCovMatrix(covMatrix);

  Double_t cp = TMath::Cos(phi);
  Double_t sp = TMath::Sin(phi);
  return
     covMatrix[0]*sp*sp  // GetCovXX
    -covMatrix[1]*2.*cp*sp  // GetCovXY
    +covMatrix[2]*cp*cp;  // GetCovYY
}

//______________________________________________________________________________
Double_t AliAODVertex::RotatedCovMatrixYZ(Double_t phi, Double_t theta) const
{
  // YZ term of covariance matrix after rotation by phi around z-axis
  // and, then, by theta around new y-axis

  if (!fCovMatrix) {
    //AliFatal("Covariance matrix not set");
    return -999.;
  }

  Double_t covMatrix[6];

  GetCovMatrix(covMatrix);

  Double_t cp = TMath::Cos(phi);
  Double_t sp = TMath::Sin(phi);
  Double_t ct = TMath::Cos(theta);
  Double_t st = TMath::Sin(theta);
  return 
     covMatrix[0]*cp*sp*st  // GetCovXX
    +covMatrix[1]*st*(sp*sp-cp*cp)  // GetCovXY
    -covMatrix[3]*sp*ct  // GetCovXZ
    -covMatrix[2]*cp*sp*st  // GetCovYY
    +covMatrix[4]*cp*ct;  // GetCovYZ
}

//______________________________________________________________________________
Double_t AliAODVertex::RotatedCovMatrixZZ(Double_t phi, Double_t theta) const
{
  // ZZ term of covariance matrix after rotation by phi around z-axis
  // and, then, by theta around new y-axis

  if (!fCovMatrix) {
    //AliFatal("Covariance matrix not set");
    return -999.;
  }

  Double_t covMatrix[6];

  GetCovMatrix(covMatrix);

  Double_t cp = TMath::Cos(phi);
  Double_t sp = TMath::Sin(phi);
  Double_t ct = TMath::Cos(theta);
  Double_t st = TMath::Sin(theta);
  return
     covMatrix[0]*cp*cp*st*st  // GetCovXX
    +covMatrix[1]*2.*cp*sp*st*st  // GetCovXY
    -covMatrix[3]*2.*cp*ct*st  // GetCovXZ
    +covMatrix[2]*sp*sp*st*st  // GetCovYY
    -covMatrix[4]*2.*sp*sp*ct*st  // GetCovYZ
    +covMatrix[5]*ct*ct;  // GetCovZZ
}

//______________________________________________________________________________
Double_t AliAODVertex::Distance2ToVertex(const AliAODVertex *vtx) const
{
  // distance in 3D to another AliAODVertex

  Double_t dx = GetX()-vtx->GetX();
  Double_t dy = GetY()-vtx->GetY();
  Double_t dz = GetZ()-vtx->GetZ();

  return dx*dx+dy*dy+dz*dz;
}

//______________________________________________________________________________
Double_t AliAODVertex::DistanceXY2ToVertex(const AliAODVertex *vtx) const
{
  // distance in XY to another AliAODVertex

  Double_t dx = GetX()-vtx->GetX();
  Double_t dy = GetY()-vtx->GetY();

  return dx*dx+dy*dy;
}

//______________________________________________________________________________
Double_t AliAODVertex::Error2DistanceToVertex(AliAODVertex *vtx) const
{
  // error on the distance in 3D to another AliAODVertex

  Double_t phi,theta;
  PhiAndThetaToVertex(vtx,phi,theta);
  // error2 due to this vertex
  Double_t error2 = RotatedCovMatrixXX(phi,theta);
  // error2 due to vtx vertex
  Double_t error2vtx = vtx->RotatedCovMatrixXX(phi,theta);

  return error2+error2vtx;
}

//______________________________________________________________________________
Double_t AliAODVertex::Error2DistanceXYToVertex(AliAODVertex *vtx) const
{
  // error on the distance in XY to another AliAODVertex

  Double_t phi,theta;
  PhiAndThetaToVertex(vtx,phi,theta);
  // error2 due to this vertex
  Double_t error2 = RotatedCovMatrixXX(phi);
  // error2 due to vtx vertex
  Double_t error2vtx = vtx->RotatedCovMatrixXX(phi);

  return error2+error2vtx;
}

//______________________________________________________________________________
template <class T, class P> 
void AliAODVertex::PhiAndThetaToVertex(AliAODVertex *vtx, P &phi, T &theta) const
{
  // rotation angles around z-axis (phi) and around new y-axis (theta)
  // with which vtx is seen (used by RotatedCovMatrix... methods)

  phi = TMath::Pi()+TMath::ATan2(-vtx->GetY()+GetY(),-vtx->GetX()+GetX());
  Double_t vtxxphi = vtx->GetX()*TMath::Cos(phi)+vtx->GetY()*TMath::Sin(phi);
  Double_t xphi = GetX()*TMath::Cos(phi)+GetY()*TMath::Sin(phi);
  theta = TMath::ATan2(vtx->GetZ()-GetZ(),vtxxphi-xphi);
}

//______________________________________________________________________________
void AliAODVertex::PrintIndices() const 
{
  // Print indices of particles originating form this vertex

  TRefArrayIter iter(&fDaughters);
  while (TObject *daugh = iter.Next()) {
    printf("Particle %p originates from this vertex.\n", static_cast<void*>(daugh));
  }
}

//______________________________________________________________________________
const char* AliAODVertex::AsString() const
{
  // Make a string describing this object
  
  TString tmp(Form("%10s pos(%7.2f,%7.2f,%7.2f)",GetTypeName((AODVtx_t)GetType()),GetX(),GetY(),GetZ()));
  
  if (GetType()==kPrimary || GetType()==kMainSPD || GetType()==kPileupSPD )
  {
    tmp += Form(" ncontrib %d chi2/ndf %4.1f",GetNContributors(),GetChi2perNDF());

  }
  
  if ( !fParent.GetObject() ) 
  {
    tmp += " no parent";
  }
  if ( fDaughters.GetEntriesFast() > 0 )
  {
    if ( fDaughters.GetEntriesFast() == 1 ) 
    {
      tmp += " origin of 1 particle";
    }
    else
    {
      tmp += Form(" origin of %2d particles",fDaughters.GetEntriesFast());
    }
  }
  
  return tmp.Data();
}

//______________________________________________________________________________
const char* AliAODVertex::GetTypeName(AODVtx_t type)
{
  // Return an ASCII version of type
  
  switch (type)
  {
    case kPrimary:
      return "primary";
      break;
    case kKink:
      return "kink";
      break;
    case kV0:
      return "v0";
      break;
    case kCascade:
      return "cascade";
      break;
    case kMainSPD:
      return "mainSPD";
      break;
    case kPileupSPD:
      return "pileupSPD";
      break;
    case kPileupTracks:
      return "pileupTRK";
      break;
    case kMainTPC:
      return "mainTPC";
      break;
    default:
      return "unknown";
      break;
  };
}

//______________________________________________________________________________
void AliAODVertex::Print(Option_t* /*option*/) const 
{
  // Print information of all data members

  printf("Vertex position:\n");
  printf("     x = %f\n", fPosition[0]);
  printf("     y = %f\n", fPosition[1]);
  printf("     z = %f\n", fPosition[2]);
  printf(" parent particle: %p\n", static_cast<void*>(fParent.GetObject()));
  printf(" origin of %d particles\n", fDaughters.GetEntriesFast());
  printf(" vertex type %d\n", fType);
  
  /*
  if (fCovMatrix) {
    printf("Covariance matrix:\n");
    printf(" %12.10f  %12.10f  %12.10f\n %12.10f  %12.10f  %12.10f\n %12.10f  %12.10f  %12.10f\n", 
	   fCovMatrix[0],
	   fCovMatrix[1],
	   fCovMatrix[3],
	   fCovMatrix[1],
	   fCovMatrix[2],
	   fCovMatrix[4],
	   fCovMatrix[3],
	   fCovMatrix[4],
	   fCovMatrix[5]); 
	   } */
  printf(" Chi^2/NDF = %f\n", fChi2perNDF);
}

