/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

/*
$Log$
*/

#include "AliCalcluster.h"
 
ClassImp(AliCalcluster) // Class implementation to enable ROOT I/O
 
AliCalcluster::AliCalcluster()
{
// Default constructer, all data is set to 0
 fCenter=0;
 fSig=0.;
 fNmods=0;
 fSig11=0.;
 fSig33=0.;
 fSig55=0.;
 fRowdisp=0.;
 fColdisp=0.;
 fNvetos=0;
 fVetos=0;
}
///////////////////////////////////////////////////////////////////////////
AliCalcluster::~AliCalcluster()
{
// Destructor to delete dynamically allocated memory
 if (fVetos)
 {
  fVetos->Delete();
  delete fVetos;
  fVetos=0;
 }
}
///////////////////////////////////////////////////////////////////////////
AliCalcluster::AliCalcluster(AliCalmodule& m)
{
// Cluster constructor with module m as center.
// Module data is only entered for a module which contains
// a signal and has not been used in a cluster yet.
// Modules situated at a detector edge are not allowed to start a cluster.

 Float_t pos[3]={0,0,0};

 if ((m.GetClusteredSignal() > 0.) && (m.GetEdgeValue() == 0))
 {
  fCenter=&m;
  m.GetPosition(pos,"sph");
  SetPosition(pos,"sph");
  fSig=m.GetClusteredSignal();
  fNmods=1;
  fSig11=m.GetClusteredSignal();
  fSig33=m.GetClusteredSignal();
  fSig55=m.GetClusteredSignal();
  fRowdisp=0.;
  fColdisp=0.;
  m.SetClusteredSignal(0.); // mark module as used in cluster
  fNvetos=0;
  fVetos=0;
 }
 else
 {
  fCenter=0;
  SetPosition(pos,"sph");
  fSig=0.;
  fNmods=0;
  fSig11=0.;
  fSig33=0.;
  fSig55=0.;
  fRowdisp=0.;
  fColdisp=0.;
  fNvetos=0;
  fVetos=0;
}
}
///////////////////////////////////////////////////////////////////////////
Int_t AliCalcluster::GetRow()
{
// Provide the row number of the cluster center
 if (fCenter)
 {
  return fCenter->GetRow();
 }
 else
 {
  return 0;
 }
}
///////////////////////////////////////////////////////////////////////////
Int_t AliCalcluster::GetColumn()
{
// Provide the column number of the cluster center
 if (fCenter)
 {
  return fCenter->GetColumn();
 }
 else
 {
  return 0;
 }
}
///////////////////////////////////////////////////////////////////////////
Float_t AliCalcluster::GetSignal(Int_t n)
{
// Provide the total signal of a module matrix of n modules around
// the cluster center.
// Example : n=9 --> total signal in the 3x3 matrix
//             1 --> signal of central module
// Note : n=0 provides the total cluster signal (Default)
 
 Float_t signal=-1;
 
 if (n == 0)  signal=fSig;
 if (n == 1)  signal=fSig11;
 if (n == 9)  signal=fSig33;
 if (n == 25) signal=fSig55;
 
 if (signal > 0.)
 {
  return signal;
 }
 else
 {
  cout << " *AliCalcluster::GetSignal* Invalid argument n = " << n << endl;
  return 0;
 }
}
///////////////////////////////////////////////////////////////////////////
Int_t AliCalcluster::GetNmodules()
{
// Provide the number of modules in the cluster
 return fNmods;
}
///////////////////////////////////////////////////////////////////////////
Float_t AliCalcluster::GetRowDispersion()
{
// Provide the normalised row dispersion of the cluster
 if (fSig > 0.)
 {
  return fRowdisp/fSig;
 }
 else
 {
  return 0.;
 }
}
///////////////////////////////////////////////////////////////////////////
Float_t AliCalcluster::GetColumnDispersion()
{
// Provide the normalised column dispersion of the cluster
 if (fSig > 0.)
 {
  return fColdisp/fSig;
 }
 else
 {
  return 0.;
 }
}
///////////////////////////////////////////////////////////////////////////
void AliCalcluster::Start(AliCalmodule& m)
{
// Reset cluster data and start with module m
// A module can only start a cluster when it contains
// a signal, has not been used in a cluster yet and is not
// situated at a detector edge

 Float_t pos[3]={0,0,0};

 if ((m.GetClusteredSignal() > 0.) && (m.GetEdgeValue() == 0))
 {
  fCenter=&m;
  m.GetPosition(pos,"sph");
  SetPosition(pos,"sph");
  fSig=m.GetSignal();
  fNmods=1;
  fSig11=m.GetSignal();
  fSig33=m.GetSignal();
  fSig55=m.GetSignal();
  fRowdisp=0.;
  fColdisp=0.;
  m.SetClusteredSignal(0.); // mark module as used in cluster
 }
 else
 {
  fCenter=0;
  SetPosition(pos,"sph");
  fSig=0.;
  fNmods=0;
  fSig11=0.;
  fSig33=0.;
  fSig55=0.;
  fRowdisp=0.;
  fColdisp=0.;
 }
}
///////////////////////////////////////////////////////////////////////////
void AliCalcluster::Add(AliCalmodule& m)
{
// Add module data to the cluster

 Float_t signal=m.GetClusteredSignal();
 
 if (signal > 0.) // only add unused modules
 {
  fSig+=signal;
  fNmods+=1;
  Int_t drow=int(fabs(double(GetRow()-m.GetRow())));       // row distance to center
  Int_t dcol=int(fabs(double(GetColumn()-m.GetColumn()))); // column distance to center
  if ((drow < 2) && (dcol < 2)) fSig33+=signal;
  if ((drow < 3) && (dcol < 3)) fSig55+=signal;
  fRowdisp+=signal*float(drow*drow);
  fColdisp+=signal*float(dcol*dcol);
  m.SetClusteredSignal(0.); // mark module as used in cluster
 }
}
///////////////////////////////////////////////////////////////////////////
void AliCalcluster::AddVetoSignal(Float_t* r,TString f,Float_t s)
{
// Associate an (extrapolated) AliSignal at location r as veto to the cluster
// Note : The default signal value (s) is 0
 if (!fVetos)
 {
  fNvetos=0;
  fVetos=new TObjArray();
 } 

 fVetos->Add(new AliSignal);
 fNvetos++;

 ((AliSignal*)fVetos->At(fNvetos-1))->SetPosition(r,f);
 ((AliSignal*)fVetos->At(fNvetos-1))->SetSignal(s);
}
///////////////////////////////////////////////////////////////////////////
Int_t AliCalcluster::GetNvetos()
{
// Provide the number of veto signals associated to the cluster
 return fNvetos;
}
///////////////////////////////////////////////////////////////////////////
AliSignal* AliCalcluster::GetVetoSignal(Int_t i)
{
// Provide access to the i-th veto signal of this cluster
// Note : The first hit corresponds to i=1

 if (i>0 && i<=fNvetos)
 {
  return (AliSignal*)fVetos->At(i-1);
 }
 else
 {
  cout << " *AliCalcluster::GetVetoSignal* Signal number " << i
       << " out of range." << endl;
  cout << " --- First signal (if any) returned." << endl;
  return (AliSignal*)fVetos->At(0);
 }
}
///////////////////////////////////////////////////////////////////////////
