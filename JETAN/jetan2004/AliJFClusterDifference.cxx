// $Id$

#include <Riostream.h>

#include <TMath.h>
#include <TParticle.h>

#include "AliJFCluster.h"
#include "AliJFClusterDifference.h"


ClassImp(AliJFClusterDifference)

AliJFClusterDifference::AliJFClusterDifference() :
                     fNLastMerge(0),fI(NULL),fJ(NULL),fDij(0)
{
}

AliJFClusterDifference::AliJFClusterDifference(const AliJFClusterDifference &copy)
{
  fNLastMerge=copy.GetNLastMerge();
  fDij=copy.GetDij();
  fI=copy.GetI();
  fJ=copy.GetJ();
}

AliJFClusterDifference::AliJFClusterDifference(Float_t dij_, AliJFCluster *i_,AliJFCluster *j_) :
                     fNLastMerge(0),fI(i_),fJ(j_),fDij(dij_)
{
}

AliJFClusterDifference::AliJFClusterDifference(AliJFCluster *i_,AliJFCluster *j_) :
  fNLastMerge(0),fI(i_),fJ(j_),fDij(0)
{
  SetValues(i_,j_);
}

Float_t AliJFClusterDifference::SetValues(AliJFCluster *i_,AliJFCluster *j_)
{
  fI=i_;
  fJ=j_;
  fDij=0;
  fNLastMerge=0;

  if(IsValidPointer()&&IsValid()){
    //cout << "valid";
    if(IsDiagonal()){
      fDij=fI->GetPt2();
      fNLastMerge=fI->GetNMerge();
      //cout << " and diagonal" << endl;
    } else {
      Float_t ret1=fI->GetPt2D();
      Float_t ret2=fJ->GetPt2D();
      if(ret1>ret2) fDij=ret2;
      else fDij=ret1;

      Float_t diff1=fI->GetY()-fJ->GetY();
      Float_t diff2=TMath::Abs(fI->GetPhi()-fJ->GetPhi());
      if(diff2>TMath::Pi()) diff2=2*TMath::Pi()-diff2;
      fDij*=(diff1*diff1+diff2*diff2);
    }  
  }

  return fDij;
}

AliJFClusterDifference& AliJFClusterDifference::operator=(const AliJFClusterDifference &copy)
{
  fDij=copy.GetDij();
  fI=copy.GetI();
  fJ=copy.GetJ();
  return *this;
}

ostream& operator<<(ostream &o, const AliJFClusterDifference &j)
{
  o << j.fDij << ": " << j.fI << " " << j.fJ << endl;

  return o;
}

