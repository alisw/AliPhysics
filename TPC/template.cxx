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

/* $Id$ */

static int PropagateTo(TVector &x, Double_t fX, Double_t xk) {
  if (TMath::Abs(x(2)*xk - x(3)) >= 0.999) {
    cerr<<"Propagation failed !\n";
    return 0;
  }

  Double_t x1=fX, x2=x1+(xk-x1), dx=x2-x1;//, y1=x(0), z1=x(1);
  Double_t c1=x(2)*x1 - x(3), r1=sqrt((1.-c1)*(1.+c1));
  Double_t c2=x(2)*x2 - x(3), r2=sqrt((1.-c2)*(1.+c2));
  
  x(0) += dx*(c1+c2)/(r1+r2);
  x(1) += dx*(c1+c2)/(c1*r2 + c2*r1)*x(4);

  return 1;
}

static int Rotate(TVector &x, Double_t fX, Double_t alpha) {
  Double_t x1=fX, y1=x(0);
  Double_t ca=cos(alpha), sa=sin(alpha);
  Double_t r1=x(2)*fX - x(3);
  
  fX = x1*ca + y1*sa;
  x(0)=-x1*sa + y1*ca;
  x(3)=x(3)*ca + (x(2)*y1 + sqrt((1.-r1)*(1.+r1)))*sa;
  
  Double_t r2=x(2)*fX - x(3);
  if (TMath::Abs(r2) >= 0.999) {
    //cerr<<"Rotation failed !\n";
    return 0;
  }
  
  Double_t y0=x(0) + sqrt((1.-r2)*(1.+r2))/x(2);
  if ((x(0)-y0)*x(2) >= 0.) {
    //cerr<<"Rotation failed !!!\n";
    return 0;
  }

  return 1;
} 

//_____________________________________________________________________________
static int templ(TVector par, Double_t x, Double_t dy, Double_t dz, 
                    const AliTPCSector *sec, int s, int rf=0) 
{
  //-----------------------------------------------------------------
  // This function tries to find a track prolongation.
  //
  // Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch
  //-----------------------------------------------------------------
  int accepted=0;
  const int ROWS_TO_SKIP=100;
  int try_again=ROWS_TO_SKIP;
  Double_t alpha=sec->GetAlpha();
  int ns=int(2*TMath::Pi()/alpha+0.5);

  Double_t xold=x;
  for (int nr=sec->GetRowNumber(x)-1; nr>=rf; nr--,xold=x) {
    //cerr<<nr<<endl;
    Double_t x=sec->GetX(nr), ymax=sec->GetMaxY(nr);
    if (!PropagateTo(par,xold,x)) return 0;

    AliTPCcluster *cl=0;
    const AliTPCRow& row=sec[s][nr];
    Double_t road=dy, y=par(0), z=par(1);

    if (road>30) {
      return 0;
    }

    if (row) {
      for (int i=row.Find(y-road); i<row; i++) {
	AliTPCcluster* c=(AliTPCcluster*)(row[i]);
	if (c->fY > y+road) break;
	if (TMath::Abs(c->fZ - z) > dz) continue;
	cl=c;       
      }
    }
    if (cl) {
      par(0)=cl->fY; par(1)=cl->fZ;
      //dy=TMath::Sqrt(cl->fSigmaY2); dz=TMath::Sqrt(cl->fSigmaZ2);
      //cerr<<nr<<' '<<cl->fTracks[0]<<' '<<cl->fTracks[1]<<' '<<cl->fTracks[2]<<endl;
      accepted++;
      try_again=ROWS_TO_SKIP;
    } else {
      if (try_again==0) break;
      if (y > ymax) {
         s = (s+1) % ns;
         if (!Rotate(par,x,alpha)) return 0;
         dy*=2;
      } else if (y <-ymax) {
         s = (s-1+ns) % ns;
         if (!Rotate(par,x,-alpha)) return 0;
         dy*=2;
      }
      try_again--;
    }
  }

  //cerr<<endl;
  return accepted;
}


