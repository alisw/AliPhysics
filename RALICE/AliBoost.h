#ifndef ALIBOOST_H
#define ALIBOOST_H
///////////////////////////////////////////////////////////////////////////
// Class AliBoost
// Perform various Lorentz transformations.
//
// Example :
// =========
//
// Float_t a[3]={0.1,0.2,0.3};
// Ali3Vector beta;
// beta.SetVector(a,"car");
//
// AliBoost b1;
// b1.SetBeta(beta);
// cout << " Boost b1 :" << endl;
// b1.Info();
//
// Float_t b[4]={14,1,2,3};
// Ali4Vector p;
// p.SetVector(b,"car");
// Ali4Vector pprim=b1.Boost(p);
// cout << " Boost b1 result p,pprim :" << endl;
// p.Info();
// pprim.Info();
//
// p=b1.Inverse(pprim);
// cout << " Inverse b1 result pprim,p :" << endl;
// pprim.Info();
// p.Info();
//
// Float_t c[4]={5,0,0,4};
// Ali4Vector q;
// q.SetVector(c,"car");
//
// AliBoost b2;
// b2.Set4Momentum(q);
// cout << " Lorbo b2 : " << endl;
// b2.Info("sph");
//
//--- NvE 14-may-1996 UU-SAP Utrecht
//--- Modified : NvE 01-apr-1999 UU-SAP Utrecht using Ali3Vector/Ali4Vector
///////////////////////////////////////////////////////////////////////////
 
#include <iostream.h>
#include <math.h>
 
#include "TObject.h"

#include "Ali4Vector.h" 

class AliBoost : public TObject
{
 public:
  AliBoost();                             // Default constructor
  ~AliBoost();                            // Default destructor
  void SetBeta(Ali3Vector b);             // Set boost parameters by beta 3-vector
  void SetGamma(Double_t g,Ali3Vector v); // Set boost parameters by gamma and direction 3-vector
  void Set4Momentum(Ali4Vector& p);       // Set boost parameters by 4-momentum
  Ali3Vector GetBetaVector();             // Provide the beta 3-vector
  Double_t GetBeta();                     // Provide norm of beta 3-vector
  Double_t GetGamma();                    // Provide gamma value
  void Info(TString f="car");             // Print boost parameter info in coord. frame f
  Ali4Vector Boost(Ali4Vector& v);        // Perform Lorentz boost on 4-vector v
  Ali4Vector Inverse(Ali4Vector& v);      // Perform inverse Lorentz boost on 4-vector v
 
 protected:
  Ali3Vector fBeta; // The beta 3-vector
  Double_t fBeta2;  // beta**2
  Double_t fGamma;  // The gamma factor
 
 ClassDef(AliBoost,1) // Class definition to enable ROOT I/O
};
#endif
