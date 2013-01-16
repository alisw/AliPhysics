/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
/* $Id$ */

#ifndef ALIFLOWMSPHISTOGRAMS_H
#define ALIFLOWMSPHISTOGRAMS_H

#include <TNamed.h>
#include <TObjArray.h>
#include <TProfile.h>
#include <TString.h>

class AliFlowMSPHistograms : public TNamed
{
public:
   AliFlowMSPHistograms();                                  // Default constructor for root IO
   AliFlowMSPHistograms(const AliFlowMSPHistograms &x);     // Copy constructor
   AliFlowMSPHistograms(const int dim, const char *name, const int nBins, const double xmin, const double xmax);  // See Init
   AliFlowMSPHistograms(const int dim, const char *name, const int nBins, const double *xBins);                   // See Init
   ~AliFlowMSPHistograms();                                 // Destructor 

   void  Init(const int dim, const char *name, const int nBins, const double xmin, const double xmax);   // Redefine all private members
   void  Init(const int dim, const char *name, const int nBins, const double *xBins);                    // Redefine all private members
   void  SetVarName(const char *name, const int i);                        // Set name of a variable
   void  SetVarName(const TString name, const int i);                      // Set name of a variable
   const char *VarName(const int i)const;                                  // Get name of a variable
   void  SetXName(const char *name){fXName=name;};                         // Set name of x-axis
   const char *XName()const{return fXName.Data();};                        // Get name of x-axis 
   void  Fill(const double x, const double *par, const double *wpar);      // Fill all profiles for bin FindBin(x)

   //Bool_t IsFolder() const { return kFALSE; };     // Tells browser that this is a folder
   //void  Print(Option_t *opt) const;               // Prints a message ...
   //void  Browse(TBrowser *b) {};                   // Allows the browser to see the contents

   int NVars()const{return fPAvg.GetEntriesFast();};              // Number of variables 
   int Nbins()const{return Avg(0)->GetNbinsX();};                 // Number of bins in x
   double XLow()const{return Avg(0)->GetBinLowEdge(1);};          // Lower end of x range
   double XHigh()const{return XLow()+Nbins()*Avg(0)->GetBinWidth(1);};  // Upper end of xrange
   int FindBin(const double x)const{return Avg(0)->FindBin(x);};  // Get bin number from x value

   double Average(const int bin, const int ivar)const;                                          // Get weighted average from the accumulated TProfile bin
   double Average(const int ivar)const{return Average(0,ivar);};                                // Weighted average summed over all bins
   double Covariance(const int bin, const int ivar, const int jvar)const;                       // Calculate covariance of two average variables
   double Covariance(const int ivar, const int jvar)const{return Covariance(0,ivar,jvar);};     // Covariance summing all X bins
   double Variance(const int bin, const int ivar)const{return Covariance(bin, ivar, ivar);};    // Calculate Variance of an average
   double Variance(const int ivar)const{return Covariance(0, ivar, ivar);};                     // Calculate Variance of bin average

   Long64_t Merge(TCollection *);                                                               // Merge all objects in the collection into this
   AliFlowMSPHistograms &operator+=(const AliFlowMSPHistograms &x);                             // Merge this with another AliFlowMSPHistograms object
   AliFlowMSPHistograms &operator=(const AliFlowMSPHistograms &x);

   TObject *Clone(const char *n=0)const;      // Create a copy of this object

private:
   TObjArray   fVarNames;        // Names of variables being correlated
   TObjArray   fPAvg;            // Averages of single variables
   TObjArray   fPProd;           // Products of each combination of two variables
   TString     fXName;           // Name of X axis for differential correlations

   inline TProfile *Avg(const int i)const{return (TProfile *)(fPAvg[i]);};                            // Index in average profile list
   inline int Idx(const int i, const int j)const{return (j<=i?((i*(i+1))>>1)+j:((j*(j+1))>>1)+i);};   // Calculate index in linear array (j<=i)
   inline TProfile *Prod(const int i, const int j)const{return (TProfile *)(fPProd[Idx(i,j)]);};      // Index in product profile list

   double Covariance(double a, double b, double ab, double wa, double wb ,double wab)const;           // Calculate weighted covariance from averages and sums of products

   void UpdateTitles(int i=-1);           // Generate the profile titles from the parameter names and X variable name

   ClassDef(AliFlowMSPHistograms,1);
};

#endif
