#ifndef ALIRSNMINIAXIS_H
#define ALIRSNMINIAXIS_H

//
// All implementations related to definition of an axis
// which is used in the output histogams.
// Simpler than TAxis, it defines an array of edges
// which is then ported to the output histogram definition.
// currently ported only in mini-package, but it could
// become a default also for general package.
//

#include "TObject.h"
#include "TArrayD.h"

class AliRsnMiniAxis : public TObject {

public:

   AliRsnMiniAxis(Int_t valID = -1)                                       : fValueID(valID), fBins(0) { }
   AliRsnMiniAxis(Int_t valID, Int_t nbins, Double_t min, Double_t max)   : fValueID(valID), fBins(0) {Set(nbins, min, max);}
   AliRsnMiniAxis(Int_t valID, Double_t min, Double_t max, Double_t step) : fValueID(valID), fBins(0) {Set(min, max, step);}
   AliRsnMiniAxis(Int_t valID, Int_t nbins, Double_t *bins)               : fValueID(valID), fBins(0) {Set(nbins, bins);}
   AliRsnMiniAxis(const AliRsnMiniAxis& copy) : TObject(copy), fValueID(copy.fValueID), fBins(copy.fBins) { }
   AliRsnMiniAxis& operator=(const AliRsnMiniAxis& copy) {fValueID = copy.fValueID; fBins = copy.fBins; return (*this);}
   
   void      SetValueID(Int_t id)     {fValueID = id;}
   Int_t     GetValueID() const       {return fValueID;}
   
   Int_t     NBins()     {return  fBins.GetSize() - 1;}
   TArrayD  *Bins()      {return &fBins;}
   Double_t *BinArray()  {return  fBins.GetArray();}
   
   void      Set(Int_t nbins, Double_t min, Double_t max);
   void      Set(Int_t nbins, Double_t *bins);
   void      Set(Double_t min, Double_t max, Double_t step);
   
private:

   Int_t    fValueID;  // index of used value in task collection
   TArrayD  fBins;     // bins
   
   ClassDef(AliRsnMiniAxis,1)
};

#endif
