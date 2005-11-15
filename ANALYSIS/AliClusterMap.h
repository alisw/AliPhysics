#ifndef ALICLUSTERMAP_H
#define ALICLUSTERMAP_H

///////////////////////////////////////////////////////////////////////////
//
// class AliClusterMap
// class that describes cluster occupation at TPC
// Each padraw has a corresponding bit in fPadRawMap
// more info: http://aliweb.cern.ch/people/skowron/analyzer/index.html
// Piotr.Skowronski@cern.ch
//
///////////////////////////////////////////////////////////////////////////


#include <TObject.h>
#include <TBits.h>

class AliTPCtrack;
class AliESDtrack;
class TBits;

class AliClusterMap: public TObject
{
  public:
   AliClusterMap();
   AliClusterMap(AliTPCtrack* track);
   AliClusterMap(AliESDtrack* track);
   virtual ~AliClusterMap(){}
   Float_t GetOverlapFactor(const AliClusterMap& clmap) const;
   Bool_t  HasClAtPadRow(Int_t i) const { return fPadRawMap.TestBitNumber(i);}
   void    Print(const Option_t * opt) const {TObject::Print(opt);}
   void    Print() const;
  protected:
  private:
   TBits    fPadRawMap;//bit vector of length 150 correspondind to total number of padraws in TPC
   static const Int_t fgkNPadRows; // Number of pad rows
   ClassDef(AliClusterMap,1)
};

#endif 
