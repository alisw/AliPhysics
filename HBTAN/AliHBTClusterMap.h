#ifndef ALIHBTCLUSTERMAP_H
#define ALIHBTCLUSTERMAP_H
//_________________________________________________
///////////////////////////////////////////////////
//
// class AliHBTClusterMap
//
// class that describes cluster occupation at TPC
// Each padraw has a corresponding bit in fPadRawMap
// 
//
// more info: http://aliweb.cern.ch/people/skowron/analyzer/index.html
// Piotr.Skowronski@cern.ch
//
///////////////////////////////////////////////////////////////////////////


#include <TObject.h>
#include <TBits.h>

class AliTPCtrack;
class AliESDtrack;
class TBits;

class AliHBTClusterMap: public TObject
{
  public:
   AliHBTClusterMap();
   AliHBTClusterMap(AliTPCtrack* track);
   AliHBTClusterMap(AliESDtrack* track);
   virtual ~AliHBTClusterMap(){}
   Float_t GetOverlapFactor(const AliHBTClusterMap& clmap) const;
   Bool_t  HasClAtPadRow(Int_t i) const { return fPadRawMap.TestBitNumber(i);}
   void    Print() const;
  protected:
  private:
   TBits    fPadRawMap;//bit vector of length 150 correspondind to total number of padraws in TPC
   static const Int_t fNPadRows;
   ClassDef(AliHBTClusterMap,1)
};

#endif 
