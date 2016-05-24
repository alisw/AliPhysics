// -*- C++ -*-
// $Id$

#ifndef _ALI_TRIGGER_STUDY_0STP_
#define _ALI_TRIGGER_STUDY_0STP_

#include <bitset>
#include <set>
#include <vector>
#include <map>

#include <TBits.h>
#include <TObject.h>

class AliTriggerStudy0STP : public TObject {
public:
  struct TableZ {
     std::bitset<20> vtx;
     std::pair<std::bitset<20>, std::bitset<20> > p;

     friend bool operator<(const TableZ& t1, const TableZ& t2) {
       return ((t1.vtx.to_ullong()      < t2.vtx.to_ullong()) ||
         (t1.p.first.to_ullong()  < t2.p.first.to_ullong()) ||
         (t1.p.second.to_ullong() < t2.p.second.to_ullong()));
     }
   } ;
   struct TablePhi {
     Int_t phi;
     std::pair<std::bitset<20>, std::bitset<40> > p;

     friend bool operator<(const TablePhi& t1, const TablePhi& t2) {
       return ((t1.phi                  < t2.phi) ||
         (t1.p.first.to_ullong()  < t2.p.first.to_ullong()) ||
         (t1.p.second.to_ullong() < t2.p.second.to_ullong()));
     }
   } ;
   struct TableDeltaPhi {
     Int_t phi1, phi2;
     std::bitset<40> p;

     friend bool operator<(const TableDeltaPhi& t1, const TableDeltaPhi& t2) {
       return ((t1.phi1          < t2.phi1) ||
         (t1.phi2          < t2.phi2) ||
         (t1.p.to_ullong() < t2.p.to_ullong()));
     }
   } ;

  AliTriggerStudy0STP(Int_t deltaPhiMin, Int_t deltaPhiMax=20)
    : TObject()
    , fOpeningAngle(9*(deltaPhiMin-3))
  {
    MakeTableZ();
    MakeTablePhi();
    MakeTableDeltaPhi(deltaPhiMin, deltaPhiMax);
  }
  virtual ~AliTriggerStudy0STP() {}

  Float_t GetOpeningAngle() const { return fOpeningAngle; }

  Bool_t          CheckForTrackletPair(const TBits* mult, std::bitset<40> &phi) const;
  std::bitset<20> CheckForVertex(const TBits* mult, const std::bitset<40> &phi, std::vector<int>& v) const;
  
protected:

  void MakeTableZ();
  void MakeTablePhi();
  void MakeTableDeltaPhi(Int_t deltaPhiMin, Int_t deltaPhiMax);

  static std::bitset<20> ExtractPhi_L0(const TBits *mult);
  static std::bitset<40> ExtractPhi_L1(const TBits *mult);
  static std::bitset<20> ExtractZ_L0(const TBits *mult, Int_t k);
  static std::bitset<20> ExtractZ_L1(const TBits *mult, Int_t k);
  
private:
  AliTriggerStudy0STP(const AliTriggerStudy0STP& );
  AliTriggerStudy0STP& operator=(const AliTriggerStudy0STP& );

  Float_t                    fOpeningAngle;        // in degrees
  std::vector<TableZ>        fLookupTableZ;        //
  std::vector<TablePhi>      fLookupTablePhi;      //
  std::vector<TableDeltaPhi> fLookupTableDeltaPhi; //

  ClassDef(AliTriggerStudy0STP, 1);
} ;

#endif // _ALI_TRIGGER_STUDY_0STP_
