// -*- C++ -*-
// $Id$

#ifndef _ALI_TRIGGER_STUDY_0STP_
#define _ALI_TRIGGER_STUDY_0STP_

#include <bitset>
#include <set>
#include <vector>
#include <map>

class TBits;
#include <TObject.h>

class AliTriggerStudy0STP : public TObject {
public:

  // OT ... online-tracklet

  typedef enum {    // online-tracklet definition kOT<inner layer>_<outer layer>
    kOT_2_3, // old 0STP is subset
    kOT_1_4, // new
    kOT_1_2  // new for test
  } OTType;

  AliTriggerStudy0STP(OTType otType,
		      Int_t deltaPhiMin)
    : TObject()
    , fOTtype(otType)
    , fDeltaPhiMin(deltaPhiMin)
    , fDeltaPhiMax(otType == kOT_2_3
		   ? 20                  // w.r.t. outer layer
		   : 10)                 // w.r.t. inner layer
    , fOpeningAngle(otType == kOT_2_3    // w.r.t. outer layer
		    ? 9*(deltaPhiMin-3)  // w.r.t. inner layer
		    : 18*(deltaPhiMin-1)) {}

  virtual ~AliTriggerStudy0STP() {}

  AliTriggerStudy0STP& InitTables() {
    for (Int_t i=fDeltaPhiMin; i<=fDeltaPhiMax; ++i)
      fMask[i] = 1;

    if (fOTtype == kOT_2_3) {
      MakeTableZ();
      MakeTablePhi();
      MakeTableDeltaPhi(fDeltaPhiMin, fDeltaPhiMax);
    }
    return *this;
  }

  Float_t GetOpeningAngle() const { return fOpeningAngle; }

  Bool_t          CheckForTrackletPair(const TBits* mult, std::bitset<40> &phi) const;
  // for now only for kOT_2_3
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

  struct TableZ {
     std::bitset<20> vtx;
     std::pair<std::bitset<20>, std::bitset<20> > p;

     friend bool operator<(const TableZ& t1, const TableZ& t2) {
       return ((t1.vtx.to_ullong()      < t2.vtx.to_ullong())     ||
	       (t1.p.first.to_ullong()  < t2.p.first.to_ullong()) ||
	       (t1.p.second.to_ullong() < t2.p.second.to_ullong()));
     }
   } ;
   struct TablePhi {
     Int_t phi;
     std::pair<std::bitset<20>, std::bitset<40> > p;

     friend bool operator<(const TablePhi& t1, const TablePhi& t2) {
       return ((t1.phi                  < t2.phi)                 ||
	       (t1.p.first.to_ullong()  < t2.p.first.to_ullong()) ||
	       (t1.p.second.to_ullong() < t2.p.second.to_ullong()));
     }
   } ;
   struct TableDeltaPhi {
     Int_t phi1, phi2;
     std::bitset<40> p;

     friend bool operator<(const TableDeltaPhi& t1, const TableDeltaPhi& t2) {
       return ((t1.phi1          < t2.phi1)         ||
	       (t1.phi2          < t2.phi2)         ||
	       (t1.p.to_ullong() < t2.p.to_ullong()));
     }
   } ;

  OTType                     fOTtype;              // type of online-tracklet used
  Float_t                    fOpeningAngle;        // in degrees
  Int_t                      fDeltaPhiMin;         //
  Int_t                      fDeltaPhiMax;         //
  std::bitset<21>            fMask;                //!
  std::vector<TableZ>        fLookupTableZ;        //!
  std::vector<TablePhi>      fLookupTablePhi;      //!
  std::vector<TableDeltaPhi> fLookupTableDeltaPhi; //!

  ClassDef(AliTriggerStudy0STP, 4);
} ;

#endif // _ALI_TRIGGER_STUDY_0STP_
