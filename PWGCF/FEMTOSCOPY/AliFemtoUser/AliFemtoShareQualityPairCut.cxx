///
/// \file AliFemtoShareQualityPairCut.cxx
///


#include "AliFemtoShareQualityPairCut.h"
#include <string>
#include <cstdio>
#include <bitset>

#ifdef __ROOT__
ClassImp(AliFemtoShareQualityPairCut)
#endif

//__________________
AliFemtoShareQualityPairCut::AliFemtoShareQualityPairCut():
  fNPairsPassed(0),
  fNPairsFailed(0),
  fShareQualityMax(1.0),
  fShareFractionMax(1.0),
  fRemoveSameLabel(false)
  , fAlternativeAlgorithm(0)
{
  // Default constructor
  // Nothing to do
}
//__________________
AliFemtoShareQualityPairCut::~AliFemtoShareQualityPairCut(){
  /* no-op */
}
AliFemtoShareQualityPairCut& AliFemtoShareQualityPairCut::operator=(const AliFemtoShareQualityPairCut& cut)
{
  if (this != &cut) {
    AliFemtoPairCut::operator=(cut);
    fNPairsPassed = 0;
    fNPairsFailed = 0;
    fShareQualityMax = cut.fShareQualityMax;
    fShareFractionMax = cut.fShareFractionMax;
    fRemoveSameLabel = cut.fRemoveSameLabel;
  }

  return *this;
}
//__________________
bool AliFemtoShareQualityPairCut::Pass(const AliFemtoPair* pair)
{
  // Check for pairs that are possibly shared/double reconstruction
  bool passes = true;

  const auto *track1 = pair->Track1()->Track(),
             *track2 = pair->Track2()->Track();

  if (fRemoveSameLabel) {
    if (abs(track1->Label()) == abs(track2->Label())) {
  //       cout << "Found a pair with same label " << pair->Track1()->Track()->Label() << endl;
  //       cout << "Quality Sharity Passed " << hsmval << " " << hsfval << " " << pair->QInv() << " " << passes << endl;
      passes = false;
    }
  }

  // only check fraction & quality if not already cut and bounds have been placed
  if (passes && (fShareFractionMax < 1.0 || fShareQualityMax < 1.0)) {

    Int_t nh = 0;
    Int_t an = 0;
    Int_t ns = 0;

    const unsigned int n_bits = track1->TPCclusters().GetNbits();

    const auto tpc_clusters_1 = track1->TPCclusters(),
               tpc_clusters_2 = track2->TPCclusters();

    const auto tpc_sharing_1 = track1->TPCsharing(),
               tpc_sharing_2 = track2->TPCsharing();

    // equivalent methods of determining the share quality and fraction
    if (fAlternativeAlgorithm == 1) {

      const auto cls1_and_cls2 = tpc_clusters_1 & tpc_clusters_2,
                 cls1_xor_cls2 = tpc_clusters_1 ^ tpc_clusters_2,
                 shr1_and_shr2 = tpc_sharing_1 & tpc_sharing_2,
                 not_sharing = ~*const_cast<TBits*>(&shr1_and_shr2);
                              // ^ weird const-cast to get around bug in ROOT5

      int cls1_xor_cls2_bits = cls1_xor_cls2.CountBits();

      nh = cls1_xor_cls2_bits + 2 * cls1_and_cls2.CountBits();
      ns = 2 * (cls1_and_cls2 & shr1_and_shr2).CountBits();
      an = cls1_xor_cls2_bits + ns / 2 - (cls1_and_cls2 & not_sharing).CountBits();

    } else if (fAlternativeAlgorithm == 2) {
      if (n_bits != 159) {
        std::cerr << "n_bits: " << n_bits << " != 159\n";
      }

      std::bitset<159> a, b, c, d;
      tpc_clusters_1.Get(reinterpret_cast<ULong64_t*>(&a));
      tpc_clusters_2.Get(reinterpret_cast<ULong64_t*>(&b));
      tpc_sharing_1.Get(reinterpret_cast<ULong64_t*>(&c));
      tpc_sharing_2.Get(reinterpret_cast<ULong64_t*>(&d));

      const auto cls1_and_cls2 = a & b,
                 cls1_xor_cls2 = a ^ b,
                 sharing = c & d,
                 not_sharing = ~sharing;

      const int cls1_xor_cls2_bits = cls1_xor_cls2.count();

      nh = cls1_xor_cls2_bits + 2 * cls1_and_cls2.count();
      ns = 2 * (cls1_and_cls2 & sharing).count();
      an = cls1_xor_cls2_bits + ns / 2 - (cls1_and_cls2 & not_sharing).count();

    } else {
      // Below is the original algorithm

    for (unsigned int imap = 0; imap < n_bits; imap++) {
        const bool cluster_bit_1 = tpc_clusters_1.TestBitNumber(imap),
                   cluster_bit_2 = tpc_clusters_2.TestBitNumber(imap);
        // If both have clusters in the same row
        if (cluster_bit_1 && cluster_bit_2) {
           // Do they share it ?
           if (tpc_sharing_1.TestBitNumber(imap) && tpc_sharing_2.TestBitNumber(imap))
           {
              //	  cout << "A shared cluster !!!" << endl;
              //	cout << "imap idx1 idx2 "
              //	     << imap << " "
              //	     << tP1idx[imap] << " " << tP2idx[imap] << endl;
              an++;
              nh+=2;
              ns+=2;
           }
           // Different hits on the same padrow
           else {
              an--;
              nh+=2;
           }
        }
        else if (cluster_bit_1 || cluster_bit_2) {
           // One track has a hit, the other does not
           an++;
           nh++;
        }
     }
    }

     Float_t hsmval = 0.0;
     Float_t hsfval = 0.0;

     if (nh > 0) {
        hsmval = an*1.0/nh;
        hsfval = ns*1.0/nh;
     }
    //  if (hsmval > -0.4) {
    //   cout << "Pair quality: " << hsmval << " " << an << " " << nh << " "
    //        << (pair->Track1()->Track()) << " "
    //        << (pair->Track2()->Track()) << endl;
    //   cout << "Bits: " << pair->Track1()->Track()->TPCclusters().GetNbits() << endl;
    //  }
    //   if (hsfval > 0.0) {
    //     cout << "Pair sharity: " << hsfval << " " << ns << " " << nh << "    " << hsmval << " " << an << " " << nh << endl;
    //   }

    if (fShareQualityMax < 1.0) {
      passes &= (hsmval < fShareQualityMax);
    }
    if (fShareFractionMax < 1.0) {
      passes &= (hsfval < fShareFractionMax);
    }
  }

  // increment appropriate counter
  // passes ? fNPairsPassed++ : ++;
  (passes ? fNPairsPassed : fNPairsFailed)++;
  return passes;
}
//__________________
AliFemtoString AliFemtoShareQualityPairCut::Report()
{
  // Prepare the report from the execution
  TString report("AliFemtoShareQuality Pair Cut - remove shared and split pairs\n");
  report += TString::Format("Number of pairs which passed:\t%ld  Number which failed:\t%ld\n", fNPairsPassed, fNPairsFailed);
  return AliFemtoString((const char *)report);
}

//__________________
TList *AliFemtoShareQualityPairCut::ListSettings()
{
  // return a list of settings in a writable form
  TList *tListSetttings = new TList();

  tListSetttings->AddVector(
    new TObjString(
      TString::Format("AliFemtoShareQualityPairCut.sharequalitymax=%f", fShareQualityMax)
    ),
    new TObjString(
      TString::Format("AliFemtoShareQualityPairCut.sharefractionmax=%f", fShareFractionMax)
    ),
    new TObjString(
      TString::Format("AliFemtoShareQualityPairCut.fRemoveSameLabel=%d", fRemoveSameLabel)
    ),
    nullptr);

  return tListSetttings;
}
