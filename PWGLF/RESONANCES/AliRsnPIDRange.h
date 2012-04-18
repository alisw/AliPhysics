#ifndef ALIRSNPIDRANGE_H
#define ALIRSNPIDRANGE_H

//
// Class for n-sigma PID cuts.
// ---
// Requires:
//
// 1) the used detector, chosen from an enumeration
// 2) the reference charged particle species, chosen from AliPID enumeration
// 3) a momentum range: outside it, the cut is never passed
//

#include <TObject.h>
#include <Rtypes.h>

class AliRsnPIDRange : public TObject {
public:

   AliRsnPIDRange(Double_t nsigma=3.0, Double_t pmin=0.0, Double_t pmax=1E20);
   AliRsnPIDRange(const AliRsnPIDRange &copy);
   AliRsnPIDRange &operator=(const AliRsnPIDRange &copy);
   virtual ~AliRsnPIDRange() { }

   Double_t PMin() const      {return fPMin;}
   Double_t PMax() const      {return fPMax;}
   Double_t NSigmaCut() const {return fNSigmaCut;}

   Bool_t IsInRange(Double_t mom)  {return (mom >= fPMin && mom <= fPMax);}
   Bool_t CutPass(Double_t nsigma) {return (nsigma <= fNSigmaCut);}

private:

   Double_t fPMin;      // lower bound of momentum range
   Double_t fPMax;      // upper bound of momentum range
   Double_t fNSigmaCut; // cut in number of sigmas

   ClassDef(AliRsnPIDRange,1)
};

#endif
