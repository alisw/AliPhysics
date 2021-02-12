#ifndef ALIVVZERO_H_
#define ALIVVZERO_H_

namespace {

    class AliVVZERO {
    public:
        Double_t _multiplicity_sum;
        Double_t GetMultiplicity(Int_t i) const
        {
            return i >= 0 && i < 64 ?
                _multiplicity_sum * (1.0 / 64.0) : 0;
        }
    };

}

#endif // ALIVVZERO_H_
