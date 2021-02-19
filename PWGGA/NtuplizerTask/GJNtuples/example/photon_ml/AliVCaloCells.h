#ifndef ALIVCALOCELLS_H_
#define ALIVCALOCELLS_H_

#include <map>

namespace {

    class AliVCaloCells {
    public:
        std::map<Int_t, Double_t> _amplitude;
        Double_t GetCellAmplitude(Int_t i)
        {
            return _amplitude[i];
        }
    };

}

#endif // ALIVCALOCELLS_H_
