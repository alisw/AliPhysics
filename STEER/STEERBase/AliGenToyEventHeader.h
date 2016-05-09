#ifndef ALIGENTOYEVENTHEADER_H
#define ALIGENTOYEVENTHEADER_H

/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliGenEventHeader.h"

class AliGenToyEventHeader : public AliGenEventHeader
{
 public:
    AliGenToyEventHeader(const char* name = "");
    AliGenToyEventHeader(const AliGenToyEventHeader &rhs);
    virtual ~AliGenToyEventHeader();

    AliGenToyEventHeader &operator=(const AliGenToyEventHeader &rhs);

    void SetCentrality(Double_t cent) { fCentrality = cent; }
    void SetValue(const std::string &key, Double_t value) { fParameters[key] = value; }

    Double_t GetCentrality() const { return fCentrality; }
    Double_t GetParameter(const std::string &key) const { return fParameters.at(key); }

protected:
    Float_t fCentrality; // centrality
    std::map<std::string, Float_t> fParameters; // additional parameters

    ClassDef(AliGenToyEventHeader, 0)
};

#endif
