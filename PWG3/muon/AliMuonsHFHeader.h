#ifndef ALIMUONSHFHEADER_H
#define ALIMUONSHFHEADER_H

#include <TNamed.h>

#include "AliVHeader.h"
#include "AliVVertex.h"

class AliMuonsHFHeader : public TNamed {
 public :

  AliMuonsHFHeader();
  ~AliMuonsHFHeader();

  void SetHeader(AliVHeader *header);
  void SetVertex(AliVVertex *vertex);

  void SetMultSingleMuon(Int_t mul) {  fMultMuon = mul; }
  void SetMultDimuon(Int_t mul)  {   fMultDimuon = mul; }
  void SetCentrality(Double_t cen) { fCentrality = cen; }

  void GetXYZ(Double_t *pos) const
    { for (Int_t i=0; i<3; i++) pos[i]=fPosition[i]; }
  Double_t GetXv() const { return fPosition[0]; }
  Double_t GetYv() const { return fPosition[1]; }
  Double_t GetZv() const { return fPosition[2]; }
  Int_t GetNContributors() const { return fNContributors; }

  ULong64_t GetTriggerMask() const { return fTriggerMask; }

  Int_t    GetMultSingleMuon() const { return fMultMuon;   }
  Int_t    GetMultDimuon()     const { return fMultDimuon; }
  Double_t GetCentrality()     const { return fCentrality; }

  Bool_t IsUnrecoVertex() const { return fUnrecoVertex; }

 private :

  ULong64_t fTriggerMask;   // trigger mask

  Double32_t fPosition[3];  // position of vtx
  Int_t fNContributors;     // number of contributor for vtx

  Int_t fMultMuon;      // event multiplicity
  Int_t fMultDimuon;
  Double_t fCentrality;     // event centrality class

  Bool_t fUnrecoVertex;

  ClassDef(AliMuonsHFHeader, 1)
};

#endif
