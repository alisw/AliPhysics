//
// Class to check the reconstruction versus the generated
// 
// TODO:
// - add a lot of histograms (much more on tracks vs mc particle) 
// - add documentation
//


#ifndef ROOT_TObject
#include "TObject.h"
#endif
#ifndef ROOT_TTree
#include "TTree.h"
#endif
#ifndef ROOT_TH3
#include "TH3.h"
#endif
#ifndef ROOT_TParticle
#include "TParticle.h"
#endif

#include "AliESDtrack.h"

class AliRecoVsGeneCheck : public TObject
{

public:
  AliRecoVsGeneCheck();
  
  void   Event(Double_t* vtx, Double_t* vtx_res, Double_t* mcvtx, Int_t n_part=-1);

  void   Track(AliESDtrack* esdTrack, TParticle* mcParticle);

  void   SetNPartAxisTitle(Char_t* t) 
    {fhVtxZResVsNPart  ->SetXTitle(t); 
    fhVtxDzNormVsNPart->SetXTitle(t);
    fhNPart->SetXTitle(t);}

  void   SaveHistograms(Char_t* dir="reco_vs_gene");

protected:

  // event specific histograms
  TH2F* fhVtzZRecoVsMC;     // z pos of vertex reco vs MC
  
  TH1F* fhVtxZRes;          // estimated z res. of vertex
  TH2F* fhVtxZResVsZ;       // estimated z res. of vertex vs z
  TH2F* fhVtxZResVsNPart;   // estimated z res. of vertex vs n part
  
  TH1F* fhVtxDzNorm;        // (z_mc - z_reco)/estimated res.
  TH2F* fhVtxDzNormVsZ;     // (z_mc - z_reco)/estimated res. vs vtx z
  TH2F* fhVtxDzNormVsNPart; // (z_mc - z_reco)/estimated res. vs n part.
  
  TH1F* fhVtxZMC;           // vtz z mc 
  TH1F* fhVtxZReco;         // vtx z reco
  
  TH1F* fhNPart;            // n charged primary particles 

  // track specific histograms
  TH3F* fhDPtVsPtVsEta;
  TH3F* fhDEtaVsPtVsEta;  

  private:
    AliRecoVsGeneCheck(const AliRecoVsGeneCheck&);
    AliRecoVsGeneCheck& operator=(const AliRecoVsGeneCheck&);

  ClassDef(AliRecoVsGeneCheck, 1)
};

