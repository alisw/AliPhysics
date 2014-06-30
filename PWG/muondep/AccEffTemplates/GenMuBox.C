#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include "TRandom.h"
#include "AliGenerator.h"
#include "AliGenBox.h"
#endif

/// Simple 3D-BOX generator for single muons
/// with a fixed fraction of + and - (50% per default)

class AliGenMuBox : public AliGenBox
{
public:
  
  AliGenMuBox(Float_t plusShare=0.50);

  virtual ~AliGenMuBox() {}
  
  void GenerateN(Int_t ntimes);

  void Generate() { GenerateN(1); }

private:
  Float_t fPlusShare; // Fraction of plus muons
  
  ClassDef(AliGenMuBox,1) // Square box random generator for muons (+ and -)
};

ClassImp(AliGenMuBox)

AliGenMuBox::AliGenMuBox(Float_t plusShare) : AliGenBox(), fPlusShare(plusShare)
{
  if ( fPlusShare <= 0.0 )
  {
    fPlusShare = 0.0;
  }
  if ( fPlusShare > 1.0 )
  {
    fPlusShare = 1.0;
  }
}

void AliGenMuBox::GenerateN(Int_t ntimes)
{
  Int_t ipart = 13;
  
  if ( fPlusShare == 1.0 )
  {
    ipart = -13;
  }
  else
  {
    Float_t x = Rndm();
  
    if ( x < fPlusShare )
    {
      ipart = -13;
    }
  }

  SetPart(ipart);
  
  AliGenBox::GenerateN(ntimes);
}

AliGenerator* GenMuBox()
{
  AliGenBox* generator = new AliGenMuBox;
  
  generator->SetNumberParticles(1);
  
  generator->SetPtRange(VAR_GENMUBOX_PTMIN,VAR_GENMUBOX_PTMAX);
  generator->SetYRange(VAR_GENMUBOX_YMIN,VAR_GENMUBOX_YMAX);
  
  generator->SetPhiRange(0., 360.);
  generator->SetTrackingFlag(1);
  
  return generator;
}
