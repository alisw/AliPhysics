#ifndef ALIPOISSONCALCULATOR_H
#define ALIPOISSONCALCULATOR_H
#include <TNamed.h>
class TH2D;
class TBrowser;


class AliPoissonCalculator : public TNamed
{
public:
  AliPoissonCalculator(); 
  AliPoissonCalculator(const char*);
  virtual ~AliPoissonCalculator();
  
  void SetEtaLumping(UShort_t n) { fEtaLumping = n; } //*MENU*
  void SetPhiLumping(UShort_t n) { fPhiLumping = n; } //*MENU*

  void Reset(const TH2D* base);
  void Fill(Double_t eta, Double_t phi, Bool_t hit, Double_t weight=1);
  void Result(TH2D* output);

  void IsFolder() const { return kTRUE; }
  void Print(const Option_t* option="") const;
  void Browse(TBrowser* b);
protected:
  AliPoissonCalculator(const AliPoissonCalculator& o);
  AliPoissonCalculator& operator=(const AliPoissonCalculator& o);

  UShort_t fEtaLumping;
  UShort_t fPhiLumping;
  TH2D*    fTotal;    // Total number of strips in a region
  TH2D*    fEmpty;    // Total number of strips in a region
  TH2D*    fBasic;    // Total number basic hits in a region

  ClassDef(AliPoissonCalculator,1) // Calculate N_ch using Poisson
};

// Local Variables:
//   mode: C++
// End:

