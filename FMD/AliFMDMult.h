#ifndef ALIFMDMULT_H
#define ALIFMDMULT_H

/* Reconstracted Particles Class: has number of reconstructed
 * particles in sectors from NumOfMinSector to NumberOfMaxSector()
 * rings from NumOfMinRing to NumOfMaxRing for each FMDvolume
 */
#ifndef ROOT_TObject
# include <TObject.h>
#endif

class AliFMDMult: public TObject
{
public:
  enum EMethod {
    kPoission, 
    kIterative, 
    kNaiive
  };
  AliFMDMult(Float_t  particles=0, UShort_t method=kNaiive);
  virtual ~AliFMDMult() {};

  Float_t         Particles() const { return fParticles; }
  UShort_t        Method()    const { return fMethod; }
  virtual Float_t Eta() const = 0;
  virtual Float_t Phi() const = 0;
  virtual void    Print(Option_t* opt="") const;
protected:
  Float_t  fParticles;       // Number of particles 
  UShort_t fMethod;          // Method use to get fParticles

  ClassDef(AliFMDMult,1)     // Base class for multiplicity data
};
#endif
//____________________________________________________________________
//
// Local Variables:
//   mode: C++
// End:
//
// EOF
//
