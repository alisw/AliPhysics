
//////////////////////////////////////////////////////////////////////////
//                                                                      //
// AliFBigBang                                                          //
//                                                                      //
// helper class to browse generated particles.                          //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <TBrowser.h>
#include <TParticle.h>
#include <TClonesArray.h>
#include "AliFast.h"
#include "AliFBigBang.h"

#include "AliFBrowsable.h"

ClassImp(AliFBigBang)



//_____________________________________________________________________________
AliFBigBang::AliFBigBang() 
                : TNamed("Histograms","Generated particles browser")
{
   fBrowsables = 0;
}

//_____________________________________________________________________________
AliFBigBang::~AliFBigBang() 
{
  if (fBrowsables) {
     fBrowsables->Delete();
     delete fBrowsables;
     fBrowsables = 0;
  }
}

//_____________________________________________________________________________
void AliFBigBang::Browse(TBrowser *b)
{
  /*

  TClonesArray *particles = mcarlo->Fruits();
  Int_t nparticles = particles->GetEntriesFast();
  TParticle *part;
  AliFBrowsable *brow;
  char name[64];
  if (!fBrowsables) fBrowsables = new TObjArray(2*nparticles);
  if (fBrowsables->GetSize() < nparticles) fBrowsables->Expand(nparticles);
  for (Int_t i=0;i<nparticles;i++) {
     part = (TParticle*)particles->UncheckedAt(i);
     if (part->GetMother(i)) continue;
     brow = GetBrowsable(i);
     sprintf(name,"%s_%d",part->GetName(),i);
     brow->SetName(name);
     brow->SetRefObject(part);
     b->Add(brow,name);
  }
  */
}

//_____________________________________________________________________________
AliFBrowsable *AliFBigBang::GetBrowsable(Int_t i)
{
  AliFBrowsable *brow = (AliFBrowsable*)fBrowsables->At(i);
  if (!brow) {
     brow = new AliFBrowsable();
     fBrowsables->AddAt(brow, i);
     brow->SetBigBang(this);
  } 
  return brow;
}







