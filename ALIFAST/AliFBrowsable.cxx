
//////////////////////////////////////////////////////////////////////////
//                                                                      //
// AliFBrowsable                                                          //
//                                                                      //
// helper class to browse generated particles.                          //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <TBrowser.h>
#include <TParticle.h>
#include <TClonesArray.h>
#include "AliFast.h"
#include "AliFBrowsable.h"
#include "AliFBigBang.h"
//#include "AliFMCMaker.h"

ClassImp(AliFBrowsable)



//_____________________________________________________________________________
AliFBrowsable::AliFBrowsable() 
{

}

//_____________________________________________________________________________
void AliFBrowsable::Browse(TBrowser *b)
{
  /*
  AliFMCMaker *mcarlo = gAliFast->MCMaker();
  TClonesArray *particles = mcarlo->Fruits();
  Int_t nparticles = particles->GetEntriesFast();
  TParticle *refpart = (TParticle*)fRefObject;
  TParticle *part;
  AliFBrowsable *brow;
  char name[64];
  Int_t iparent;

  for (Int_t i=0;i<nparticles;i++) {
     part = (TParticle*)particles->UncheckedAt(i);
     iparent = part->GetMother(0);
     if (!iparent) continue;
     if (particles->UncheckedAt(iparent-1) != refpart) continue;
     brow = fBigBang->GetBrowsable(i);
     sprintf(name,"%s_%d",part->GetName(),i);
     brow->SetName(name);
     brow->SetRefObject(part);
     b->Add(brow,name);
  }
  */
}















