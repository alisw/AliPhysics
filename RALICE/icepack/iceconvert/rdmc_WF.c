
/******** implements functions for the WF structure *********/

#include <string.h>
#include <stdlib.h>

#include "rdmc.h"


void rdmc_init_WF(waveform *wf){
  wf->id = wf->om = wf->ch = wf->pa = RDMC_NA ;
  wf->ndigi   = 0;
  wf->t_start = wf->t_bin = 0.;
  wf->baseline =   RDMC_WF_BASELINE_NA ;
  wf->digi=NULL;
}

void rdmc_clear_WF(waveform *wf){

  rdmc_free_WF(wf);
  rdmc_init_WF(wf);
}


void rdmc_free_WF(waveform *wf){

  if (wf->digi !=NULL) free(wf->digi);
    wf->ndigi=0;
    wf->digi=NULL;
}



/****************************************************************************
 * copies the waveform *in to *out by overwriting.
 * If in == out, nothing is done.
 ****************************************************************************/
void rdmc_cp_WF(waveform *ou, waveform *in)
{
 if ((in == NULL) || (ou == NULL) || (in == ou))
    return;

  memcpy(ou, in, sizeof(waveform));

  if (ou->ndigi > 0){

    ou->digi = malloc(in->ndigi*sizeof(float));

    if (!(ou->digi)){
      ou->ndigi = 0;
    } else {
       memcpy(ou->digi, in->digi, in->ndigi*sizeof(float));
    }
  }

} /* rdmc_cp_WF() */
