
/* implement functions special for the f2k 1.2 format */
#include <stdlib.h>

#include "rdmc.h"
#include "amanda.h"
#include "f2k.h"

const f2000_line_t WF_line_2004x1 = 
  {"WF ", WF_LINE, COMP_STRINGWISE, rdmc_amanda_WF_2004x1 };



const f2000_event_t f2000_mevt_2004x1 =
{ RDMC_EVENT_MUON,
  { &(EM_line_1x1),
    NULL
  },
  {
    &(TR_line_1x1),
    &(HT_line_1x2),  &(CH_line_1x1),    
    &(WF_line_2004x1),
    &(STATUS_line_1x1),   &(US_line_1x1),
    &(FIT_block_1x1),  &(FRESULT_line_1x1),
    &(TRIG_block_1x1), &(USES_line_1x1),
    &(COMMENT_line_1x1)
    ,  NULL
  },
  { &(EE_line_1x1)
    , NULL},
  rdmc_mevt_f2k_1x1,
  rdmc_wevt_f2k_1x1_2
};


const f2000_event_t  * f2000_events_2004x1[] 
  = { 
    &f2000_mevt_2004x1,
    NULL 
  };


int rdmc_amanda_WF_2004x1( mcfile *fp , array *a, mevt *e, void *tmp)
{
  rdmc_f2k_buffer_t *f2k_buff = fp->info.f2000.f2k_buffer;
  char *s=f2k_buff->line_pt[f2k_buff->iline];
  char **args=NULL;
  int nargs=0;
  int i;
  waveform wf;
  float channel;


  /* start parsing*/
  args = rdmc_f2k_tokenize(s,&nargs);                   //check the parsing string

  if ( nargs < 8) return RDMC_ILF;                      // check # 

  rdmc_init_WF(&wf);  
  channel    = rdmc_amanda_strtof(args[1], RDMC_NA);   // fill waveform struct 
  wf.om      = rdmc_nint(channel);
  wf.ch      = (channel - wf.om) * 100.;
  wf.id      = rdmc_amanda_strtoi(args[2] ,RDMC_NA);
  wf.pa      = rdmc_amanda_strtoi(args[3] ,RDMC_PARENT_NA);
  wf.ndigi   = atoi(args[4]);
  wf.t_start = rdmc_amanda_strtof(args[5], RDMC_NA);

  wf.t_bin   = rdmc_amanda_strtof(args[6], RDMC_NA);
  wf.baseline = rdmc_amanda_strtof(args[7], RDMC_WF_BASELINE_NA);

  if ( nargs < 8 + wf.ndigi)  return RDMC_ILF;          // check # 

  if( wf.ndigi <= 0 || wf.ndigi > RDMC_MAXDIGI )
    return RDMC_ILF;
      
  wf.digi   = malloc(wf.ndigi*sizeof(float));           // fill digit. of waveform
  for (i = 0; i < wf.ndigi; i++)
    wf.digi[i] = rdmc_amanda_strtof(args[8+i], RDMC_TDC_NA);


  rdmc_add_WF(e,&wf,e->nwf);                            // add waveform to event
  rdmc_clear_WF(&wf);

  return RDMC_IO_OK;
} /* rdmc_amanda_WF() */

// pd - end


/****************************************************************************
 ********************************** E O F ***********************************
 ****************************************************************************/



