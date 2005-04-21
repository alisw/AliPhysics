
/* implements the array_calib_t structure */
#include <string.h>

#include "rdmc.h"

void rdmc_init_array_calib(array_calib_t *c){
    /* tdc */
    c->t_0 = 0.0;
    c->beta_t = 0.0;
    c->alpha_t = 0.0;
    /* adc */
    c->ped = 0.0;
    c->beta_a = 1.0;
    c->kappa = 0.0;
    /* tot */
    c->ped_tot = 0.0;
    c->beta_tot = 1.0;
    c->kappa_tot = 0.0;

    /* flag */
    c->flag = 0; /* this was for longtime 1 but is now change to 0 
			    because cause this is not valid but just the 
			    default init */ 
}

void rdmc_clear_array_calib(array_calib_t *c){
  rdmc_free_array_calib(c);
  rdmc_init_array_calib(c);
}

void rdmc_free_array_calib(array_calib_t *c){
  /* nothing */
}

void rdmc_init_array_calib_stat(array_calib_stat_t *cs){
  cs->geo=0;
  cs->adc=0;
  cs->tdc=0;
  cs->tot=0;
  cs->utc=0;
}

void rdmc_clear_array_calib_stat(array_calib_stat_t *cs){
  rdmc_free_array_calib_stat(cs);
  rdmc_init_array_calib_stat(cs);
}

void rdmc_free_array_calib_stat(array_calib_stat_t *cs){
  /* nothing */
}

void rdmc_init_array_calib_utc(array_calib_utc_t *cu){
  strcpy(cu->utc_src,"?");
  cu->secs = 0;
  cu->nsecs = 0;
}

void rdmc_clear_array_calib_utc(array_calib_utc_t *cu){
  rdmc_free_array_calib_utc(cu);
  rdmc_init_array_calib_utc(cu);
}

void rdmc_free_array_calib_utc(array_calib_utc_t *cu){
  /* nothing */
}


