
/*
 *  functions for the array_hdef structure 
 */

#include <stdlib.h>
#include <string.h>

#include "rdmc.h"

static int rdmc_exists_hdef_tag(const array_hdef_t *p, 
			       int ndef, const char *tag);



/* Get a unique id for a new defpar object */
int rdmc_unique_hdef_id(const array_hdef_t *p, int ndef) 
{
  int id;
  for (id = ndef; 1; id++)
    if (rdmc_get_hdef_id(p, ndef, id) == RDMC_NA )
      return (id);
}


/* get a unique tag for a new defpar object. It is the tag_root itself
 * as long it is unique, or the tag_root with a number at the end */
int rdmc_unique_hdef_tag(char *tag, const array_hdef_t *p
			, int ndef, const char *tag_root){
  int i;
  
  /* First check if we can use the orginal tag */
  if (!rdmc_exists_hdef_tag(p, ndef, tag_root)) {
    strcpy(tag, tag_root);
    return 0;
  }
  /* now try to construct a new one with "Tag" + "number" */
  for (i = 1; 1; i++) {
    sprintf(tag, "%s%i", tag_root, i);
    if (!rdmc_exists_hdef_tag(p, ndef, tag)) 
      return 0;
  }
}


static int rdmc_exists_hdef_tag(const array_hdef_t *p, int ndef, const char *tag){
  int i;

  for (i = 0; i < ndef; i++){
    if (rdmc_get_hdef_tag(p, ndef, tag) != RDMC_NA) {
      return 1;
    }
  }
  return 0;
}

/* add a new initilized header def */
int  rdmc_add_array_hdef(array_hdef_t **list,
			 int *count, array_hdef_t *new, int ipos){
  int i;

  if( ( list == NULL )
      || (new == NULL )
      || (count == NULL )
      || (*count < 0)
      || (ipos < 0)
      || (ipos > *count)
      )
    return 1;

  (*count)++;
  *list = (array_hdef_t *) 
    realloc(*list,sizeof(array_hdef_t)*(*count));

  for (i = (*count - 1) ; i > ipos ; i-- ){
    rdmc_cp_hdef( &((*list)[i]), &((*list)[i-1]) );
  }
  rdmc_cp_hdef( &((*list)[ipos]), new );

  return 0;
}


/* remove  header def */
int  rdmc_del_array_hdef(array_hdef_t **list, 
			       int *count, int ipos){
  int i;
  if( (list == NULL )
      || (*list == NULL )
      || (*count <= 0)
      || (ipos < 0)
      || (ipos >= *count)
      )
    return 1;
  for (i = ipos ; i < (*count-1) ; i++ ){
    rdmc_cp_hdef( &((*list)[i]), &((*list)[i+1]) );
  }
  (*count)--;
  if (*count > 0)
    *list = (array_hdef_t *) 
      realloc(*list,sizeof(array_hdef_t)*(*count));
  else {
    free(*list);
    *list=NULL;
  }
  return 0;
}



void  rdmc_init_array_hdef(array_hdef_t *def){/*init the structures */
  int i;
  def->id = RDMC_NA;
  strcpy(def->tag,"unknown");
  def->nwords=0;
  for (i =0 ; i < RDMC_MAXTOKEN_PER_LINE ; i++){
    def->words[i][0] = '\0';
  }
  def->npars=0;
  for (i =0 ; i < RDMC_MAXTOKEN_PER_LINE ; i++){
    def->pars[i][0] = '\0';
  }
} /* init_array_header_def */

void  rdmc_clear_array_hdef(array_hdef_t *def){/*init the structures */
  rdmc_free_array_hdef(def);
  rdmc_init_array_hdef(def);
} /* clear_array_header_def */

void  rdmc_free_array_hdef(array_hdef_t *def){/*init the structures */
} /* free_array_header_def */


void  rdmc_cp_hdef(array_hdef_t *out, const array_hdef_t *in){
  memcpy(out,in,sizeof(array_hdef_t));
} /* init_array_header_def */

/* search an array of ndef array_header_def_t elements */
/* find the number of a header definition */
/* returns the number  0..(ndef-1) or RDMC_NA */
int  rdmc_get_hdef_tag(const array_hdef_t *hd, int ndef, const char *tag){
  int i;
  for (i=0 ;  i < ndef ; i++ ){
#if 0
    rdmc_msgprintf("%s:%s",hd[i].tag,tag);
#endif
    if (strcmp(hd[i].tag,tag) == 0)
      return i;
  }
  return RDMC_NA;
}

int rdmc_get_hdef_id(const array_hdef_t *hd, int ndef,int id){ /*according to an id */
  int i;
  for (i=0 ;  i < ndef ; i++ ){
    if (hd[i].id == id)
      return i;
  }
  return RDMC_NA;
}

/* returns 1 if differetn  */
int rdmc_comp_array_hdef(const array_hdef_t *d1, const array_hdef_t *d2){
  int i;
  if(d1->id != d2->id)
    return 1;
  if(strcmp(d1->tag,d2->tag))
    return 1;
  if(d1->nwords != d2->nwords)
    return 1;
  for ( i =0 ; i < d1->nwords ; i++){
    if(strcmp(d1->words[i],d2->words[i]))
      return 1;
  }
  if(d1->npars != d2->npars)
    return 1;
  for ( i =0 ; i < d1->npars ; i++){
    if(strcmp(d1->pars[i],d2->pars[i]))
      return 1;
  }
  return 0;
}

/* returns the index of the string token in the list of hdef parameters */
/* if the token is not found RDMC_NA is returned */
int rdmc_token_in_hdef(const array_hdef_t *p, const char *token){
  register int i;
  for (i=0 ; i < p->nwords ; i++ ){
    if (! strcmp(p->words[i],token) )
      return i;
  }
  return RDMC_NA;
}
