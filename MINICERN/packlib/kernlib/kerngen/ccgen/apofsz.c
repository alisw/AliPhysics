/*
 * Get the file size when using Apogee compiler
 * on a Sun Solaris box. 
 */
#include <sys/stat.h>
#include <stdlib.h>

int
apofsz_(s,nbytes,slen)
char *s;
int *nbytes;
int slen;
{
  int i;
  char *ss = (char *)malloc(slen + 1);
  struct stat sbuf;
  for (i = 0; i < slen; i++) {
    if (s[i] == ' ') break;
    else ss[i] = s[i];
  }
  ss[i] = '\0';
  i = stat(ss,&sbuf);
  free(ss);
  *nbytes = (int)sbuf.st_size;
  return i;
}
