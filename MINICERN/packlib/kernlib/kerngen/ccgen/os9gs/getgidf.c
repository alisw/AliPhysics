/*
 * $Id$
 *
 * $Log$
 * Revision 1.1.1.1  1996/02/15 17:49:31  mclareni
 * Kernlib
 *
 */
/*>    ROUTINE GETGIDF
  CERN PROGLIB# Z265    GETGIDF         .VERSION KERNOS9  1.01  940727
  ORIG. 01/04/93, JS
  Fortran interface routine to _os_id (). 
  Replaces call to getgid for OS-9 systems.
  MOD.  27/07/94, MM
*/
#include <process.h>
#include <types.h>

void getgidf_(int *gid)
{
      process_id      proc_id;
      u_int16         priority;
      u_int16         age;
      int32           schedule;
      u_int16         group;
      u_int16         user;

      error_code      err;

      err = _os_id (&proc_id, &priority, &age, &schedule, &group, &user);
      *gid = group;
      return;
}
/*> END <----------------------------------------------------------*/
