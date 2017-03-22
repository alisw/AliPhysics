  #define _GNU_SOURCE
    #include <stdio.h>
    #include <stdlib.h>
    #include <fenv.h>

    void enable_exceptions_()
    {
         int retval;
         /* feenableexcept returns the previous exceptions that were  enabled
            on success, otherwise it returns -1
         */
         retval=feenableexcept( FE_DIVBYZERO | FE_INVALID |  FE_OVERFLOW | FE_UNDERFLOW );
         if ( retval == -1 )
         {
             fprintf(stderr, "Warning: call to feenableexcept() failed \n");
         }
    }

    /* This second routine is for Fortran compilers such as g77 and  pathf90
       which follow the f2c name mangling style
    */
    void enable_exceptions__()
    {
        enable_exceptions_();
    }
