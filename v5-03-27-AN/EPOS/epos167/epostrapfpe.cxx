//
// from: http://www.fortran-2000.com/ArnaudRecipes/CompilerTricks.html
//

// Starting with glibc 2.2, the following C99-style (but glibc specific)
// code is preferred.

// #define _GNU_SOURCE 1
// #include <fenv.h>
//  static void __attribute__ ((constructor)) trapfpe(void)
//  {
//    /* Enable some exceptions.  At startup all exceptions are masked. */
//    feenableexcept(FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
//  }



// Previous versions of glibc require the following code. 
//                             trapfpe.c

#include <fpu_control.h>
static void __attribute__ ((constructor)) trapfpe(void)
{
 fpu_control_t cw = 
  _FPU_DEFAULT & ~(_FPU_MASK_IM | _FPU_MASK_ZM | _FPU_MASK_OM);
  _FPU_SETCW(cw);
  /* On x86, this expands to: */
  /* unsigned int cw = 0x037f & ~(0x01 | 0x04 | 0x08); */
  /* __asm__ ("fldcw %0" : : "m" (*&cw));              */ 
}
