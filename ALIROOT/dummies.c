#ifdef __linux
#include <fpu_control.h>
void __attribute__ ((constructor))
trapfpe () {
  (void) __setfpucw (_FPU_DEFAULT &
                     ~(_FPU_MASK_IM | _FPU_MASK_ZM | _FPU_MASK_OM));
}
void MAIN__()  {}
#endif
void izrtoc_() {}
void igmess_() {}
void igloc2_() {}
void igpxmp_() {}
void izitoc_() {}
