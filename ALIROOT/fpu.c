#ifdef __linux
#include <fpu_control.h>
void __attribute__ ((constructor))
     trapfpe () {
  fpu_control_t cw = _FPU_DEFAULT & ~(_FPU_MASK_IM | _FPU_MASK_ZM |
				      _FPU_MASK_OM);
  _FPU_SETCW(cw);
}
#endif
