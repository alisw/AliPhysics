#if defined(CERNLIB_WINNT)
  #define gcaddc GCADDC
  #define gcaddf GCADDF
  #define gcaddd GCADDD
  #define gcaddi GCADDI
  #define gcaddl GCADDL
  #define type_of_call _stdcall
#else
  #define gcaddc gcaddc_
  #define gcaddf gcaddf_
  #define gcaddd gcaddd_
  #define gcaddi gcaddi_
  #define gcaddl gcaddl_
  #define type_of_call
#endif

extern "C" char* type_of_call gcaddc(char *arg)
{
  return arg;
}
extern "C" double* type_of_call gcaddd(double *arg)
{
  return arg;
}
extern "C" int*  type_of_call gcaddi(int  *arg)
{
  return arg;
}
extern "C" float* type_of_call gcaddf(float *arg)
{
  return arg;
}
extern "C" int* type_of_call gcaddl(int *arg)
{
  return arg;
}
