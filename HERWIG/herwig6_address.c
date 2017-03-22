#  define herwig6_addressc herwig6_addressc_
#  define herwig6_addressf herwig6_addressf_
#  define herwig6_addressi herwig6_addressi_
#  define herwig6_addressd herwig6_addressd_
#  define herwig6_addressl herwig6_addressl_
#  define herwig6_addressdc herwig6_addressdc_
#  define type_of_call

struct dbcomplex {double r, i;};

char* type_of_call herwig6_addressc(char *arg)
{
  return arg;
}
int*  type_of_call herwig6_addressi(int  *arg)
{
  return arg;
}
float* type_of_call herwig6_addressf(float *arg)
{
  return arg;
}
double* type_of_call herwig6_addressd(double *arg)
{
  return arg;
}

int* type_of_call herwig6_addressl(int *arg)
{
  return arg;
}

struct dbcomplex* type_of_call herwig6_addressdc(struct dbcomplex *arg)
{
  return arg;
}
