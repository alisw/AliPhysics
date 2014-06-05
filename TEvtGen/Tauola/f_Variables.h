#ifndef _f_TauolaVariables_h_included_
#define _f_TauolaVariables_h_included_

/**
 * This file contains definitions of tauola fortran routines and common
 * blocks so they can be access by the C++ code.
 *
 * @author Nadia Davidson
 * @date June 17 2008
 */

namespace Tauolapp
{

extern "C" {
  extern struct { //particle masses
    float amtau;
    float amnuta;
    float amell;
    float amnue;
    float ammu;
    float amnumu;
    float ampiz;
    float ampi;
    float amro;
    float gamro;
    float ama1;
    float gama1;
    float amk;
    float amkz;
    float amkst;
    float gamkst;
  } parmas_;

  extern struct {
    int jak1;
    int jak2;
    int jakp;
    int jakm;
    int ktom;
  } jaki_;


  extern struct {
    double xk0dec;
    int   itdkrc;
  } taurad_;

  extern struct {
    float gamprt[30];
    int   jlist[30];
    int   nchan;
  } taubra_;

  extern struct {
    float bra1,brk0,brk0b,brks;
  } taukle_;

  //extern float amas4_(float*);
  //extern void bostr3_(float*, float*, float*);
  extern void filhep_(int * N, int * IST, int * ID,
                      int * JMO1, int * JMO2, int * JDA1, int * JDA2, 
                      float P4[4], float * PINV, bool * PHFLAG);

  extern void luhepc_(float flag=2);
  extern void lulist_(float flag=2);
  
  // Initialization of RChL currents
  // Dummy for default CLEO installation
  extern void inirchl_(int *flag);

}

} // namespace Tauolapp
#endif
