/*******************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: EvtGenBase
 *    File: $Id: EvtBlattWeisskopf.hh,v 1.3 2004/09/14 03:24:11 gowdy Exp $
 *  Author: Alexei Dvoretskii, dvoretsk@slac.stanford.edu, 2001-2002
 *
 * Copyright (C) 2002 Caltech
 *******************************************************************************/

// Blatt-Weisskopf penetration form factor for a resonance R->AB.
// Taken from CLEO preprint 00-23 (hep-ex/0011065)

#ifndef EVT_BLATT_WEISSKOPF_HH
#define EVT_BLATT_WEISSKOPF_HH

class EvtBlattWeisskopf {


public:

  EvtBlattWeisskopf(int LL, double R, double p0);
  EvtBlattWeisskopf(const EvtBlattWeisskopf&);
  ~EvtBlattWeisskopf();

  double operator()(double p) const;

private:

  int    _LL;   // angular momentum of daughters
  double _radial;    // resonance radial parameter
  double _p0;

  double _F0;   // formula evaluated at _p0
  double compute(double p) const;

};

#endif


