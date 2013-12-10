#ifndef FJ_includes_H
#define FJ_includes_H

// $Id$

#if !defined(__CINT__) && !defined(__MAKECINT__)
#if defined __GNUC__
#pragma GCC system_header
#endif
#include <fastjet/config.h>
#include <fastjet/PseudoJet.hh>
#include <fastjet/JetDefinition.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/ClusterSequenceArea.hh>
#include <fastjet/AreaDefinition.hh>
#include <fastjet/SISConePlugin.hh>
#ifdef FASTJET_VERSION
#include <fastjet/tools/JetMedianBackgroundEstimator.hh>
#endif
#endif
#endif
