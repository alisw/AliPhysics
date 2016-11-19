#ifndef __CINT__
#include "AliAnalysisTaskEmcalJetCDF.h"
#endif

/// \file AddTaskEmcalJetCDF.C
/// \brief Adds a AliAnalysisTaskEmcalJetCDF analysis task and coresponding containers by calling NS_AliAnalysisTaskEmcalJetCDF::AddTaskEmcalJetCDF(...)
///
/// Analysis of jets (all and leading) distribution of pt and multiplicity, R distribution, N{70,75,80,85,90} and Pt{70,75,80,85,90}
///
/// \author Adrian SEVCENCO <Adrian.Sevcenco@cern.ch>, Institute of Space Science, Romania
/// \date Mar 19, 2015

// definition in AliAnalysisTaskEmcalJetCDF.h - namespace NS_AliAnalysisTaskEmcalJetCDF

/// Add a AliAnalysisTaskEmcalJetCDF task - detailed signature
/// \param ntracks name of tracks collection
/// \param nclusters name of clusters collection
/// \param ncells name of EMCAL cell collection
/// \param tag tag name of analysis task
/// \return AliAnalysisTaskEmcalJetCDF* task
AliAnalysisTaskEmcalJetCDF* AddTaskEmcalJetCDF ( const char* ntracks = "usedefault", const char* nclusters = "usedefault", const char* ncells = "usedefault", const char* tag = "CDF" )
{ return NS_AliAnalysisTaskEmcalJetCDF::AddTaskEmcalJetCDF ( ntracks, nclusters, ncells, tag); }

// kate: indent-mode none; indent-width 2; replace-tabs on;
