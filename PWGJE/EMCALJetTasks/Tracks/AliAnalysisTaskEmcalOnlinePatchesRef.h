/************************************************************************************
 * Copyright (C) 2015, Copyright Holders of the ALICE Collaboration                 *
 * All rights reserved.                                                             *
 *                                                                                  *
 * Redistribution and use in source and binary forms, with or without               *
 * modification, are permitted provided that the following conditions are met:      *
 *     * Redistributions of source code must retain the above copyright             *
 *       notice, this list of conditions and the following disclaimer.              *
 *     * Redistributions in binary form must reproduce the above copyright          *
 *       notice, this list of conditions and the following disclaimer in the        *
 *       documentation and/or other materials provided with the distribution.       *
 *     * Neither the name of the <organization> nor the                             *
 *       names of its contributors may be used to endorse or promote products       *
 *       derived from this software without specific prior written permission.      *
 *                                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND  *
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED    *
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE           *
 * DISCLAIMED. IN NO EVENT SHALL ALICE COLLABORATION BE LIABLE FOR ANY              *
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES       *
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;     *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND      *
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT       *
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS    *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                     *
 ************************************************************************************/
#ifndef ALIANALYSISTASKEMCALONLINEPATCHESREF_H
#define ALIANALYSISTASKEMCALONLINEPATCHESREF_H

#include <vector>
#include <TString.h>

#include <AliAnalysisTaskEmcal.h>
#include <AliEMCALTriggerDataGrid.h>
#include <AliCutValueRange.h>

class THistManager;
class AliOADBContainer;
class AliEMCALRecoUtils;
class AliEMCALTriggerPatchInfo;

namespace PWGJE {

namespace EMCALJetTasks {

class AliAnalysisTaskEmcalOnlinePatchesRef: public AliAnalysisTaskEmcal {
public:
  AliAnalysisTaskEmcalOnlinePatchesRef();
  AliAnalysisTaskEmcalOnlinePatchesRef(const char *name);
  virtual ~AliAnalysisTaskEmcalOnlinePatchesRef();

  void SetOnlineTriggerClass(const char *triggerclass) { fOnlineTriggerClass = triggerclass; };
  void SetCellTimeCut(double min, double max) { fCellTimeCut.SetLimits(min, max); }

  void DefineMaskedFastorOADB(const char *oadbname) { fNameMaskedFastorOADB = oadbname; }
  void DefineMaskedCellOADB(const char *oadbname) { fNameMaskedCellOADB = oadbname; }

  static AliAnalysisTaskEmcalOnlinePatchesRef *AddTaskOnlinePatchesRef(const char *name, const char *suffix);

protected:
  virtual void UserCreateOutputObjects();
  virtual void UserExecOnce();
  virtual bool Run();
  virtual bool IsTriggerSelected();
  virtual void RunChanged(int newrun);
  
  void LoadCellEnergies();
  void LoadFastorEnergies();

  bool IsCellMasked(int absCellID) const;
  bool IsFastORMasked(int absFastORID) const;
  bool IsPhosHole(int icol, int row) const;

  bool SelectPatch(const AliEMCALTriggerPatchInfo &patch) const;
  std::pair<int,int> GetNumberOfFastors(const AliEMCALTriggerPatchInfo &patch) const;
  void MarkFastorsContributing(const AliEMCALTriggerPatchInfo &patch) const;

private:
  AliAnalysisTaskEmcalOnlinePatchesRef(const AliAnalysisTaskEmcalOnlinePatchesRef &);
  AliAnalysisTaskEmcalOnlinePatchesRef &operator=(const AliAnalysisTaskEmcalOnlinePatchesRef &);

  THistManager                            *fHistos;               //!<! Histogram manager
  TString                                 fOnlineTriggerClass;    ///< Name of the trigger class
  AliEMCALTriggerDataGrid<double>         *fFastOREnergy;         //!<! Grid for energies calculated from FastOR ADCs
  AliEMCALTriggerDataGrid<double>         *fFEEnergy;             //!<! Grid for energies calculated from FEE Energies
  AliEMCALTriggerDataGrid<int>            *fInOnlinePatch;        //!<! Grid indicating FastORs contributing to an online trigger patch
  AliEMCALTriggerDataGrid<int>            *fMaskedCellsFastor;    //!<! Number of masked cells in Fastor
  std::vector<int>                        fMaskedFastors;         ///< List of masked fastors
  AliCutValueRange<double>                fCellTimeCut;  ///< Cell time cut
  TString                                 fNameMaskedFastorOADB;  ///< Name of the OADB container with masked fastors
  TString                                 fNameMaskedCellOADB;    ///< Name of the OADB container with masked cells
  AliOADBContainer                        *fMaskedFastorOADB;     //!<! OADB container with masked fastors
  AliOADBContainer                        *fMaskedCellOADB;       //!<! OADB container with masked cells
  AliEMCALRecoUtils                       *fRecoUtils;            //!<! Reco utils for bad channel handling

  ClassDef(AliAnalysisTaskEmcalOnlinePatchesRef, 1);
};

} /* namespace EMCALJetTasks */

} /* namespace PWGJE */

#endif /* ALIANALYSISTASKEMCALONLINEPATCHESREF_H */
