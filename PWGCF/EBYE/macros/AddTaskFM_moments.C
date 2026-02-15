#include "AliAnalysisManager.h"
#include "AliAnalysisTaskFM_moments.h"

AliAnalysisTaskFM_moments *AddTaskFM_moments(Int_t _YEAR_ = 2015, Bool_t _IS_MC_ = kFALSE, const char *suffix = "") {
  // Initialize variables
  Int_t _n_pt_bins = 1;
  Double_t _pt_bins[10] = {0.4, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  Int_t _max_m_bins = 82;
  Int_t _bField = 0;
  Int_t _cent_min = 0, _cent_max = 5;
  Bool_t _pileup_rejection = kTRUE, _two_track_QA = kFALSE, _mixed_event = kFALSE;
  Double_t _eta_min = -0.8, _eta_max = 0.8;
  Double_t _vtx_x_max = 0.3, _vtx_y_max = 0.4, _vtx_z_max = 10.0;
  if (_YEAR_ == 2010) _pileup_rejection = kFALSE; // pileup rejection not available for 2010 data

  // all the redundant variables:: not being used now
  Int_t fBit = 128;
  Double_t fDeta = 0.02;
  Double_t fDphi = 0.02;
  Bool_t fTwoTrack = kFALSE;
  Double_t fSharedFraction = 0.05;
  Bool_t fSharity = kFALSE;
  Bool_t fRejectElectrons = kFALSE;
  Int_t fDCAXYCut = -1; // 0 for loose, 1 for tight, -1 for no cut
  Double_t fDCAZRangeMax = 0.0;
  Double_t fITSClusterCut = 0.0;
  Double_t fTPCClusterCut = 50.0;
  Double_t fnTPCrossedrows = 0.0;
  Double_t fSharedCls = 0.0;
  Double_t fSharedRows = 0.0;
  Double_t fFindableCls = 0.0;
  Bool_t defSharedFrac = kFALSE;
  Bool_t fSelfAffAnalysis = kFALSE;
  ////////////
  
  std::cout << "\033[31m == :: == task marooz == :: == " << " with the options: " << "\n"
    << " number of pt bins: " << _n_pt_bins << " :: "<< _pt_bins[0] << "-" << _pt_bins[1] << "-" << _pt_bins[2] << "," << _pt_bins[3] << "-" << _pt_bins[4] << "," <<  _pt_bins[5] << "-" << _pt_bins[6] << "," << _pt_bins[7] << "-" << _pt_bins[8] << "," << _pt_bins[9] << "\n"
    << " number of max M bins: " << _max_m_bins << "\n"
    << " pileup: " << _pileup_rejection << " :: eta range: " << _eta_min << "-" << _eta_max << "\n"
    << " centrality range: " << _cent_min << "-" << _cent_max << "\033[0m" << std::endl;
    
  std::vector<Double_t> _pt_vector;
  for (Double_t _pt : _pt_bins) if (_pt != 0.0) _pt_vector.push_back(_pt);
  if (_pt_vector.size() / 2 != _n_pt_bins) {
    std::cout << "\033[31m" << " == :: == task marooz == :: == " << " ERROR: pt bins not set correctly " << "\033[0m" << std::endl;
    return 0x0;
  }
  TArrayD _pt_array(_n_pt_bins * 2, _pt_vector.data());

  Int_t _m_bins[52], _n_bins[52];
  for (size_t i = 0; i < 52; i++) _m_bins[i] = (_max_m_bins == 123) ? 3 * (i + 2) : 2 * (i + 2);
  TArrayI _m_array(52, _m_bins), _n_array(52, _n_bins);

  AliAnalysisManager *_mgr = AliAnalysisManager::GetAnalysisManager();
  if (!_mgr) return 0x0;

  AliAnalysisTaskFM_moments *_task = new AliAnalysisTaskFM_moments(suffix);
  if (!_task) return 0x0;
  std::cout << "\033[31m == :: == task marooz == :: == " << " created task and adding PHYSICS SELECTION " << "\033[0m" << std::endl;

  if (_YEAR_ == 2010) _task->SelectCollisionCandidates(AliVEvent::kMB); // select minimum bias events for LHC10h
  if (_YEAR_ == 2015) _task->SelectCollisionCandidates(AliVEvent::kINT7); // select minimum bias events for LHC15o
  if (_YEAR_ == 2018) _task->SelectCollisionCandidates(AliVEvent::kINT7);

  TString _GEN_ = "Hijing";

  _task->SetIsMC(_IS_MC_, _GEN_);
  _task->SetYear(_YEAR_);
  _task->SetPSbinning(_m_array, _n_array, _max_m_bins);
  _task->SetRejectPileup(_pileup_rejection);
  _task->SetEta(_eta_min, _eta_max);
  _task->SetCentLim(_cent_min, _cent_max);
  _task->SetVtxCut(_vtx_x_max, _vtx_y_max, _vtx_z_max);
  _task->SetPtArray(_pt_array, _n_pt_bins);
  _task->SetBfield(_bField);
  _task->SetTwoTrackCuts(fDeta, fDphi, fTwoTrack, _two_track_QA);
  _task->SetMixedEvent(_mixed_event);

  // all the redundant initializations:: not being used now
  _task->Setfbit(fBit);
  _task->SetSharingFraction(fSharedFraction, fSharity);
  _task->SetRejectElectrons(fRejectElectrons);
  _task->SetDCAXYRangeMax(fDCAXYCut);
  _task->SetDCAZRangeMax(fDCAZRangeMax);     // 1
  _task->SetITSClusterCut(fITSClusterCut);   // chi2 per ITS 36
  _task->SetTPCClusterCut(fTPCClusterCut);   // chi2 per TPC 4
  _task->SetnTPCrossedrows(fnTPCrossedrows); // 70
  _task->SetSelfAffine(fSelfAffAnalysis);
  _task->SetSharedCls(fSharedCls, fSharedRows, fFindableCls, defSharedFrac);
  ////////////

  _mgr->AddTask(_task);

  const Int_t _NUM_CON_ = 6; // Number of conditions (e.g., 1 or 2)
  const Int_t _N_OUTPUTS_ = _NUM_CON_ * (_n_pt_bins + 1) + (_IS_MC_ ? _n_pt_bins : 0);
  AliAnalysisDataContainer *_output_list[_N_OUTPUTS_];
  
  TString str;
  TString _con_str[_NUM_CON_] = { "fb128", "fb768", "fb768&&rchi2>80", "fb768&&shclsncls<0.3&&shclsncrows<0.25&&nfindableclsncls>0.8", "fb768&&rchi2>80&&shclsncls<0.3&&shclsncrows<0.25&&nfindableclsncls>0.8", "fb768&&pTDCAxy&&rchi2>80&&shclsncls<0.3&&shclsncrows<0.25&&nfindableclsncls>0.8" };

  for (Int_t i = 0; i < 3; i++) std::cout << "\033[31m =================================================== \033[0m" << std::endl;
  
  std::cout << "\033[31m == :: == task marooz == :: == " << " number of output containers: " << _N_OUTPUTS_ << "\n"
            << " number of conditions: " << _NUM_CON_ << "\033[0m" << std::endl;
  
  for (Int_t i = 0; i < _NUM_CON_; i++) std::cout << "\033[31m == :: == task marooz == :: == " << " condition: " << _con_str[i] << "\033[0m" << std::endl;
 
  for (Int_t i = 0; i < 3; i++) std::cout << "\033[31m =================================================== \033[0m" << std::endl;
  
  for (Int_t i = 0; i < _NUM_CON_; i++) {
    if (_IS_MC_) str = Form("list_%dMC_%s_%dpts_%d%dcent_%s", _YEAR_, _con_str[i].Data(), _n_pt_bins, _cent_min, _cent_max, suffix);
    else
    str = Form("list_%dDATA_%s_%dpts_%d%dcent_%s", _YEAR_, _con_str[i].Data(), _n_pt_bins, _cent_min, _cent_max, suffix);
    _output_list[i] = _mgr->CreateContainer(str, TList::Class(), AliAnalysisManager::kOutputContainer, AliAnalysisManager::GetCommonFileName());

    std::cout << "\033[31m == :: == task marooz == :: == " << " created output container: " << str << "\033[0m" << std::endl;
  }
  
  char *name = "ntplist";
  if (_IS_MC_) name = "ntplistRe";
  for (Int_t i = 0; i < _NUM_CON_; i++) {
    for (Int_t nbins = 0; nbins < _n_pt_bins; nbins++) {
      _output_list[_NUM_CON_ + i * _n_pt_bins + nbins] = _mgr->CreateContainer(Form("%s%d%s_cond%d",name, nbins + 1, suffix, i + 1), TList::Class(), AliAnalysisManager::kOutputContainer, AliAnalysisManager::GetCommonFileName());
      std::cout << "\033[31m == :: == task marooz == :: == " << " created output container: " << Form("%s%d%s_cond%d",name, nbins + 1, suffix, i + 1) << "\033[0m" << std::endl;
    } // loop over pt bins
  } // loop over conditions
  // Place ntplistGen (if MC) at the end of the array
  if (_IS_MC_) {
    for (Int_t i = 0; i < _n_pt_bins; i++) {
      _output_list[_NUM_CON_ * (_n_pt_bins + 1) + i] = _mgr->CreateContainer(Form("ntplistGen%d%s", i + 1, suffix), TList::Class(), AliAnalysisManager::kOutputContainer, AliAnalysisManager::GetCommonFileName());
      std::cout << "\033[31m == :: == task marooz == :: == " << " created output container: " << Form("ntplistGen%d%s", i + 1, suffix) << "\033[0m" << std::endl;
    } // loop over pt bins
  } // end of MC condition

  _mgr->ConnectInput(_task, 0, _mgr->GetCommonInputContainer());
  for (Int_t i = 0; i < _N_OUTPUTS_; i++) _mgr->ConnectOutput(_task, i + 1, _output_list[i]);
  std::cout << "\033[31m == :: == task marooz == :: == " << " connected i/p o/p  \033[0m" << std::endl;

  return _task;
}
