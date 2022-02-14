#include "AliTimeRangeMasking.h"

AliOADBContainer fContainer("TimeRangeMasking");

using time_type = ULong64_t;
using bitmap_type = UShort_t;
using TimeRangeMask = AliTimeRangeMask<time_type, bitmap_type>;

using mask_reasons_type = std::vector<TimeRangeMask::MaskReason>;

using range_mask_type = std::tuple<time_type, time_type, mask_reasons_type>;

void AddMasking(int run, std::string pass, std::vector<range_mask_type> rangeMasks);


void SetupTimeRangeMasking(const std::string outFile = "$ALICE_PHYSICS_SRC/OADB/COMMON/PHYSICSSELECTION/data/TimeRangeMasking.root")
//void SetupTimeRangeMasking(const std::string outFile = "/tmp/COMMON/PHYSICSSELECTION/data/TimeRangeMasking.root")
{
  // LHC18r
  AddMasking(296749, "pass1", {
      { 628880853571, 684122084392, {TimeRangeMask::kBadTPCPID} }, // hMIP_gid_296749_C06
      });
  AddMasking(296750, "pass1", {
      { 687976579337, 734727672930, {TimeRangeMask::kBadTPCPID} }, // hMIP_gid_296750_C06
      });
  AddMasking(296849, "pass1", {
      { 483958451208, 745436019833, {TimeRangeMask::kBadTPCPID} }, // hMIP_gid_296849_A05
      //{ 611163754863, 745436019833, {TimeRangeMask::kBadTPCPID} }, // hMIP_gid_296849_C06
      });
  AddMasking(296890, "pass1", {
      { 498344702183, 582200169595, {TimeRangeMask::kBadTPCPID} }, // hMIP_gid_296890_A05
      });
  AddMasking(297029, "pass1", {
      { 436213319386, 494857925357, {TimeRangeMask::kBadTPCPID} }, // hMIP_gid_297029_A05
      });
  AddMasking(297194, "pass1", {
      { 508681589448, 648815037601, {TimeRangeMask::kBadTPCPID} }, // hMIP_gid_297194_A05
      });
  AddMasking(297219, "pass1", {
      { 678000000000, 690000000000, {TimeRangeMask::kBadTPCPID} }, // hMIP_gid_297194_C06
      });
  AddMasking(297481, "pass1", {
      { 495708154506, 698497854077, {TimeRangeMask::kBadTPCPID} }, // hMIP_gid_297481_A05
      //{ 622854077253, 643615879828, {TimeRangeMask::kBadTPCPID} }, // hMIP_gid_297481_C06
      });

  AddMasking(296749, "pass3", {
      { 628880853571, 684122084392, {TimeRangeMask::kBadTPCPID} }, // hMIP_gid_296749_C06
      });
  AddMasking(296750, "pass3", {
      { 687976579337, 734727672930, {TimeRangeMask::kBadTPCPID} }, // hMIP_gid_296750_C06
      });
  AddMasking(296849, "pass3", {
      { 483958451208, 745436019833, {TimeRangeMask::kBadTPCPID} }, // hMIP_gid_296849_A05
      //{ 611163754863, 745436019833, {TimeRangeMask::kBadTPCPID} }, // hMIP_gid_296849_C06
      });
  AddMasking(296890, "pass3", {
      { 498344702183, 582200169595, {TimeRangeMask::kBadTPCPID} }, // hMIP_gid_296890_A05
      });
  AddMasking(297029, "pass3", {
      { 436213319386, 494857925357, {TimeRangeMask::kBadTPCPID} }, // hMIP_gid_297029_A05
      });
  AddMasking(297194, "pass3", {
      { 508681589448, 648815037601, {TimeRangeMask::kBadTPCPID} }, // hMIP_gid_297194_A05
      });
  AddMasking(297219, "pass3", {
      { 678000000000, 690000000000, {TimeRangeMask::kBadTPCPID} }, // hMIP_gid_297194_C06
      });
  AddMasking(297481, "pass3", {
      { 495708154506, 698497854077, {TimeRangeMask::kBadTPCPID} }, // hMIP_gid_297481_A05
      //{ 622854077253, 643615879828, {TimeRangeMask::kBadTPCPID} }, // hMIP_gid_297481_C06
      });
  fContainer.WriteToFile(outFile.data());
}

bitmap_type getMask(mask_reasons_type reasons)
{
  bitmap_type mask{};

  for (const auto& reason : reasons) {
    mask |= BIT(reason);
  }

  return mask;
}

void AddMasking(int run, std::string pass, std::vector<range_mask_type> rangeMasks)
{
  auto timeRangeMask = new AliTimeRangeMasking<time_type, bitmap_type>;
  for (auto& rangeMask : rangeMasks) {
    const auto& start = std::get<0>(rangeMask);
    const auto& stop  = std::get<1>(rangeMask);
    const auto& masks = std::get<2>(rangeMask);
    const auto mask = getMask(masks);
    timeRangeMask->AddTimeRangeMask(start, stop, mask);
  }

  std::cout << "Setting up time range mask for run " << run << " pass " << pass<< "\n";
  timeRangeMask->Print();
  std::cout << "\n";

  fContainer.AppendObject(timeRangeMask, run, run, pass.data());
}
