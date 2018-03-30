#include "AliMagFast.h"
#include "AliMagFastDip5k.h"
typedef unsigned short ushort;
const float Infinity = INFINITY;
AliMagFast::SegmentEnd_t dip5k_z_slices[] = {
};
AliMagFast::SegmentSearch_t dip5k_z_segments[] = {
};
AliMagFast::SegmentSearch_t dip5k_z = {-1, 0.0, -1760.0, dip5k_z_slices, dip5k_z_segments};
AliMagFast::ChebFormula_t dip5k_params[] = {
};
