/*
 * AliDPG - ALICE Experiment Data Preparation Group
 * Utility functions
 *
 */

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/

Int_t
RunToYear(Int_t run)
{
  if      (run < 105524) return 2009;
  else if (run < 139699) return 2010;
  else if (run < 170724) return 2011;
  else if (run < 194309) return 2012;
  else if (run < 200008) return 2013;
  else if (run < 208365) return 2014;
  else if (run < 247171) return 2015;
  else                   return 2016;
}
