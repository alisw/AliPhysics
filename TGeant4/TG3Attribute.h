// $Id$
// Category: visualization

#ifndef TG3_ATTRIBUTE_H
#define TG3_ATTRIBUTE_H

enum TG3Attribute
{
  kWORK,	// Set the volume active/inactive for tracking (not used!)
  kSEEN,	// Set visibility : 0-invisible , 1-visible,
		//		   -1-volume and daughters invisible,
		//		   -2-volume visible but daughters invisible.
  kLSTY,	// Set line style : 1-unbroken(default)
		//		    2-dashed
		//		    3-dotted
		//   	negative values do the same for daughters.
  kLWID,	// Set line width : 
  kCOLO,	// Set colour :	   1-7 -G3 base colours(default=1)
		//	1-black	   n=7+m, m=1,9, grey with increasing luminosity,
		//	2-red	   n=17+m, m=1,25,
		//	3-green	   n=67+m, m=1,25,
		//	4-blue     n=118+m, m=1,25,
		//	5-yellow   n=42+m, m=1,25,
		//	6-violet   n=142+m, m=1,25,
		//	7-turquoise	n=92+m, m=1,25,
		//   	negative values do the same for daughters.
  kFILL,	// Set fill style: 0- forces drawing style to wireframe(default)
		//		   1- forces solid drawing style.
		//   	negative values do the same for daughters.
  kSET,	        // Set number associated to volume name (not used!)
  kDET,	        // Set detector number associated to volume name (not used!)
  kDTYP,	// Set detector type (not used!)
  kUNKNOWN
};

#endif //TG3_ATTRIBUTE_H
