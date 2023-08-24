#pragma once

#include <cstdlib>
#include <cmath>
#include <cstring>

namespace GCLF
{

#ifdef BYTE_ORDER
# undef WORDS_BIG_ENDIAN
# undef WORDS_LITTLE_ENDIAN
# if BYTE_ORDER == BIG_ENDIAN
#  define WORDS_BIG_ENDIAN 1
# endif
# if BYTE_ORDER == LITTLE_ENDIAN
#  define WORDS_LITTLE_ENDIAN 1
# endif
#endif

union ieee754_double
{
	double d;
	/* This is the IEEE 754 double-precision format. */
	struct
	{
#ifdef WORDS_BIGENDIAN
		unsigned int negative : 1;
		unsigned int exponent : 11;
		/* Together these comprise the mantissa.  */
		unsigned int mantissa0 : 20;
		unsigned int mantissa1 : 32;
#else
		/* Together these comprise the mantissa.  */
		unsigned int mantissa1 : 32;
		unsigned int mantissa0 : 20;
		unsigned int exponent : 11;
		unsigned int negative : 1;
#endif
	} ieee;
};

union ieee754_float
{
	float f;
	/* This is the IEEE 754 float-precision format. */
	struct
	{
#ifdef WORDS_BIGENDIAN
		unsigned int negative : 1;
		unsigned int exponent : 8;
		unsigned int mantissa : 23;
#else
		unsigned int mantissa : 23;
		unsigned int exponent : 8;
		unsigned int negative : 1;
#endif
	} ieee;
};

inline ieee754_double to_i754d(double d)
{
	ieee754_double i754d;
	i754d.d = d;
	return i754d;
}

inline double epsilon_diff(ieee754_double i754d)
{
	ieee754_double new_i754d = i754d;
	new_i754d.ieee.mantissa1 ^= 0x1u;		// reverse lowest bit
	return fabs(i754d.d - new_i754d.d);
}

inline double epsilon_diff(double d)
{
	ieee754_double  i754d = to_i754d(d);
	return epsilon_diff(i754d);
}
}// namespace GCLF