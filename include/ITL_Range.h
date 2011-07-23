/*
 * ITL_Range.h
 *
 *  Created on: Jul 21, 2011
 *      Author: leeten
 */

#ifndef ITL_RANGE_H_
#define ITL_RANGE_H_

struct CRange{
	double dMin;
	double dMax;

	CRange()
	{
		dMin = -HUGE_VAL;
		dMax = +HUGE_VAL;
	}

	void _Set(double dMin, double dMax)
	{
		this->dMin = dMin;
		this->dMax = dMax;
	}

	double DClamp(double d) const
	{
		return min(max(d, dMin), dMax);
	}
};
#endif /* ITL_RANGE_H_ */
