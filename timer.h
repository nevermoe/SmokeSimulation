#ifndef _TIMER_H
#define _TIMER_H

#include "core.h"

class Timer {
public:
	Timer();
	void StartTimer();
	double StopTimer();
private:
	double _startTime;
	double _stopTime;
	double _elapsedTime;
};

#endif
