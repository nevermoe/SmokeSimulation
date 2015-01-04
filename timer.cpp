#include "core.h"
#include "timer.h"


Timer::Timer()
{
	_startTime = _stopTime = _elapsedTime = 0;
}

void Timer::StartTimer()
{
	_startTime = glfwGetTime();
}

double Timer::StopTimer()
{
	_stopTime = glfwGetTime();
	return (_stopTime - _startTime);
}
