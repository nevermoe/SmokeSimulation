#include "core.h"
#include "timer.h"


Timer::Timer()
{
	startTime_ = stopTime_ = elapsedTime_ = 0;
}

void Timer::StartTimer()
{
	startTime_ = glfwGetTime();
}

double Timer::StopTimer()
{
	stopTime_ = glfwGetTime();
	return (stopTime_ - startTime_);
}
