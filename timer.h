#include "core.h"

class Timer {
public:
	Timer();
	void StartTimer();
	double StopTimer();
private:
	double startTime_;
	double stopTime_;
	double elapsedTime_;
};
