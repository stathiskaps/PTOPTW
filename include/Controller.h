#include <chrono>
#include <csignal>
#include <thread>
#include "Divider.h"
#include "Processor.h"
#include "TaskManager.h"
#include "Unifier.h"

#ifndef CONTROLLER_H
#define CONTROLLER_H

class Controller {
private:
	int mOpenTime, mCloseTime, mWalks;
	int m_tmax;
	std::vector<TA*> mAttractions;
	std::vector<Point> mPoints;
	TA* m_depot;
protected:
	
public:
	Controller();
	Controller(int openTime, int closeTime, std::vector<TA*> attractions,  std::vector<Point> points, int walks);
	~Controller();
	void exec();
};

#endif