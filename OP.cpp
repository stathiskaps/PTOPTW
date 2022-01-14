#include "OP.h"

OP::OP() {
	mDepot = nullptr;
	mStartTime = OPEN_DAY_TIME;
	mEndTime = MAX_MINUTES;
}

//Constructor uses the received unvisited TAs, it doesn't create new ones
OP::OP(std::vector<TA*> unvisitedVec, std::vector<std::vector<double>> ttMatrix, TA* depot, double startTime, double endTime) {
	ListTA unvisited = ListTA(unvisitedVec);
	//mDepot = (unvisited.grabUsingPos(0));	//remove first node from unvisited and store it as depot
	mTtMatrix = ttMatrix;
	mDepot = depot;
	mStartTime = startTime;
	mEndTime = endTime;
	mProcessSolution = Solution(depot, unvisited, mStartTime, mEndTime);
}

OP::OP(std::vector<TA*> unvisitedVec, Walk walk, std::vector<std::vector<double>> ttMatrix, double startTime, double endTime) {
	ListTA unvisited = ListTA(unvisitedVec);
	//mDepot = (unvisited.grabUsingPos(0));	//remove first node from unvisited and store it as depot
	mTtMatrix = ttMatrix;
	mStartTime = startTime;
	mEndTime = endTime;
	mProcessSolution = Solution(walk, unvisited, mStartTime, mEndTime);
}

OP::~OP() {

}

std::tuple<int, double, int, int> OP::getBestPos(TA* ta) {

	double minShift = DBL_MAX;
	int length = mProcessSolution.mWalk.getLength();
	double shift;
	int pos = 1;
	std::tuple<int, double, int, int> res = std::make_tuple(-1, -1, -1, -1);
	
	
	TA* temp = mProcessSolution.mWalk.first();

	TA* curr = temp->next;

	Sight* s = nullptr;
	Route* r = nullptr;

	int tempArrPointId, tempDepPointId;

	if (length < 2) {
		std::cerr << "Invalid length of route" << std::endl;
		return res;
	}
	else {

		if ((s = dynamic_cast<Sight*>(ta))) {
			while (curr != nullptr) {
				shift = mTtMatrix[temp->depPointId][s->point.id] + s->visitDuration + mTtMatrix[s->point.id][curr->arrPointId] - mTtMatrix[temp->depPointId][curr->arrPointId];

				if (shift <= curr->maxShift) {	//check if insertion is possible
					if (shift < minShift) {
						minShift = shift;
						std::get<0>(res) = pos;
						std::get<1>(res) = minShift;
						std::get<2>(res) = s->point.id;
						std::get<3>(res) = s->point.id;
					}
				}

				pos++;
				curr = curr->next;
				temp = temp->next;
			}
		}
		else if ((r = dynamic_cast<Route*>(ta))) {

			while (curr != nullptr) {
				//We got two points to check
				int shift1 = mTtMatrix[temp->depPointId][r->edges[0].id] + r->visitDuration + mTtMatrix[r->edges[1].id][curr->arrPointId] - mTtMatrix[temp->depPointId][curr->arrPointId];

				int shift2 = mTtMatrix[temp->depPointId][r->edges[1].id] + r->visitDuration + mTtMatrix[r->edges[0].id][curr->arrPointId] - mTtMatrix[temp->depPointId][curr->arrPointId];

				if (shift1 < shift2) {
					tempArrPointId = r->edges[0].id;
					tempDepPointId = r->edges[1].id;
					shift = shift1;
				}
				else {
					tempArrPointId = r->edges[1].id;
					tempDepPointId = r->edges[0].id;
					shift = shift2;
				}

				if (shift <= curr->maxShift) {	//check if insertion is possible
					if (shift < minShift) {
						minShift = shift;
						std::get<0>(res) = pos;
						std::get<1>(res) = minShift;
						std::get<2>(res) = tempArrPointId;
						std::get<3>(res) = tempDepPointId;
					}
				}

				pos++;
				curr = curr->next;
				temp = temp->next;
			}
		}

	}
	return res;
}

bool OP::validate() {
	std::cout << "Validating route:\t"; mBestSolution.mWalk.print();
	std::string msg = "";
	TA* curr = mBestSolution.mWalk.first();
	if (curr == nullptr) {
		std::cout << "List is empty" << std::endl;
		return false;
	}
	while (curr != nullptr) {
		if (curr->startOfVisitTime != curr->arrTime) {
			std::cout << "StartOfVisitTime(" << curr->startOfVisitTime << ") != arrTime(" << curr->arrTime << ")" << std::endl;
			return false;
		}

		if (curr->depTime != curr->startOfVisitTime + curr->visitDuration) {
			std::cout << "depTime(" << curr->depTime << ") != startOfVisitTime(" << curr->startOfVisitTime << ") + visitDuration(" << curr->visitDuration << ")" << std::endl;
			return false;
		}

		if (curr->next != nullptr) {
			if (curr->next->arrTime != curr->depTime + mTtMatrix[curr->depPointId][curr->next->arrPointId]) {
				std::cout << "next->arrTime(" << curr->next->arrTime << ") != depTime(" << curr->depTime << ") + timeTravel[" << curr->depPointId << "][" << curr->next->arrPointId << "](" << mTtMatrix[curr->depPointId][curr->next->arrPointId] << ")" << std::endl;
				return false;
			}
		}
		curr = curr->next;
	}
	return true;
}

bool OP::updateTimes(int startIndex, bool smart) {
	//update times first
	TA* curr = mProcessSolution.mWalk.get(startIndex);

	if (curr == mProcessSolution.mWalk.first()) {
		return false;
	}
	while (curr != nullptr) {
		curr->arrTime = curr->prev->depTime + mTtMatrix[curr->prev->depPointId][curr->arrPointId];
		curr->startOfVisitTime = curr->arrTime;
		curr->depTime = curr->startOfVisitTime + curr->visitDuration;
		curr = curr->next;
	}

	//update max shifts
	curr = mProcessSolution.mWalk.last();
	curr->maxShift = curr->timeWindow.closeTime - curr->depTime;
	curr = curr->prev;
	while (curr != nullptr) {
		if (curr->next != nullptr) {
			curr->maxShift = curr->next->maxShift;
		}
		curr = curr->prev;
	}

	return true;
}

void OP::print(std::string msg) {
	std::cout << "=======================================================" << std::endl;
	std::cout << msg << std::endl;
	std::cout << "processSolution: ";
	mProcessSolution.print();
	std::cout << "bestSolution: ";
	mBestSolution.print();
	std::cout << "startTime: " << mStartTime << std::endl;
	std::cout << "endTime:" << mEndTime << std::endl;
	std::cout << "=======================================================" << std::endl;
}

std::tuple<double, double, double> OP::calcTimeEventCut(){
    double dist;
    double min = DBL_MAX;
    double timePoint = -1;
	TA* curr = mProcessSolution.mUnvisited.first();

	std::vector<double> timeWindowEvents;


	while(curr != nullptr){
		timeWindowEvents.push_back(curr->timeWindow.openTime);
		timeWindowEvents.push_back(curr->timeWindow.closeTime);
		curr = curr->next;
	}

	sort(timeWindowEvents.begin(), timeWindowEvents.end());
	timeWindowEvents.erase(unique(timeWindowEvents.begin(), timeWindowEvents.end()), timeWindowEvents.end());

	double average = reduce(timeWindowEvents.begin(), timeWindowEvents.end()) / timeWindowEvents.size();

	curr = mProcessSolution.mUnvisited.first();
	while(curr != nullptr){
 		dist = abs(curr->timeWindow.openTime - average);
        if(dist < min){
            min = dist;
            timePoint = curr->timeWindow.openTime;
        }
        dist = abs(curr->timeWindow.closeTime - average);
        if(dist < min){
            min = dist;
            timePoint = curr->timeWindow.closeTime;
        }
		curr = curr->next;
	}

    return {timeWindowEvents.front(), timePoint, timeWindowEvents.back()};
}

void OP::setBucketActivityDurations() {
	auto [minEvent, avgEvent, maxEvent] = calcTimeEventCut();

	TA* curr = mProcessSolution.mWalk.first();
	while (curr != nullptr) {
		curr->bucketActivities[0].duration = avgEvent - curr->timeWindow.openTime;
		curr->bucketActivities[1].duration = curr->timeWindow.closeTime - avgEvent;
		curr = curr->next;
	}

}


void OP::Local(double avgEvent) {

	TA* curr = mProcessSolution.mUnvisited.first();
	TA* temp;
	
	size_t walk_length = mProcessSolution.mWalk.getLength();
	if (walk_length == 2) {
		ListTA UnvisitedA, UnvisitedB;
		
		while (curr != nullptr) {
			temp = curr;
			curr = curr->next;
			if (curr->bucketActivities[0].duration >= curr->bucketActivities[1].duration) {
				UnvisitedA.push(temp);
			}
			else {
				UnvisitedB.push(temp);
			}
			
		}



	}
	else if (walk_length > 2) {

	}
	else {
		std::cerr << "Invalid walk length" << std::endl;
		std::exit(1);
	}

}

////appends list2 after list1 (list1 + vector2)
//ListTA operator+(ListTA& l1, ListTA& l2) {
//	ListTA l3;
//
//	TA* curr = l1.first();
//	while(curr != nullptr){
//		l3.pushNew(curr);
//		curr = curr->next;
//	}
//
//	curr = l2.first();
//	while(curr != nullptr){
//		l3.pushNew(curr);
//		curr = curr->next;
//	}
//
//	return l3;
//}

bool OP::compareByProfit(const TA a, const TA b) {
    return a.profit < b.profit;
}

int OP::LocalSearch(double avgPoint) {
	double minShift;
	int pos, bestPos;
	double maxRation = DBL_MIN, ratio;
	TA* curr = mProcessSolution.mUnvisited.first();

	std::vector<TA*> UnvisitedA, UnvisitedB;

	while(curr != nullptr){
		double leftDuration = avgPoint - curr->timeWindow.openTime;
        double rightDuration = curr->timeWindow.closeTime - avgPoint;

        if(leftDuration > rightDuration){
			TA* temp = curr;
            UnvisitedA.push_back(temp);
        } else {
			TA* temp = curr;
            UnvisitedB.push_back(temp);
        }
		curr = curr->next;
	}

	size_t walkLength = mProcessSolution.mWalk.getLength();

	if (walkLength == 2) {
		OPTW subproblemA = OPTW(UnvisitedA, mTtMatrix, mDepot, mStartTime, avgPoint);
		Walk walkA = subproblemA.Construct();
		TA* lastNode = walkA.get(walkA.getLength() - 2);
		Walk temp = walkA.copyPart(walkA.getLength() - 2, walkA.getLength() - 1);
		temp.last()->reset();
		temp.last()->timeWindow.openTime = avgPoint;
		temp.last()->timeWindow.closeTime = mEndTime;
		temp.last()->maxShift = mEndTime - avgPoint;
		/*temp.first()->reset();*/

		
		OPTW subproblemB = OPTW(UnvisitedB, temp, mTtMatrix, avgPoint, mEndTime);

		//HERE after the construct, walkB remains with the startings node, i don't think it should. The time space betweeen the two nodes is too small
		Walk walkB = subproblemB.Construct();

		walkA.removeLast();
		walkB.removeFirst();
		//walkB.removeLast();

		//mProcessSolution.mUnvisited = subproblemA.mProcessSolution.mUnvisited + subproblemB.mProcessSolution.mUnvisited;
		//mProcessSolution.mWalk = walkA + walkB;
		mProcessSolution.mScore = mProcessSolution.mWalk.collectProfit();

	}
	else if (walkLength > 2) {
		
		curr = mProcessSolution.mWalk.first();
		auto [wa, wb] = mProcessSolution.mWalk.splitOnDepTime(avgPoint);

		size_t wa_length = wa.getLength(), wb_length = wb.getLength();

		Walk walkA, walkB;

		OPTW subproblemA = OPTW(UnvisitedA, wa, mTtMatrix, mStartTime, avgPoint);
		OPTW subproblemB = OPTW(UnvisitedB, wb, mTtMatrix, avgPoint, mEndTime);

		if (wa_length == 1) {
			walkA = wa;
		}
		else if (wa_length > 1) {
			TA* last = subproblemA.mProcessSolution.mWalk.last();
			last->timeWindow.closeTime = avgPoint;
			last->maxShift = last->timeWindow.closeTime - last->timeWindow.openTime;
			walkA = subproblemA.Construct();
		}

		if (wb_length == 1) {
			walkB = wb;
		}
		else if (wb_length > 1) {
			//subproblemB.mProcessSolution.mWalk.first()->arrTime = walkA.last()->depTime + mTtMatrix[walkA.last()->depPointId][subproblemB.mProcessSolution.mWalk.first()->arrPointId];
			walkB = subproblemB.Construct();
		}

		//walkB.removeFirst();
		//mProcessSolution.mWalk = walkA + walkB;
		updateTimes(1, true);
		//mProcessSolution.mUnvisited = subproblemA.mProcessSolution.mUnvisited + subproblemB.mProcessSolution.mUnvisited;
		mProcessSolution.mScore = mProcessSolution.mWalk.collectProfit();

	}
	else {
		std::cerr << "Invalid walk length" << std::endl;
		std::exit(1);
	}

	if (!mProcessSolution.mWalk.isValid(true, mEndTime, mTtMatrix)) {
		std::cout << "Process Solution Walk is not valid" << std::endl;
		updateTimes(0, true);
		if (!mProcessSolution.mWalk.isValid(true, mEndTime, mTtMatrix)) {
			std::cout << "Process Solution Walk is still not valid LAWL" << std::endl;
			std::exit(1);
		}
	}

	return mProcessSolution.mScore;

}

Walk OP::Construct(){
	double minShift;
	int pos, bestPos = DEFAULT_POS;
	double maxRatio = DBL_MIN, ratio;
	std::string nodeId;
	TA* nodeToInsert = nullptr;
	TA* curr;

	int bestArrPointId = DEFAULT_POINT_ID, bestDepPointId = DEFAULT_POINT_ID;
	int arrPointId = DEFAULT_POINT_ID, depPointId = DEFAULT_POINT_ID;

	while (true)
	{
		maxRatio = DBL_MIN;
		nodeId = DEFAULT_TA_ID;
		curr = mProcessSolution.mUnvisited.first();


		while (curr != nullptr) {
			minShift = DBL_MAX;

			std::tie(pos, minShift, arrPointId, depPointId) = getBestPos(curr);

			if (pos == -1) { //insertion is not feasible
				curr = curr->next;
				continue;
			}

			ratio = pow(curr->profit, 2) / minShift;	//Ratioi = (Si)^2/Shifti

			if (maxRatio < ratio) //check each ratio
			{
				maxRatio = ratio;
				nodeId = curr->id;
				bestPos = pos;
				bestArrPointId = arrPointId;
				bestDepPointId = depPointId;
			}

			curr = curr->next;
		}

		if (nodeId == DEFAULT_TA_ID) {
			break;
		}
		
		nodeToInsert = mProcessSolution.mUnvisited.grabUsingID(nodeId);
		if(nodeToInsert != nullptr){
			mProcessSolution.mWalk.insertAt(nodeToInsert, bestPos, bestArrPointId, bestDepPointId, mTtMatrix); //add it to route
			updateTimes(bestPos, true);
		} else {
			std::cerr << "Didn't find node " << nodeId << std::endl;
			std::exit(1);
		}
	}
	return mProcessSolution.mWalk;
}

int OP::Insert() {
	double minShift;
	int pos, bestPos = DEFAULT_POS;
	double maxRatio = DBL_MIN, ratio;
	std::string nodeId;
	//List<TA> unvisited = sol.mUnvisited.copy();
	TA* nodeToInsert = nullptr;
	TA* curr;

	int bestArrPointId = DEFAULT_POINT_ID, bestDepPointId = DEFAULT_POINT_ID;
	int arrPointId = DEFAULT_POINT_ID, depPointId = DEFAULT_POINT_ID;

	while (true)
	{
		maxRatio = DBL_MIN;
		nodeId = DEFAULT_TA_ID;
		curr = mProcessSolution.mUnvisited.first();
		//numOfUnvisitedNodes = ProcessSolution.UnvisitedPOIs.GetNumOfNodes();

		while (curr != nullptr) {
			minShift = DBL_MAX;
			std::tie(pos, minShift, arrPointId, depPointId) = getBestPos(curr);

			if (pos == -1) { //insertion is not feasible
				curr = curr->next;
				continue;
			}

			ratio = pow(curr->profit, 2) / minShift;	//Ratioi = (Si)^2/Shifti

			if (maxRatio < ratio) //check each ratio
			{
				maxRatio = ratio;
				nodeId = curr->id;
				bestPos = pos;
				bestArrPointId = arrPointId;
				bestDepPointId = depPointId;
			}

			curr = curr->next;
		}

		//if we have reached a point where there is no node which can be inserted, then we stop the search
		if (nodeId == DEFAULT_TA_ID) {
			break;
		}
		
		nodeToInsert = mProcessSolution.mUnvisited.grabUsingID(nodeId);
		if(nodeToInsert != nullptr){
			mProcessSolution.mWalk.insertAt(nodeToInsert, bestPos, bestArrPointId, bestDepPointId, mTtMatrix); //add it to route
			updateTimes(bestPos, true);
		} else {
			std::cerr << "Didn't find node " << nodeId << std::endl;
			std::exit(1);
		}
		

	}

	return mProcessSolution.mWalk.collectProfit();
	//mProcessSolution.SetRoute(route);
	//mProcessSolution.SetScore(route.collectProfit());

}

void OP::Shake(int S, int R, int numOfPois) {
	if(numOfPois == 0){
		return;
	} else if (numOfPois > 0) {
		List<TA> nodes = mProcessSolution.mWalk.grabPart(S, std::min(S+R-1, numOfPois+1));
		mProcessSolution.mUnvisited.append(nodes);
		updateTimes(1, false);
	} else {
		std::cerr << "Route has invalid length" << std::endl;
		exit(1);
	}
}

void OP::SaveSolution(Solution sol) {
	mBestSolution.mUnvisited = sol.mUnvisited.copy();
	mBestSolution.mWalk = sol.mWalk.copy();
	mBestSolution.mScore = sol.mScore;
}

Solution OP::solve() {

	int counter = 0;
	int S = 1, R = 1;
	int timesNotImproved = 0;
	TaskManager* taskManager = TaskManager::GetInstance();
	int bestScore = INT_MIN;

	setBucketActivityDurations();

	auto [minPoint, avgPoint, maxPoint] = calcTimeEventCut();

	while (timesNotImproved < MAX_TIMES_NOT_IMPROVED) {
		counter++;
		std::cout << "this is revision: " << counter << std::endl;

		int score = LocalSearch(avgPoint);

		//int score = Insert();

		if (score > bestScore) {
			bestScore = score;
			mBestSolution.reset();
			mBestSolution.copy(mProcessSolution);
			mBestSolution.SetScore(score);
			R = 1;
			timesNotImproved = 0;

		}
		else {
			timesNotImproved++;
		}

		int numOfPois = mProcessSolution.mWalk.getLength() - 2;
		if (S > numOfPois) {
			S = 1;
			// S = S - numOfPois;
			// if (S <= 0)
			// 	S = 1;
		}
		if (R > numOfPois/2){
			R = 1;
		}

		Shake(S, R, numOfPois);

		S += R;
		R += 1;
		

	}

	bool valid = validate();
	std::string msg = valid ? "yes" : "no";
	std::cout << "valid solution? " << msg << std::endl;

	return mBestSolution;
}


bool OPTW::validate() {
	std::cout << "Validating route:\t"; mBestSolution.mWalk.print();
	std::string msg = "";

	TA* curr = mBestSolution.mWalk.first();
	if (curr == nullptr) {
		std::cout << "List is empty" << std::endl;
		return false;
	}
	while (curr != nullptr) {
		if (curr->startOfVisitTime != curr->arrTime + curr->waitDuration) {
			std::cout << "StartOfVisitTime(" << curr->startOfVisitTime << ") != arrTime(" << curr->arrTime << ")" << std::endl;
			return false;
		}

		if (curr->depTime != curr->startOfVisitTime + curr->visitDuration) {
			std::cout << "depTime(" << curr->depTime << ") != startOfVisitTime(" << curr->startOfVisitTime << ") + visitDuration(" << curr->visitDuration << ")" << std::endl;
			return false;
		}

		if (curr->depTime > curr->timeWindow.closeTime) {
			std::cout << "depTime(" << curr->depTime << ") > closeTime(" << curr->timeWindow.closeTime << ")" << std::endl;
			return false;
		}

		if (curr->next != nullptr) {
			if (curr->next->arrTime != curr->depTime + mTtMatrix[curr->depPointId][curr->next->arrPointId]) {
				std::cout << "next->arrTime(" << curr->next->arrTime << ") != depTime(" << curr->depTime << ") + timeTravel[" << curr->depPointId << "][" << curr->next->arrPointId << "](" << mTtMatrix[curr->depPointId][curr->next->arrPointId] << ")" << std::endl;
				return false;
			}
		}

		curr = curr->next;
	}

	return true;
}


bool TOPTW::validate() {
	std::cout << "Validating route:\t"; mBestSolution.mWalk.print();
	std::string msg = "";

	TA* curr = mBestSolution.mWalk.first();
	if (curr == nullptr) {
		std::cout << "List is empty" << std::endl;
		return false;
	}
	while (curr != nullptr) {
		if (curr->startOfVisitTime != curr->arrTime + curr->waitDuration) {
			std::cout << "StartOfVisitTime(" << curr->startOfVisitTime << ") != arrTime(" << curr->arrTime << ")" << std::endl;
			return false;
		}

		if (curr->depTime != curr->startOfVisitTime + curr->visitDuration) {
			std::cout << "depTime(" << curr->depTime << ") != startOfVisitTime(" << curr->startOfVisitTime << ") + visitDuration(" << curr->visitDuration << ")" << std::endl;
			return false;
		}

		if (curr->depTime > curr->timeWindow.closeTime) {
			std::cout << "depTime(" << curr->depTime << ") > closeTime(" << curr->timeWindow.closeTime << ")" << std::endl;
			return false;
		}

		if (curr->next != nullptr) {
			if (curr->next->arrTime != curr->depTime + mTtMatrix[curr->depPointId][curr->next->arrPointId]) {
				std::cout << "next->arrTime(" << curr->next->arrTime << ") != depTime(" << curr->depTime << ") + timeTravel[" << curr->depPointId << "][" << curr->next->arrPointId << "](" << mTtMatrix[curr->depPointId][curr->next->arrPointId] << ")" << std::endl;
				return false;
			}
		}

		curr = curr->next;
	}

	return true;
}

std::tuple<int, double, int, int> OPTW::getBestPos(TA* ta)
{
	double minShift = DBL_MAX;
	int length = mProcessSolution.mWalk.getLength();
	double arrTime, waitDuration, startOfVisitTime, startOfVisitTime1, startOfVisitTime2, depTime, depTime1, depTime2, shift;
	int pos = 1;
	std::tuple<int, double, int, int> res = std::make_tuple(-1, -1, -1, -1);

	TA* temp = mProcessSolution.mWalk.first();
	TA* curr = temp->next;

	Sight* s = nullptr;
	Route* r = nullptr;

	int tempArrPointId, tempDepPointId;

	if (length < 2) {
		std::cerr << "Invalid length of route" << std::endl;
		return res;
	}
	else {
		if ((s = dynamic_cast<Sight*>(ta))) {

			while (curr != nullptr) {
				arrTime = temp->depTime + mTtMatrix[temp->depPointId][s->point.id];
				waitDuration = std::max(0.0, s->timeWindow.openTime - arrTime);
				startOfVisitTime = arrTime + waitDuration;
				depTime = startOfVisitTime + s->visitDuration;
				shift = mTtMatrix[temp->depPointId][s->point.id] + waitDuration + s->visitDuration + mTtMatrix[s->point.id][curr->arrPointId] - mTtMatrix[temp->depPointId][curr->arrPointId];
				
				if (depTime <= s->timeWindow.closeTime && shift <= curr->maxShift) {
					if (shift < minShift) {
						minShift = shift;
						std::get<0>(res) = pos;
						std::get<1>(res) = minShift;
						std::get<2>(res) = s->point.id;
						std::get<3>(res) = s->point.id;
					}
				}

				pos++;
				curr = curr->next;
				temp = temp->next;
			}

			
		}
		else if ((r = dynamic_cast<Route*>(ta))) {

			while (curr != nullptr) {

				double shift1 = -1, shift2 = -1;
				//We got two points to check
				arrTime = temp->depTime + mTtMatrix[temp->depPointId][r->edges[0].id];
				waitDuration = std::max(0.0, r->timeWindow.openTime - arrTime);
				startOfVisitTime1 = arrTime + waitDuration;
				depTime1 = startOfVisitTime1 + r->visitDuration;
				shift1 = mTtMatrix[temp->depPointId][r->edges[0].id] + waitDuration + r->visitDuration + mTtMatrix[r->edges[1].id][curr->arrPointId] - mTtMatrix[temp->depPointId][curr->arrPointId];

				arrTime = temp->depTime + mTtMatrix[temp->depPointId][r->edges[1].id];
				waitDuration = std::max(0.0, r->timeWindow.openTime - arrTime);
				startOfVisitTime2 = arrTime + waitDuration;
				depTime2 = startOfVisitTime2 + r->visitDuration;
				shift2 = mTtMatrix[temp->depPointId][r->edges[1].id] + waitDuration + r->visitDuration + mTtMatrix[r->edges[0].id][curr->arrPointId] - mTtMatrix[temp->depPointId][curr->arrPointId];

				if (shift1 < shift2) {
					tempArrPointId = r->edges[0].id;
					tempDepPointId = r->edges[1].id;
					startOfVisitTime = startOfVisitTime1;
					depTime = depTime1;
					shift = shift1;
				}
				else {
					tempArrPointId = r->edges[1].id;
					tempDepPointId = r->edges[0].id;
					startOfVisitTime = startOfVisitTime2;
					depTime = depTime2;
					shift = shift2;
				}

				if (depTime <= r->timeWindow.closeTime && shift <= curr->maxShift && shift < minShift) {	//check if insertion is possible
						minShift = shift;
						std::get<0>(res) = pos;
						std::get<1>(res) = minShift;
						std::get<2>(res) = tempArrPointId;
						std::get<3>(res) = tempDepPointId;
				}

				pos++;
				curr = curr->next;
				temp = temp->next;
			}
		
		}
	}
	return res;
}

std::tuple<int, double, int, int> TOPTW::getBestPos(TA* ta) // explicitly name what type to specialize
{
	double minShift = DBL_MAX;
	int length = mProcessSolution.mWalk.getLength();
	double arrTime, waitDuration, startOfVisitTime, depTime, shift;
	int pos = 1;
	std::tuple<int, double, int, int> res = std::make_tuple(-1, -1, -1, -1);

	TA* temp = mProcessSolution.mWalk.first();
	TA* curr = temp->next;

	Sight* s = nullptr;
	Route* r = nullptr;

	if (length < 2) {
		std::cerr << "Invalid length of route" << std::endl;
		return res;
	}
	else {
		if ((s = dynamic_cast<Sight*>(ta))) {
			while (curr != nullptr) {
				arrTime = temp->depTime + mTtMatrix[temp->depPointId][s->point.id];
				waitDuration = std::max(0.0, s->timeWindow.openTime - arrTime);
				startOfVisitTime = arrTime + waitDuration;
				depTime = startOfVisitTime + s->visitDuration;
				shift = mTtMatrix[temp->depPointId][s->point.id] + waitDuration + s->visitDuration + mTtMatrix[s->point.id][curr->arrPointId] - mTtMatrix[temp->depPointId][curr->arrPointId];

				if (depTime <= s->timeWindow.closeTime && shift <= curr->maxShift) {	//check if insertion is possible
					if (shift < minShift) {
						minShift = shift;
						std::get<0>(res) = pos;
						std::get<1>(res) = minShift;
						std::get<2>(res) = s->point.id;
						std::get<3>(res) = s->point.id;
					}
				}

				pos++;
				curr = curr->next;
				temp = temp->next;
			}
		}
		else if ((r = dynamic_cast<Route*>(ta))) {
			while (curr != nullptr) {

				pos++;
				curr = curr->next;
				temp = temp->next;
			}
		}
	}
	return res;
}

bool OPTW::updateTimes(int startIndex, bool smart) {
	//update times first
	TA* curr = mProcessSolution.mWalk.get(startIndex);

	if (curr == mProcessSolution.mWalk.first()) {
		return false;
	}
	while (curr != nullptr) {

		curr->arrTime = curr->prev->depTime + mTtMatrix[curr->prev->depPointId][curr->arrPointId];
		curr->waitDuration = std::max(0.0, curr->timeWindow.openTime - curr->arrTime);
		curr->startOfVisitTime = curr->arrTime + curr->waitDuration;
		curr->depTime = curr->startOfVisitTime + curr->visitDuration;

		/*if(curr->data.waitDuration > 0){
			break;
		}*/
		curr = curr->next;

	}

	//update max shifts
	curr = mProcessSolution.mWalk.last();
	curr->maxShift = curr->timeWindow.closeTime - curr->depTime;
	curr = curr->prev;
	while (curr != nullptr) {
		if (curr->next != nullptr) {
			curr->maxShift = std::min(curr->timeWindow.closeTime - curr->depTime, curr->next->waitDuration + curr->next->maxShift);
		}
		else {
			curr->maxShift = curr->timeWindow.closeTime - curr->depTime;
		}
		curr = curr->prev;
	}

	return true;

}

bool TOPTW::updateTimes(int startIndex, bool smart) {
	//update times first
	TA* curr = mProcessSolution.mWalk.get(startIndex);

	if (curr == mProcessSolution.mWalk.first()) {
		return false;
	}
	while (curr != nullptr) {

		curr->arrTime = curr->prev->depTime + mTtMatrix[curr->prev->depPointId][curr->arrPointId];
		curr->waitDuration = std::max(0.0, curr->timeWindow.openTime - curr->arrTime);
		curr->startOfVisitTime = curr->arrTime + curr->waitDuration;
		curr->depTime = curr->startOfVisitTime + curr->visitDuration;

		/*if(curr->data.waitDuration > 0){
			break;
		}*/
		curr = curr->next;

	}

	//update max shifts
	curr = mProcessSolution.mWalk.last();
	curr->maxShift = curr->timeWindow.closeTime - curr->depTime;
	curr = curr->prev;
	while (curr != nullptr) {
		if (curr->next != nullptr) {
			curr->maxShift = std::min(curr->timeWindow.closeTime - curr->depTime, curr->next->waitDuration + curr->next->maxShift);
		}
		else {
			curr->maxShift = curr->timeWindow.closeTime - curr->depTime;
		}
		curr = curr->prev;
	}

	return true;

}
