#include "ILS.h"

Solution ILS::LocalSearch(Solution solution, double avgPoint) {

	TA* curr = solution.mUnvisited.first();
	TA* temp;

	size_t walk_length = solution.mWalk.getLength();
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

	return solution;

}

Solution ILS::Shake(Solution solution, int S, int R, int numOfPois, std::vector<std::vector<double>> ttMatrix) {
	if (numOfPois == 0) {
		return;
	}
	else if (numOfPois > 0) {
		List<TA> nodes = solution.mWalk.grabPart(S, std::min(S + R - 1, numOfPois + 1));
		solution.mUnvisited.append(nodes);
		solution = updateTimes(solution, 1, false, ttMatrix);
	}
	else {
		std::cerr << "Route has invalid length" << std::endl;
		exit(1);
	}
}

std::tuple<double, double, double> ILS::calcTimeEventCut(ListTA& unvisited) {
	double dist;
	double min = DBL_MAX;
	double timePoint = -1;
	TA* curr = unvisited.first();

	std::vector<double> timeWindowEvents;


	while (curr != nullptr) {
		timeWindowEvents.push_back(curr->timeWindow.openTime);
		timeWindowEvents.push_back(curr->timeWindow.closeTime);
		curr = curr->next;
	}

	sort(timeWindowEvents.begin(), timeWindowEvents.end());
	timeWindowEvents.erase(unique(timeWindowEvents.begin(), timeWindowEvents.end()), timeWindowEvents.end());

	double average = reduce(timeWindowEvents.begin(), timeWindowEvents.end()) / timeWindowEvents.size();

	curr = unvisited.first();
	while (curr != nullptr) {
		dist = abs(curr->timeWindow.openTime - average);
		if (dist < min) {
			min = dist;
			timePoint = curr->timeWindow.openTime;
		}
		dist = abs(curr->timeWindow.closeTime - average);
		if (dist < min) {
			min = dist;
			timePoint = curr->timeWindow.closeTime;
		}
		curr = curr->next;
	}

	return { timeWindowEvents.front(), timePoint, timeWindowEvents.back() };
}


ListTA ILS::setBucketActivityDurations(ListTA unvisited) {
	auto [minEvent, avgEvent, maxEvent] = calcTimeEventCut(unvisited);

	TA* curr = unvisited.first();
	while (curr != nullptr) {
		curr->bucketActivities[0].duration = avgEvent - curr->timeWindow.openTime;
		curr->bucketActivities[1].duration = curr->timeWindow.closeTime - avgEvent;
		curr = curr->next;
	}
	return unvisited;

}


Solution ILS::Solve(ListTA& unvisited, TA* start, TA* end, std::vector<std::vector<double>> ttMatrix) {

	Walk walk;
	walk.push(start);
	walk.push(end);

	Solution processSolution = Solution(walk, unvisited, start->timeWindow.openTime, end->timeWindow.closeTime);
	Solution bestSolution = Solution();

	int counter = 0;
	int S = 1, R = 1;
	int timesNotImproved = 0;
	TaskManager* taskManager = TaskManager::GetInstance();
	int bestScore = INT_MIN;

	processSolution = setBucketActivityDurations(processSolution.mUnvisited);

	auto [minEvent, avgEvent, maxEvent] = calcTimeEventCut(walk);

	while (timesNotImproved < MAX_TIMES_NOT_IMPROVED) {
		counter++;
		std::cout << "this is revision: " << counter << std::endl;

		processSolution = LocalSearch(processSolution, avgEvent);

		int score = walk.collectProfit();

		//int score = Insert();

		if (score > bestScore) {
			bestScore = score;
			bestSolution.reset();
			bestSolution.copy(processSolution);
			bestSolution.SetScore(score);
			R = 1;
			timesNotImproved = 0;

		}
		else {
			timesNotImproved++;
		}

		int numOfPois = processSolution.mWalk.getLength() - 2;
		if (S > numOfPois) {
			S = 1;
			// S = S - numOfPois;
			// if (S <= 0)
			// 	S = 1;
		}
		if (R > numOfPois / 2) {
			R = 1;
		}

		processSolution = Shake(processSolution, S, R, numOfPois, ttMatrix);

		S += R;
		R += 1;


	}

	return bestSolution;


	
}


int ILS::Insert(ListTA& walk, ListTA& unvisited, std::vector<std::vector<double>> ttMatrix) {
	double minShift;
	int pos, bestPos = DEFAULT_POS;
	double maxRatio = DBL_MIN, ratio;
	std::string nodeId;
	TA* nodeToInsert = nullptr;
	TA* curr;

	int bestArrPointId = DEFAULT_POINT_ID, bestDepPointId = DEFAULT_POINT_ID;
	int arrPointId = DEFAULT_POINT_ID, depPointId = DEFAULT_POINT_ID;

	while (true) {
		maxRatio = DBL_MIN;
		nodeId = DEFAULT_TA_ID;
		curr = unvisited.first();
		//numOfUnvisitedNodes = ProcessSolution.UnvisitedPOIs.GetNumOfNodes();

		while (curr != nullptr) {
			minShift = DBL_MAX;
			std::tie(pos, minShift, arrPointId, depPointId) = getBestPos(curr, walk, ttMatrix);

			if (pos == -1) { //insertion is not feasible
				curr = curr->next;
				continue;
			}

			ratio = pow(curr->profit, 2) / minShift;	//Ratioi = (Si)^2/Shifti

			if (maxRatio < ratio) //check each ratio
			{
				maxRatio = ratio;
				nodeToInsert = curr;
				bestPos = pos;
				bestArrPointId = arrPointId;
				bestDepPointId = depPointId;
			}

			curr = curr->next;
		}

		//if we have reached a point where there is no node which can be inserted, then we stop the search
		if (nodeToInsert == nullptr) {
			break;
		}

		if (nodeToInsert != nullptr) {
			walk.insertAt(nodeToInsert, bestPos, bestArrPointId, bestDepPointId, ttMatrix); //add it to route
			updateTimes(walk, bestPos, true, ttMatrix);
		}
		else {
			std::cerr << "Didn't find node " << nodeId << std::endl;
			std::exit(1);
		}
	}
}

std::tuple<int, double, int, int> ILS::getBestPos(TA* ta, ListTA& walk, std::vector<std::vector<double>> ttMatrix) {
	double minShift = DBL_MAX;
	int length = walk.getLength();
	double shift;
	int pos = 1;
	std::tuple<int, double, int, int> res = std::make_tuple(-1, -1, -1, -1);

	TA* temp = walk.first();

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
				shift = ttMatrix[temp->depPointId][s->point.id] + s->visitDuration + ttMatrix[s->point.id][curr->arrPointId] - ttMatrix[temp->depPointId][curr->arrPointId];

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
				int shift1 = ttMatrix[temp->depPointId][r->edges[0].id] + r->visitDuration + ttMatrix[r->edges[1].id][curr->arrPointId] - ttMatrix[temp->depPointId][curr->arrPointId];

				int shift2 = ttMatrix[temp->depPointId][r->edges[1].id] + r->visitDuration + ttMatrix[r->edges[0].id][curr->arrPointId] - ttMatrix[temp->depPointId][curr->arrPointId];

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


Solution ILS::updateTimes(Solution solution, int startIndex, bool smart, std::vector<std::vector<double>> ttMatrix) {
	//update times first
	TA* curr = solution.mWalk.get(startIndex);

	if (curr == solution.mWalk.first()) {
		std::cerr << "i should think what happens here" << std::endl;
		std::exit(1);
	}
	while (curr != nullptr) {
		curr->arrTime = curr->prev->depTime + ttMatrix[curr->prev->depPointId][curr->arrPointId];
		curr->startOfVisitTime = curr->arrTime;
		curr->depTime = curr->startOfVisitTime + curr->visitDuration;
		curr = curr->next;
	}

	//update max shifts
	curr = solution.mWalk.last();
	curr->maxShift = curr->timeWindow.closeTime - curr->depTime;
	curr = curr->prev;
	while (curr != nullptr) {
		if (curr->next != nullptr) {
			curr->maxShift = curr->next->maxShift;
		}
		curr = curr->prev;
	}

	return solution;
}

std::tuple<bool, std::string> ILS::validate(ListTA& walk, std::vector<std::vector<double>> ttMatrix) {
	std::cout << "Validating route:\t"; walk.print();
	std::string error = "";
	bool valid = true;
	TA* curr = walk.first();
	if (curr == nullptr) {
		return { false, "List is empty" };
	}
	while (curr != nullptr) {
		if (curr->startOfVisitTime != curr->arrTime) {
			error = "StartOfVisitTime(" + std::to_string(curr->startOfVisitTime) + ") != arrTime(" + std::to_string(curr->arrTime) + ")";
			valid = false;
			break;
		}

		if (curr->depTime != curr->startOfVisitTime + curr->visitDuration) {
			error = "depTime(" + std::to_string(curr->depTime) + ") != startOfVisitTime(" + std::to_string(curr->startOfVisitTime) + ") + visitDuration(" + std::to_string(curr->visitDuration) + ")";
			valid = false;
			break;
		}

		if (curr->next != nullptr) {
			if (curr->next->arrTime != curr->depTime + ttMatrix[curr->depPointId][curr->next->arrPointId]) {
				error = "next->arrTime(" + std::to_string(curr->next->arrTime) + ") != depTime(" + std::to_string(curr->depTime) + ") + timeTravel[" + std::to_string(curr->depPointId) + "][" + std::to_string(curr->next->arrPointId) + "](" + std::to_string(ttMatrix[curr->depPointId][curr->next->arrPointId]) + ")";
				valid = false;
				break;
			}
		}	
		curr = curr->next;
	}
	return {valid, error};
}

std::tuple<int, double, int, int> ILS_OPTW::getBestPos(TA* ta, ListTA& walk, std::vector<std::vector<double>> ttMatrix) {
	double minShift = DBL_MAX;
	int length = walk.getLength();
	double arrTime, waitDuration, startOfVisitTime, startOfVisitTime1, startOfVisitTime2, depTime, depTime1, depTime2, shift;
	int pos = 1;
	std::tuple<int, double, int, int> res = std::make_tuple(-1, -1, -1, -1);

	TA* temp = walk.first();
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
				arrTime = temp->depTime + ttMatrix[temp->depPointId][s->point.id];
				waitDuration = std::max(0.0, s->timeWindow.openTime - arrTime);
				startOfVisitTime = arrTime + waitDuration;
				depTime = startOfVisitTime + s->visitDuration;
				shift = ttMatrix[temp->depPointId][s->point.id] + waitDuration + s->visitDuration + ttMatrix[s->point.id][curr->arrPointId] - ttMatrix[temp->depPointId][curr->arrPointId];

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
				arrTime = temp->depTime + ttMatrix[temp->depPointId][r->edges[0].id];
				waitDuration = std::max(0.0, r->timeWindow.openTime - arrTime);
				startOfVisitTime1 = arrTime + waitDuration;
				depTime1 = startOfVisitTime1 + r->visitDuration;
				shift1 = ttMatrix[temp->depPointId][r->edges[0].id] + waitDuration + r->visitDuration + ttMatrix[r->edges[1].id][curr->arrPointId] - ttMatrix[temp->depPointId][curr->arrPointId];

				arrTime = temp->depTime + ttMatrix[temp->depPointId][r->edges[1].id];
				waitDuration = std::max(0.0, r->timeWindow.openTime - arrTime);
				startOfVisitTime2 = arrTime + waitDuration;
				depTime2 = startOfVisitTime2 + r->visitDuration;
				shift2 = ttMatrix[temp->depPointId][r->edges[1].id] + waitDuration + r->visitDuration + ttMatrix[r->edges[0].id][curr->arrPointId] - ttMatrix[temp->depPointId][curr->arrPointId];

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


Solution ILS_OPTW::updateTimes(Solution solution, int startIndex, bool smart, std::vector<std::vector<double>> ttMatrix) {
	//update times first
	TA* curr = solution.mWalk.get(startIndex);

	if (curr == solution.mWalk.first()) {
		return solution;
	}
	while (curr != nullptr) {

		curr->arrTime = curr->prev->depTime + ttMatrix[curr->prev->depPointId][curr->arrPointId];
		curr->waitDuration = std::max(0.0, curr->timeWindow.openTime - curr->arrTime);
		curr->startOfVisitTime = curr->arrTime + curr->waitDuration;
		curr->depTime = curr->startOfVisitTime + curr->visitDuration;

		/*if(curr->data.waitDuration > 0){
			break;
		}*/
		curr = curr->next;

	}

	//update max shifts
	curr = solution.mWalk.last();
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

	return solution;
}

std::tuple<bool, std::string> ILS_OPTW::validate(ListTA& walk, std::vector<std::vector<double>> ttMatrix) {
	std::cout << "Validating route:\t"; walk.print();
	std::string error = "";
	bool valid = true;
	TA* curr = walk.first();
	if (curr == nullptr) {
		return { false, "List is empty" };
	}
	while (curr != nullptr) {
		if (curr->startOfVisitTime != curr->arrTime) {
			error = "StartOfVisitTime(" + std::to_string(curr->startOfVisitTime) + ") != arrTime(" + std::to_string(curr->arrTime) + ")";
			valid = false;
			break;
		}

		if (curr->depTime != curr->startOfVisitTime + curr->visitDuration) {
			error = "depTime(" + std::to_string(curr->depTime) + ") != startOfVisitTime(" + std::to_string(curr->startOfVisitTime) + ") + visitDuration(" + std::to_string(curr->visitDuration) + ")";
			valid = false;
			break;
		}

		if (curr->depTime > curr->timeWindow.closeTime) {
			error = "depTime(" + std::to_string(curr->depTime) + ") > closeTime(" + std::to_string(curr->timeWindow.closeTime) + ")";
			valid = false;
			break;
		}

		if (curr->next != nullptr) {
			if (curr->next->arrTime != curr->depTime + ttMatrix[curr->depPointId][curr->next->arrPointId]) {
				error = "next->arrTime(" + std::to_string(curr->next->arrTime) + ") != depTime(" + std::to_string(curr->depTime) + ") + timeTravel[" + std::to_string(curr->depPointId) + "][" + std::to_string(curr->next->arrPointId) + "](" + std::to_string(ttMatrix[curr->depPointId][curr->next->arrPointId]) + ")";
				valid = false;
				break;
			}
		}
		curr = curr->next;
	}
	return { valid, error };
}

std::tuple<int, double, int, int> ILS_TOPTW::getBestPos(TA* ta, ListTA& walk, std::vector<std::vector<double>> ttMatrix) {
	double minShift = DBL_MAX;
	int length = walk.getLength();
	double arrTime, waitDuration, startOfVisitTime, depTime, shift;
	int pos = 1;
	std::tuple<int, double, int, int> res = std::make_tuple(-1, -1, -1, -1);

	TA* temp = walk.first();
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
				arrTime = temp->depTime + ttMatrix[temp->depPointId][s->point.id];
				waitDuration = std::max(0.0, s->timeWindow.openTime - arrTime);
				startOfVisitTime = arrTime + waitDuration;
				depTime = startOfVisitTime + s->visitDuration;
				shift = ttMatrix[temp->depPointId][s->point.id] + waitDuration + s->visitDuration + ttMatrix[s->point.id][curr->arrPointId] - ttMatrix[temp->depPointId][curr->arrPointId];

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


Solution ILS_TOPTW::updateTimes(Solution solution, int startIndex, bool smart, std::vector<std::vector<double>> ttMatrix) {
	//update times first
	TA* curr = solution.mWalk.get(startIndex);

	if (curr == solution.mWalk.first()) {
		return solution;
	}
	while (curr != nullptr) {

		curr->arrTime = curr->prev->depTime + ttMatrix[curr->prev->depPointId][curr->arrPointId];
		curr->waitDuration = std::max(0.0, curr->timeWindow.openTime - curr->arrTime);
		curr->startOfVisitTime = curr->arrTime + curr->waitDuration;
		curr->depTime = curr->startOfVisitTime + curr->visitDuration;

		/*if(curr->data.waitDuration > 0){
			break;
		}*/
		curr = curr->next;

	}

	//update max shifts
	curr = solution.mWalk.last();
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

	return solution;
}

std::tuple<bool, std::string> ILS_TOPTW::validate(ListTA& walk, std::vector<std::vector<double>> ttMatrix) {
	std::cout << "Validating route:\t"; walk.print();
	std::string error = "";
	bool valid = true;
	TA* curr = walk.first();
	if (curr == nullptr) {
		return { false, "List is empty" };
	}
	while (curr != nullptr) {
		if (curr->startOfVisitTime != curr->arrTime) {
			error = "StartOfVisitTime(" + std::to_string(curr->startOfVisitTime) + ") != arrTime(" + std::to_string(curr->arrTime) + ")";
			valid = false;
			break;
		}

		if (curr->depTime != curr->startOfVisitTime + curr->visitDuration) {
			error = "depTime(" + std::to_string(curr->depTime) + ") != startOfVisitTime(" + std::to_string(curr->startOfVisitTime) + ") + visitDuration(" + std::to_string(curr->visitDuration) + ")";
			valid = false;
			break;
		}

		if (curr->depTime > curr->timeWindow.closeTime) {
			error = "depTime(" + std::to_string(curr->depTime) + ") > closeTime(" + std::to_string(curr->timeWindow.closeTime) + ")";
			valid = false;
			break;
		}

		if (curr->next != nullptr) {
			if (curr->next->arrTime != curr->depTime + ttMatrix[curr->depPointId][curr->next->arrPointId]) {
				error = "next->arrTime(" + std::to_string(curr->next->arrTime) + ") != depTime(" + std::to_string(curr->depTime) + ") + timeTravel[" + std::to_string(curr->depPointId) + "][" + std::to_string(curr->next->arrPointId) + "](" + std::to_string(ttMatrix[curr->depPointId][curr->next->arrPointId]) + ")";
				valid = false;
				break;
			}
		}
		curr = curr->next;
	}
	return { valid, error };
}