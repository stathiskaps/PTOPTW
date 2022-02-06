#include "ILS.h"

ILS::ILS() {}

ILS::ILS(int bucketsNum) : mBucketsNum(bucketsNum) {}

ILS::~ILS() {}

//appends list2 after list1 (list1 + list2)
ListTA operator+(ListTA & l1, ListTA & l2) {
	ListTA l3;

	TA* curr = l1.first();
	while (curr != nullptr) {
		l3.pushClone(curr);
		curr = curr->next;
	}

	curr = l2.first();
	while (curr != nullptr) {
		l3.pushClone(curr);
		curr = curr->next;
	}

	return l3;
}

struct {

};

bool areNeighbors(TA* ta1, TA* ta2, double maxtime, std::vector<std::vector<double>> ttMatrix) {
	bool neighbors = true;
	//two nodes are neighbors if satisfy the following requirements
	//1st requirement: travel time between them is less than n
	if (ttMatrix[ta1->point.id][ta2->point.id] > maxtime) {
		neighbors = false;
	}

	//2nd requirement: they have at least m shared minutes in their time windows 
	//TimeWindow common_time_window = TimeWindow{ std::max(p->timeWindow.openTime, q->timeWindow.openTime), std::min(p->timeWindow.closeTime, q->timeWindow.closeTime) };
	//if (common_time_window.length() < m) {
	//	neighbors = false;
	//}
	
	// [3,9] [5,7] [6, 8]
	// [2,4] [7,9]
	// [2,6] [4,8] [4, 6]

	return neighbors;
}

std::tuple<bool, TimeWindow> commonTimeWindow(std::vector<TA*> nodes) {
	double maxOpenTime = DBL_MIN;
	double minCloseTime = DBL_MAX;
	for (auto& n : nodes) {
		if (n->timeWindow.openTime > maxOpenTime) {
			maxOpenTime = n->timeWindow.openTime;
		}
		if (n->timeWindow.closeTime < minCloseTime) {
			minCloseTime = n->timeWindow.closeTime;
		}
	}
	if (maxOpenTime == DBL_MIN || minCloseTime == DBL_MAX || minCloseTime-maxOpenTime <= 0) {
		return { false, TimeWindow{0, 0} };
	}
	return { true, TimeWindow{ maxOpenTime, minCloseTime } };
}

std::vector<TA*> rangeQuery(TA* Q, ListTA nodes, double radius, double minCommonTimeWindow, std::vector<std::vector<double>> ttMatrix) {
	TA* curr = nodes.first();
	std::vector<TA*> neighbors;

	while (curr != nullptr) {
		if (areNeighbors(curr, Q, radius, ttMatrix)) {
			neighbors.push_back(curr);
		}
		curr = curr->next;
	}

	//auto [common, window] = commonTimeWindow(neighbors);
	//if (!common) {
	//	neighbors.clear();
	//}
	return neighbors;
}

bool compareById(const TA* a, const TA* b)
{
	return a->id < b->id;
}

bool sameId(const TA* a, const TA* b)
{
	return a->id == b->id;
}

void dbScan(ListTA list, std::vector<std::vector<double>> ttMatrix) {
	
	TA* curr = list.first();
	const int radius = 10;
	const int common_minutes = 200;
	const int minPoints = 3;
	
	int C = 0;
	while (curr != nullptr) {
		if (curr->cluster != UNDEFINED) {
			curr = curr->next;
			continue;
		}
		std::vector<TA*> neighbors = rangeQuery(curr, list, radius, common_minutes, ttMatrix);
		if (neighbors.size() < minPoints) {
			curr->cluster = NOISE;
			curr = curr->next;
			continue;
		}
		C++;
		curr->cluster = C;
		size_t old_size = neighbors.size();
		size_t erased = std::erase_if(neighbors, [curr](TA* t) {return t->id == curr->id; });
		if (erased != 1) {
			throw std::out_of_range("node wasn't found to erase");
		}

		size_t neighborsSize = neighbors.size();
		for (size_t i = 0; i < neighborsSize; ++i) {
			if (neighbors[i]->cluster == NOISE) {
				neighbors[i]->cluster = LEAF;
			}

			if (neighbors[i]->cluster != UNDEFINED) {
				continue;
			}

			neighbors[i]->cluster = C;
			std::vector<TA*> nb = rangeQuery(neighbors[i], list, radius, common_minutes, ttMatrix);
			if (nb.size() >= minPoints) {
				neighbors.insert(neighbors.end(), nb.begin(), nb.end());
				/*std::sort(neighbors.begin(), neighbors.end(), compareById);
				neighbors.erase(unique(neighbors.begin(), neighbors.end(), sameId), neighbors.end());*/
				neighborsSize += nb.size();
			}
		}


		curr = curr->next;
	}
	
}

void ILS::LocalSearch(std::vector<Solution>& solutions, std::vector<double> cuts, OP& op) {

	size_t solutionsSize = solutions.size();

	for (int i = 0; i < solutionsSize; ++i) {

		Walk walk;
		double startTime = 0, endTime = cuts[i];

		int currentWalkLength = solutions[i].mWalk.getLength();
		TA *startDepot, *endDepot;
		if (currentWalkLength == 0) {
			//TODO: check if there is a case where previous walk is empty
			if (i == 0) {
				startDepot = op.mStartDepot->clone();
			}
			else {
				startDepot = solutions[i - 1].mWalk.last()->clone();
			}
			solutions[i].mWalk.pushFront(startDepot);
		}
		else {
			if (i > 0) {
				startDepot = solutions[i - 1].mWalk.last()->clone();
				solutions[i].mWalk.pushFront(startDepot);
			}
			else {
				if (solutions[i].mWalk.first()->id != op.mStartDepot->id) {
					solutions[i].mWalk.pushFront(op.mStartDepot->clone());
				}
			}
			//if i == 0 keep the same startNode
		}

		if (i == solutionsSize - 1) {
			endDepot = op.mEndDepot->clone();
			endDepot->maxShift = endDepot->timeWindow.closeTime - endDepot->depTime;
			if (solutions[i].mWalk.last()->id != endDepot->id) {
				solutions[i].mWalk.pushBack(endDepot);
			}
			
		}
		else {
			Point cnext;
			int nextWalkLength = solutions[i + 1].mWalk.getLength();
			if (nextWalkLength == 0) {
				cnext = solutions[i + 1].mUnvisited.getWeightedCentroid();
			}
			else {
				int n = std::min(3, nextWalkLength);
				cnext = solutions[i + 1].mUnvisited.copyPart(0,n-1).getWeightedCentroid();
			}
			op.AddPointToGraph(cnext);
			//op.PrintTravelTimes("New travel times");
			endDepot = new TA(cnext); //todo: delete endDepot
			endDepot->timeWindow.closeTime = cuts[i];
			//endDepot->maxShift = endDepot->timeWindow.closeTime - endDepot->depTime;
			solutions[i].mWalk.pushBack(endDepot);
		}

		updateTimes(solutions[i], 0, false, op.mTravelTimes);


		construct(solutions[i], op.mTravelTimes);
		if (i > 0) {
			solutions[i].mWalk.trimLeft(1);
		}

		if (i < solutionsSize - 1) {
			solutions[i].mWalk.trimRight(1);
		}
		

		


		//TA* curr = solutions[i].mUnvisited.first();
		//TA* temp;
		//TA* depot = solutions[i].mWalk.first();

		//ListTA unvisitedA, unvisitedB;
		//int unvisitedACount = 0, unvisitedBCount = 0;
		//while (curr != nullptr) {
		//	temp = curr;
		//	curr = curr->next;
		//	if (temp->bucketActivities[0].duration >= temp->bucketActivities[1].duration) {
		//		unvisitedA.pushClone(temp);
		//		unvisitedACount++;
		//	}
		//	else {
		//		unvisitedB.pushClone(temp);
		//		unvisitedBCount;
		//	}
		//}
		//size_t walk_length = solutions[i].mWalk.getLength();
		//if (walk_length == 2) {

		//	Point cb = unvisitedB.getWeightedCentroid();
		//	cb.id = op.mTravelTimes.size();
		//	op.AddPointToGraph(cb);

		//	//construct depends completely on first node's departure time and on last node's close time
		//	TA* temp = new TA(cb);
		//	Solution solA = Solution(depot, temp, unvisitedA, depot->timeWindow.openTime, cuts[i]);
		//	delete(temp);
		//	solA = construct(solA, op.mTravelTimes);

		//	//get last two elements
		//	Walk walkB = solA.mWalk.copyPart(-2, -1);
		//	Solution solB = Solution(walkB.prelast(), depot, unvisitedB, cuts[i], solutions[i].mWalk.last()->timeWindow.closeTime);
		//	solB = construct(solB, op.mTravelTimes);


		//	solA.mWalk.trimRight(1);
		//	solB.mWalk.trimLeft(1);
		//	//walkB.removeLast();

		//	solutions[i].mUnvisited = solA.mUnvisited + solB.mUnvisited;
		//	solutions[i].mWalk = solA.mWalk + solB.mWalk;

		//}
		//else if (walk_length > 2) {
		//	curr = solutions[i].mWalk.first();
		//	auto [wa, wb] = solutions[i].mWalk.splitOnDepTime(cuts[i]);

		//	size_t wa_length = wa.getLength(), wb_length = wb.getLength();

		//	Solution solA, solB;

		//	if (wa_length <= 1 || unvisitedACount == 0) {
		//		solA = Solution(wa, unvisitedA);
		//	}
		//	else {
		//		wa.last()->timeWindow.closeTime = cuts[i];
		//		wa.last()->maxShift = cuts[i] - wa.last()->depTime;
		//		//we changed the maxshift of the last node, we need to change maxShift of previous nodes as well
		//		updateMaxShifts(wa, op.mTravelTimes);

		//		solA = Solution(wa, unvisitedA);
		//		solA = construct(solA, op.mTravelTimes);
		//	}

		//	if (wb_length <= 1 || unvisitedBCount == 0) {
		//		solB = Solution(wb, unvisitedB);
		//	}
		//	else {
		//		solB = Solution(wb, unvisitedB);
		//		solB = construct(solB, op.mTravelTimes);
		//	}

		//	//walkB.removeFirst();

		//	solutions[i].mWalk = solA.mWalk + solB.mWalk;
		//	solutions[i].mUnvisited = solA.mUnvisited + solB.mUnvisited;

		//	//because we may added a TA to walkA, some departure times may moved forward, which affects time of walkB
		//	solutions[i] = updateTimes(solutions[i], 1, true, op.mTravelTimes);
		//}
		//else {
		//	std::cerr << "Invalid walk length" << std::endl;
		//	std::exit(1);
		//}

		//auto [valid, error] = validate(solutions[i].mWalk, op.mTravelTimes);
		//if (!valid) {
		//	std::cout << error << std::endl;
		//	std::exit(1);
		//}
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


ListTA ILS::setBucketActivityDurations(ListTA& unvisited, double avgEvent) {

	TA* curr = unvisited.first();
	while (curr != nullptr) {
		curr->bucketActivities[0].duration = avgEvent - curr->timeWindow.openTime;
		curr->bucketActivities[1].duration = curr->timeWindow.closeTime - avgEvent;
		curr = curr->next;
	}
	return unvisited;

}

bool compareTimeWindowCenter(const TA* left, const TA* right) {
	return (left->timeWindow.openTime + left->timeWindow.closeTime) / 2 < (right->timeWindow.openTime + right->timeWindow.closeTime) / 2;
}

std::vector<std::vector<TA*>> ILS::getBuckets(std::vector<TA*> nodes, int m) {
	
	std::vector<std::vector<TA*>> buckets;
	int offset = 0, remainder = nodes.size() % m, bucketSize = nodes.size() / m;

	std::sort(nodes.begin(), nodes.end(), &compareTimeWindowCenter);

	for (int i = 0; i < m; ++i) {
		std::vector<TA*> slice(nodes.begin() + offset, nodes.begin() + offset + bucketSize + remainder);
		offset += bucketSize + remainder;
		if (remainder > 0) {
			remainder--;
		}
		buckets.push_back(slice);
	}

	return buckets;
}

std::vector<double> ILS::getTimeCuts(std::vector<std::vector<TA*>> buckets) {
	std::vector<double> cuts;

	size_t bucketsSize = buckets.size();
	for (int i = 0; i < bucketsSize - 1; ++i) {
		double average = (buckets[i].back()->timeWindow.openTime + 
			buckets[i].back()->timeWindow.closeTime + 
			buckets[i + 1].front()->timeWindow.openTime + 
			buckets[i + 1].front()->timeWindow.closeTime) / 4;
		cuts.push_back(average);
	}
	return cuts;
}

int ILS::collectScore(std::vector<Solution> solutions){
	int totalScore = 0;
	for (auto& sol : solutions) {
		totalScore += sol.getScore();
	}
	return totalScore;
}

std::tuple<int, int> ILS::getMinMaxLength(std::vector<Solution> solutions) {
	int min = INT_MAX, max = INT_MIN;

	for (auto& sol : solutions) {
		int walkLength = sol.mWalk.getLength();
		if (walkLength > max) {
			max = walkLength;
		}

		if (walkLength < min) {
			min = walkLength;
		}
	}

	return { min, max };
}

Solution ILS::connectSolutions(std::vector<Solution> solutions) {
	Walk walk;
	ListTA unvisited;
	for (auto& s : solutions) {
		walk.append(s.mWalk.copy());
		unvisited.append(s.mUnvisited.copy());
	}
	return Solution(walk, unvisited);
}

void ILS::Shake(std::vector<Solution>& solutions, std::vector<ShakeParameters> params, OP& op) {
	size_t solutionsSize = solutions.size();
	for (int i = 0; i < solutionsSize; ++i) {
		int walkLength = solutions[i].mWalk.getLength();

		//todo: check what happens when s > length - 1
		ListTA nodes = solutions[i].mWalk.grabPart(params[i].S, std::min(params[i].S + params[i].R - 1, walkLength - 1));
		nodes.removeByIds({ op.mStartDepot->id, op.mEndDepot->id });
		solutions[i].mUnvisited.append(nodes);
		updateTimes(solutions[i], 1, false, op.mTravelTimes);
	}
}

Solution ILS::Solve(OP& op) {

	
	auto started = std::chrono::high_resolution_clock::now();

	ListTA unvisited = ListTA(op.mAttractions);
	
	//initialize num different buckets, in which we will run the problems. We want to swap nodes inside these buckets
	auto bins = getBuckets(op.mAttractions, mBucketsNum);
	auto cuts = getTimeCuts(bins);
	cuts.push_back(op.mEndDepot->timeWindow.closeTime);

	//intialize solutions
	std::vector<Solution> processSolutions, best;
	std::vector<ShakeParameters> parameters;
	for (auto& bin : bins) {
		processSolutions.push_back(Solution(bin));
		parameters.push_back(ShakeParameters());
	}

	int solutionsSize = processSolutions.size();



	Walk walk;
	walk.pushClone(op.mStartDepot);
	walk.pushClone(op.mEndDepot);

	Solution processSolution = Solution(walk, unvisited, op.mStartDepot->timeWindow.openTime, op.mEndDepot->timeWindow.closeTime);
	Solution bestSolution = Solution();

	int counter = 0;
	int S = 1, R = 1;
	int timesNotImproved = 0;
	int bestScore = INT_MIN;

	auto [minEvent, avgEvent, maxEvent] = calcTimeEventCut(processSolution.mUnvisited);
	processSolution.mUnvisited = setBucketActivityDurations(processSolution.mUnvisited, avgEvent);

	dbScan(processSolution.mUnvisited, op.mTravelTimes);

	//processSolution.mUnvisited.foreach([](TA* ta) {std::cout << "id: " << ta->id << ", cluster: " << ta->cluster << std::endl; });

	while (timesNotImproved < MAX_TIMES_NOT_IMPROVED) {
		counter++;
		LocalSearch(processSolutions, cuts, op);
		int score = collectScore(processSolutions);

		if (score > bestScore) {
			bestScore = score;
			bestSolution.reset();
			bestSolution = connectSolutions(processSolutions);
			for (auto& params : parameters) {
				params.R = 1;
			}
			timesNotImproved = 0;

		}
		else {
			timesNotImproved++;
		}

		auto [min, max] = getMinMaxLength(processSolutions);

		for (int i = 0; i < solutionsSize; ++i) {
			int walkLength = processSolutions[i].mWalk.getLength();
			if (parameters[i].S > walkLength - 1) {
				parameters[i].S = 0;
			}
			if (parameters[i].R > walkLength / 2) {
				parameters[i].R = 1;
			}
		}

		Shake(processSolutions, parameters, op);

		for (auto& params : parameters) {
			params.S += 1;
			params.R += 1;
		}


	}
	auto done = std::chrono::high_resolution_clock::now();
	int finalScore = bestSolution.getScore();
	std::cout.clear();
	std::cout << "final score: " << finalScore << std::endl;
	std::cout << "time taken: " << std::chrono::duration_cast<std::chrono::milliseconds>(done - started).count() << "ms" << std::endl;

	return bestSolution;


	
}


//construct depends completely on first node's departure time and on last node's close time
void ILS::construct(Solution& sol, std::vector<std::vector<double>>& ttMatrix) {
	double minShift;
	int pos, bestPos = DEFAULT_POS;
	double maxRatio = DBL_MIN, ratio;
	TA* nodeToInsert = nullptr;
	TA* curr;

	int bestArrPointId = DEFAULT_POINT_ID, bestDepPointId = DEFAULT_POINT_ID;
	int arrPointId = DEFAULT_POINT_ID, depPointId = DEFAULT_POINT_ID;

	while (true) {
		maxRatio = DBL_MIN;
		nodeToInsert = nullptr;
		curr = sol.mUnvisited.first();
		//numOfUnvisitedNodes = ProcessSolution.UnvisitedPOIs.GetNumOfNodes();

		while (curr != nullptr) {
			minShift = DBL_MAX;
			std::tie(pos, minShift, arrPointId, depPointId) = getBestPos(curr, sol.mWalk, ttMatrix);

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

		nodeToInsert = sol.mUnvisited.grab(nodeToInsert);
		sol.mWalk.insertAt(nodeToInsert, bestPos, bestArrPointId, bestDepPointId, ttMatrix); //add it to route
		updateTimes(sol, bestPos, true, ttMatrix);
	}
}

//void ILS::correct(Solution& solution, double startLimit, double endLimit) {
//	
//}



std::tuple<int, double, int, int> ILS::getBestPos(TA* ta, ListTA walk, std::vector<std::vector<double>> ttMatrix) {
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


void ILS::updateTimes(Solution& solution, int startIndex, bool smart, std::vector<std::vector<double>>& ttMatrix) {
	//update times first
	TA* curr = solution.mWalk.get(startIndex);
	if (curr == nullptr) {
		return;
	}

	if (startIndex == 0) {
		curr->startOfVisitTime = curr->arrTime + curr->waitDuration;
		curr->depTime = curr->startOfVisitTime + curr->visitDuration;
		curr = curr->next;
	}

	while (curr != nullptr) {
		curr->arrTime = curr->prev->depTime + ttMatrix[curr->prev->depPointId][curr->arrPointId];
		curr->startOfVisitTime = curr->arrTime;
		curr->depTime = curr->startOfVisitTime + curr->visitDuration;
		curr = curr->next;
	}

	//update max shifts
	updateMaxShifts(solution.mWalk, ttMatrix);

}

void ILS::updateMaxShifts(Walk& walk, std::vector<std::vector<double>>& ttMatrix) {
	TA* curr = walk.last();
	if (curr == nullptr) return;
	curr->maxShift = curr->timeWindow.closeTime - curr->depTime;
	curr = curr->prev;
	while (curr != nullptr) {
		curr->maxShift = curr->next->maxShift;
		curr = curr->prev;
	}

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
		if (curr->startOfVisitTime != curr->arrTime + curr->waitDuration) {
			error = "StartOfVisitTime(" + std::to_string(curr->startOfVisitTime) + ") != arrTime(" + std::to_string(curr->arrTime) + ")" + " + waitDuration(" + std::to_string(curr->waitDuration) + ")";
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

std::tuple<int, double, int, int> ILS_OPTW::getBestPos(TA* ta, ListTA walk, std::vector<std::vector<double>> ttMatrix) {
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


void ILS_OPTW::updateTimes(Solution& solution, int startIndex, bool smart, std::vector<std::vector<double>>& ttMatrix) {
	//update times first
	TA* curr = solution.mWalk.get(startIndex);
	if (curr == nullptr) {
		return;
	}

	if (startIndex == 0) {
		curr->waitDuration = std::max(0.0, curr->timeWindow.openTime - curr->arrTime);
		curr->startOfVisitTime = curr->arrTime + curr->waitDuration;
		curr->depTime = curr->startOfVisitTime + curr->visitDuration;
		curr = curr->next;
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
	updateMaxShifts(solution.mWalk, ttMatrix);

}

void ILS_OPTW::updateMaxShifts(Walk& walk, std::vector<std::vector<double>>& ttMatrix) {
	TA* curr = walk.last();
	if (curr == nullptr) return;
	curr->maxShift = curr->timeWindow.closeTime - curr->depTime;
	curr = curr->prev;
	while (curr != nullptr) {
		curr->maxShift = std::min(curr->timeWindow.closeTime - curr->depTime, curr->next->waitDuration + curr->next->maxShift);
		curr = curr->prev;
	}
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
		if (curr->startOfVisitTime != curr->arrTime + curr->waitDuration) {
			error = "StartOfVisitTime(" + std::to_string(curr->startOfVisitTime) + ") != arrTime(" + std::to_string(curr->arrTime) + ")" + " + waitDuration(" + std::to_string(curr->waitDuration) + ")";
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

std::tuple<int, double, int, int> ILS_TOPTW::getBestPos(TA* ta, ListTA walk, std::vector<std::vector<double>> ttMatrix) {
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


void ILS_TOPTW::updateTimes(Solution& solution, int startIndex, bool smart, std::vector<std::vector<double>>& ttMatrix) {
	//update times first
	TA* curr = solution.mWalk.get(startIndex);
	if (curr == nullptr) {
		return;
	}

	if (startIndex == 0) {
		curr->waitDuration = std::max(0.0, curr->timeWindow.openTime - curr->arrTime);
		curr->startOfVisitTime = curr->arrTime + curr->waitDuration;
		curr->depTime = curr->startOfVisitTime + curr->visitDuration;
		curr = curr->next;
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

	updateMaxShifts(solution.mWalk, ttMatrix);
	
}

void ILS_TOPTW::updateMaxShifts(Walk& walk, std::vector<std::vector<double>>& ttMatrix) {
	TA* curr = walk.last();
	if (curr == nullptr) return;
	curr->maxShift = curr->timeWindow.closeTime - curr->depTime;
	curr = curr->prev;
	while (curr != nullptr) {
		curr->maxShift = std::min(curr->timeWindow.closeTime - curr->depTime, curr->next->waitDuration + curr->next->maxShift);
		curr = curr->prev;
	}

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
		if (curr->startOfVisitTime != curr->arrTime + curr->waitDuration) {
			error = "StartOfVisitTime(" + std::to_string(curr->startOfVisitTime) + ") != arrTime(" + std::to_string(curr->arrTime) + ")" + " + waitDuration(" + std::to_string(curr->waitDuration) + ")";
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

