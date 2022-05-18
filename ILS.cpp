#include "ILS.h"

ILS::ILS() : mBucketsNum(0) {}

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

	if (neighbors.size() >= 3) {
		Divider divider;
		std::vector<TA*> group;
		group = divider.findBestAttractionCombination(neighbors);
		if (group.size() > 1) {
			for (auto& ta : group) {
				std::cout << ta->id << std::endl;
			}
		}
	}

	//Step 2: filter neighbors based on their common time window

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

//In each revision we want to do two things:
//1) swap unvisited nodes between buckets
//2) find small clusters, construct small solutions, use them as route

Solution ILS::Preprocess(OP& op) {

	Walk walk;
	walk.pushClone(op.mStartDepot);
	walk.pushClone(op.mEndDepot);
	ListTA unvisited = ListTA(op.mAttractions);
	Solution processSolution = Solution(walk, unvisited, op.mStartDepot->timeWindow.openTime, op.mEndDepot->timeWindow.closeTime);
	Solution bestSolution = Solution();
	int times_not_improved = 0;

	{
		auto started = std::chrono::high_resolution_clock::now();
		dbScan(processSolution.mUnvisited, op.mTravelTimes);
		auto done = std::chrono::high_resolution_clock::now();
		std::cout << "dbscan time taken: " << std::chrono::duration_cast<std::chrono::milliseconds>(done - started).count() << "ms" << std::endl;
	}
	

	processSolution.mUnvisited.foreach([](TA* ta) {std::cout << "id: " << ta->id << ", cluster: " << ta->cluster << std::endl; });

	int maxCluster = INT_MIN;
	processSolution.mUnvisited.foreach([&maxCluster](TA* curr) {
		if (curr->cluster > maxCluster) { maxCluster = curr->cluster; }
	});

	std::vector<ListTA> groups;
	Divider divider;
	for (int i = 1; i <= maxCluster; ++i) {
		ListTA group = processSolution.mUnvisited.map([i](TA* ta) { return ta->cluster == i; });
		if (group.size() > 2) {
			std::vector<TA*> combination = divider.findBestAttractionCombination(group.toVecPtr());
			Solution sol = Solution(Walk(), combination);
			ILS ils = ILS();
			//ils.construct(sol, op.mTravelTimes);

			Solution solution = Solution(Walk(), group);
		}
		
	}

	//while (times_not_improved < MAX_TIMES_NOT_IMPROVED){
	//	//Step 1: choose small clusters
	//}
	


	return bestSolution;
}

void ILS::SolveNew(OP& op) {
	std::cout.clear();

	//todo: for pr13 and buckets_num 4 got invalid bucket error
	const int num_locations = op.mAttractions.size();
	const int max_to_remove = num_locations / (3 * op.m_walks_num);
	const int buckets_num = 4;

	CustomList<TA> unvisited;
	for (auto ta : op.mAttractions) {
		unvisited.push_back(*ta);
	}

	std::vector<double> cuts = Preprocessing(op.mAttractions, buckets_num, op.mEndDepot->timeWindow.closeTime);
	std::map<std::string, Activity> registry = initializeRegistry(unvisited, cuts);
	CustomSolution process_solution{ *op.mStartDepot, *op.mEndDepot, unvisited, op.mStartDepot->timeWindow.openTime, op.mEndDepot->timeWindow.closeTime, op.m_walks_num }, best_solution{};

	int counter{}, S{ 1 }, R{ 1 }, times_not_improved{ 0 }, best_score{ INT_MIN };

	while (times_not_improved < MAX_TIMES_NOT_IMPROVED) {
		counter++;

		std::vector<CustomSolution> process_solutions = splitSolution(process_solution, cuts, registry);
		SplitSearch(process_solutions, cuts, op, registry);
		process_solution = connectSolutions(process_solutions, op.m_walks_num);
		int score = process_solution.getScores();
		//construct(process_solution, op.mTravelTimes);
		//int score = process_solution.getScores();
		if (score > best_score) {
			best_score = score;
			//best_solution.reset();
			best_solution = process_solution;
			R = 1;
			times_not_improved = 0;
		}
		else {
			times_not_improved++;
		}

		Shake(process_solution, S, R, op, max_to_remove);

	}
	validate(best_solution.m_walks, op.mTravelTimes);
	std::cout << "Best score: " << best_score << std::endl;
	std::cout << "Visits: " << best_solution.getVisits() << std::endl;
	std::cout << std::endl;
}

std::vector<double> ILS::Preprocessing(std::vector<TA*> unvisited, int bins_num, double day_close_time) {
	std::vector<double> cuts;
	cuts.reserve(bins_num + 1);
	std::vector<Bin> bins;
	bins.reserve(bins_num);
	uint64_t remainder{ unvisited.size() % bins_num }, bucketSize{ unvisited.size() / bins_num };

	std::sort(unvisited.begin(), unvisited.end(), [](const TA* left,const TA* right) -> bool {
		return left->timeWindow.openTime + left->timeWindow.closeTime < right->timeWindow.openTime + right->timeWindow.closeTime;
	});

	for (int i = 0, offset = 0; i < bins_num; ++i) {
		if (remainder > 0) {
			std::vector<TA*> slice(unvisited.begin() + offset, unvisited.begin() + offset + bucketSize + 1);
			Bin b;
			b.unvisited.reserve(slice.size());
			for (auto& p : slice) b.unvisited.push_back(*p);
			bins.push_back(b);
			offset += bucketSize + 1;
			remainder--;
		}
		else {
			std::vector<TA*> slice(unvisited.begin() + offset, unvisited.begin() + offset + bucketSize);
			Bin b;
			b.unvisited.reserve(slice.size());
			for (auto& p : slice) b.unvisited.push_back(*p);
			offset += bucketSize;
			bins.push_back(b);
		}
	}

	cuts.push_back(0);
	for (std::vector<Bin>::iterator it = bins.begin() + 1; it != bins.end(); ++it) {
		const TA &ta = it->unvisited.front();
		cuts.push_back((ta.timeWindow.openTime + ta.timeWindow.closeTime) / 2 - 1);
	}
	cuts.push_back(day_close_time);
	

	for (std::vector<TA*>::iterator it = unvisited.begin(); it != unvisited.end(); ++it) {
		for (std::vector<double>::iterator left = cuts.begin(), right = cuts.begin() + 1; right != cuts.end(); ++left, ++right) {
			double duration = std::min((*it)->timeWindow.closeTime, * right) - std::max((*it)->timeWindow.openTime, *left);
			if (duration < 0) duration = 0;
			ActivityInBucket activity{ .inBucket = 1, .inWalk = 1, .duration = duration };
			(*it)->stats.buckets.push_back(activity);
		}
	}

	return cuts;
}

std::vector<CustomSolution> ILS::splitSolution(CustomSolution& sol, const std::vector<double>& cuts, std::map<std::string, Activity>& registry) {
	std::vector<CustomSolution> solutions(cuts.size() - 1, CustomSolution());
	
	std::cout << "Splitting unvisited" << std::endl;
	//Split Unvisited
	for (auto& ta : sol.m_unvisited) {
		const double total_duration{ ta.timeWindow.closeTime - ta.timeWindow.openTime };
		double max_score{ -1 };
		std::vector<Bucket>::iterator best_b{ registry[ta.id].buckets.end()};
		for (std::vector<Bucket>::iterator b = registry[ta.id].buckets.begin(); b != registry[ta.id].buckets.end(); ++b) {
			double score = (b->duration / total_duration) * (b->inSolution / b->inBucket);
			if (score > max_score) {
				max_score = score;
				best_b = b;
			}
		}
		if (best_b == registry[ta.id].buckets.end()) {
			std::cerr << "invalid bucket" << std::endl;
			std::exit(1);
		}

		const int index = best_b - registry[ta.id].buckets.begin();
		solutions[index].m_unvisited.push_back(ta);
		best_b->inBucket++;
	}

	sol.m_unvisited.clear();

	//Split Walks
	for (auto& s : solutions) {
		for (size_t i = 0; i < sol.m_walks.size(); ++i) {
			s.m_walks.push_back(CustomList<TA>());
		}
	}
	
	for (Walks::iterator walk_it = sol.m_walks.begin(); walk_it != sol.m_walks.end(); ++walk_it) {
		if (walk_it->size() == 2) continue;
		for (CustomList<TA>::iterator ta_it = walk_it->begin(); ta_it != walk_it->end(); ++ta_it) {
			for (std::vector<double>::const_iterator left = cuts.begin(), right = cuts.begin() + 1; right != cuts.end(); ++left, ++right) {
				if (ta_it.iter->data.depTime >= *left && ta_it.iter->data.depTime < *right) {
					int64_t sol_index{ left - cuts.begin() }, walk_index{walk_it - sol.m_walks.begin()};
					solutions[sol_index].m_walks[walk_index].push_back(ta_it.iter->data);
					break;
				}
			}
		}
		walk_it->clear();
	}
	
	return solutions;
}

CustomSolution ILS::connectSolutions(std::vector<CustomSolution>& sols, const size_t walks_num) {
	CustomSolution solution;
	//Connect unvisited
	for (auto& s : sols) {
		solution.m_unvisited.append(s.m_unvisited);
		s.m_unvisited.clear();
	}
	for (size_t i = 0; i < walks_num; ++i) {
		CustomList<TA> walk;
		for (auto& s : sols) {
			walk.append(s.m_walks[i]);
			s.m_walks[i].clear();
		}
		solution.m_walks.push_back(walk);
	}

	//registry.clear();

	return solution;
}

int ILS::collectScores(std::vector<CustomSolution> sols) {
	int total_score{};
	for (auto& sol : sols) {
		total_score += sol.getScores();
	}
	return total_score;
}

std::map<std::string, Activity> ILS::initializeRegistry(const CustomList<TA>& unvisited, const std::vector<double>& cuts) {
	std::map<std::string, Activity> reg;
	for (CustomList<TA>::iterator ta_it = unvisited.begin(); ta_it != unvisited.end(); ++ta_it) {
		double total_duration = ta_it.iter->data.timeWindow.closeTime - ta_it.iter->data.timeWindow.openTime;
		reg[ta_it.iter->data.id] = Activity{ .duration = ta_it.iter->data.timeWindow.closeTime - ta_it.iter->data.timeWindow.openTime, .buckets = std::vector<Bucket>() };
		for (std::vector<double>::const_iterator left = cuts.begin(), right = cuts.begin() + 1; right != cuts.end(); ++left, ++right) {
			double duration = std::min(ta_it.iter->data.timeWindow.closeTime, *right) - std::max(ta_it.iter->data.timeWindow.openTime, *left);
			if (duration < 0) duration = 0;
			Bucket b{ .inBucket = 1, .inSolution = 1, .duration = duration, .ratio = duration/total_duration };
			reg[ta_it.iter->data.id].buckets.push_back(b);
		}
	}
	return reg;
}

void ILS::Shake(CustomSolution& sol, int& S, int& R, OP& op, const int& max_to_remove) {
	int minWalkSize = sol.getMinWalkSize();
	//std::cout << "maxToRemove: " << max_to_remove << "\t";
	//std::cout << "minWalkSize: " << minWalkSize << "\t";
	if (S >= minWalkSize) {
		S -= minWalkSize;
		R = 1;
		if (S < 1) S = 1;
	}

	if (R == max_to_remove) {
		R = 1;
	}
	//std::cout << "S: " << S << "\t";
	//std::cout << "R: " << R << "\t" << std::endl;

	for (auto& walk : sol.m_walks) {

		int start_pos{ S }, end_pos{ std::min(S + R, (int)walk.size() - 1) };
		CustomList<TA>::iterator start_it{ walk.begin() + start_pos }, end_it{ walk.begin() + end_pos };
		CustomList<TA> part(start_it, end_it);
		CustomList<TA>::iterator next = walk.erase(start_it, end_it);
		sol.m_unvisited.append(part);
		part.clear();
		updateTimes(walk, next, false, op.mTravelTimes);
	}

	S += R;
	R++;
}

int ILS::collectProfit(const CustomList<TA>::iterator& first, const CustomList<TA>::iterator& last) const {
	int sum{};
	for (CustomList<TA>::iterator it = first; it != last; ++it) { sum += it.iter->data.profit; }
	return sum;
}

Point ILS::getWeightedCentroid(const CustomList<TA>::iterator& first, const CustomList<TA>::iterator& last) {

	int total_profit{ collectProfit(first, last) };
	Position pos;
	double profit_ratio{};
	for (auto it = first; it != last; ++it) {
		profit_ratio = it.iter->data.profit / (double)total_profit;
		pos.lat += it.iter->data.point.pos.lat * profit_ratio;
		pos.lon += it.iter->data.point.pos.lon * profit_ratio;
	}
	return Point{ DEFAULT_POINT_ID, pos };
	
}

void ILS::SplitSearch(std::vector<CustomSolution>& solutions, const std::vector<double>& cuts, OP& op, std::map<std::string, Activity>& registry) {

	const int min_size = 3;

	for (size_t i = 0; i < solutions.size(); ++i) {
		if (i == 0) { //first solution
			for (size_t j = 0; j < solutions[i].m_walks.size(); ++j) { //foreach walk
				if (solutions[i].m_walks[j].empty()) { // if empty
					solutions[i].m_walks[j].push_front(*op.mStartDepot); //start point
				}

				Point cnext;
				if (solutions[i + 1].m_walks[j].size() > min_size) {
					const CustomList<TA>& next_solution_walk = solutions[i + 1].m_walks[j];
					cnext = getWeightedCentroid(next_solution_walk.at(0), next_solution_walk.at(min_size));
				}
				else {
					cnext = getWeightedCentroid(solutions[i + 1].m_unvisited.begin(), solutions[i + 1].m_unvisited.end());
				}
				op.AddPointToGraph(cnext);
				TA endDepot = TA(DEPOT_ID, cnext); //todo: delete endDepot
				endDepot.timeWindow.closeTime = cuts[i + 1];
				//endDepot->maxShift = endDepot->timeWindow.closeTime - endDepot->depTime;
				solutions[i].m_walks[j].push_back(endDepot);
				
			}

		} 
		else if (i == solutions.size() - 1) {
			for (size_t j = 0; j < solutions[i].m_walks.size(); ++j) { //foreach walk
				if (solutions[i].m_walks[j].empty()) { // if empty
					solutions[i].m_walks[j].push_back(*op.mEndDepot); //start point
				}

				int prev_index = i-1;
				while (solutions[prev_index].m_walks[j].empty()) prev_index--;

				TA prev_ta = solutions[prev_index].m_walks[j].back();
				solutions[i].m_walks[j].push_front(prev_ta);
			}
		}
		else {
			for (size_t j = 0; j < solutions[i].m_walks.size(); ++j) { //foreach walk
				//start depot
				int prev_index = i - 1;
				while (solutions[prev_index].m_walks[j].empty()) prev_index--;

				TA prev_ta = solutions[prev_index].m_walks[j].back();
				solutions[i].m_walks[j].push_front(prev_ta);

				//end depot
				Point cnext;
				if (solutions[i + 1].m_walks[j].size() > min_size) {
					const CustomList<TA>& next_solution_walk = solutions[i + 1].m_walks[j];
					cnext = getWeightedCentroid(next_solution_walk.at(0), next_solution_walk.at(min_size));
				}
				else {
					cnext = getWeightedCentroid(solutions[i + 1].m_unvisited.begin(), solutions[i + 1].m_unvisited.end());
				}
				op.AddPointToGraph(cnext);
				TA endDepot = TA(DEPOT_ID, cnext); //todo: delete endDepot
				endDepot.timeWindow.closeTime = cuts[i + 1];
				solutions[i].m_walks[j].push_back(endDepot);
			}
		}

		for (size_t j = 0; j < solutions[i].m_walks.size(); ++j) { //foreach walk
			updateTimes(solutions[i].m_walks[j], solutions[i].m_walks[j].begin(), false, op.mTravelTimes);
		}
		
		std::vector<std::string> ids = construct(solutions[i], op.mTravelTimes);
		for (auto& id : ids) {
			registry[id].buckets[i].inSolution++;
		}

		if (i > 0) {
			for (auto& walk : solutions[i].m_walks) {
				if(!walk.empty())
					walk.pop_front();
			}
				
		}

		if (i < solutions.size() - 1) {
			for (auto& walk : solutions[i].m_walks) {
				if (!walk.empty()) {
					walk.pop_back();
				}
			} 
		}
	}
}


std::vector<std::string> ILS::construct(CustomSolution& sol, const Vector2D<double>& travel_times) {
	
	std::vector<std::string> inserted_ids;

	if (sol.m_unvisited.empty()) {
		return inserted_ids;
	}

	double min_shift{}, max_ratio{ DBL_MIN }, ratio{};
	int best_arr_point_id{ DEFAULT_POINT_ID }, best_dep_point_id{ DEFAULT_POINT_ID }, 
		arr_point_id{ DEFAULT_POINT_ID }, dep_point_id{ DEFAULT_POINT_ID };
	Walks::iterator best_walk_it;
	CustomList<TA>::iterator pos{ }, best_pos{ };

	CustomList<TA>::iterator insert_it{ sol.m_unvisited.end() }, curr{}, inserted_it{};

	while (true) {
		max_ratio = DBL_MIN;
		insert_it = sol.m_unvisited.end();
		curr = sol.m_unvisited.begin();

		while (curr != sol.m_unvisited.end()) {
			min_shift = DBL_MAX;
			auto [walk_it, pos, min_shift, arr_point_id, dep_point_id] = getBestPos(curr.iter->data, sol.m_walks, travel_times);
			//default value of walk_it is sol.m_walks.end()
			//default value of pos is sol.m_walks.front().end()
			//if these iterators have these values then insertion is infeasible

			if (walk_it == sol.m_walks.end()) {
				curr++;
				continue;
			}

			ratio = pow(curr.iter->data.profit, 2) / min_shift;	//Ratioi = (Si)^2/Shifti

			if (max_ratio < ratio) //check each ratio
			{
				max_ratio = ratio;
				insert_it = curr;
				best_pos = pos;
				best_walk_it = walk_it;
				best_arr_point_id = arr_point_id;
				best_dep_point_id = dep_point_id;
			}
			curr++;
		}

		if (insert_it == sol.m_unvisited.end()) {
			break;
		}


		inserted_it = best_walk_it->insert(best_pos, insert_it.iter->data);
		inserted_ids.push_back(inserted_it.iter->data.id);
		inserted_it.iter->data.arrPointId = best_arr_point_id;
		inserted_it.iter->data.depPointId = best_dep_point_id;
		sol.m_unvisited.erase(insert_it);
		updateTimes(*best_walk_it, inserted_it, true, travel_times);
	}
	return inserted_ids;
}

std::tuple<Walks::iterator, CustomList<TA>::iterator, double, int, int> ILS::getBestPos(const TA& ta, Walks& walks, const Vector2D<double>& travel_times) {

	Walks::iterator best_walk{ walks.end() };
	CustomList<TA>::iterator best_pos{ walks.front().end() };
	int arr_point_id{ DEFAULT_POINT_ID }, dep_point_id{ DEFAULT_POINT_ID };
	double min_shift{ DBL_MAX };


	for (Walks::iterator walk_it = walks.begin(); walk_it != walks.end(); ++walk_it) {
		if (walk_it->size() < 2) {
			std::cerr << "invalid length of route" << std::endl;
			std::exit(1);
		}

		double shift{};
		CustomList<TA>::iterator left{ walk_it->begin() }, right{ walk_it->begin() + 1 };
		int temp_arr_point_id{ DEFAULT_POINT_ID }, temp_dep_point_id{ DEFAULT_POINT_ID };

		while (right != walk_it->end()) {
			shift = travel_times[left.iter->data.depPointId][ta.point.id]
				+ ta.visitDuration
				+ travel_times[ta.point.id][right.iter->data.arrPointId]
				- travel_times[left.iter->data.depPointId][right.iter->data.arrPointId];

			if (shift <= right.iter->data.maxShift) {	//check if insertion is possible
				if (shift < min_shift) {
					best_walk = walk_it;
					best_pos = right;
					min_shift = shift;
					arr_point_id = ta.point.id;
					dep_point_id = ta.point.id;
				}
			}

			left++; right++;
		}
	}


	
	return { best_walk, best_pos, min_shift, arr_point_id, dep_point_id };
}

std::tuple<CustomList<TA>::iterator, double, int, int> ILS::getBestPos(const TA& ta, const CustomList<TA>& walk, const Vector2D<double>& travel_times) {

	if (walk.size() < 2) {
		std::cerr << "invalid length of route" << std::endl;
		std::exit(1);
	}

	const TA* ptr{ &ta };
	double min_shift{ DBL_MAX }, shift{};
	CustomList<TA>::iterator best_pos{ walk.end() };

	CustomList<TA>::iterator left{ walk.begin() }, right{ walk.begin() + 1};
	
	int arr_point_id{ DEFAULT_POINT_ID }, dep_point_id{ DEFAULT_POINT_ID },
		temp_arr_point_id{ DEFAULT_POINT_ID }, temp_dep_point_id{ DEFAULT_POINT_ID };

	while (right != walk.end()) {
		shift = travel_times[left.iter->data.depPointId][ta.point.id]
			+ ta.visitDuration
			+ travel_times[ta.point.id][right.iter->data.arrPointId]
			- travel_times[left.iter->data.depPointId][right.iter->data.arrPointId];

		if (shift <= right.iter->data.maxShift) {	//check if insertion is possible
			if (shift < min_shift) {
				best_pos = right;
				min_shift = shift;
				arr_point_id = ta.point.id;
				dep_point_id = ta.point.id;
			}
		}

		left++; right++;
	}

	return { best_pos, min_shift, arr_point_id, dep_point_id };

}

void ILS::updateTimes(CustomList<TA>& walk, const CustomList<TA>::iterator& start_pos, const bool smart, const Vector2D<double>& travel_times) {
	//update times first
	CustomList<TA>::iterator it{ start_pos };

	while (it != walk.end()) {
		it.iter->data.arrTime = it.iter->previous->data.depTime + travel_times[it.iter->previous->data.depPointId][it.iter->data.arrPointId];
		it.iter->data.startOfVisitTime = it.iter->data.arrTime;
		it.iter->data.depTime = it.iter->data.startOfVisitTime + it.iter->data.visitDuration;
		++it;
	}

	updateMaxShifts(walk, travel_times);
}

void ILS::updateMaxShifts(const CustomList<TA>& li, const Vector2D<double>& travel_times) {
	CustomList<TA>::iterator it{ li.end() - 1 };
	it.iter->data.maxShift = it.iter->data.timeWindow.closeTime - it.iter->data.depTime;
	it--;
	do {
		it.iter->data.maxShift = (it + 1).iter->data.maxShift;
	} while (it-- != li.begin());
}

void ILS::print(const CustomList<TA>& li) {
	for (auto it : li) {
		std::cout << it.id << " ";
	}
	std::cout << std::endl;
}

void ILS::validate(const std::vector<CustomList<TA>>& walks, const Vector2D<double>& travel_times) {
	for (auto& walk : walks) {
		validate(walk, travel_times);
	}
}

void ILS::validate(const CustomList<TA>& walk, const Vector2D<double>& travel_times) {
	std::cout << "validating walk: "; print(walk);
	std::string msg{};
	bool valid{ true };
	CustomList<TA>::iterator next{};
	for (CustomList<TA>::iterator it = { walk.begin() }; it != walk.end(); ++it) {
		if (it.iter->data.startOfVisitTime != it.iter->data.arrTime + it.iter->data.waitDuration) {
			msg = "StartOfVisitTime(" 
				+ std::to_string(it.iter->data.startOfVisitTime) 
				+ ") != arrTime(" 
				+ std::to_string(it.iter->data.arrTime) 
				+ ")" 
				+ " + waitDuration(" 
				+ std::to_string(it.iter->data.waitDuration) 
				+ ")";
			valid = false;
			break;
		}

		if (it.iter->data.depTime != it.iter->data.startOfVisitTime + it.iter->data.visitDuration) {
			msg = "depTime(" 
				+ std::to_string(it.iter->data.depTime) 
				+ ") != startOfVisitTime(" 
				+ std::to_string(it.iter->data.startOfVisitTime) 
				+ ") + visitDuration(" 
				+ std::to_string(it.iter->data.visitDuration) 
				+ ")";
			valid = false;
			break;
		}

		next = it + 1;
		if (next != walk.end()) {
			if (next.iter->data.arrTime != it.iter->data.depTime + travel_times[it.iter->data.depPointId][next.iter->data.arrPointId]) {
				msg = "next->arrTime(" 
					+ std::to_string(next.iter->data.arrTime) 
					+ ") != depTime(" 
					+ std::to_string(it.iter->data.depTime) 
					+ ") + timeTravel[" 
					+ std::to_string(it.iter->data.depPointId) 
					+ "][" 
					+ std::to_string(next.iter->data.arrPointId) 
					+ "](" 
					+ std::to_string(travel_times[it.iter->data.depPointId][next.iter->data.arrPointId]) 
					+ ")";
				valid = false;
				break;
			}
		}
	}

	std::cout << (valid ? "walk is valid" : msg) << std::endl;
	if (!valid) {
		std::cerr << "walk is invalid: exiting.. " << std::endl; 
		std::exit(1);
	}
}

std::tuple<Walks::iterator, CustomList<TA>::iterator, double, int, int> ILS_TOPTW::getBestPos(const TA& ta, Walks& walks, const Vector2D<double>& travel_times) {

	Walks::iterator best_walk{ walks.end() };
	CustomList<TA>::iterator best_pos{ walks.front().end()};
	int arr_point_id{ DEFAULT_POINT_ID }, dep_point_id{ DEFAULT_POINT_ID };
	double min_shift{ DBL_MAX };

	for (Walks::iterator walk_it = walks.begin(); walk_it != walks.end(); ++walk_it) {
		if (walk_it->size() < 2) {
			std::cerr << "invalid length of route" << std::endl;
			std::exit(1);
		}

		double arr_time{}, wait_dur{}, start_of_visit_time{}, dep_time{}, shift{};
		CustomList<TA>::iterator left{ walk_it->begin() }, right{ walk_it->begin() + 1 };
		int temp_arr_point_id{ DEFAULT_POINT_ID }, temp_dep_point_id{ DEFAULT_POINT_ID };

		while (right != walk_it->end()) {
			arr_time = left.iter->data.depTime + travel_times[left.iter->data.depPointId][ta.point.id];
			wait_dur = std::max(0.0, ta.timeWindow.openTime - arr_time);
			start_of_visit_time = arr_time + wait_dur;
			dep_time = start_of_visit_time + ta.visitDuration;
			shift = travel_times[left.iter->data.depPointId][ta.point.id]
				+ wait_dur
				+ ta.visitDuration
				+ travel_times[ta.point.id][right.iter->data.arrPointId]
				- travel_times[left.iter->data.depPointId][right.iter->data.arrPointId];

			if (dep_time <= ta.timeWindow.closeTime && shift <= right.iter->data.maxShift) {	//check if insertion is possible
				if (shift < min_shift) {
					best_walk = walk_it; 
					best_pos = right;
					min_shift = shift;
					arr_point_id = ta.point.id;
					dep_point_id = ta.point.id;
				}
			}

			left++; right++;
		}
	}

	return { best_walk, best_pos, min_shift, arr_point_id, dep_point_id };
}

std::tuple<CustomList<TA>::iterator, double, int, int> ILS_TOPTW::getBestPos(const TA& ta, const CustomList<TA>& walk, const Vector2D<double>& travel_times) {

	if (walk.size() < 2) {
		std::cerr << "invalid length of route" << std::endl;
		std::exit(1);
	}

	double min_shift{ DBL_MAX }, arr_time{}, wait_dur{}, start_of_visit_time{}, dep_time{}, shift{};
	CustomList<TA>::iterator best_pos{ walk.end() }, left{ walk.begin() }, right{ walk.begin() + 1 };

	int arr_point_id{ DEFAULT_POINT_ID }, dep_point_id{ DEFAULT_POINT_ID },
		temp_arr_point_id{ DEFAULT_POINT_ID }, temp_dep_point_id{ DEFAULT_POINT_ID };

	while (right != walk.end()) {
		arr_time = left.iter->data.depTime + travel_times[left.iter->data.depPointId][ta.point.id];
		wait_dur = std::max(0.0, ta.timeWindow.openTime - arr_time);
		start_of_visit_time = arr_time + wait_dur;
		dep_time = start_of_visit_time + ta.visitDuration;
		shift = travel_times[left.iter->data.depPointId][ta.point.id]
			+ wait_dur
			+ ta.visitDuration
			+ travel_times[ta.point.id][right.iter->data.arrPointId]
			- travel_times[left.iter->data.depPointId][right.iter->data.arrPointId];

		if (dep_time <= ta.timeWindow.closeTime && shift <= right.iter->data.maxShift) {	//check if insertion is possible
			if (shift < min_shift) {
				best_pos = right;
				min_shift = shift;
				arr_point_id = ta.point.id;
				dep_point_id = ta.point.id;
			}
		}

		left++; right++;
	}

	return { best_pos, min_shift, arr_point_id, dep_point_id };
}

void ILS_TOPTW::updateTimes(CustomList<TA>& walk, const CustomList<TA>::iterator& start_pos, const bool smart, const Vector2D<double>& travel_times) {
	//update times first
	CustomList<TA>::iterator it{ start_pos };

	if (it == walk.begin()) {
		it.iter->data.waitDuration = std::max(0.0, it.iter->data.timeWindow.openTime - it.iter->data.arrTime);
		it.iter->data.startOfVisitTime = it.iter->data.arrTime + it.iter->data.waitDuration;
		it.iter->data.depTime = it.iter->data.startOfVisitTime + it.iter->data.visitDuration;
		it++;
	}
	
	while (it != walk.end()) {
		it.iter->data.arrTime = it.iter->previous->data.depTime + travel_times[it.iter->previous->data.depPointId][it.iter->data.arrPointId];
		it.iter->data.waitDuration = std::max(0.0, it.iter->data.timeWindow.openTime - it.iter->data.arrTime);
		it.iter->data.startOfVisitTime = it.iter->data.arrTime + it.iter->data.waitDuration;
		it.iter->data.depTime = it.iter->data.startOfVisitTime + it.iter->data.visitDuration;
		it++;
	}

	updateMaxShifts(walk, travel_times);
}

void ILS_TOPTW::updateMaxShifts(const CustomList<TA>& li, const Vector2D<double>& travel_times) {
	CustomList<TA>::iterator it{ li.end() - 1 };
	it.iter->data.maxShift = it.iter->data.timeWindow.closeTime - it.iter->data.depTime;
	it--;
	do {
		it.iter->data.maxShift = std::min(it.iter->data.timeWindow.closeTime - it.iter->data.depTime, (it + 1).iter->data.waitDuration + (it + 1).iter->data.maxShift);
	} while (it-- != li.begin());
}

void ILS_TOPTW::validate(const CustomList<CustomList<TA>>& walks, const Vector2D<double>& travel_times) {
	for (auto& walk : walks) {
		validate(walk, travel_times);
	}
}

/// <summary>
/// Validating walk:
/// In order of a walk to be valid it needs to fullfill the following requirements:
/// 1)the departure time at each sight should be less than its closing time. 
/// 2)the closing time of the last node, indicates the time budget of the problem so the total walk time should not surpass it.
/// 3)all the variable times/durations should be correct (arrTime, waitDuration, depTime) based on the travelling times between the attractions
/// </summary>
/// <param name="walk">Holds the walk that gets examined</param>
/// <param name="ttMatrix">Holds the travel times between the points of the problem</param>
void ILS_TOPTW::validate(const CustomList<TA>& walk, const Vector2D<double>& travel_times) {
	std::cout << "validating walk: "; print(walk);
	std::string msg{};
	bool valid{ true };
	CustomList<TA>::iterator next{};
	for (CustomList<TA>::iterator it = { walk.begin() }; it != walk.end(); ++it) {
		if (it.iter->data.startOfVisitTime != it.iter->data.arrTime + it.iter->data.waitDuration) {
			msg = "StartOfVisitTime("
				+ std::to_string(it.iter->data.startOfVisitTime)
				+ ") != arrTime("
				+ std::to_string(it.iter->data.arrTime)
				+ ")"
				+ " + waitDuration("
				+ std::to_string(it.iter->data.waitDuration)
				+ ")";
			valid = false;
			break;
		}

		if (it.iter->data.depTime != it.iter->data.startOfVisitTime + it.iter->data.visitDuration) {
			msg = "depTime("
				+ std::to_string(it.iter->data.depTime)
				+ ") != startOfVisitTime("
				+ std::to_string(it.iter->data.startOfVisitTime)
				+ ") + visitDuration("
				+ std::to_string(it.iter->data.visitDuration)
				+ ")";
			valid = false;
			break;
		}

		if (it.iter->data.depTime > it.iter->data.timeWindow.closeTime) {
			msg = "depTime("
				+ std::to_string(it.iter->data.depTime)
				+ ") > closeTime("
				+ std::to_string(it.iter->data.timeWindow.closeTime)
				+ ")";
			valid = false;
			break;
		}

		next = it + 1;
		if (next != walk.end()) {
			if (next.iter->data.arrTime != it.iter->data.depTime + travel_times[it.iter->data.depPointId][next.iter->data.arrPointId]) {
				msg = "next->arrTime("
					+ std::to_string(next.iter->data.arrTime)
					+ ") != depTime("
					+ std::to_string(it.iter->data.depTime)
					+ ") + timeTravel["
					+ std::to_string(it.iter->data.depPointId)
					+ "]["
					+ std::to_string(next.iter->data.arrPointId)
					+ "]("
					+ std::to_string(travel_times[it.iter->data.depPointId][next.iter->data.arrPointId])
					+ ")";
				valid = false;
				break;
			}
		}
	}

	std::cout << (valid ? "walk is valid" : msg) << std::endl;
	if (!valid) {
		std::cerr << "walk is invalid: exiting.. " << std::endl;
		std::exit(1);
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


ListTA ILS::setBucketActivityDurations(ListTA& unvisited, double avgEvent, std::vector<double> cuts) {

	TA* curr = unvisited.first();
	while (curr != nullptr) {
		for (size_t i = 0; i < cuts.size() - 1; ++i) {
			//Bucket bucket{ .duration = avgEvent };
			//curr->metrics.bucketActivities.push_back(activity);
		}
		// for(auto& c: cuts){
		// 	BucketActivity activity{.duration = avgEvent}
		// }
		// BucketActivity activity{.duration = avgEvent}
		// curr->metrics.bucketActivities.push_back(std::abs())
		// curr->bucketActivities[0].duration = avgEvent - curr->timeWindow.openTime;
		// curr->bucketActivities[1].duration = curr->timeWindow.closeTime - avgEvent;
		curr = curr->next;
	}
	return unvisited;

}

//std::vector<std::vector<TA*>> ILS::getBuckets(std::vector<TA*> nodes, int m) {
//
//	std::vector<std::vector<TA*>> buckets;
//	int offset = 0, remainder = nodes.size() % m, bucketSize = nodes.size() / m;
//
//	std::sort(nodes.begin(), nodes.end(), &compareTimeWindowCenter);
//
//	for (int i = 0; i < m; ++i) {
//		if (remainder > 0) {
//			std::vector<TA*> slice(nodes.begin() + offset, nodes.begin() + offset + bucketSize + 1);
//			buckets.push_back(slice);
//			remainder--;
//		}
//		else {
//			std::vector<TA*> slice(nodes.begin() + offset, nodes.begin() + offset + bucketSize);
//			buckets.push_back(slice);
//		}
//	}
//
//	return buckets;
//}

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

int ILS::collectScore(std::vector<Solution> solutions) {
	int totalScore = 0;
	for (auto& sol : solutions) {
		totalScore += sol.getScore();
	}
	return totalScore;
}

std::tuple<int, int> ILS::getMinMaxLength(std::vector<Solution> solutions) {
	int min = INT_MAX, max = INT_MIN;

	for (auto& sol : solutions) {
		int walkLength = sol.mWalk.size();
		if (walkLength > max) {
			max = walkLength;
		}

		if (walkLength < min) {
			min = walkLength;
		}
	}

	return { min, max };
}
