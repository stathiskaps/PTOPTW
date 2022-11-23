#include "ILS.h"

ILS::ILS() : mBucketsNum(0) {}

ILS::ILS(int bucketsNum, int intervalsNum) : mBucketsNum(bucketsNum), mIntervalsNum(intervalsNum) {}

ILS::~ILS() {}

void ILS::dbScan(OP& op){
	typedef ClusteringDatum<2, double, 1, float, uint16_t> CDat_t;

    constexpr size_t MaxElementsInANode = 6; // 16, 32, 128, 256, ... ?
    typedef boost::geometry::index::rstar<MaxElementsInANode> RTreeParameter_t;
    typedef boost::geometry::index::rtree<CDat_t,RTreeParameter_t> RTree_t;
    typedef boost::geometry::model::box<CDat_t> Box_t;

    RTree_t rtree;
	double mean_travel_time = 0;
	int n = 0;
	std::vector<std::vector<int>>::const_iterator row; 
    std::vector<int>::const_iterator col; 

	for (std::vector<std::vector<double>>::iterator row = op.mTravelTimes.begin(); row != op.mTravelTimes.end(); ++row){
		for(std::vector<double>::iterator col = row->begin(); col != row->end(); ++col){
			mean_travel_time += *col;
			n++;
		}
	}
	mean_travel_time = mean_travel_time / n;
	std::cout << "mean_travel_time: " << mean_travel_time << std::endl;

	for(std::vector<TA*>::iterator ta_it = op.mAttractions.begin(); ta_it != op.mAttractions.end(); ++ta_it){
		rtree.insert(CDat_t( { (*ta_it)->point.pos.lat, (*ta_it)->point.pos.lon }, {0.0f} ));
	}

    //const size_t MinPts = CDat_t::SpatialDimensionCount_*2;
    //const SpatialQueryTechnique UsersSpatialQueryTechnique = SpatialQueryTechnique::UseWithin;
    const double Eps = 10.0;
    
    DBSCAN<RTree_t,CDat_t>(rtree,Eps);

	//Print out the points with cluster info.
    if(true){
        constexpr auto RTreeSpatialQueryGetAll = [](const CDat_t &) -> bool { return true; };
        RTree_t::const_query_iterator it;
        it = rtree.qbegin(boost::geometry::index::satisfies( RTreeSpatialQueryGetAll ));
        for( ; it != rtree.qend(); ++it){
            std::cout << "  Point: " << boost::geometry::wkt(*it)
                      << "\t\t Attribute[0]: " << it->Attributes[0]
                      << "\t\t ClusterID: " << it->CID.ToText()
                      << std::endl;
        }

		std::cout << "finished" << std::endl;
    }

	if(true){
        constexpr auto RTreeSpatialQueryGetAll = [](const CDat_t &) -> bool { return true; };
        std::ofstream svg("Visualized.svg");
        boost::geometry::svg_mapper<CDat_t> mapper(svg, 1280, 1024);

        //Add the items to the map. The virtual bounds will be computed to accomodate.
        // Also keep a record of the distinct clusters encountered.
        std::set<typename CDat_t::ClusterIDType_> RawCIDs;

        RTree_t::const_query_iterator it;
        it = rtree.qbegin(boost::geometry::index::satisfies( RTreeSpatialQueryGetAll ));
        for( ; it != rtree.qend(); ++it){
            mapper.add(*it);
            RawCIDs.insert( it->CID.Raw );
        }
        std::cout << RawCIDs.size() << " distinct ClusterIDs encountered." << std::endl;

        //Create a mapping between the ClusterIDs and a pseudo-random RGB colour.
        size_t ColourSeed = 9137;
        std::mt19937 re(ColourSeed);
        std::uniform_real_distribution<> rdC(0.0, 1.0);
        std::uniform_int_distribution<> rdA( 50, 210);
        std::uniform_int_distribution<> rdB( 20, 125);
        std::map<typename CDat_t::ClusterIDType_, std::string> ColoursA, ColoursB;
        for(const auto RawCID : RawCIDs){
            ColoursA[RawCID] = std::to_string((rdC(re) > 0.33) ? rdA(re) : 230) + std::string(",")
                             + std::to_string((rdC(re) > 0.33) ? rdA(re) : 230) + std::string(",")
                             + std::to_string((rdC(re) > 0.33) ? rdA(re) : 230);
            ColoursB[RawCID] = std::to_string((rdC(re) > 0.33) ? rdB(re) :  10) + std::string(",")
                             + std::to_string((rdC(re) > 0.33) ? rdB(re) :  10) + std::string(",")
                             + std::to_string((rdC(re) > 0.33) ? rdB(re) :  10);
        }

        //Actually draw the items.
        it = rtree.qbegin(boost::geometry::index::satisfies( RTreeSpatialQueryGetAll ));
        for( ; it != rtree.qend(); ++it){
            std::stringstream ss;
            ss << "fill-opacity:0.80; "
               << "fill:rgb(" << ColoursA[it->CID.Raw] << "); "
               << "stroke-opacity:0.90; "
               << "stroke:rgb(" << ColoursB[it->CID.Raw] << "); "
               << "stroke-width:1";
            //mapper.map(*it, "fill-opacity:0.75; fill:rgb(75,100,0); stroke:rgb(30,40,0); stroke-width:2", 5);
            mapper.map(*it, ss.str(), 6);
        }
    }
}

std::map<std::string, std::vector<double>> ILS::getActivities(List<TA>& unvisited, std::vector<ILS::Interval> intervals){
	std::map<std::string, std::vector<double>> activities;
	for(auto& n : unvisited) {
		std::vector<double> durations;
		for(auto& i : intervals) {
			double active_dur = std::min(i.end_time, n.timeWindow.closeTime) - std::max(i.start_time, n.timeWindow.openTime);
			if (active_dur < 0) active_dur = 0;
			durations.push_back(active_dur/(n.timeWindow.closeTime-n.timeWindow.openTime));
		}
		activities.insert({n.id, durations});
	}
	return activities;
}

std::vector<ILS::Interval> ILS::getIntervals(std::vector<TA*> unvisited, int intervals_num, double day_start_time, double day_close_time){
	std::vector<ILS::Interval> intervals;
	double diff = day_close_time - day_start_time;
	double dur = diff/intervals_num;

	double start = day_start_time;
	double end = day_start_time + dur;

	for(size_t i = 0; i < intervals_num; ++i){
		intervals.push_back(ILS::Interval{start, end});
		start+=dur;
		end+=dur;
	}

	return intervals;
}

std::map<std::string, std::vector<ILS::Usage>> ILS::initRegistry(List<TA>& unvisited, std::vector<ILS::Interval> intervals) {
	std::map<std::string, std::vector<ILS::Usage>> reg;

	for(auto& n : unvisited){
		std::vector<ILS::Usage> usages;
		for(auto &i : intervals){
			usages.push_back(Usage());
		}
		reg.insert({n.id, usages});
	}
	return reg;
}

void ILS::printSolution(const std::string tag, const Solution& sol){
	std::cout << "===================================================" << std::endl;
	std::cout << tag << std::endl;
	sol.m_unvisited.print("Unvisited");
	for(std::vector<List<TA>>::const_iterator walk_it = sol.m_walks.begin(); walk_it != sol.m_walks.end(); ++walk_it){
		const int index = walk_it - sol.m_walks.begin();
		walk_it->print("Walk " + std::to_string(index));
	}
	std::cout << "===================================================" << std::endl;
}


void ILS::printSolutions(const std::string tag, const std::vector<Solution>& sols){
	std::cout << tag << std::endl;
	for(std::vector<Solution>::const_iterator sol_it = sols.begin(); sol_it != sols.end(); ++sol_it){
		const int index = sol_it - sols.begin();
		const std::string tag = "Solution [" + std::to_string(index) + "]: ";
		printSolution(tag, *sol_it);
		std::cout << std::endl;
	}
}

void ILS::InitSolutions(std::vector<Solution>& sols, const OP& op){
	//go to each walk of last solution and add an endDepot
	for (auto &walk : sols.front().m_walks) {
		if (walk.front().id != START_DEPOT_ID) {
			walk.push_front(*op.mStartDepot);
		}
	}

	//go to each walk of last solution and add an endDepot
	for (auto &walk : sols.back().m_walks) {
		if (walk.back().id != END_DEPOT_ID) {
			walk.push_back(*op.mEndDepot);
		}
	}
}


void ILS::SolveNew(OP& op) {

	// std::cout.setstate(std::ios_base::failbit);
	// dbScan(op);
	// return;

	auto start = std::chrono::steady_clock::now();

	std::cout << "Problem will be divided in " << mIntervalsNum << " intervals." << std::endl;

	//todo: for pr13 and buckets_num 4 got invalid bucket error
	const int num_locations = op.mAttractions.size();
	const int max_to_remove = num_locations / (3 * op.m_walks_num);

	List<TA> unvisited;
	for (auto ta : op.mAttractions) {
		unvisited.push_back(*ta);
	}

	std::vector<ILS::Interval> intervals = getIntervals(op.mAttractions, mIntervalsNum, op.mStartDepot->timeWindow.openTime, op.mEndDepot->timeWindow.closeTime);
	const std::map<std::string, std::vector<double>> activities = getActivities(unvisited, intervals);
	std::map<std::string, std::vector<ILS::Usage>> reg = initRegistry(unvisited, intervals); 

	int counter{}, S{ 1 }, R{ 1 }, times_not_improved{ 0 }, best_score{ INT_MIN };
	Solution best_solution;

	std::vector<ILS::SR> shake_settings;
	for(size_t i = 0; i < mIntervalsNum; ++i){
		shake_settings.push_back(ILS::SR{1, 1});
	}

	std::vector<Solution> proc_solutions(mIntervalsNum, Solution()), best_solutions;
	//Initialize solutions
	for (auto& s : proc_solutions) {
		for (size_t i = 0; i < op.m_walks_num; ++i) {
			s.m_walks.push_back(List<TA>());
		}
	}
	InitSolutions(proc_solutions, op);
	List<TA> pool = std::move(unvisited);


	while (times_not_improved < MAX_TIMES_NOT_IMPROVED) {
		counter++;

		// std::cout << "Revision " << counter << std::endl;
		
		auto split_search_start = std::chrono::steady_clock::now();
		splitUnvisitedList(proc_solutions, pool, mIntervalsNum, reg, activities);
		SplitSearch(proc_solutions, intervals, op, reg);
		auto split_search_end = std::chrono::steady_clock::now();
		auto split_search_diff = split_search_end - split_search_start;

		// std::cout << "Split search execution time (rev " << std::to_string(counter) << "): " << std::chrono::duration <double, std::milli> (split_search_diff).count() << " ms" << std::endl;

		int score = collectScores(proc_solutions);

		if (score > best_score) {
			best_score = score;
			best_solutions = proc_solutions;
			for(auto& sr : shake_settings){
				sr.R = 1;
			}
			times_not_improved = 0;
		}
		else {
			times_not_improved++;
		}

		int removed_counter = SplitShake(proc_solutions, shake_settings, op, max_to_remove);
		// std::cout << "Removed " << removed_counter << " during split shake" << std::endl;
		gatherUnvisited(proc_solutions, pool);

	}
	auto end = std::chrono::steady_clock::now();
	auto diff = end - start;

	std::cout.clear();
	best_solution = connectSolutions(proc_solutions, op.m_walks_num);
	for(auto& walk : best_solution.m_walks){
		updateTimes(walk, walk.begin(), -1, false, op.mTravelTimes);
	}
	validate(best_solution.m_walks, op.mTravelTimes);
	std::cout << "Best score: " << best_score << std::endl;
	std::cout << "Visits: " << best_solution.getVisits() << std::endl;
	std::cout << "Execution time: " << std::chrono::duration <double, std::milli> (diff).count() << " ms" << std::endl;
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
			for (auto& p : slice) b.unvisited.push_back(*p);
			bins.push_back(b);
			offset += bucketSize + 1;
			remainder--;
		}
		else {
			std::vector<TA*> slice(unvisited.begin() + offset, unvisited.begin() + offset + bucketSize);
			Bin b;
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

void ILS::gatherUnvisited(std::vector<Solution>& solutions, List<TA>& pool){
	pool.clear();
	for(auto& sol : solutions){
		pool.append(sol.m_unvisited);
		sol.m_unvisited.clear();
	}
}

std::vector<List<TA>> ILS::splitUnvisitedList(std::vector<Solution>& solutions, List<TA>& pool, int intervals_num, std::map<std::string, std::vector<ILS::Usage>>& reg, std::map<std::string, std::vector<double>> activities) {

	if(solutions.size() != intervals_num){
		throw std::runtime_error("solutions size is indifferent from intervals number");
	}


	std::vector<List<TA>> lists(intervals_num, List<TA>());
	//Split Unvisited
	for (auto& ta : pool) {
		const double total_duration{ ta.timeWindow.closeTime - ta.timeWindow.openTime };
		double max_score{ -1 };

		if(reg[ta.id].size() != activities[ta.id].size()){
			throw std::runtime_error("vectors are not of equal size");
		}

		double best_score = DBL_MIN;
		std::vector<ILS::Usage>::iterator best_it{ reg[ta.id].end() };
		std::vector<ILS::Usage>::iterator usage_it;
		std::vector<double>::iterator duratio_it;
		for(usage_it = reg[ta.id].begin(), duratio_it = activities[ta.id].begin(); usage_it != reg[ta.id].end() ; ++usage_it, ++duratio_it){
			double score = *duratio_it * (usage_it->solved / static_cast<double>(usage_it->imported));
			if (score>best_score){
				best_score = score;
				best_it = usage_it;
			}
		}
		if(best_it == reg[ta.id].end()){
			throw std::runtime_error("didn't find a good interval for node " + ta.id);
		}

		const int index = best_it - reg[ta.id].begin();
		solutions.at(index).m_unvisited.push_back(ta);
		reg[ta.id].at(index).imported += 1;
	}

	//check if i can implement emplace back instead of creating new object
	pool.clear();

	return lists;
}

Solution ILS::connectSolutions(std::vector<Solution>& sols, const size_t walks_num) {
	Solution solution;
	//Connect unvisited
	for (auto& s : sols) {
		gatherUnvisited(sols, solution.m_unvisited);
	}
	for (size_t i = 0; i < walks_num; ++i) {
		List<TA> walk;
		for (auto& s : sols) {
			walk.append(s.m_walks[i]);
			s.m_walks[i].clear();
		}
		solution.m_walks.push_back(walk);
	}

	return solution;
}

int ILS::collectScores(std::vector<Solution> sols) {
	int total_score{};
	for (auto& sol : sols) {
		total_score += sol.getScores();
	}
	return total_score;
}

void ILS::PrepareForShake(std::vector<Solution>& sols){
	for(std::vector<Solution>::iterator sol_it = sols.begin(); sol_it != sols.end(); ++sol_it){
		if(sol_it == sols.begin()){
			for(std::vector<List<TA>>::iterator walk_it = sol_it->m_walks.begin(); walk_it != sol_it->m_walks.end(); ++walk_it){
				// if(walk_it->empty()) continue;
				walk_it->push_back(walk_it->back());
				walk_it->back().id = DUMMY_ID;;
			}
		} else if (sol_it == sols.end() - 1){
			for(std::vector<List<TA>>::iterator walk_it = sol_it->m_walks.begin(); walk_it != sol_it->m_walks.end(); ++walk_it){
				// if(walk_it->empty()) continue;
				walk_it->push_front(walk_it->front());
				walk_it->front().id = DUMMY_ID;
			}
		} else {
			for(std::vector<List<TA>>::iterator walk_it = sol_it->m_walks.begin(); walk_it != sol_it->m_walks.end(); ++walk_it){
				// if(walk_it->empty()) continue;
				walk_it->push_back(walk_it->back());
				walk_it->back().id = DUMMY_ID;;

				walk_it->push_front(walk_it->front());
				walk_it->front().id = DUMMY_ID;
			}
		}
	}
}

void ILS::RemoveDummyNodes(std::vector<Solution>& sols){
	for(std::vector<Solution>::iterator sol_it = sols.begin(); sol_it != sols.end(); ++sol_it){
		for(std::vector<List<TA>>::iterator walk_it = sol_it->m_walks.begin(); walk_it != sol_it->m_walks.end(); ++walk_it){
			for(List<TA>::iterator ta_it = walk_it->begin(); ta_it != walk_it->end();){
				if(ta_it.iter->data.id == DUMMY_ID){
					ta_it = walk_it->erase(ta_it);
				} else {
					ta_it++;
				}
			}
		}
	}
}

int ILS::SplitShake(std::vector<Solution>& sols, std::vector<ILS::SR>& shake_settings, OP& op, const int& max_to_remove){
	int removed_counter = 0;
	PrepareForShake(sols);
	for(auto sol_it = sols.begin(); sol_it != sols.end(); ++sol_it){
		const int sol_index = sol_it - sols.begin();
		int removed = Shake(*sol_it, shake_settings[sol_index].S, shake_settings[sol_index].R, op, max_to_remove);
		removed_counter += removed;
	}
	RemoveDummyNodes(sols);
	return removed_counter;
}

int ILS::Shake(Solution& sol, int& S, int& R, OP& op, const int& max_to_remove) {
	int minWalkSize = sol.getMinWalkSize();
	int removed_counter = 0;
	//std::cout << "maxToRemove: " << max_to_remove << "\t";
	//std::cout << "minWalkSize: " << minWalkSize << "\t";
	if (S >= minWalkSize - 2) {
		S = 1;
		R = 1;
	}

	// if (S < 1) S = 1;

	if (R == max_to_remove) {
		R = 1;
	}
	//std::cout << "S: " << S << "\t";
	//std::cout << "R: " << R << "\t" << std::endl;

	

	for (auto& walk : sol.m_walks) {
		if(walk.size() == 2) continue;
		int start_pos{ S }, end_pos{ std::min(S + R, (int)walk.size() - 1) };
		removed_counter += end_pos - start_pos;
		List<TA>::iterator start_it{ walk.begin() + start_pos }, end_it{ walk.begin() + end_pos };
		List<TA> part(start_it, end_it);
		List<TA>::iterator next = walk.erase(start_it, end_it);
		sol.m_unvisited.append(part);
		part.clear();
		updateTimes(walk, next, -1, false, op.mTravelTimes);
	}
	S += R;
	R++;

	return removed_counter;
}

int ILS::collectProfit(const List<TA>::iterator& first, const List<TA>::iterator& last) const {
	int sum{};
	for (List<TA>::iterator it = first; it != last; ++it) { sum += it.iter->data.profit; }
	return sum;
}

Point ILS::getWeightedCentroid(const List<TA>::iterator& first, const List<TA>::iterator& last) {

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

bool ILS::hasWeightedCentroid(const Solution& sol, const int walk_index, const int min_walk_size) {
	return sol.m_walks[walk_index].size() >= min_walk_size || !sol.m_unvisited.empty();
}

// double ILS::getShift(const List<TA>::iterator left, const List<TA>::iterator ins, const List<TA>::iterator right, const Vector2D<double>& travel_times){
// 	double start_time = left.iter->data
// 	double shift = travel_times[left.iter->data.depPointId][ins.point.id]
// 				+ ins.waitDuration
// 				+ ins.visitDuration
// 				+ travel_times[ins.point.id][right.iter->data.arrPointId]
// 				- travel_times[left.iter->data.depPointId][right.iter->data.arrPointId];

// 	return shift;
// }

// bool ILS::TrimLeft(TA ta, List<TA> walk){
// 	List<TA>::iterator curr = walk.begin();
// 	for(List<TA>::iterator curr = walk.begin(); curr!=walk)
// }

void ILS::AddStartDepots(std::vector<Solution>& solutions, const std::vector<ILS::Interval>& intervals, const int i, const OP& op){
	for (size_t j = 0; j < solutions[i].m_walks.size(); ++j) {
		int prev_sol_index = i-1;
		
		if(first_solution){
			solutions[i].m_walks[j].push_front(*op.mStartDepot);
			solutions[i].m_walks[j].front().depTime = intervals[i].start_time;
		} else {
			while(solutions[prev_sol_index].m_walks[j].empty() && prev_sol_index >= 0){
				prev_sol_index--;
				continue;
			}

			if(prev_sol_index == -1){
				throw std::runtime_error("didn't find a valid previous walk to get a startDepot");
			}

			List<TA>::iterator prev_it = solutions[prev_sol_index].m_walks[j].end() - 1;

			if(!solutions[i].m_walks[j].empty()) { 
				List<TA>::iterator curr = solutions[i].m_walks[j].begin();
				while(curr != solutions[i].m_walks[j].end()){ // left trim walk
					auto [valid, _] = insertionIsValid(prev_it.iter->data, curr.iter->data, op.mTravelTimes);
					if(!valid){
						solutions[i].m_unvisited.push_back(curr.iter->data);
						curr = solutions[i].m_walks[j].erase(curr);
						updateTimes(solutions[i].m_walks[j], curr, intervals[i].start_time, false, op.mTravelTimes);
					} else {
						break;
					}
				}
			}

			solutions[i].m_walks[j].push_front(prev_it.iter->data);
		}

	}
}

void ILS::AddEndDepots(std::vector<Solution>& solutions, const std::vector<ILS::Interval>& intervals, const int i, OP& op){
	const int min_size = 3;
	for (size_t j = 0; j < solutions[i].m_walks.size(); ++j) {
		//add endpoint
		int next_sol_index = i + 1;
		while (next_sol_index < solutions.size() && !hasWeightedCentroid(solutions[next_sol_index], j, min_size)) {
			next_sol_index++;
		}

		if (last_solution || next_sol_index == solutions.size()) {
			if(solutions[i].m_walks[j].back().id != op.mEndDepot->id){
				solutions[i].m_walks[j].push_back(*op.mEndDepot);
				solutions[i].m_walks[j].back().timeWindow.closeTime = intervals[i].end_time;
			}
		}
		else {
			Point cnext;
			if (solutions[next_sol_index].m_walks[j].size() > min_size) {
				const List<TA>& next_solution_walk = solutions[next_sol_index].m_walks[j];
				cnext = getWeightedCentroid(next_solution_walk.at(0), next_solution_walk.at(min_size));
			}
			else {
				cnext = getWeightedCentroid(solutions[next_sol_index].m_unvisited.begin(), solutions[next_sol_index].m_unvisited.end());
			}
			op.AddPointToGraph(cnext);
			TA endDepot = TA(CNEXT_ID, cnext); //todo: delete endDepot
			endDepot.timeWindow.closeTime = intervals[i].end_time;
			solutions[i].m_walks[j].push_back(endDepot);
		}
	}
}

void ILS::SplitSearch(std::vector<Solution>& solutions, const std::vector<ILS::Interval>& intervals, OP& op, std::map<std::string, std::vector<ILS::Usage>>& reg) {

	const int min_size = 3;
	const size_t solutions_size = solutions.size();

	for (size_t i = 0; i < solutions.size(); ++i) {

		if (solutions[i].m_unvisited.empty()) continue;

		if (first_solution) {
			AddEndDepots(solutions, intervals, i, op);
		} 
		else if (last_solution) {
			AddStartDepots(solutions, intervals, i, op);
		}
		else if (middle_solution){
			AddStartDepots(solutions, intervals, i, op);
			AddEndDepots(solutions, intervals, i, op);
		}

		for (size_t j = 0; j < solutions[i].m_walks.size(); ++j) { //foreach walk
			updateTimes(solutions[i].m_walks[j], solutions[i].m_walks[j].begin(), -1, false, op.mTravelTimes);
		}
		
		std::vector<std::string> ids = construct(solutions[i], op.mTravelTimes);
		for (auto& id : ids) {
			reg[id].at(i).solved++;
			// registry[id].buckets[i].inSolution++;
		}

		if (i > 0) {
			for(std::vector<List<TA>>::iterator walk_it = solutions[i].m_walks.begin(); walk_it != solutions[i].m_walks.end(); ++walk_it){
				walk_it->pop_front();
			}
		}

		if (i < solutions.size() - 1) {
			for (auto& walk : solutions[i].m_walks) {
				walk.pop_back();
			} 
		}
	}
}


std::vector<std::string> ILS::construct(Solution& sol, const Vector2D<double>& travel_times) {
	
	std::vector<std::string> inserted_ids;

	if (sol.m_unvisited.empty()) {
		return inserted_ids;
	}

	double min_shift{}, max_ratio{ DBL_MIN }, ratio{};
	int best_arr_point_id{ DEFAULT_POINT_ID }, best_dep_point_id{ DEFAULT_POINT_ID }, 
		arr_point_id{ DEFAULT_POINT_ID }, dep_point_id{ DEFAULT_POINT_ID };
	Walks::iterator best_walk_it;
	List<TA>::iterator pos{ }, best_pos{ };

	List<TA>::iterator insert_it{ sol.m_unvisited.end() }, curr{}, inserted_it{};

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
		inserted_it.iter->data.walk = best_walk_it - sol.m_walks.begin();
		inserted_ids.push_back(inserted_it.iter->data.id);
		inserted_it.iter->data.arrPointId = best_arr_point_id;
		inserted_it.iter->data.depPointId = best_dep_point_id;
		sol.m_unvisited.erase(insert_it);
		updateTimes(*best_walk_it, inserted_it, -1, true, travel_times);
	}
	return inserted_ids;
}

std::tuple<bool, double> ILS::insertionIsValid(const TA& left, const TA& ta, const TA& right, const Vector2D<double>& travel_times){
	double shift = travel_times[left.depPointId][ta.point.id]
		+ ta.visitDuration
		+ travel_times[ta.point.id][right.arrPointId]
		- travel_times[left.depPointId][right.arrPointId];
	return { shift <= right.maxShift, shift };
}

std::tuple<bool, double> ILS::insertionIsValid(const TA& ta, const TA& right, const Vector2D<double>& travel_times){
	double shift = ta.visitDuration
		+ travel_times[ta.point.id][right.arrPointId];
	return { shift <= right.maxShift, shift };
}

std::tuple<Walks::iterator, List<TA>::iterator, double, int, int> ILS::getBestPos(const TA& ta, Walks& walks, const Vector2D<double>& travel_times) {

	Walks::iterator best_walk{ walks.end() };
	List<TA>::iterator best_pos{ walks.front().end() };
	int arr_point_id{ DEFAULT_POINT_ID }, dep_point_id{ DEFAULT_POINT_ID };
	double min_shift{ DBL_MAX };


	for (Walks::iterator walk_it = walks.begin(); walk_it != walks.end(); ++walk_it) {
		if (walk_it->size() < 2) {
			std::cerr << "invalid length of route" << std::endl;
			std::exit(1);
		}

		double shift{};
		List<TA>::iterator left{ walk_it->begin() }, right{ walk_it->begin() + 1 };
		int temp_arr_point_id{ DEFAULT_POINT_ID }, temp_dep_point_id{ DEFAULT_POINT_ID };

		while (right != walk_it->end()) {
			auto [valid, shift] = insertionIsValid(left.iter->data, ta, right.iter->data, travel_times);
			if(valid){
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

void ILS::updateTimes(List<TA>& walk, const List<TA>::iterator& start_pos, const double start_arr_time, const bool smart, const Vector2D<double>& travel_times) {
	//update times first
	List<TA>::iterator it{ start_pos };

	while (it != walk.end()) {
		it.iter->data.arrTime = it.iter->previous->data.depTime + travel_times[it.iter->previous->data.depPointId][it.iter->data.arrPointId];
		it.iter->data.startOfVisitTime = it.iter->data.arrTime;
		it.iter->data.depTime = it.iter->data.startOfVisitTime + it.iter->data.visitDuration;
		++it;
	}

	updateMaxShifts(walk, travel_times);
}

void ILS::updateMaxShifts(const List<TA>& li, const Vector2D<double>& travel_times) {
	List<TA>::iterator it{ li.end() - 1 };
	it.iter->data.maxShift = it.iter->data.timeWindow.closeTime - it.iter->data.depTime;
	it--;
	do {
		it.iter->data.maxShift = (it + 1).iter->data.maxShift;
	} while (it-- != li.begin());
}

void ILS::validate(const std::vector<Solution>& sols, const Vector2D<double>& travel_times){
	for(auto sol_it = sols.begin(); sol_it != sols.end(); ++sol_it){
		const int index = sol_it - sols.begin();
		std::cout << "Checking solution [" << index << "]: ";
		validate(sol_it->m_walks, travel_times);
	}
}

void ILS::validate(const std::vector<List<TA>>& walks, const Vector2D<double>& travel_times) {
	for (auto& walk : walks) {
		validate(walk, travel_times);
	}
}

void ILS::validate(const List<TA>& walk, const Vector2D<double>& travel_times) {
	walk.print("validating walk");
	std::string msg{};
	bool valid{ true };
	List<TA>::iterator next{};
	for (List<TA>::iterator it = { walk.begin() }; it != walk.end(); ++it) {
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

std::tuple<bool, double> ILS_TOPTW::insertionIsValid(const TA& left, const TA& ta, const TA& right, const Vector2D<double>& travel_times){
	double arr_time{}, wait_dur{}, start_of_visit_time{}, dep_time{}, shift{};

	arr_time = left.depTime + travel_times[left.depPointId][ta.point.id];
	wait_dur = std::max(0.0, ta.timeWindow.openTime - arr_time);
	start_of_visit_time = arr_time + wait_dur;
	dep_time = start_of_visit_time + ta.visitDuration;
	shift = travel_times[left.depPointId][ta.point.id]
		+ wait_dur
		+ ta.visitDuration
		+ travel_times[ta.point.id][right.arrPointId]
		- travel_times[left.depPointId][right.arrPointId];

	return { dep_time <= ta.timeWindow.closeTime && shift <= right.maxShift, shift };
}

std::tuple<bool, double> ILS_TOPTW::insertionIsValid(const TA& ta, const TA& right, const Vector2D<double>& travel_times){
	double arr_time{}, wait_dur{}, start_of_visit_time{}, dep_time{}, shift{};
	shift = ta.visitDuration
		+ travel_times[ta.point.id][right.arrPointId];

	return { shift <= right.maxShift, shift };
}

std::tuple<Walks::iterator, List<TA>::iterator, double, int, int> ILS_TOPTW::getBestPos(const TA& ta, Walks& walks, const Vector2D<double>& travel_times) {

	Walks::iterator best_walk{ walks.end() };
	List<TA>::iterator best_pos{ walks.front().end()};
	int arr_point_id{ DEFAULT_POINT_ID }, dep_point_id{ DEFAULT_POINT_ID };
	double min_shift{ DBL_MAX };

	for (Walks::iterator walk_it = walks.begin(); walk_it != walks.end(); ++walk_it) {
		if (walk_it->size() < 2) {
			std::cerr << "invalid length of route" << std::endl;
			std::exit(1);
		}

		double arr_time{}, wait_dur{}, start_of_visit_time{}, dep_time{}, shift{};
		List<TA>::iterator left{ walk_it->begin() }, right{ walk_it->begin() + 1 };
		int temp_arr_point_id{ DEFAULT_POINT_ID }, temp_dep_point_id{ DEFAULT_POINT_ID };

		while (right != walk_it->end()) {
			auto [valid, shift] = insertionIsValid(left.iter->data, ta, right.iter->data, travel_times);
			if (valid) {	//check if insertion is possible
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

void ILS_TOPTW::updateTimes(List<TA>& walk, const List<TA>::iterator& start_pos, const double start_arr_time, const bool smart, const Vector2D<double>& travel_times) {
	//update times first
	List<TA>::iterator it{ start_pos };
	if (it == walk.begin()) {
		if(start_arr_time >= 0){ it.iter->data.arrTime = start_arr_time; }
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

void ILS_TOPTW::updateMaxShifts(const List<TA>& li, const Vector2D<double>& travel_times) {
	List<TA>::iterator it{ li.end() - 1 };
	it.iter->data.maxShift = it.iter->data.timeWindow.closeTime - it.iter->data.depTime;
	it--;
	do {
		it.iter->data.maxShift = std::min(it.iter->data.timeWindow.closeTime - it.iter->data.depTime, (it + 1).iter->data.waitDuration + (it + 1).iter->data.maxShift);
	} while (it-- != li.begin());
}

void ILS_TOPTW::validate(const List<List<TA>>& walks, const Vector2D<double>& travel_times) {
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
void ILS_TOPTW::validate(const List<TA>& walk, const Vector2D<double>& travel_times) {
	walk.print("validating walk");
	std::string msg{};
	bool valid{ true };
	List<TA>::iterator next{};
	for (List<TA>::iterator it = { walk.begin() }; it != walk.end(); ++it) {
		if (it.iter->data.startOfVisitTime != it.iter->data.arrTime + it.iter->data.waitDuration) {
			msg = "["+it.iter->data.id+"]StartOfVisitTime("
				+ std::to_string(it.iter->data.startOfVisitTime)
				+ ") != ["+it.iter->data.id+"]arrTime("
				+ std::to_string(it.iter->data.arrTime)
				+ ")"
				+ "+ ["+it.iter->data.id+"]waitDuration("
				+ std::to_string(it.iter->data.waitDuration)
				+ ")";
			valid = false;
			break;
		}

		if (it.iter->data.depTime != it.iter->data.startOfVisitTime + it.iter->data.visitDuration) {
			msg = "["+it.iter->data.id+"]depTime("
				+ std::to_string(it.iter->data.depTime)
				+ ") != "+it.iter->data.id+"startOfVisitTime("
				+ std::to_string(it.iter->data.startOfVisitTime)
				+ ") + ["+it.iter->data.id+"]visitDuration("
				+ std::to_string(it.iter->data.visitDuration)
				+ ")";
			valid = false;
			break;
		}

		if (it.iter->data.depTime > it.iter->data.timeWindow.closeTime) {
			msg = "["+it.iter->data.id+"]depTime("
				+ std::to_string(it.iter->data.depTime)
				+ ") > ["+it.iter->data.id+"]closeTime("
				+ std::to_string(it.iter->data.timeWindow.closeTime)
				+ ")";
			valid = false;
			break;
		}

		next = it + 1;
		if (next != walk.end()) {
			if (next.iter->data.arrTime != it.iter->data.depTime + travel_times[it.iter->data.depPointId][next.iter->data.arrPointId]) {
				msg = "["+it.iter->data.id+"]next->arrTime("
					+ std::to_string(next.iter->data.arrTime)
					+ ") != ["+it.iter->data.id+"]depTime("
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
