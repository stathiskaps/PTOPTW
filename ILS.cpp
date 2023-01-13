#include "ILS.h"

ILS::ILS() {}

ILS::ILS(int intervalsNum) : mIntervalsNum(intervalsNum) {
	metrics.local_search = std::vector<double>(intervalsNum, 0);
	metrics.split_unvisited = 0.0;
	metrics.shake = 0.0;
	metrics.final_pos = 0;
	metrics.middle_pos = 0;
}

ILS::~ILS() {}

std::map<std::string, std::vector<double>> ILS::getActivities(List<TA>& unvisited, std::vector<TimeWindow> intervals){
	std::map<std::string, std::vector<double>> activities;
	for(auto& n : unvisited) {
		std::vector<double> durations;
		for(auto& i : intervals) {
			double active_dur = std::min(i.closeTime, n.timeWindow.closeTime) - std::max(i.openTime, n.timeWindow.openTime);
			if (active_dur < 0) active_dur = 0;
			durations.push_back(active_dur/(n.timeWindow.closeTime-n.timeWindow.openTime));
		}
		activities.insert({n.id, durations});
	}
	return activities;
}

std::vector<TimeWindow> ILS::getIntervals(std::vector<TA> unvisited, int intervals_num, double day_start_time, double day_close_time) {

	std::vector<TimeWindow> intervals;

	if(intervals_num == 1) {
		intervals.push_back(TimeWindow{day_start_time, day_close_time});
		return intervals;
	}

	const size_t unvisited_size = unvisited.size();
	double diff = day_close_time - day_start_time;
	double dur = diff/intervals_num;

	double a = 0.2;

	double score, best_score = DBL_MAX;

	uint32_t not_improved = 0, max_times_not_improved = 20;

	std::vector<double> time_cuts, best_time_cuts;
	double time_cut = day_start_time;
	time_cuts.reserve(intervals_num+1);

	//Initialize time cuts
	time_cuts.push_back(day_start_time);
	for(size_t i = 0; i < intervals_num; ++i){
		time_cuts.push_back(time_cuts.back() + dur);
	}

	//Acceptance criterion
	while(true) {

		std::vector<Bin> bins;
		bins.reserve(intervals_num);
		for(uint32_t i = 0; i < intervals_num; ++i) {
			bins.push_back(Bin{i, i+1});
		}

		for(auto& ta: unvisited) {
			double max_activity = 0;
			std::vector<Bin>::iterator best_bin;
			for(std::vector<Bin>::iterator it = bins.begin(); it != bins.end(); ++it) {
				const double left_time_cut = time_cuts[it->left_cut_index];
				const double right_time_cut = time_cuts[it->right_cut_index];
				const double activity = std::min(right_time_cut, ta.timeWindow.closeTime) - std::max(left_time_cut, ta.timeWindow.openTime);

				if(activity > max_activity) {
					best_bin = it;
					max_activity = activity;
				}
			}
			best_bin->count++;
		}
		std::vector<Bin>::iterator most_used_bin = std::max_element(bins.begin(), bins.end(), [](Bin a, Bin b) {
			return a.count < b.count;
		});

		std::vector<Bin>::iterator least_used_bin= std::min_element(bins.begin(), bins.end(), [](Bin a, Bin b) {
			return a.count < b.count;
		});

		score = most_used_bin->count - least_used_bin->count;

		if(score < best_score) {
			best_score = score; //minimization
			best_time_cuts = time_cuts;
			not_improved = 0;
			
			//Improve intervals
			const double bin_duration = time_cuts[most_used_bin->right_cut_index] - time_cuts[most_used_bin->left_cut_index];

			double reduce_left = 0, reduce_right = 0;
			const double reduce_amount = bin_duration * a;

			//TODO:: check at the start of the function if there is only one bin
			if(most_used_bin > bins.begin() && most_used_bin < bins.end()){
				std::vector<Bin>::iterator left_bin = most_used_bin-1;
				std::vector<Bin>::iterator right_bin = most_used_bin+1;

				if(left_bin->count == right_bin->count){
					reduce_left = reduce_amount * 50/100;
					reduce_right = reduce_left;
				} else {
					if(left_bin->count < right_bin->count) {
						reduce_left = reduce_amount - (left_bin->count/right_bin->count)*reduce_amount;
						reduce_right = reduce_amount - reduce_left;
					} else {
						reduce_right = reduce_amount - (right_bin->count/left_bin->count)*reduce_amount;
						reduce_left = reduce_amount - reduce_right;
					}
				}
			} else {
				if(most_used_bin == bins.begin()){
					reduce_right = reduce_amount;
				} else {
					reduce_left = reduce_amount;
				}
			}

			time_cuts[most_used_bin->left_cut_index]+=reduce_left;
			time_cuts[most_used_bin->right_cut_index]-=reduce_right;

		} else {
			break;
		}
	}

	for(size_t i = 0; i < mIntervalsNum; ++i){
		intervals.push_back(TimeWindow{best_time_cuts[i], best_time_cuts[i+1]});
	}

	return intervals;
}

std::map<std::string, std::vector<ILS::Usage>> ILS::initRegistry(List<TA>& unvisited, std::vector<TimeWindow> intervals) {
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

void ILS::InitSolutions(std::vector<Solution>& sols, const std::vector<TimeWindow> intervals,  const OP& op){

	std::vector<Solution>::iterator first_sol = sols.begin(), last_sol = sols.end() - 1;

	for(auto& walk : first_sol->m_walks) {
		walk.push_front(*op.mStartDepot);
		walk.front().timeWindow = TimeWindow{intervals[0].openTime, intervals[0].closeTime};
		updateTimes(walk, walk.begin(), false, op.mTravelTimes, intervals[0]);
	}

	for(auto& walk : last_sol->m_walks) {
		walk.push_back(*op.mEndDepot);
		walk.back().timeWindow = TimeWindow{intervals[intervals.size()-1].openTime, intervals[intervals.size()-1].closeTime};
		updateTimes(walk, walk.begin(), false, op.mTravelTimes, intervals[intervals.size()-1]);
	}

}

ILS* g_CurrentInstance;

extern "C"
void drawCallback(){
	g_CurrentInstance->draw();
}

void ILS::setupDrawCallback(){
	::g_CurrentInstance = this;
	::glutDisplayFunc(::drawCallback);
}

void ILS::Solve(OP& op) {

	// std::cout.setstate(std::ios_base::failbit);
	
	// myInit();
	// // Set the OpenGL display mode
	// glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
	// // Set the initial window size
	// glutInitWindowSize(WINDOW_WIDTH, WINDOW_HEIGHT);
	// // Create the window with the given title
	// glutCreateWindow("Best Solution");

	// //Set the display callback function
	// setupDrawCallback();
	// glutReshapeFunc(onResize);

	// // Enter the GLUT main loop
	// std::thread glutThread(glutMainLoop);

	auto start = std::chrono::high_resolution_clock::now();
	const int num_locations = op.mAttractions.size();
	const int max_to_remove = num_locations / (3 * op.m_walks_num);

	std::vector<TA> unvisitedVec;
	List<TA> unvisited;
	for (auto ta : op.mAttractions) {
		unvisited.push_back(*ta);
		unvisitedVec.push_back(*ta);
	}

	std::vector<TimeWindow> intervals = getIntervals(unvisitedVec, mIntervalsNum, op.mStartDepot->timeWindow.openTime, op.mEndDepot->timeWindow.closeTime);
	const std::map<std::string, std::vector<double>> activities = getActivities(unvisited, intervals);
	std::map<std::string, std::vector<ILS::Usage>> reg = initRegistry(unvisited, intervals); 

	int counter{}, times_not_improved{ 0 }, best_score{ INT_MIN };

	std::vector<ILS::SR> shake_settings;
	for(size_t i = 0; i < mIntervalsNum; ++i){
		shake_settings.push_back(ILS::SR{1, 1});
	}

	std::vector<Solution> proc_solutions(mIntervalsNum, Solution());
	//Initialize solutions
	for (auto& s : proc_solutions) {
		for (size_t i = 0; i < op.m_walks_num; ++i) {
			s.m_walks.push_back(List<TA>());
		}
	}
	InitSolutions(proc_solutions, intervals, op);
	List<TA> pool = std::move(unvisited);

 
	while (times_not_improved < MAX_TIMES_NOT_IMPROVED) {
		counter++;

		std::cout << "Revision " << counter << std::endl;
		
		splitUnvisitedList(proc_solutions, pool, mIntervalsNum, reg, activities);
		SplitSearch(proc_solutions, intervals, op, reg);
		int score = collectScores(proc_solutions);

		// for(size_t i = 0; i < proc_solutions.size(); ++i){
		// 	proc_solutions[i].print("Solution "+std::to_string(i));
		// }

		if (score > best_score) {
			best_score = score;
			best_solutions = proc_solutions;
			for(auto& sr : shake_settings){
				sr.R = 1;
			}
			times_not_improved = 0;
			// glutPostRedisplay();
		}
		else {
			times_not_improved++;
		}

		int removed_counter = SplitShake(proc_solutions, shake_settings, op, max_to_remove, intervals);
		gatherUnvisited(proc_solutions, pool);

	}
	auto end = std::chrono::high_resolution_clock::now();
	auto elapsed_time = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();

	std::cout.clear();
	best_solution = connectSolutions(best_solutions, op.m_walks_num);
	for(auto walk_it = best_solution.m_walks.begin(); walk_it != best_solution.m_walks.end(); ++walk_it){
		const size_t index = walk_it - best_solution.m_walks.begin();
		updateTimes(*walk_it, walk_it->begin(), false, op.mTravelTimes, op.mTimeWindow);
	}
	validate(best_solution.m_walks, op.mTravelTimes, false);
	std::cout << "Best score: " << best_score << std::endl;
	std::cout << "Visits: " << best_solution.getVisits() << std::endl;

	best_solution.print("Best solution");

	std::cout << std::endl;
	std::cout << "-------------Metrics-----------------" << std::endl;
	std::cout << "Split Local Search: " << std::endl;
	for(size_t i = 0; i < metrics.local_search.size(); ++i){
		std::cout << "Solution " << i << ": " << metrics.local_search[i] << " seconds" << std::endl;
	}

	std::cout << std::endl;
	std::cout << "Split Unvisited List: " << metrics.split_unvisited << " seconds\n" << std::endl;

	std::cout << "Shake: " << metrics.shake << " seconds\n" << std::endl;

	std::cout << "Total Times middle position was best: " << metrics.middle_pos << "\n" << std::endl;
	std::cout << "Total Times final position was best: " << metrics.final_pos << "\n" << std::endl;

	std::cout << "Total execution time: " << elapsed_time << " seconds" << std::endl;
	std::cout << "-------------------------------------" << std::endl;
	std::cout << std::endl;

	std::cout << best_score << "\t" << elapsed_time << "\t" << best_solution.getVisits() << std::endl;
	std::cout << std::endl;

	// Wait for the GLUT thread to finish
	// glutMainLoop();
    // glutThread.join();
}

void ILS::gatherUnvisited(std::vector<Solution>& solutions, List<TA>& pool){
	pool.clear();
	for(auto& sol : solutions){
		pool.append(sol.m_unvisited);
		sol.m_unvisited.clear();
	}
	std::cout << std::endl;
}

std::vector<List<TA>> ILS::splitUnvisitedList(std::vector<Solution>& solutions, List<TA>& pool, int intervals_num, std::map<std::string, std::vector<ILS::Usage>>& reg, std::map<std::string, std::vector<double>> activities) {

	auto start = std::chrono::high_resolution_clock::now();
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

	auto end = std::chrono::high_resolution_clock::now();
	auto elapsed_time = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
	metrics.split_unvisited += elapsed_time;

	return lists;
}

Solution ILS::connectSolutions(std::vector<Solution>& sols, const size_t walks_num) {
	Solution solution;
	//Connect unvisited

	gatherUnvisited(sols, solution.m_unvisited);
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
		if(sol_it == sols.begin()){ //first sol
			for(std::vector<List<TA>>::iterator walk_it = sol_it->m_walks.begin(); walk_it != sol_it->m_walks.end(); ++walk_it){
				walk_it->push_back(walk_it->back());
				walk_it->back().id = DUMMY_ID;;
			}
		} 
		else if (sol_it == sols.end() - 1){ //last sol
			for(std::vector<List<TA>>::iterator walk_it = sol_it->m_walks.begin(); walk_it != sol_it->m_walks.end(); ++walk_it){
				walk_it->push_front(walk_it->front());
				walk_it->front().id = DUMMY_ID;
			}
		} 
		else { //middle sol
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

int ILS::SplitShake(std::vector<Solution>& sols, std::vector<ILS::SR>& shake_settings, OP& op, const int& max_to_remove, const std::vector<TimeWindow> intervals){
	auto start = std::chrono::high_resolution_clock::now();
	int removed_counter = 0;
	PrepareForShake(sols);
	for(auto sol_it = sols.begin(); sol_it != sols.end(); ++sol_it){
		const int sol_index = sol_it - sols.begin();
		int removed = Shake(*sol_it, shake_settings[sol_index].S, shake_settings[sol_index].R, op, max_to_remove, intervals[sol_index]);
		removed_counter += removed;
	}
	RemoveDummyNodes(sols);
	auto end = std::chrono::high_resolution_clock::now();
	auto elapsed_time = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
	metrics.shake += elapsed_time;
	return removed_counter;
}

int ILS::Shake(Solution& sol, int& S, int& R, OP& op, const int& max_to_remove, const TimeWindow time_budget) {
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

	

	for (auto walk_it = sol.m_walks.begin(); walk_it != sol.m_walks.end(); ++walk_it) {
		const size_t index = walk_it - sol.m_walks.begin();
		if(walk_it->size() == 2) continue;
		int start_pos{ S }, end_pos{ std::min(S + R, (int)walk_it->size() - 1) };
		removed_counter += end_pos - start_pos;
		List<TA>::iterator start_it{ walk_it->begin() + start_pos }, end_it{ walk_it->begin() + end_pos };
		List<TA> part(start_it, end_it);
		List<TA>::iterator next = walk_it->erase(start_it, end_it);
		sol.m_unvisited.append(part);
		part.clear();
		updateTimes(*walk_it, next, false, op.mTravelTimes, time_budget);
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

Point ILS::getWeightedCentroid(const List<TA>::iterator& first, const List<TA>::iterator& last, const int pointId) {

	int total_profit{ collectProfit(first, last) };
	Position pos;
	double profit_ratio{};
	for (auto it = first; it != last; ++it) {
		profit_ratio = it.iter->data.profit / (double)total_profit;
		pos.lat += it.iter->data.point.pos.lat * profit_ratio;
		pos.lon += it.iter->data.point.pos.lon * profit_ratio;
	}
	return Point{ pointId, pos };
	
}

bool ILS::hasWeightedCentroid(const Solution& sol, const int walk_index, const int min_walk_size) {
	return sol.m_walks[walk_index].size() >= min_walk_size || !sol.m_unvisited.empty();
}

std::tuple<bool, double> ILS::CandidateStartDepotIsValid(const List<TA>& walk, const TA& ta, const double start_time, const Vector2D<double>& travel_times){
	auto first_it = walk.begin();
	double wait_dur = std::max(0.0, ta.timeWindow.openTime - start_time);
	double start_of_visit_time = start_time + wait_dur;
	double dep_time = start_of_visit_time + ta.visitDuration;
	double shift = wait_dur
		+ ta.visitDuration
		+ travel_times[ta.point.id][first_it.iter->data.arrPointId];

	return { dep_time <= ta.timeWindow.closeTime && shift <= first_it.iter->data.maxShift, shift };
}

void ILS::AddStartDepots(std::vector<Solution>& solutions, const std::vector<TimeWindow>& intervals, const int i, const OP& op){

	TimeWindow timeBudget = TimeWindow{intervals[i].openTime, intervals[i].closeTime};
	for (size_t j = 0; j < solutions[i].m_walks.size(); ++j) {
		int prev_sol_index = i-1;
		
		if(first_solution){ //we have already added a start depot at first solution's walks
			return;
		}

		while(solutions[prev_sol_index].m_walks[j].empty() && prev_sol_index >= 0){
			prev_sol_index--;
			continue;
		}

		if(prev_sol_index == -1){
			throw std::runtime_error("didn't find a valid previous walk to get a startDepot");
		}

		TA &prev_last = solutions[prev_sol_index].m_walks[j].back();
		TA candidate = prev_last;
		candidate.neutralize(DUMMY_START_DEPOT);
		candidate.timeWindow = intervals[i];

		//left trim invalid nodes
		if(!solutions[i].m_walks[j].empty()) { 
			List<TA>::iterator curr = solutions[i].m_walks[j].begin();
			while(curr != solutions[i].m_walks[j].end()){
				auto [valid, _] = CandidateStartDepotIsValid(solutions[i].m_walks[j], candidate, intervals[i].openTime, op.mTravelTimes);
				if(!valid){
					solutions[i].m_unvisited.push_back(curr.iter->data);
					curr = solutions[i].m_walks[j].erase(curr);
					updateTimes(solutions[i].m_walks[j], curr, false, op.mTravelTimes, intervals[i]);
				} else {
					break;
				}
			}
		}

		solutions[i].m_walks[j].push_front(candidate);
		updateTimes(solutions[i].m_walks[j], solutions[i].m_walks[j].begin(), false, op.mTravelTimes, intervals[i]);

	}
}

std::tuple<bool, double> ILS::CandidateEndDepotIsValid(const List<TA>& walk, const TA candidateEndDepot, const TimeWindow timebudget){
	double shift = candidateEndDepot.waitDuration + candidateEndDepot.visitDuration;
	auto last_it = walk.end()-1;
	if(walk.empty()){
		return { timebudget.openTime + shift <= timebudget.closeTime, shift };
	}
	const double distance = last_it.iter->data.point.euclidean_distance(candidateEndDepot.point);
	shift += distance;
	return { last_it.iter->data.depTime + shift <= timebudget.closeTime, shift };
}

void ILS::AddEndDepots(std::vector<Solution>& solutions, const std::vector<TimeWindow>& intervals, const int i, OP& op){
	const int min_size = 3;
	TimeWindow timeBudget = TimeWindow{intervals[i].openTime, intervals[i].closeTime};
	for (size_t j = 0; j < solutions[i].m_walks.size(); ++j) {
		//add endpoint
		Point candidateEndPoint, finalEndPoint;
		TA last = solutions[i].m_walks[j].back(), endDepot;
		std::string endPointId;

		int next_sol_index = i + 1;
		while (next_sol_index < solutions.size() && !hasWeightedCentroid(solutions[next_sol_index], j, min_size)) {
			next_sol_index++;
		}

		if (last_solution || next_sol_index == solutions.size()) {
			if(solutions[i].m_walks[j].back().id != op.mEndDepot->id){
				candidateEndPoint = op.mEndDepot->point;
			}
		}
		else {
			if (solutions[next_sol_index].m_walks[j].size() > min_size) {
				const List<TA>& next_solution_walk = solutions[next_sol_index].m_walks[j];
				candidateEndPoint = getWeightedCentroid(next_solution_walk.at(0), next_solution_walk.at(min_size), DEFAULT_POINT_ID);
			}
			else {
				candidateEndPoint = getWeightedCentroid(solutions[next_sol_index].m_unvisited.begin(), solutions[next_sol_index].m_unvisited.end(), DEFAULT_POINT_ID);
			}
		}
		
		//We should check if candidateEndPoint can be added at the end of the walk as endDepot
		//If we use function insertionAfterIsValid, we need first to add the point to the graph and calculate
		//the travel times for and to each other point. But we only need to know the distance between walk's last point
		//and candidateEndPoint to figure out if the insertion is feasible. So we won't add the candidate end point in the graph.
		const TA candidateEndDepot = TA(CANDIDATE_END_DEPOT_ID, candidateEndPoint);
		auto [valid, _] = CandidateEndDepotIsValid(solutions[i].m_walks[j], candidateEndDepot, timeBudget);
		if (!valid){
			finalEndPoint = last.point.findPointWithDistance(candidateEndPoint, intervals[i].closeTime - last.depTime - 1);
		} else {
			finalEndPoint = candidateEndPoint;
		}

		op.AddPointToGraph(finalEndPoint); //finalEndPoint gets an id here
		endDepot = TA(DUMMY_END_DEPOT_ID, finalEndPoint); //endDepot arrPointId and dePointId get the finalPoint.id value
		endDepot.timeWindow = timeBudget;
		solutions[i].m_walks[j].push_back(endDepot);
		updateTimes(solutions[i].m_walks[j], solutions[i].m_walks[j].end() - 1, false, op.mTravelTimes, intervals[i]);

	}
}

std::vector<Point> ILS::getTargets(const std::vector<Solution>& solutions, const int i, const OP& op){
	if(last_solution){
		return std::vector<Point>(solutions[i].m_walks.size(), Point()); //DEFAULT_POINT_ID
	}

	const int min_size = 3;
	std::vector<Point> targets;
	for(size_t j = 0; j < solutions[i].m_walks.size(); ++j){
		int next_sol_index = i + 1;
		while (next_sol_index < solutions.size() && !hasWeightedCentroid(solutions[next_sol_index], j, min_size)) {
			next_sol_index++;
		}

		if(next_sol_index == solutions.size()){
			Point p = op.mEndDepot->point;
			p.id = TARGET_POINT_ID;
			targets.push_back(p);
			continue;
		}

		Point target;
		const List<TA>& next_solution_walk = solutions[next_sol_index].m_walks[j];
		if (next_solution_walk.size() > min_size) {
			const List<TA>& next_solution_walk = solutions[next_sol_index].m_walks[j];
			target = getWeightedCentroid(next_solution_walk.at(0), next_solution_walk.at(min_size), TARGET_POINT_ID);
		} else {
			target = getWeightedCentroid(solutions[next_sol_index].m_unvisited.begin(), solutions[next_sol_index].m_unvisited.end(), TARGET_POINT_ID);
		}
		targets.push_back(target);
	}
	return targets;

}

void ILS::SplitSearch(std::vector<Solution>& solutions, const std::vector<TimeWindow>& intervals, OP& op, std::map<std::string, std::vector<ILS::Usage>>& reg) {

	const int min_size = 3;
	const int solutions_size = solutions.size();

	for (int i = 0; i < solutions.size(); ++i) {

		auto start = std::chrono::high_resolution_clock::now();

		if (solutions[i].m_unvisited.empty()) continue;

		if(i > 0) {
			AddStartDepots(solutions, intervals, i, op);
		}

		

		std::vector<Point> targets = getTargets(solutions, i, op);
		std::vector<std::string> ids = construct(solutions[i], op.mTravelTimes, targets, intervals[i], !(last_solution));

		for (auto& id : ids) {
			reg[id].at(i).solved++;
		}

		if (i > 0) {

			for(std::vector<List<TA>>::iterator walk_it = solutions[i].m_walks.begin(); walk_it != solutions[i].m_walks.end(); ++walk_it){
				walk_it->pop_front();
				updateTimes(*walk_it, walk_it->begin(), false, op.mTravelTimes, intervals[i]);
			}

			//Make a local search betweem current solution and previous
			for (size_t j = 0; j < solutions[i].m_walks.size(); ++j) {
				if(solutions[i-1].m_walks[j].size() > 2 && solutions[i].m_walks[j].size() > 2) {
					Solution sol;
					Walk walk;
					TA& first = solutions[i-1].m_walks[j].back();
					TA& last = solutions[i].m_walks[j].front();
					TimeWindow interval = {first.depTime, last.depTime + last.maxShift};
					walk.push_front(first);
					walk.push_back(last);
					updateTimes(walk, walk.begin(), false, op.mTravelTimes, interval);
					sol.m_walks.push_back(walk);

					sol.m_unvisited.append(solutions[i-1].m_unvisited);
					sol.m_unvisited.append(solutions[i].m_unvisited);
					
					std::vector<Point> fake_targets = std::vector<Point>(solutions[i].m_walks.size(), Point());

					std::vector<std::string> inserted_ids = construct(sol, op.mTravelTimes, fake_targets, interval, false);

					bool changed_previous = false, changed_current = false;

					List<TA>::iterator it;
					for(it = sol.m_walks[j].begin()+1; it != sol.m_walks[j].end()-1; ++it) {
						if(it.iter->data.depTime > intervals[i].openTime){
							break;
						}
					}
					std::pair<List<TA>, List<TA>> walks = sol.m_walks[j].split(it);
					
					const TA& temp = walks.second.front();
					if(!walks.second.empty() && temp.arrTime < intervals[i].openTime){
						const double dep_time = temp.dep_time(intervals[i].openTime);
						if(temp.maxShift < intervals[i].openTime - temp.arrTime || dep_time > temp.timeWindow.closeTime){
							walks.second.trim_left(1);
						} else {
							std::cout << "Will shift walk by " << intervals[i].openTime - temp.arrTime << " time units to the right" << std::endl;
						}
					}

					for(List<TA>::iterator it = walks.first.begin() + 1; it != walks.first.end(); ++it){
						solutions[i-1].m_walks[j].push_back(*it);
						if(!solutions[i-1].m_unvisited.eraseId(it.iter->data.id) && !solutions[i].m_unvisited.eraseId(it.iter->data.id)) {
							throw std::runtime_error("unvisited node "+it.iter->data.id+" wasn't found");
						}
						reg[it.iter->data.id].at(i-1).solved++;
						changed_previous = true;
					}

					for(List<TA>::iterator it = walks.second.end() - 2; it != walks.second.end(); --it){
						solutions[i].m_walks[j].push_front(*it);
						if(!solutions[i-1].m_unvisited.eraseId(it.iter->data.id) && !solutions[i].m_unvisited.eraseId(it.iter->data.id)) {
							throw std::runtime_error("unvisited node "+it.iter->data.id+" wasn't found");
						}
						reg[it.iter->data.id].at(i).solved++;
						changed_current = true;
					}

					if(changed_previous){
						updateTimes(solutions[i-1].m_walks[j], solutions[i-1].m_walks[j].begin(), false, op.mTravelTimes, intervals[i-1]);
					}

					if(changed_current) {
						updateTimes(solutions[i].m_walks[j], solutions[i].m_walks[j].begin(), false, op.mTravelTimes, intervals[i]);
					}

					// std::cout << std::endl;

				}
			}
		}

		validate(solutions[i].m_walks, op.mTravelTimes, false);

		auto end = std::chrono::high_resolution_clock::now();
		auto elapsed_time = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
		metrics.local_search[i] += elapsed_time;

	}
}

void ILS::drawSolutions(const std::vector<Solution>& solutions){
	Bounds bounds;
	const double radius = 2;
	const std::string filename = "Solutions.svg";
	std::vector<std::string> colors = {"black", "steelblue", "firebrick", "darkgray", "midnightblue"};
	int sol_counter{};
	for(auto& sol : solutions){
		for(List<TA>::iterator ta_it = sol.m_unvisited.begin(); ta_it != sol.m_unvisited.end(); ++ta_it){
			if(ta_it.iter->data.point.pos.lat < bounds.minLat){
				bounds.minLat = ta_it.iter->data.point.pos.lat;
			}
			if(ta_it.iter->data.point.pos.lat > bounds.maxLat){
				bounds.maxLat = ta_it.iter->data.point.pos.lat;
			}
			if(ta_it.iter->data.point.pos.lon < bounds.minLon){
				bounds.minLon = ta_it.iter->data.point.pos.lon;
			}
			if(ta_it.iter->data.point.pos.lon > bounds.maxLon){
				bounds.maxLon = ta_it.iter->data.point.pos.lon;
			}
		}

		for(auto& walk : sol.m_walks){
			for(List<TA>::iterator ta_it = walk.begin(); ta_it != walk.end(); ++ta_it){
				if(ta_it.iter->data.point.pos.lat < bounds.minLat){
					bounds.minLat = ta_it.iter->data.point.pos.lat;
				}
				if(ta_it.iter->data.point.pos.lat > bounds.maxLat){
					bounds.maxLat = ta_it.iter->data.point.pos.lat;
				}
				if(ta_it.iter->data.point.pos.lon < bounds.minLon){
					bounds.minLon = ta_it.iter->data.point.pos.lon;
				}
				if(ta_it.iter->data.point.pos.lon > bounds.maxLon){
					bounds.maxLon = ta_it.iter->data.point.pos.lon;
				}
			}
		}
	}

	// Open an output stream to write the SVG file
    std::ofstream out(filename);

	// Write the SVG file header
    out << "<svg viewBox=\"" << bounds.minLat - radius << " " << bounds.minLon - radius << " " <<
	bounds.maxLat - bounds.minLat + radius*2 << " " << bounds.maxLon - bounds.minLon + radius*2 <<
	"\" xmlns=\"http://www.w3.org/2000/svg\">" << std::endl;

	// Write the routes as lines in the SVG file
	for(auto sol_it = solutions.begin(); sol_it != solutions.end(); ++sol_it){
		const int index = sol_it - solutions.begin();
		for(auto& walk : sol_it->m_walks){
			List<TA>::iterator left, right;
			for(left = walk.begin(), right = walk.begin() + 1; right != walk.end(); ++left, ++right){
				out << "<line x1=\"" << left.iter->data.point.pos.lat << "\" y1=\"" << left.iter->data.point.pos.lon <<
				"\" x2=\"" << right.iter->data.point.pos.lat << "\" y2=\"" << right.iter->data.point.pos.lon << "\" stroke=\" "<< colors[index] << "\" stroke-width=\"0.5\" />" << std::endl;
			}
		}
	}

	//Write nodes
	for(auto sol_it = solutions.begin(); sol_it != solutions.end(); ++sol_it){
		const int index = sol_it - solutions.begin();
		for(List<TA>::iterator ta_it = sol_it->m_unvisited.begin(); ta_it != sol_it->m_unvisited.end(); ++ta_it){
        	out << "<circle cx=\"" << ta_it.iter->data.point.pos.lat << "\" cy=\"" << ta_it.iter->data.point.pos.lon  << "\" fill=\"" << colors[index] << "\" r=\"" << radius << "\" />" << std::endl;
		}
		for(auto& walk : sol_it->m_walks){
			for(List<TA>::iterator ta_it = walk.begin(); ta_it != walk.end(); ++ta_it){
				out << "<circle cx=\"" << ta_it.iter->data.point.pos.lat << "\" cy=\"" << ta_it.iter->data.point.pos.lon  << "\" fill=\"" << colors[index] << "\" r=\"" << radius << "\" />" << std::endl;
			}
		}
	}

	//Write texts
	for(auto sol_it = solutions.begin(); sol_it != solutions.end(); ++sol_it){
		const int index = sol_it - solutions.begin();
		for(List<TA>::iterator ta_it = sol_it->m_unvisited.begin(); ta_it != sol_it->m_unvisited.end(); ++ta_it){
        	out << "<text x=\""<< ta_it.iter->data.point.pos.lat << "\" y=\""<< ta_it.iter->data.point.pos.lon << "\" text-anchor=\"middle\" font-size=\"2px\" fill=\"white\" alignment-baseline=\"middle\">" << ta_it.iter->data.point.id  << "</text>" << std::endl;
    	}
		for(auto& walk : sol_it->m_walks){
			for(List<TA>::iterator ta_it = walk.begin(); ta_it != walk.end(); ++ta_it){
				out << "<text x=\""<< ta_it.iter->data.point.pos.lat << "\" y=\""<< ta_it.iter->data.point.pos.lon << "\" text-anchor=\"middle\" font-size=\"2px\" fill=\"white\" alignment-baseline=\"middle\">" << ta_it.iter->data.point.id  << "</text>" << std::endl;
			}
		}
	}

	// Write the SVG file footer
    out << "</svg>" << std::endl;

    // Close the output stream
    out.close();
}


std::vector<std::string> ILS::construct(Solution& sol, const Vector2D<double>& travel_times, const std::vector<Point> targets, 
	const TimeWindow time_budget, const bool open) {
	
	std::vector<std::string> inserted_ids;
	if (sol.m_unvisited.empty()) {
		return inserted_ids;
	}

	double min_shift{}, max_score{ DBL_MIN }, score{};
	int best_arr_point_id{ DEFAULT_POINT_ID }, best_dep_point_id{ DEFAULT_POINT_ID }, 
		arr_point_id{ DEFAULT_POINT_ID }, dep_point_id{ DEFAULT_POINT_ID };
	Walks::iterator best_walk_it;
	List<TA>::iterator pos{ }, best_pos{ };

	List<TA>::iterator insert_it{ sol.m_unvisited.end() }, curr{}, inserted_it{};

	//calculate the distances from all unvisited nodes to all target points
	std::map<int, std::vector<double>> distances;
	for(List<TA>::iterator ta_it = sol.m_unvisited.begin(); ta_it != sol.m_unvisited.end(); ++ta_it){
		distances[ta_it.iter->data.point.id] = std::vector<double>();
		distances[ta_it.iter->data.point.id].reserve(targets.size());
		for(const auto &t : targets){
			if(t.id == DEFAULT_POINT_ID){
				distances[ta_it.iter->data.point.id].push_back(0);
			} else {
				distances[ta_it.iter->data.point.id].push_back(ta_it.iter->data.point.euclidean_distance(t));
			}
			
		}
	}

	while (true) {
		max_score = DBL_MIN;
		insert_it = sol.m_unvisited.end();
		curr = sol.m_unvisited.begin();

		while (curr != sol.m_unvisited.end()) {
			min_shift = DBL_MAX;
			std::vector<double> dists = distances[curr.iter->data.point.id];
			auto [walk_it, pos, min_shift, arr_point_id, dep_point_id] = getBestPos(curr.iter->data, sol.m_walks, travel_times, dists, time_budget, open);
			//default value of walk_it is sol.m_walks.end()
			//default value of pos is sol.m_walks.front().end()
			//if these iterators have these values then insertion is infeasible

			if (walk_it == sol.m_walks.end()) {
				curr++;
				continue;
			}

			const int walk_index = walk_it - sol.m_walks.begin();
			const double distance_to_target = distances[curr.iter->data.depPointId].at(walk_index);
			score = pow(curr.iter->data.profit, 2) / min_shift;	//score = profit^2/ minShifti

			//check each score 
			if (max_score < score){ 
				max_score = score;
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
		updateTimes(*best_walk_it, inserted_it, true, travel_times, time_budget);
	}
	return inserted_ids;
}

std::tuple<bool, double> ILS::insertionBetweenIsValid(const TA& left, const TA& ta, const TA& right, 
	const Vector2D<double>& travel_times){
	double shift = travel_times[left.depPointId][ta.point.id]
		+ ta.visitDuration
		+ travel_times[ta.point.id][right.arrPointId]
		- travel_times[left.depPointId][right.arrPointId];
	return { shift <= right.maxShift, shift };
}

//Check if inserting ta before node K is valid
std::tuple<bool, double> ILS::insertionBeforeIsValid(const TA& ta, const TA& right, const double start_time, const Vector2D<double>& travel_times){
	double shift = ta.visitDuration
		+ travel_times[ta.point.id][right.arrPointId];
	return { shift <= right.maxShift, shift };
}

//Check if inserting ta after node K is valid. Close time of route is needed.
std::tuple<bool, double> ILS::insertionAfterIsValid(const TA& left, const TA& ta, const double close_time, const Vector2D<double>& travel_times){
	double shift = travel_times[left.depPointId][ta.point.id]
		+ ta.visitDuration;

	double arr_time = left.depTime + travel_times[left.depPointId][ta.point.id];
	double dep_time = arr_time + ta.visitDuration;

	return { dep_time < close_time, shift };
}

std::tuple<Walks::iterator, List<TA>::iterator, double, int, int> ILS::getBestPos(const TA& ta, Walks& walks, 
	const Vector2D<double>& travel_times, const std::vector<double> distances, const TimeWindow time_budget, const bool open) {

	Walks::iterator best_walk{ walks.end() };
	List<TA>::iterator best_pos{ walks.front().end() };
	int arr_point_id{ DEFAULT_POINT_ID }, dep_point_id{ DEFAULT_POINT_ID };
	double min_shift{ DBL_MAX };


	for (Walks::iterator walk_it = walks.begin(); walk_it != walks.end(); ++walk_it) {
		if (walk_it->empty()) {
			std::cerr << "walk is empty, minimum size: 1" << std::endl;
			std::exit(1);
		}

		double shift{};
		List<TA>::iterator left{ walk_it->begin() }, right{ walk_it->begin() + 1 };
		int temp_arr_point_id{ DEFAULT_POINT_ID }, temp_dep_point_id{ DEFAULT_POINT_ID };

		while (right != walk_it->end()) {
			auto [valid, shift] = insertionBetweenIsValid(left.iter->data, ta, right.iter->data, travel_times);
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

void ILS::updateTimes(List<TA>& walk, const List<TA>::iterator& start_pos, 
	const bool smart, const Vector2D<double>& travel_times, const TimeWindow time_budget) {
	//update times first
	List<TA>::iterator it{ start_pos };

	while (it != walk.end()) {
		it.iter->data.arrTime = it.iter->previous->data.depTime + travel_times[it.iter->previous->data.depPointId][it.iter->data.arrPointId];
		it.iter->data.startOfVisitTime = it.iter->data.arrTime;
		it.iter->data.depTime = it.iter->data.startOfVisitTime + it.iter->data.visitDuration;
		++it;
	}

	updateMaxShifts(walk, travel_times, time_budget);
}

void ILS::updateMaxShifts(const List<TA>& li, const Vector2D<double>& travel_times, const TimeWindow time_budget) {
	List<TA>::iterator it{ li.end() - 1 };
	it.iter->data.maxShift = std::min(it.iter->data.timeWindow.closeTime - it.iter->data.depTime, time_budget.closeTime - it.iter->data.depTime);
	if(it.iter->data.maxShift < 0) {
		throw std::runtime_error("invalid maxshift number");
	}
	it--;
	while(it != li.end()){
		it.iter->data.maxShift = (it + 1).iter->data.maxShift;
		it--;
	}
}

void ILS::validate(const std::vector<Solution>& sols, const Vector2D<double>& travel_times, const bool verbose){
	for(auto sol_it = sols.begin(); sol_it != sols.end(); ++sol_it){
		const int index = sol_it - sols.begin();
		std::cout << "Checking solution [" << index << "]: ";
		validate(sol_it->m_walks, travel_times, verbose);
	}
}

void ILS::validate(const std::vector<List<TA>>& walks, const Vector2D<double>& travel_times, const bool verbose) {
	for (auto& walk : walks) {
		validate(walk, travel_times, verbose);
	}
}

void ILS::validate(const List<TA>& walk, const Vector2D<double>& travel_times, const bool verbose) {
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

	if(verbose) {
		std::cout << (valid ? "walk is valid" : msg) << std::endl;
	}
	
	if (!valid) {
		std::cerr << "walk is invalid: exiting.. " << std::endl; 
		std::exit(1);
	}
}

std::tuple<bool, double> ILS_TOPTW::insertionBetweenIsValid(const TA& left, const TA& ta, const TA& right, 
	const Vector2D<double>& travel_times){

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

//Check if inserting ta before route is valid. The start_time of the route is needed.
std::tuple<bool, double> ILS_TOPTW::insertionBeforeIsValid(const TA& ta, const TA& right, const double start_time, const Vector2D<double>& travel_times){
	double wait_dur = std::max(0.0, ta.timeWindow.openTime - start_time);
	double start_of_visit_time = start_time + wait_dur;
	double dep_time = start_of_visit_time + ta.visitDuration;
	double shift = wait_dur
		+ ta.visitDuration
		+ travel_times[ta.point.id][right.arrPointId];

	return { dep_time <= ta.timeWindow.closeTime && shift <= right.maxShift, shift };
}

//Check if inserting ta after left is valid. The close time of route is needed.
std::tuple<bool, double> ILS_TOPTW::insertionAfterIsValid(const TA& left, const TA& ta, const double close_time, const Vector2D<double>& travel_times){
	double arr_time = left.depTime + travel_times[left.depPointId][ta.point.id];
	double wait_duration = std::max(0.0, ta.timeWindow.openTime - arr_time);
	double start_of_visit_time = arr_time + wait_duration;
	double dep_time = start_of_visit_time + ta.visitDuration;

	double shift = travel_times[left.depPointId][ta.point.id]
	+ ta.waitDuration
	+ ta.visitDuration;

	return { dep_time <= close_time, shift };
}


std::tuple<Walks::iterator, List<TA>::iterator, double, int, int> ILS_TOPTW::getBestPos(const TA& ta, 
	Walks& walks, const Vector2D<double>& travel_times, const std::vector<double> distances, const TimeWindow time_budget, const bool open) {

	Walks::iterator best_walk{ walks.end() };
	List<TA>::iterator best_pos{ walks.front().end()};
	int arr_point_id{ DEFAULT_POINT_ID }, dep_point_id{ DEFAULT_POINT_ID };
	double min_shift{ DBL_MAX };

	for (Walks::iterator walk_it = walks.begin(); walk_it != walks.end(); ++walk_it) {
		if (walk_it->empty()) {
			std::cerr << "walk is empty, minimum size: 1" << std::endl;
			std::exit(1);
		}

		const size_t walk_index = walk_it - walks.begin();
		const double distance_to_next_target = distances[walk_index];

		double arr_time{}, wait_dur{}, start_of_visit_time{}, dep_time{}, shift{};

		const size_t possible_positions { walk_it->size() };
		int temp_arr_point_id{ DEFAULT_POINT_ID }, temp_dep_point_id{ DEFAULT_POINT_ID }, pos_counter{};
		int insertion_pos_index = 1;

		List<TA>::iterator left{ walk_it->begin() }, right{ left + 1 };
		while(right != walk_it->end()) {
			arr_time = left.iter->data.depTime + travel_times[left.iter->data.depPointId][ta.point.id];
			wait_dur = std::max(0.0, ta.timeWindow.openTime - arr_time);
			start_of_visit_time = arr_time + wait_dur;
			dep_time = start_of_visit_time + ta.visitDuration;
			shift = travel_times[left.iter->data.depPointId][ta.point.id]
				+ wait_dur
				+ ta.visitDuration
				+ travel_times[ta.point.id][right.iter->data.arrPointId]
				- travel_times[left.iter->data.depPointId][right.iter->data.arrPointId]
				+ distance_to_next_target * (insertion_pos_index/possible_positions);

			if(dep_time <= ta.timeWindow.closeTime && shift <= right.iter->data.maxShift) {
				if (shift < min_shift) {
					best_walk = walk_it; 
					best_pos = right;
					min_shift = shift;
					arr_point_id = ta.point.id;
					dep_point_id = ta.point.id;
				}
			}
			left++; right++; insertion_pos_index++;
		}

		if(open){
			//here we check the insertion at the end of the walk
			arr_time = left.iter->data.depTime + travel_times[left.iter->data.depPointId][ta.point.id];
			wait_dur = std::max(0.0, ta.timeWindow.openTime - arr_time);
			start_of_visit_time = arr_time + wait_dur;
			dep_time = start_of_visit_time + ta.visitDuration;
			shift = travel_times[left.iter->data.depPointId][ta.point.id]
				+ wait_dur
				+ ta.visitDuration
				+ distance_to_next_target * (insertion_pos_index/possible_positions);


			if(dep_time <= ta.timeWindow.closeTime && dep_time <= time_budget.closeTime) {
				if (shift < min_shift) {
					metrics.final_pos++;
					best_walk = walk_it; 
					best_pos = right;
					min_shift = shift;
					arr_point_id = ta.point.id;
					dep_point_id = ta.point.id;
				} else {
					metrics.middle_pos++;
				}
			}
		}
		
	}

	return { best_walk, best_pos, min_shift, arr_point_id, dep_point_id };
}

void ILS_TOPTW::updateTimes(List<TA>& walk, const List<TA>::iterator& start_pos, 
	const bool smart, const Vector2D<double>& travel_times, const TimeWindow time_budget) {
	//update times first
	List<TA>::iterator it{ start_pos };
	if (it == walk.begin()) {
		it.iter->data.arrTime = time_budget.openTime;
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

	updateMaxShifts(walk, travel_times, time_budget);
}

void ILS_TOPTW::updateMaxShifts(const List<TA>& li, const Vector2D<double>& travel_times, const TimeWindow time_budget) {
	List<TA>::iterator it{ li.end() - 1 };
	it.iter->data.maxShift = std::min(it.iter->data.timeWindow.closeTime - it.iter->data.depTime, time_budget.closeTime - it.iter->data.depTime);
	it--;
	while(it != li.end()){
		it.iter->data.maxShift = std::min(it.iter->data.timeWindow.closeTime - it.iter->data.depTime, (it + 1).iter->data.waitDuration + (it + 1).iter->data.maxShift);
		it--;
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
void ILS_TOPTW::validate(const List<TA>& walk, const Vector2D<double>& travel_times, const bool verbose) {
	if(verbose){
		walk.print("validating walk");
	}
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
	
	if (!valid) {
		std::cerr << "walk is invalid: exiting.. " << std::endl;
		std::exit(1);
	}
}

void ILS::drawSolution(const Solution& sol){

	const float factor = 2.0f;

	//Write vertices
	for(List<TA>::iterator ta_it = sol.m_unvisited.begin(); ta_it != sol.m_unvisited.end(); ++ta_it) {
		glPointSize(20.0f);
		glColor3f(1.0f, 1.0f, 1.0f);
		glBegin(GL_POINTS);
		glVertex2d(ta_it.iter->data.point.pos.lat*factor, ta_it.iter->data.point.pos.lon*factor);
		glEnd();

		glPointSize(18.0f);
		glColor3ub(170, 170, 170); //grey
		glBegin(GL_POINTS);
		glVertex2d(ta_it.iter->data.point.pos.lat*factor, ta_it.iter->data.point.pos.lon*factor);
		glEnd();
	}
	
	for(std::vector<Walk>::const_iterator walk_it = sol.m_walks.begin(); walk_it != sol.m_walks.end(); ++walk_it) {
		for(List<TA>::iterator ta_it = walk_it->begin(); ta_it != walk_it->end(); ++ta_it) {
			glPointSize(20.0f);
			glColor3f(1.0f, 1.0f, 1.0f);
			glBegin(GL_POINTS);
			glVertex2d(ta_it.iter->data.point.pos.lat*factor, ta_it.iter->data.point.pos.lon*factor);
			glEnd();

			glPointSize(18.0f);
			glColor3ub(200, 103, 51); //orange
			glBegin(GL_POINTS);
			glVertex2d(ta_it.iter->data.point.pos.lat*factor, ta_it.iter->data.point.pos.lon*factor);
			glEnd();
		}
	}

	//Write routes
	// glColor3ub(200, 200, 200);
	glLineWidth(2.0);
	for(std::vector<Walk>::const_iterator walk_it = sol.m_walks.begin(); walk_it != sol.m_walks.end(); ++walk_it) {
		glBegin(GL_LINE_LOOP);
		for(List<TA>::iterator ta_it = walk_it->begin(); ta_it != walk_it->end(); ++ta_it) {
			glVertex2d(ta_it.iter->data.point.pos.lat*factor, ta_it.iter->data.point.pos.lon*factor);
		}
		glEnd();
	}


	//Write labels
	glColor3ub(255, 255, 255); //white

	for(List<TA>::iterator ta_it = sol.m_unvisited.begin(); ta_it != sol.m_unvisited.end(); ++ta_it) {
		double lat = ta_it.iter->data.point.pos.lat*factor;
		double lon = ta_it.iter->data.point.pos.lon*factor;
		std::string id = std::to_string(ta_it.iter->data.point.id);
		int length = id.size();


		glRasterPos2d(lat-length, lon-1);
		// Use a loop to draw each character of the text string
		for (const auto& c : id) {
			glutBitmapCharacter(GLUT_BITMAP_8_BY_13, c);
		}
	}

	for(std::vector<Walk>::const_iterator walk_it = sol.m_walks.begin(); walk_it != sol.m_walks.end(); ++walk_it) {
		for(List<TA>::iterator ta_it = walk_it->begin(); ta_it != walk_it->end(); ++ta_it) {
			double lat = ta_it.iter->data.point.pos.lat*factor;
			double lon = ta_it.iter->data.point.pos.lon*factor;
			std::string id = std::to_string(ta_it.iter->data.point.id);
			int length = id.size();


			glRasterPos2d(lat-length, lon-1);
			// Use a loop to draw each character of the text string
			for (const auto& c : id) {
				glutBitmapCharacter(GLUT_BITMAP_8_BY_13, c);
			}
		}
	}

}

void ILS::draw(){

	// Clear the screen to black
	glClearColor(0.0, 0.0, 0.0 , 1.0);
	glClear(GL_COLOR_BUFFER_BIT);
	glLoadIdentity();

	const float pointSize = 20.0f;
	
	// Use glPointSize() to set the size of the points
	glPointSize(pointSize);
	// Set the drawing color to orange
	glColor3ub(200, 103, 51);

	drawSolution(best_solution);

	// for(auto sol_it = best_solutions.begin(); sol_it != best_solutions.end(); ++sol_it){
	// 	drawSolution(*sol_it);
	// }

	// Flush the OpenGL buffers to the screen
	glFlush();
	std::cout << std::endl;
}
