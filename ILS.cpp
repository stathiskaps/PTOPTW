#include "ILS.h"

ILS::ILS() {}

ILS::ILS(int intervalsNum, std::string instance, Configuration conf) : mIntervalsNum(intervalsNum), mInstance(instance), mConf(conf) {
	metrics.local_search = std::vector<double>(intervalsNum, 0);
	metrics.split_unvisited = 0.0;
	metrics.shake = 0.0;
	metrics.final_pos = 0;
	metrics.middle_pos = 0;
	metrics.validation_time = 0.0;
	metrics.second_phase_counter = 0;
	metrics.second_phase_window_sum = 0;
	metrics.second_phase_improved = 0;
	metrics.comparisons = 0;
}

ILS::~ILS() {
	best_solutions.clear();
}

std::map<std::string, std::vector<double>> ILS::getActivities(List<TA>& unvisited, std::vector<TimeWindow> intervals){
	std::map<std::string, std::vector<double>> activities;
	for(auto& n : unvisited) {
		std::vector<double> durations;
		durations.reserve(intervals.size());
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

	uint8_t not_improved = 0;
	const uint8_t max_times_not_improved = 25;

	std::vector<double> time_cuts, best_time_cuts;
	double time_cut = day_start_time;
	time_cuts.reserve(intervals_num+1);

	//Initialize time cuts
	time_cuts.push_back(day_start_time);
	for(size_t i = 0; i < intervals_num; ++i){
		time_cuts.push_back(time_cuts.back() + dur);
	}

	//Acceptance criterion
	while(not_improved < max_times_not_improved) {

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
		} else {
			not_improved++;
		}
		
		//Improve intervals
		const double bin_duration = time_cuts[most_used_bin->right_cut_index] - time_cuts[most_used_bin->left_cut_index];

		double reduce_left = 0, reduce_right = 0;
		const double reduce_amount = bin_duration * a;

		if(most_used_bin > bins.begin() && most_used_bin < bins.end()-1){
			std::vector<Bin>::iterator left_bin = most_used_bin-1;
			std::vector<Bin>::iterator right_bin = most_used_bin+1;

			if(right_bin->count == 0){
				if(left_bin->count == right_bin->count){
					reduce_left = reduce_amount * 1/2;
					reduce_right = reduce_left;
				} else {
					reduce_right = reduce_amount;
					reduce_left = 0;
				}
			} else {
				const double ratio = left_bin->count/right_bin->count;
				reduce_left = std::ceil(ratio * reduce_amount / (ratio + 1)) ;
				reduce_right = reduce_amount - reduce_left;
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
		usages.reserve(intervals.size());
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
	sol.m_unvisited.print("Unvisited", false);
	for(std::vector<List<TA>>::const_iterator walk_it = sol.m_walks.begin(); walk_it != sol.m_walks.end(); ++walk_it){
		const int index = walk_it - sol.m_walks.begin();
		walk_it->print("Walk " + std::to_string(index), false);
	}
	std::cout << "===================================================" << std::endl;
}


void ILS::printSolutions(const std::string tag, const std::vector<Solution>& sols){
	std::cout << tag << std::endl;
	for(std::vector<Solution>::const_iterator sol_it = sols.begin(); sol_it != sols.end(); ++sol_it){
		const int index = sol_it - sols.begin();
		const std::string tag = "Solution [" + std::to_string(index) + "]: ";
		printSolution(tag, *sol_it);
	}
}

void ILS::InitSolutions(std::vector<Solution>& sols, const std::vector<TimeWindow> intervals,  const OP& op){

	std::vector<Solution>::iterator first_sol = sols.begin(), last_sol = sols.end() - 1;

	for(auto& walk : first_sol->m_walks) {
		walk.push_front(op.mStartDepot);
		walk.front().timeWindow = TimeWindow{intervals[0].openTime, intervals[0].closeTime};
		updateTimes(walk, walk.begin(), false, op.mTravelTimes, intervals[0]);
	}

	// for(auto& walk : last_sol->m_walks) {
	// 	walk.push_back(op.mEndDepot);
	// 	walk.back().timeWindow = TimeWindow{intervals[intervals.size()-1].openTime, intervals[intervals.size()-1].closeTime};
	// 	updateTimes(walk, walk.begin(), false, op.mTravelTimes, intervals[intervals.size()-1]);
	// }

}

void ILS::connectAndValidateSolutions(const std::vector<Solution>& solutions, const OP& op){
	Solution temp;
	for(auto& sol : solutions){
		temp.m_unvisited.append(sol.m_unvisited);
	}

	for (size_t i = 0; i < op.m_walks_num; ++i) {
		List<TA> walk;
		for (auto& s : solutions) {
			walk.append(s.m_walks[i]);
		}
		updateTimes(walk, walk.begin(), false, op.mTravelTimes, op.mTimeWindow);
		temp.m_walks.push_back(walk);
	}

	validateDirectedSolution(temp, op, false);
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

size_t ILS::countNodes(const std::vector<Solution>& sols){
	size_t counter {};
	for(auto& sol : sols){
		for(auto& walk : sol.m_walks){
			counter += walk.size();
		}
		counter += sol.m_unvisited.size();
	}
	return counter;
}

int ILS::Solve(OP& op) {

	// std::cout.setstate(std::ios_base::failbit);
	
	// Graphics::myInit();
	// // Set the OpenGL display mode
	// glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
	// // Set the initial window size
	// glutInitWindowSize(WINDOW_WIDTH, WINDOW_HEIGHT);
	// // Create the window with the given title
	// glutCreateWindow("Best Solution");
	// glutMouseFunc(Graphics::mouseButton);
	// //Set the display callback function
	// setupDrawCallback();
	// glutReshapeFunc(Graphics::onResize);

	// // Enter the GLUT main loop
	// std::thread glutThread(glutMainLoop);

	std::cout << mInstance << " -m " << op.m_walks_num << " -s " << mIntervalsNum << std::endl;
	auto start = std::chrono::high_resolution_clock::now();

	List<TA> unvisited(op.mAttractions);
	const TimeWindow time_budget{op.mStartDepot.timeWindow.openTime, op.mEndDepot.timeWindow.closeTime};

	std::vector<TimeWindow> intervals = getIntervals(op.mAttractions, mIntervalsNum, op.mStartDepot.timeWindow.openTime, op.mEndDepot.timeWindow.closeTime);
	const std::map<std::string, std::vector<double>> activities = getActivities(unvisited, intervals);
	std::map<std::string, std::vector<ILS::Usage>> reg = initRegistry(unvisited, intervals); 

	int times_not_improved{ 0 }, best_score{ INT_MIN }, bestCounter{};
	size_t counter{};
	const size_t split_every_n_iterations = 10;

	std::vector<ILS::SR> shake_settings;
	for(size_t i = 0; i < mIntervalsNum; ++i){
		shake_settings.push_back(ILS::SR{1, 1});
	}

	std::vector<Solution> process_solutions(mIntervalsNum, Solution());
	//Initialize solutions
	for (auto& s : process_solutions) {
		for (size_t i = 0; i < op.m_walks_num; ++i) {
			s.m_walks.push_back(List<TA>());
		}
	}
	InitSolutions(process_solutions, intervals, op);
	List<TA> pool = std::move(unvisited);
	const size_t split_iteration = 10;

	std::function<bool()> stopping_criterion = [&]() { return times_not_improved > MAX_TIMES_NOT_IMPROVED; };
	if(mConf.time_limited_execution){
		stopping_criterion = [&]() {
			auto end = std::chrono::high_resolution_clock::now();
			auto time_elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);
			return std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count() >= mConf.execution_time_limit;
    	};
	}
 
	while (!stopping_criterion()) {

		if(pool.size() > op.mAttractions.size()){
			throw std::runtime_error("Something wrong with the pool");
		}

		splitUnvisitedList(process_solutions, pool, mIntervalsNum, reg, activities);
		connectAndValidateSolutions(process_solutions, op);
		SplitSearch(process_solutions, pool, intervals, op, reg);
		checkSolutions(process_solutions, intervals, op);
		connectAndValidateSolutions(process_solutions, op);
		int score = collectScores(process_solutions);

		if (score > best_score) {
			bestCounter = counter;
			best_score = score;
			best_solutions = process_solutions;
			validateTimes(best_solutions, op.mTravelTimes, false);
			for(auto& sr : shake_settings){
				sr.R = 1;
			}
			times_not_improved = 0;
		}
		else {
			times_not_improved++;
		}

		int removed_counter = SplitShake(process_solutions, shake_settings, op, intervals);
		// pool.clear();
		gatherUnvisited(process_solutions, pool);

		counter++;

	}
	auto end = std::chrono::high_resolution_clock::now();
	metrics.total_execution_time = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
	metrics.total_execution_time = std::round(metrics.total_execution_time * 1000) / 1000;

	std::cout.clear();
	best_solution = connectSolutions(best_solutions, op.m_walks_num);
	for(auto walk_it = best_solution.m_walks.begin(); walk_it != best_solution.m_walks.end(); ++walk_it){
		const size_t index = walk_it - best_solution.m_walks.begin();
		walk_it->push_back(op.mEndDepot);
		updateTimes(*walk_it, walk_it->begin(), false, op.mTravelTimes, op.mTimeWindow);
	}
	
	validateDirectedSolution(best_solution, op, false);
	std::cout << "Best counter: " << bestCounter << std::endl;
	std::cout << "Best score: " << best_score << std::endl;
	std::cout << "Execution time(s): " << metrics.total_execution_time << std::endl;
	std::cout << "Visits: " << best_solution.getVisits(2) << std::endl << std::endl;
	
	if(mConf.write_output) {
		std::ofstream outputFile;
		outputFile.open("output.txt", std::ios::out | std::ios::app);
		outputFile << mInstance << " -m " << op.m_walks_num << " -s " << mIntervalsNum << ":\t" << 
			best_score << "\t" << metrics.total_execution_time << "\t" << best_solution.getVisits(2) << "\t" << std::endl;
		outputFile.close();
	}

	best_solution.output();

	// Wait for the GLUT thread to finish
	// glutMainLoop();

	return best_score;


}

void ILS::printMetrics(){
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
	std::cout << "Total execution time: " << metrics.total_execution_time << " seconds" << std::endl;
	std::cout << "Total comparisons made: " << metrics.comparisons << std::endl;
	std::cout << "-------------------------------------" << std::endl;
	std::cout << std::endl;
}

void ILS::gatherUnvisited(std::vector<Solution>& solutions, List<TA>& pool){
	for(auto& sol : solutions){
		pool.append(std::move(sol.m_unvisited));
	}
}

std::vector<List<TA>> ILS::splitUnvisitedList(std::vector<Solution>& solutions, List<TA>& pool, int intervals_num, std::map<std::string, std::vector<ILS::Usage>>& reg, std::map<std::string, std::vector<double>> activities) {

	auto start = std::chrono::high_resolution_clock::now();
	if(solutions.size() != intervals_num){
		throw std::runtime_error("solutions size is indifferent from intervals number");
	}


	std::vector<List<TA>> lists(intervals_num, List<TA>());
	//Split Unvisited
	for (List<TA>::iterator ta_it = pool.begin(), next_ta; ta_it != pool.end(); ta_it = next_ta) {
		next_ta = ta_it.next();
		const double total_duration{ ta_it.iter->data.timeWindow.closeTime - ta_it.iter->data.timeWindow.openTime };
		double max_score{ -1 };

		if(reg[ta_it.iter->data.id].size() != activities[ta_it.iter->data.id].size()){
			throw std::runtime_error("vectors are not of equal size");
		}

		double best_score = DBL_MIN;
		std::vector<ILS::Usage>::iterator best_it{ reg[ta_it.iter->data.id].end() };
		std::vector<ILS::Usage>::iterator usage_it;
		std::vector<double>::iterator duratio_it;
		for(usage_it = reg[ta_it.iter->data.id].begin(), duratio_it = activities[ta_it.iter->data.id].begin(); usage_it != reg[ta_it.iter->data.id].end() ; ++usage_it, ++duratio_it){
			double score = *duratio_it * (usage_it->solved / static_cast<double>(usage_it->imported));
			// double score = *duratio_it;
			if (score>best_score){
				best_score = score;
				best_it = usage_it;
			}
		}
		if(best_it == reg[ta_it.iter->data.id].end()){
			throw std::runtime_error("didn't find a good interval for node " + ta_it.iter->data.id);
		}

		const int index = best_it - reg[ta_it.iter->data.id].begin();
		solutions.at(index).m_unvisited.emplace_back(pool, ta_it);
		
		reg[ta_it.iter->data.id].at(index).imported += 1;
	}

	auto end = std::chrono::high_resolution_clock::now();
	auto elapsed_time = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
	metrics.split_unvisited += elapsed_time;

	return lists;
}

Solution ILS::connectSolutions(std::vector<Solution>& sols, const size_t walks_num) {
	Solution solution;
	gatherUnvisited(sols, solution.m_unvisited);
	for (size_t i = 0; i < walks_num; ++i) {
		List<TA> walk;
		for (auto& s : sols) {
			walk.append(std::move(s.m_walks[i]));
		}
		solution.m_walks.push_back(walk);
	}

	return solution;
}

int ILS::collectScores(const std::vector<Solution>& sols) const {
	int total_score{};
	for (auto& sol : sols) {
		total_score += sol.getScores();
	}
	return total_score;
}

void ILS::PrepareForShake(std::vector<Solution>& sols){
	for(std::vector<Solution>::iterator sol_it = sols.begin(); sol_it != sols.end(); ++sol_it){
		if(sol_it != sols.begin()){
			for(std::vector<List<TA>>::iterator walk_it = sol_it->m_walks.begin(); walk_it != sol_it->m_walks.end(); ++walk_it){
				walk_it->push_front(walk_it->front());
				walk_it->front().id = DUMMY_ID;
			}
		}
		// if(sol_it != sols.end() - 1) {

		// }
		for(std::vector<List<TA>>::iterator walk_it = sol_it->m_walks.begin(); walk_it != sol_it->m_walks.end(); ++walk_it){
			walk_it->push_back(walk_it->back());
			walk_it->back().id = DUMMY_ID;;
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

int ILS::SplitShake(std::vector<Solution>& sols, std::vector<ILS::SR>& shake_settings, OP& op, const std::vector<TimeWindow> intervals){
	auto start = std::chrono::high_resolution_clock::now();
	int removed_counter = 0;
	PrepareForShake(sols);
	for(auto sol_it = sols.begin(); sol_it != sols.end(); ++sol_it){
		const int sol_index = sol_it - sols.begin();
		int removed = Shake(*sol_it, shake_settings[sol_index].S, shake_settings[sol_index].R, op, intervals[sol_index]);
		removed_counter += removed;
	}
	RemoveDummyNodes(sols);
	auto end = std::chrono::high_resolution_clock::now();
	auto elapsed_time = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
	metrics.shake += elapsed_time;
	return removed_counter;
}

int ILS::Shake(Solution& sol, int& S, int& R, OP& op, const TimeWindow time_budget) {
	int minWalkSize = sol.getMinWalkSize();
	int num_locations = 0;
	for(auto& walk : sol.m_walks){
		num_locations += walk.size();
	}
	const int max_to_remove = num_locations / (2 * sol.m_walks.size());

	int removed_counter = 0;
	if (S >= minWalkSize - 2) {
		S = 1;
		R = 1;
	}

	if (R == max_to_remove) {
		R = 1;
	}

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
	R += 1;

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

std::tuple<bool, double> ILS::CandidateStartDepotIsValid(const List<TA>& walk, const TA& candidate, const double start_time, const Vector2D<double>& travel_times){
	auto first_it = walk.begin();
	double wait_dur = std::max(0.0, candidate.timeWindow.openTime - start_time);
	double start_of_visit_time = start_time + wait_dur;
	double dep_time = start_of_visit_time + candidate.visitDuration;
	double shift = wait_dur
		+ candidate.visitDuration
		+ travel_times[candidate.point.id][first_it.iter->data.arrPointId];

	return { dep_time <= candidate.timeWindow.closeTime && shift <= first_it.iter->data.maxShift, shift };
}


void ILS::RemoveUnfeasibleVisits(std::vector<Solution>& solutions, const int i, const size_t j){
	for(List<TA>::iterator it = solutions[i].m_walks[j].begin(); it != solutions[i].m_walks[j].end(); ++it){
		if(it.iter->data.maxShift < 0){
			solutions[i].m_unvisited.push_back(it.iter->data);
			solutions[i].m_walks[j].erase(it);
		}
	}
}

//i is the solution index
//j is the walk index
TA ILS::getValidPreviousTA(std::vector<Solution>& solutions, const int i, const size_t j) {
	int prev_sol_index = i-1;
	while(solutions[prev_sol_index].m_walks[j].empty() && prev_sol_index >= 0){
		prev_sol_index--;
		continue;
	}

	if(prev_sol_index == -1){
		throw std::runtime_error("didn't find a valid previous walk to get a startDepot");
	}
	TA ta(solutions[prev_sol_index].m_walks[j].back());
	return ta;
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
		TA dummy = prev_last;
		dummy.neutralize(DUMMY_START_DEPOT);
		dummy.timeWindow = intervals[i];

		//left trim invalid nodes
		if(!solutions[i].m_walks[j].empty()) { 
			List<TA>::iterator curr = solutions[i].m_walks[j].begin();
			while(curr != solutions[i].m_walks[j].end()){
				auto [valid, _] = CandidateStartDepotIsValid(solutions[i].m_walks[j], dummy, intervals[i].openTime, op.mTravelTimes);
				if(!valid){ 
					solutions[i].m_unvisited.push_back(curr.iter->data);
					curr = solutions[i].m_walks[j].erase(curr);
					updateTimes(solutions[i].m_walks[j], curr, false, op.mTravelTimes, intervals[i]);
				} else {
					break;
				}
			}
		}
		solutions[i].m_walks[j].push_front(dummy);
		updateTimes(solutions[i].m_walks[j], solutions[i].m_walks[j].begin(), false, op.mTravelTimes, intervals[i]);
	}
}

std::vector<Point> ILS::getTargets(const std::vector<Solution>& solutions, const int i, const OP& op){
	if(last_solution){
		return std::vector<Point>(solutions[i].m_walks.size(), op.mEndDepot.point); //END_DEPOT
	}

	const int min_size = 3;
	std::vector<Point> targets;
	for(size_t j = 0; j < solutions[i].m_walks.size(); ++j){
		int next_sol_index = i + 1;
		while (next_sol_index < solutions.size() && !hasWeightedCentroid(solutions[next_sol_index], j, min_size)) {
			next_sol_index++;
		}

		if(next_sol_index == solutions.size()){
			Point p = op.mEndDepot.point;
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

void ILS::SplitSearch(std::vector<Solution>& solutions, List<TA>& pool, const std::vector<TimeWindow>& intervals, OP& op, std::map<std::string, std::vector<ILS::Usage>>& reg) {

	const int min_size = 3;
	const int solutions_size = solutions.size();
	int counter1{}, counter2{}, tempScore1{}, tempScore2{};

	for (size_t i = 0; i < solutions.size(); ++i) {

		auto start = std::chrono::high_resolution_clock::now();

		if (solutions[i].m_unvisited.empty()) continue;

		std::vector<bool> added_start_depot;
		if(i > 0) {
			AddStartDepots(solutions, intervals, i, op);
		}

		std::vector<Point> targets = getTargets(solutions, i, op);
		std::vector<std::string> ids = construct(solutions[i], op.mTravelTimes, targets, intervals[i], true);
		counter1 += ids.size();

		for (auto& id : ids) {
			reg[id].at(i).solved++;
		}

		if (i > 0) {
			
			//Make a local search betweem current solution and previous
			for (size_t j = 0; j < solutions[i].m_walks.size(); ++j) {
				solutions[i].m_walks[j].pop_front();

				//TODO: check if these comments are needed
				// if(!solutions[i].m_walks[j].empty()){
				// 	const TA prev_last = getValidPreviousTA(solutions, i, j);
				// 	TA &curr_first = solutions[i].m_walks[j].front();
				// 	double start_time = std::max(prev_last.depTime+op.mTravelTimes[prev_last.point.id][curr_first.point.id], intervals[i].openTime);
				// 	updateTimes(solutions[i].m_walks[j], solutions[i].m_walks[j].begin(), false, op.mTravelTimes, TimeWindow{start_time, intervals[i].closeTime});
				// }
				
			}
		}

		validateTimes(solutions[i], op.mTravelTimes, false);

		auto end = std::chrono::high_resolution_clock::now();
		auto elapsed_time = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
		metrics.local_search[i] += elapsed_time;

	}
}

std::vector<std::string> ILS::construct(Solution& sol, const Vector2D<double>& travel_times, const std::vector<Point> targets, 
	const TimeWindow time_budget, const bool open) {
	
	std::vector<std::string> inserted_ids;
	if (sol.m_unvisited.empty()) {
		return inserted_ids;
	}

	double best_pos_score{}, max_score{ DBL_MIN }, score{};
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
			best_pos_score = DBL_MAX;
			std::vector<double> dists = distances[curr.iter->data.point.id];
			auto [walk_it, pos, best_pos_score, arr_point_id, dep_point_id] = getBestPos(curr.iter->data, sol.m_walks, travel_times, dists, time_budget, open);
			//default value of walk_it is sol.m_walks.end()
			//default value of pos is sol.m_walks.front().end()
			//if these iterators have these values then insertion is infeasible

			if (walk_it == sol.m_walks.end()) {
				curr++;
				continue;
			}

			const int walk_index = walk_it - sol.m_walks.begin();
			const double distance_to_target = distances[curr.iter->data.depPointId].at(walk_index);
			const size_t best_pos_index = pos - walk_it->begin();
			const size_t total_pos = walk_it->size();
			score = pow(curr.iter->data.profit, 2) / best_pos_score;	//score = profit^2/ minShifti

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


std::tuple<Walks::iterator, List<TA>::iterator, double, int, int> ILS::getBestPos(const TA& ta, 
	Walks& walks, const Vector2D<double>& travel_times, const std::vector<double> distances, const TimeWindow time_budget, const bool open) {

	Walks::iterator best_walk{ walks.end() };
	List<TA>::iterator best_pos{ walks.front().end()};
	int arr_point_id{ DEFAULT_POINT_ID }, dep_point_id{ DEFAULT_POINT_ID };
	double min_score{ DBL_MAX };

	for (Walks::iterator walk_it = walks.begin(); walk_it != walks.end(); ++walk_it) {
		if (walk_it->empty()) {
			std::cerr << "walk is empty, minimum size: 1" << std::endl;
			std::exit(1);
		}

		const size_t walk_index = walk_it - walks.begin();
		const double distance_to_next_target = distances[walk_index];

		double arr_time{}, wait_dur{}, start_of_visit_time{}, dep_time{}, shift{}, pos_score{};

		const size_t possible_positions { walk_it->size() };
		int temp_arr_point_id{ DEFAULT_POINT_ID }, temp_dep_point_id{ DEFAULT_POINT_ID }, pos_counter{};
		int insertion_pos_index = 1;

		List<TA>::iterator left{ walk_it->begin() }, right{ left + 1 };
		while(right != walk_it->end()) {
			metrics.comparisons++;
			arr_time = left.iter->data.depTime + travel_times[left.iter->data.depPointId][ta.point.id];
			wait_dur = std::max(0.0, ta.timeWindow.openTime - arr_time);
			start_of_visit_time = arr_time + wait_dur;
			dep_time = start_of_visit_time + ta.visitDuration;
			shift = travel_times[left.iter->data.depPointId][ta.point.id]
				+ wait_dur
				+ ta.visitDuration
				+ travel_times[ta.point.id][right.iter->data.arrPointId]
				- travel_times[left.iter->data.depPointId][right.iter->data.arrPointId];
			pos_score = shift + distance_to_next_target * (insertion_pos_index/possible_positions);
			if(dep_time <= ta.timeWindow.closeTime && shift <= right.iter->data.maxShift + right.iter->data.waitDuration) {
				if (pos_score < min_score) {
					best_walk = walk_it;
					best_pos = right;
					min_score = pos_score;
					arr_point_id = ta.point.id;
					dep_point_id = ta.point.id;
				}
			}
			left++; right++; insertion_pos_index++;
		}

		if(open){
			metrics.comparisons++;
			//here we check the insertion at the end of the walk
			arr_time = left.iter->data.depTime + travel_times[left.iter->data.depPointId][ta.point.id];
			wait_dur = std::max(0.0, ta.timeWindow.openTime - arr_time);
			start_of_visit_time = arr_time + wait_dur;
			dep_time = start_of_visit_time + ta.visitDuration;
			shift = travel_times[left.iter->data.depPointId][ta.point.id]
				+ wait_dur
				+ ta.visitDuration;
			pos_score = shift + distance_to_next_target * (insertion_pos_index/possible_positions);

			if(dep_time <= ta.timeWindow.closeTime && dep_time <= time_budget.closeTime) {
				if (pos_score < min_score) {
					metrics.final_pos++;
					best_walk = walk_it; 
					best_pos = right;
					min_score = pos_score;
					arr_point_id = ta.point.id;
					dep_point_id = ta.point.id;
				} else {
					metrics.middle_pos++;
				}
			}
		}
	}

	return { best_walk, best_pos, min_score, arr_point_id, dep_point_id };
}

void ILS::updateBasicTimes(List<TA>::iterator& it, const Vector2D<double>& travel_times, const TimeWindow time_budget){
	it.iter->data.arrTime = it.prev().iter->data.depTime + travel_times[it.prev().iter->data.depPointId][it.iter->data.arrPointId];
	it.iter->data.waitDuration = std::max(0.0, it.iter->data.timeWindow.openTime - it.iter->data.arrTime);
	it.iter->data.startOfVisitTime = it.iter->data.arrTime + it.iter->data.waitDuration;
	it.iter->data.depTime = it.iter->data.startOfVisitTime + it.iter->data.visitDuration;
}

bool ILS::updateTimes(List<TA>& walk, const List<TA>::iterator& start_pos, 
	const bool smart, const Vector2D<double>& travel_times, const TimeWindow time_budget) {
	if(walk.empty()) return true;
	//update times first
	List<TA>::iterator it{ start_pos };
	while (it != walk.end()) {
		if(it == walk.begin()) {
			it.iter->data.arrTime = time_budget.openTime;
			it.iter->data.waitDuration = std::max(0.0, it.iter->data.timeWindow.openTime - it.iter->data.arrTime);
			it.iter->data.startOfVisitTime = it.iter->data.arrTime + it.iter->data.waitDuration;
			it.iter->data.depTime = it.iter->data.startOfVisitTime + it.iter->data.visitDuration;
			if(it.next() != walk.end()){
				it.iter->data.shift = it.iter->data.waitDuration
				+ it.iter->data.visitDuration
				+ travel_times[it.iter->data.depPointId][it.next().iter->data.arrPointId];
			} else {
				it.iter->data.shift = it.iter->data.waitDuration + it.iter->data.visitDuration;
			}
			
		} else if(it == walk.end()-1){
			updateBasicTimes(it, travel_times, time_budget);
			it.iter->data.shift = travel_times[it.prev().iter->data.depPointId][it.iter->data.arrPointId]
				+ it.iter->data.waitDuration
				+ it.iter->data.visitDuration;
		} else {
			updateBasicTimes(it, travel_times, time_budget);
			it.iter->data.shift = travel_times[it.prev().iter->data.depPointId][it.iter->data.arrPointId]
				+ it.iter->data.waitDuration
				+ it.iter->data.visitDuration
				+ travel_times[it.iter->data.depPointId][it.next().iter->data.arrPointId];
		}
		it++;
	}

	updateMaxShifts(walk, travel_times, time_budget);
	return walk.back().maxShift >= 0;
}

void ILS::updateMaxShifts(const List<TA>& li, const Vector2D<double>& travel_times, const TimeWindow time_budget) {
	List<TA>::iterator it{ li.end() - 1 };
	//TODO: write in documentation that maxShift of last node depends not only to its window but also to time budget
	//cause now the last node isn't the endDepot of the problem
	it.iter->data.maxShift = std::min(it.iter->data.timeWindow.closeTime - it.iter->data.depTime, time_budget.closeTime - it.iter->data.depTime);
	it--;
	while(it != li.end()){
		it.iter->data.maxShift = std::min(it.iter->data.timeWindow.closeTime - it.iter->data.depTime, (it + 1).iter->data.waitDuration + (it + 1).iter->data.maxShift);
		it--;
	}
}

bool compareByShift(const TA &a, const TA &b) {
    return a.shift < b.shift;
}

std::vector<std::string> ILS::fixWalk(List<TA>& walk, const OP& op, TimeWindow time_budget){
	if(walk.size() < 2){
		throw std::runtime_error("Invalid walk size: " + std::to_string(walk.size()));
	}
	std::vector<std::string> remove_nodes;
	while(!updateTimes(walk, walk.begin(), false, op.mTravelTimes, time_budget)){
		double minScore = DBL_MAX;
		List<TA>::iterator remove_it = walk.end();
		for(List<TA>::iterator it = walk.begin()+1; it != walk.end()-1; ++it){
			double score = pow(it.iter->data.profit, 2) / it.iter->data.shift;
			if(score < minScore){
				remove_it = it;
				minScore = score;
			}
		}

		if(remove_it != walk.end()){
			remove_nodes.push_back(remove_it.iter->data.id);
			walk.erase(remove_it);
		}
	}
	return remove_nodes;
}

void ILS::checkSolutions(std::vector<Solution>& sols,  const std::vector<TimeWindow>& intervals, const OP& op){
	TimeWindow custom_budget = {op.mTimeWindow.openTime, op.mTimeWindow.closeTime/2}; //TODO: remove that
	for(size_t j = 0; j < op.m_walks_num; ++j){
		List<TA> walk;
		for(auto& sol: sols){
			walk.append(sol.m_walks[j]);
		}
		walk.push_back(op.mEndDepot);
		std::vector<std::string> unvisit_ids = fixWalk(walk, op, op.mTimeWindow);
		if(unvisit_ids.empty()) {
			continue;
		}
		std::sort(unvisit_ids.begin(), unvisit_ids.end());

		for(size_t i = 0; i < sols.size(); ++i){
			bool modified = false;
			for(List<TA>::iterator ta_it = sols[i].m_walks[j].begin(); ta_it != sols[i].m_walks[j].end(); ){
				auto it = std::lower_bound(unvisit_ids.begin(), unvisit_ids.end(), ta_it.iter->data.id);
				if(*it == ta_it.iter->data.id){
					sols[i].m_unvisited.push_back(ta_it.iter->data);
					ta_it = sols[i].m_walks[j].erase(ta_it);
					unvisit_ids.erase(it);
					modified = true;
					if(unvisit_ids.empty()) {
						updateTimes(sols[i].m_walks[j], sols[i].m_walks[j].begin(), false, op.mTravelTimes, intervals[i]);
						goto done;
					}
				} else {
					++ta_it;
				}
			}
			if(modified){
				updateTimes(sols[i].m_walks[j], sols[i].m_walks[j].begin(), false, op.mTravelTimes, intervals[i]);
			}
		}

		done:
		continue;
	}
	
}

void ILS::validateDirectedSolution(const Solution& solution, const OP& op, const bool verbose){
	std::vector<std::string> v1;
	for(auto& ta : solution.m_unvisited){
		v1.push_back(ta.id);
	}

	for(auto& walk : solution.m_walks){
		for(auto& ta : walk){
			v1.push_back(ta.id);
		}
	}

	v1.erase(std::remove(v1.begin(), v1.end(), op.mStartDepot.id), v1.end());
	v1.erase(std::remove(v1.begin(), v1.end(), op.mEndDepot.id), v1.end());

	std::vector<std::string> v2;
	for(auto& ta : op.mAttractions){
		v2.push_back(ta.id);
	}

	std::sort(v1.begin(), v1.end());
	std::sort(v2.begin(), v2.end());

	if(v1.size() != v2.size()){
		throw std::runtime_error("Solution has different size of nodes than the problem");
	}

	for(auto& walk : solution.m_walks){
		validateTimes(walk, op.mTravelTimes, verbose);
	}
}

void ILS::validateTimes(const Solution& sol, const Vector2D<double>& travel_times, const bool verbose) {
	for(auto walk_it = sol.m_walks.begin(); walk_it != sol.m_walks.end(); ++walk_it){
		validateTimes(*walk_it, travel_times, verbose);
	}
}

void ILS::validateTimes(const std::vector<Solution>& sols, const Vector2D<double>& travel_times, const bool verbose){
	for(auto sol_it = sols.begin(); sol_it != sols.end(); ++sol_it){
		const int index = sol_it - sols.begin();
		validateTimes(*sol_it, travel_times, verbose);
	}
}

void ILS::validateTimes(const List<TA>& walk, const Vector2D<double>& travel_times, const bool verbose) {
	if(verbose){
		walk.print("validating walk", verbose);
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
		std::cerr << msg << std::endl;
		std::cerr << "walk is invalid: exiting.. " << std::endl;
		std::exit(1);
	}

}

void ILS::drawSolution(const Solution& sol){

	const float factor = 2.0f;

	//Draw vertices
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
			const size_t i = walk_it - sol.m_walks.begin();
			glPointSize(20.0f);
			glColor3f(1.0f, 1.0f, 1.0f);
			glBegin(GL_POINTS);
			glVertex2d(ta_it.iter->data.point.pos.lat*factor, ta_it.iter->data.point.pos.lon*factor);
			glEnd();

			glPointSize(18.0f);
			glColor3ub(Graphics::colors[i].r, Graphics::colors[i].g, Graphics::colors[i].b); //orange
			glBegin(GL_POINTS);
			glVertex2d(ta_it.iter->data.point.pos.lat*factor, ta_it.iter->data.point.pos.lon*factor);
			glEnd();
		}
	}

	//Draw routes
	glLineWidth(2.0);
	for(std::vector<Walk>::const_iterator walk_it = sol.m_walks.begin(); walk_it != sol.m_walks.end(); ++walk_it) {
		const size_t i = walk_it - sol.m_walks.begin();
		glColor3ub(Graphics::colors[i].r, Graphics::colors[i].g, Graphics::colors[i].b);
		glBegin(GL_LINE_LOOP);
		for(List<TA>::iterator ta_it = walk_it->begin(); ta_it != walk_it->end(); ++ta_it) {
			glVertex2d(ta_it.iter->data.point.pos.lat*factor, ta_it.iter->data.point.pos.lon*factor);
		}
		glEnd();
	}


	//Draw labels
	glColor3ub(255, 255, 255); //white

	for(List<TA>::iterator ta_it = sol.m_unvisited.begin(); ta_it != sol.m_unvisited.end(); ++ta_it) {
		double lat = ta_it.iter->data.point.pos.lat*factor;
		double lon = ta_it.iter->data.point.pos.lon*factor;
		std::string id = std::to_string(ta_it.iter->data.point.id);
		int length = id.size();

		glRasterPos2d(lat-length-1, lon-1);
		// Use a loop to draw each character of the text string
		for (const auto& c : id) {
			glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_10, c);
		}
	}

	for(std::vector<Walk>::const_iterator walk_it = sol.m_walks.begin(); walk_it != sol.m_walks.end(); ++walk_it) {
		for(List<TA>::iterator ta_it = walk_it->begin(); ta_it != walk_it->end(); ++ta_it) {
			double lat = ta_it.iter->data.point.pos.lat*factor;
			double lon = ta_it.iter->data.point.pos.lon*factor;
			std::string id = std::to_string(ta_it.iter->data.point.id);
			int length = id.size();


			glRasterPos2d(lat-length-1, lon-1);
			// Use a loop to draw each character of the text string
			for (const auto& c : id) {
				glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_10, c);
			}
		}
	}

}

void ILS::draw(){

	Graphics::myInit();

	// Clear the screen to black
	glClearColor(0.0, 0.0, 0.0 , 1.0);
	glClear(GL_COLOR_BUFFER_BIT);
	glLoadIdentity();

	glScalef(Graphics::zoom, Graphics::zoom, 1.0);

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
