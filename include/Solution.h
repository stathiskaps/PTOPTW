#pragma once
#include <stdexcept>
#include <climits>
#include "List.h"
#include "Custom.h"
#include <fstream>
#include "nlohmann/json.hpp"

using json = nlohmann::json;

class Solution {

	struct Bounds{
		double minLat, minLon, maxLat, maxLon;
	};

	friend class ILS;
	friend class ILS_OPTW;
	friend class ILS_TOPTW;
	List<TA> m_unvisited;
	std::vector<List<TA>> m_walks;
public:
	Solution();
	~Solution();
	Solution(List<TA>);
	Solution(TA, TA, List<TA>, double, double, int); 
	Solution(TA, TA, List<TA>, size_t);

	int getScores() const;
	size_t getVisits(size_t);
	int getMinWalkSize() const;
	void draw(std::string);
	void print(std::string, const bool verbose);
	void jsonOutput();


};