#pragma once
#include <stdexcept>
#include <climits>
#include "List.h"
#include "Custom.h"
#include <fstream>


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
	Solution(List<TA> unvisited);
	Solution(TA start, TA end, List<TA> unvisited, double startTime, double endTime, int walksNum); 

	int getScores();
	int getVisits();
	int getMinWalkSize() const;
	void draw(std::string, std::string);
	void print(std::string, const bool verbose);



};