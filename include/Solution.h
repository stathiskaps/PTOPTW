#pragma once
#include <stdexcept>
#include <climits>
#include "List.h"
#include <fstream>


class Solution {

	struct Bounds{
		double minLat, minLon, maxLat, maxLon;
	};

	friend class ILS;
	friend class ILS_OPTW;
	friend class ILS_TOPTW;
	std::list<TA> m_unvisited;
	std::vector<std::list<TA>> m_walks;
public:
	Solution();
	~Solution();
	Solution(std::list<TA>);
	Solution(TA, TA, std::list<TA>, double, double, int); 
	Solution(TA, TA, std::list<TA>, size_t);

	int getScores() const;
	int getVisits();
	int getMinWalkSize() const;
	void draw(std::string, std::string);
	void print(std::string, const bool verbose);



};