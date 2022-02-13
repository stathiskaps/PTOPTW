#include <iostream>
#include <vector>
#include <algorithm>
#include <set>
#include <map>
#include <climits>
#include "List.h"


#ifndef DIVIDER_H
#define DIVIDER_H

enum class timeWindowEventType { none, open, close };

struct TimeWindowEvent {
	int id;
	std::string nodeId;
	double time;
	timeWindowEventType type;
	TimeWindowEvent() :id(0), nodeId(0), time(0), type(timeWindowEventType::none) {};
	TimeWindowEvent(int pId, std::string pNodeId, double pTime, timeWindowEventType pType) :
		id(pId),
		nodeId(pNodeId),
		time(pTime),
		type(pType) {};

	void print() {
		std::string typeStr;
		switch (type) {
		case timeWindowEventType::open:
			typeStr = "open";
			break;
		case timeWindowEventType::close:
			typeStr = "close";
			break;
		default:
			typeStr = "none";
		}
		std::cout << id
			<< " " << nodeId
			<< " " << time
			<< " " << typeStr
			<< std::endl;
	}
};

struct TimeEdge {
	int id;
	double time;
	std::vector<TimeWindowEvent> events;
	std::vector<std::string> activeNodes, openingNodes, closingNodes;
	TimeEdge() : id(0), time(0.0) {}

	void print() {
		std::cout << "TimeEdge id: " << id
			<< "\tTime:  " << time
			<< std::endl;

		std::cout << "/=========Time Window Events=========/" << std::endl;
		for (auto& e : events) {
			e.print();
		}
		std::cout << "=====================================\n" << std::endl;

		std::cout << "/=========Opening Nodes=========/" << std::endl;
		std::cout << "Number: " << openingNodes.size() << std::endl;
		for (auto& n : openingNodes) {
			std::cout << n << "\t";
		}
		std::cout << "\n=====================================\n" << std::endl;

		std::cout << "/=========Closing Nodes=========/" << std::endl;
		std::cout << "Number: " << closingNodes.size() << std::endl;
		for (auto& n : closingNodes) {
			std::cout << n << "\t";
		}
		std::cout << "\n=====================================\n" << std::endl;
	}
};

class Divider {
private:
	int mWalksNum;
	std::vector<TA*> mAttractions;
public:
	Divider();
	Divider(std::vector<TA*>, int);
	~Divider();

	std::vector<std::vector<Cluster>> exec();

	std::vector<Cluster> kmeans(std::vector<TA*>);

	std::vector<std::string> getActiveNodes(std::vector<TimeEdge>&, int, int);
	double calcScore(int, int, double, double);
	double calcScore2(int, int, double, double);
	std::vector<Cluster> intervals(std::vector<TA*>);
	std::vector<TimeEdge> getTimeEdges(const std::vector<TA*>&);
	std::tuple<int, int> getBestInterval(std::vector<TimeEdge>&, std::map<std::string, int>&, int);
	static bool compareInterval(TimeWindowEvent, TimeWindowEvent);
	std::vector<TA*> findBestAttractionCombination(std::vector<TA*>);
};

#endif