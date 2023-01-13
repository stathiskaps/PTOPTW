#pragma once

#include <iostream>
#include <algorithm>
#include <tuple>
#include <limits>
#include <cfloat>
#include <type_traits> // enable_if, conjuction
#include <cstdarg>
#include <functional>
#include "Definitions.h"

enum class touristAttractionType {none, sight, route};
enum class sightCategory {none, depot, museum, agora};

struct Position {
	double lat, lon;

	Position() {
		lat = 0;
		lon = 0;
	}
	Position(double pLat, double pLon) {
		lat = pLat;
		lon = pLon;
	}
	~Position() {}
};

struct Point {
	int id;
	Position pos;

	Point() {
		id = DEFAULT_POINT_ID;
	}

	Point(int pId){
		id = pId;
	}

	Point(int pId, Position pPos) {
		id = pId;
		pos = pPos;
	}

	Point(int pId, double pLat, double pLon) {
		id = pId;
		pos = Position(pLat, pLon);
	}

	Point(double pLat, double pLon) {
		id = DEFAULT_POINT_ID;
		pos = Position(pLat, pLon);
	}

	~Point() {}

	double euclidean_distance(Point p) {
		double distance = sqrt(pow(p.pos.lat - this->pos.lat, 2) + pow(p.pos.lon - this->pos.lon, 2));;
		return distance;
	}

	//Returns a point towards point p2, with distance dt from current point
	Point findPointWithDistance(Point p2, double dt){
		double d = euclidean_distance(p2);
		double t = dt/d;
		return Point(((1-t)*this->pos.lat + t*p2.pos.lat), ((1-t)*this->pos.lon + t*p2.pos.lon));
	}
};

struct TimeWindow {
	double openTime;
	double closeTime;

	double duration() const {
		return closeTime - openTime;
	}

	double center() const {
		return (openTime + closeTime) / 2;
	}
};

typedef struct TouristAttraction {
	std::string id;
	Point point;
	TimeWindow timeWindow;
	touristAttractionType type;
	sightCategory category;
	double visitDuration;
	int profit;
	int walk;
	double arrTime, waitDuration, startOfVisitTime, depTime, shift, maxShift;
	int route, cluster;
	double minDist;  // default infinite dist to nearest cluster
	int arrPointId, depPointId;


	//default constructor
	TouristAttraction()
		: id(DEFAULT_TA_ID),
		point(Point()),
		visitDuration(.0),
		profit(0),
		timeWindow(TimeWindow{ .0, .0 }),
		arrTime(.0),
		waitDuration(.0),
		startOfVisitTime(.0),
		depTime(.0),
		shift(.0),
		maxShift(.0),
		route(-1),
		walk(UNDEFINED),
		cluster(UNDEFINED),
		minDist(DBL_MAX),
		arrPointId(DEFAULT_POINT_ID),
		depPointId(DEFAULT_POINT_ID),
		type(touristAttractionType::sight),
		category(sightCategory::none)
	{
	}

	TouristAttraction(const TouristAttraction &other)
		: id(other.id),
		point(other.point),
		visitDuration(other.visitDuration),
		profit(other.profit),
		timeWindow(other.timeWindow),
		arrTime(other.arrTime),
		waitDuration(other.waitDuration),
		startOfVisitTime(other.startOfVisitTime),
		depTime(other.depTime),
		shift(other.shift),
		maxShift(other.maxShift),
		route(other.route),
		walk(other.walk),
		cluster(other.cluster),
		minDist(other.minDist),
		arrPointId(other.arrPointId),
		depPointId(other.depPointId),
		type(other.type),
		category(other.category)
	{
	}

	TouristAttraction(std::string pId, Point p)
		: id(pId),
		point(p),
		visitDuration(.0),
		profit(0),
		timeWindow(TimeWindow{ .0, .0 }),
		arrTime(.0),
		waitDuration(.0),
		startOfVisitTime(.0),
		depTime(.0),
		shift(.0),
		maxShift(.0),
		route(-1),
		walk(UNDEFINED),
		cluster(UNDEFINED),
		minDist(DBL_MAX),
		arrPointId(p.id),
		depPointId(p.id),
		type(touristAttractionType::sight),
		category(sightCategory::none)
	{
	}

	TouristAttraction(std::string pId) : TouristAttraction(pId, Point()) {}

	TouristAttraction(Point p) : TouristAttraction(DEPOT_ID, p) {}

	TouristAttraction(std::string pId, Point p, double pVisitDuration, int pProfit, double pOpenTime, double pCloseTime)
		: id(pId),
		point(p),
		visitDuration(pVisitDuration),
		profit(pProfit),
		timeWindow(TimeWindow{ pOpenTime, pCloseTime }),
		arrTime(.0),
		waitDuration(.0),
		startOfVisitTime(.0),
		depTime(.0),
		shift(.0),
		maxShift(.0),
		route(-1),
		walk(UNDEFINED),
		cluster(UNDEFINED),
		minDist(DBL_MAX),
		arrPointId(p.id),
		depPointId(p.id),
		type(touristAttractionType::sight),
		category(sightCategory::none)
	{
	}

	//constructor for depots
	TouristAttraction(Point p, double pDepTime, double pOpenTime, double pCloseTime)
		: id(DEPOT_ID),
		point(p),
		visitDuration(0),
		profit(0),
		timeWindow(TimeWindow{ pOpenTime, pCloseTime }),
		arrTime(.0),
		waitDuration(.0),
		startOfVisitTime(.0),
		depTime(pDepTime),
		shift(.0),
		maxShift(.0),
		route(-1),
		walk(UNDEFINED),
		cluster(UNDEFINED),
		minDist(DBL_MAX),
		arrPointId(p.id),
		depPointId(p.id),
		type(touristAttractionType::sight),
		category(sightCategory::none)
	{
	}

	void neutralize(std::string pId){
		id = pId;
		visitDuration = 0;
	}

	double dep_time(const double arr_time) const {
		const double wait_duration = std::max(0.0, this->timeWindow.openTime - arr_time);
		const double start_of_visit_time = arr_time + wait_duration;
		const double dep_time = start_of_visit_time + this->visitDuration;
		return dep_time;
	}

	double euclidean_distance(TouristAttraction ta) { 
		return (ta.point.pos.lat - point.pos.lat) * (ta.point.pos.lat - point.pos.lat) + (ta.point.pos.lon - point.pos.lon) * (ta.point.pos.lon - point.pos.lon);
	}

	virtual void print() {
		std::cout << "default ta print" << std::endl;
	}

	virtual void printMin(){
		std::cout << "default ta printMin" << std::endl;
	}

	virtual TouristAttraction* clone() const {
		return new TouristAttraction(*this);
	};

	void reset() {
		arrTime = .0;
		waitDuration = .0;
		startOfVisitTime = .0;
		depTime = .0;
		shift = .0;
		maxShift = .0;
	}
	

	virtual ~TouristAttraction() {} //virtual destructor to ensure our subclasses are correctly dealocated
} TA;

struct Sight : TA {
	using TouristAttraction::TouristAttraction;
	

	void print() override {
		std::cout << "id: " << id
			<< "\n point.id: " << point.id
			<< "\n visitDuration: " << visitDuration
			<< "\n profit: " << profit
			<< "\n openTime: " << timeWindow.openTime
			<< "\n closeTime: " << timeWindow.closeTime
			<< "\n arrTime: " << arrTime
			<< "\n waitDuration: " << waitDuration
			<< "\n startOfVisitTime: " << startOfVisitTime
			<< "\n depTime: " << depTime
			<< "\n shift: " << shift
			<< "\n maxShift: " << maxShift
			<< "\n route: " << route
			<< "\n cluster: " << cluster
			<< "\n minDist: " << minDist
			<< std::endl;
	}

	TouristAttraction* clone() const override {
		return new Sight(*this);
	}

};

struct Route : TA {

	Point edges[2];

	Route(std::string pId, Point p, Point p1, Point p2, double pVisitDuration, int pProfit, double pOpenTime, double pCloseTime) :
		TouristAttraction(pId, p, pVisitDuration, pProfit, pOpenTime, pCloseTime), edges{p1, p2}
	{}

	void print() override {
		std::cout << "id: " << id
			<< "\n point.id: " << point.id
			<< "\n edge1.id " << edges[0].id
			<< "\n edge2.id " << edges[1].id
			<< "\n visitDuration: " << visitDuration
			<< "\n profit: " << profit
			<< "\n openTime: " << timeWindow.openTime
			<< "\n closeTime: " << timeWindow.closeTime
			<< "\n arrTime: " << arrTime
			<< "\n waitDuration: " << waitDuration
			<< "\n startOfVisitTime: " << startOfVisitTime
			<< "\n depTime: " << depTime
			<< "\n shift: " << shift
			<< "\n maxShift: " << maxShift
			<< "\n route: " << route
			<< "\n cluster: " << cluster
			<< "\n minDist: " << minDist
			<< std::endl;
	}

	TouristAttraction* clone() const override {
		return new Route(*this);
	}

};

struct Cluster {
	
	int id, openTime, closeTime;
	Point centroid;
	std::vector<TA*> nodes;

	Cluster() : id(DEFAULT_CLUSTER_ID), openTime(-1), closeTime(-1), nodes(){}

	Cluster(int openTime, int closeTime, std::vector<TA*> attractions, double** ttMatrix) :
		id(DEFAULT_CLUSTER_ID), openTime(openTime), closeTime(closeTime), nodes(attractions) {}

	void print() {
		std::cout << 
			"Cluster: " << id << 
			" - openTme: " << openTime <<
			" - closeTime: " << closeTime << std::endl << "Sights: " << std::endl;
		for (auto p = nodes.begin(); p != nodes.end(); ++p) {
			(*p)->print();
		}
	}

	void setCentroid(int id) {
		double totalLat = .0;
		double totalLon = .0;
		size_t nodessize = nodes.size();
		for (auto& node : nodes) {
			totalLat += node->point.pos.lat;
			totalLon += node->point.pos.lon;
		}
		centroid = Point(id, totalLat / nodessize, totalLon / nodessize);
	}
};

