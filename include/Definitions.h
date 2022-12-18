#include <random>
#include <sstream>
#include <iostream>
#include <ctime>
#include <list>

#define DEFAULT_GENERIC_CLUSTER_ID -10
#define DEFAULT_DISTANCE_CLUSTER_ID -20
#define DEFAULT_TIME_WINDOW_CLUSTER_ID -30
#define MAX_TIMES_NOT_IMPROVED 150
#define KMEANS_EPOCHS 50

//#define DEFAULT_ID -1
//#define DEPOT_ID  -2

#define DEFAULT_CENTROID_ID "c"
#define DEFAULT_TA_ID "default"
#define TARGET_TA_ID "target"
#define NULL_TA_ID "null"
#define DEPOT_ID "depot"
#define START_DEPOT_ID "start_depot"
#define DUMMY_START_DEPOT "dummy_start_depot"
#define END_DEPOT_ID "end_depot"
#define CANDIDATE_END_DEPOT_ID "candidate_end_depot"
#define DUMMY_END_DEPOT_ID "dummy_end_depot"
#define DUMMY_ID "dummy"
#define CNEXT_ID "cnext"
#define NEAREST_NEXT "nearest_next"
#define DEFAULT_POINT_ID -1
#define TARGET_POINT_ID -3
#define DEFAULT_CLUSTER_ID -2
#define DEFAULT_POS -1
#define DEFAULT_WALK_ID -1

#define LEAF -1
#define NOISE -2
#define UNDEFINED -3

#define MAX_MINUTES 1440
#define OPEN_DAY_TIME 0
#define CLOSE_DAY_TIME 1260

bool contains(std::list<std::string>&, std::string);
bool contains(const std::vector<std::string>&, std::string);
double GetEuclideanDistance(int x1, int y1, int x2, int y2);

#ifndef ID_H
#define ID_H

namespace id {
	static std::random_device              rd;
	static std::mt19937                    gen(rd());
	static std::uniform_int_distribution<> dis(0, 20000);

	std::string generate();
}

#endif
