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
#define DEFAULT_TA_ID ""
#define DEFAULT_DEPOT_ID "depot"
#define DEFAULT_POINT_ID -1
#define DEFAULT_CLUSTER_ID -2
#define DEFAULT_POS -1

#define LEAF -1
#define NOISE -2
#define UNDEFINED -3


#define MAX_MINUTES 1440
#define OPEN_DAY_TIME 0
#define CLOSE_DAY_TIME 1260

#ifndef S_H
#define S_H

class S {
public:
	static int ta_id;
};

#endif
//namespace uuid {
//	static std::random_device              rd;
//	static std::mt19937                    gen(rd());
//	static std::uniform_int_distribution<> dis(0, 9);
//	static std::uniform_int_distribution<> dis2(8, 11);
//
//	//std::string generate_uuid_v3() {
//	//	std::stringstream ss;
//	//	int i;
//	//	ss << std::hex;
//	//	for (i = 0; i < 8; i++) {
//	//		ss << dis(gen);
//	//	}
//	//	ss << "-";
//	//	for (i = 0; i < 4; i++) {
//	//		ss << dis(gen);
//	//	}
//	//	ss << "-4";
//	//	for (i = 0; i < 3; i++) {
//	//		ss << dis(gen);
//	//	}
//	//	ss << "-";
//	//	ss << dis2(gen);
//	//	for (i = 0; i < 3; i++) {
//	//		ss << dis(gen);
//	//	}
//	//	ss << "-";
//	//	for (i = 0; i < 12; i++) {
//	//		ss << dis(gen);
//	//	};
//	//	return ss.str();
//	//}
//}

bool contains(std::list<std::string>&, std::string);
bool contains(const std::vector<std::string>&, std::string);
#ifndef ID_H
#define ID_H

namespace id {
	static std::random_device              rd;
	static std::mt19937                    gen(rd());
	static std::uniform_int_distribution<> dis(0, 20000);

	std::string generate();
}

#endif
