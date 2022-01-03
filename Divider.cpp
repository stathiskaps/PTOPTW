#include "Divider.h"

Divider::Divider() {
	mWalksNum = 0;
	std::cout << "Constructed divider component" << std::endl;
}

Divider::Divider(std::vector<TA*> attractions, int walks) {
	mAttractions = attractions;
	mWalksNum = walks;
	std::cout << "Constructed divider component" << std::endl;
}

Divider::~Divider() {

}

std::vector<std::string> operator||(const std::vector<std::string>& v1, const std::vector<std::string>& v2) {

	std::set<std::string> all;

	for (auto i : v1) {
		all.insert(i);
	}
	for (auto i : v2) {
		all.insert(i);
	}

	std::vector<std::string> v3(all.begin(), all.end());
	return v3;
}

std::vector<std::string> operator&&(const std::vector<std::string>& v1, const std::vector<std::string>& v2) {
	std::vector<std::string> v3;
	std::set_intersection(v1.begin(), v1.end(),
		v2.begin(), v2.end(),
		back_inserter(v3));
	return v3;
}

//appends vector2 after vector1 (vector1 + vector2)
std::vector<std::string> operator+(std::vector<std::string> v1, std::vector<std::string> v2) {
	v1.insert(v1.end(), v2.begin(), v2.end());
	std::sort(v1.begin(), v1.end());
	return v1;
}

//returns difference of two vectors
std::vector<std::string> operator-(const std::vector<std::string>& v1, const std::vector<std::string>& v2) {
	std::vector<std::string> v3;
	std::set_difference(v1.begin(), v1.end(), v2.begin(), v2.end(),
		std::inserter(v3, v3.end()));
	return v3;
}

//removes v2 elements from v1 (not working)
void operator-=(const std::vector<std::string>& v1, const std::vector<std::string>& v2) {
	std::vector<std::string> v3;
	std::set_difference(v1.begin(), v1.end(), v2.begin(), v2.end(),
		std::inserter(v3, v3.end()));
}

std::vector<std::vector<Cluster>> Divider::exec() {
	std::vector<std::vector<Cluster>> clusters;
	std::vector<Cluster> distance_clusters = kmeans(mAttractions);
	for (auto& dc : distance_clusters) {
		std::vector<Cluster> time_window_clusters = intervals(dc.nodes);
		clusters.push_back(time_window_clusters);
	}
	return clusters;
}

//uses attractions of this cluster instance, to generate distance clusters
std::vector<Cluster> Divider::kmeans(std::vector<TA*> attractions) {

	std::vector<Cluster> clusters(mWalksNum); 
	std::vector<TA*> centroids;
	int attractionsSize = attractions.size();

	//Initialize centroids randomly
	srand(time(0));  // need to set the random seed
	for (int i = 0; i < mWalksNum; ++i) {
		centroids.push_back(attractions.at(rand() % attractionsSize));
	}

	for (std::vector<TA*>::iterator c = centroids.begin(); c != centroids.end(); ++c) {
		// quick hack to get cluster index
		int clusterId = c - std::begin(centroids);
		for (std::vector<TA*>::iterator it = attractions.begin(); it != attractions.end(); ++it) {

			TA* p = *it;
			double dist = (*c)->euclidean_distance(*p);
			if (dist < p->minDist) {

				p->minDist = dist;
				p->cluster = clusterId;

			}
			*it = p;
		}
	}

	for (int i = 0; i < KMEANS_EPOCHS; i++) {

		std::vector<int> nAttractions;
		std::vector<double> sumLat, sumLon;

		// Initialise with zeroes
		for (int j = 0; j < mWalksNum; ++j) {
			nAttractions.push_back(0);
			sumLat.push_back(0.0);
			sumLon.push_back(0.0);
		}

		// Iterate over Sights to append data to centroids
		for (std::vector<TA*>::iterator it = attractions.begin(); it != attractions.end(); ++it) {
			int clusterId = (*it)->cluster;
			nAttractions[clusterId] += 1;
			sumLat[clusterId] += (*it)->point.pos.lat;
			sumLon[clusterId] += (*it)->point.pos.lon;

			(*it)->minDist = DBL_MAX;  // reset distance
		}

		// Compute the new centroids
		for (std::vector<TA*>::iterator c = centroids.begin(); c != centroids.end(); ++c) {
			int clusterId = c - centroids.begin();
			(*c)->point.pos.lat = sumLat[clusterId] / nAttractions[clusterId];
			(*c)->point.pos.lon = sumLon[clusterId] / nAttractions[clusterId];
		}

		for (std::vector<TA*>::iterator c = centroids.begin(); c != centroids.end(); ++c) {
			// quick hack to get cluster index
			int clusterId = c - centroids.begin();

			for (std::vector<TA*>::iterator it = attractions.begin(); it != attractions.end(); ++it) {

				TA* p = *it;
				double dist = (*c)->euclidean_distance(*p);
				if (dist < p->minDist) {
					p->minDist = dist;
					p->cluster = clusterId;
				}
				*it = p;
			}
		}
	}

	// Iterate over Sights and push them to clusters
	for (std::vector<TA*>::iterator it = attractions.begin(); it != attractions.end(); ++it) {
		clusters.at((*it)->cluster).nodes.push_back(*it);
	}

	for (std::vector<Cluster>::iterator c = clusters.begin(); c != clusters.end(); ++c) {
		c->id = c - clusters.begin();
	}

	return clusters;
}

std::vector<std::string> Divider::getActiveNodes(std::vector<TimeEdge>& edges, int leftEdge, int rightEdge) {
	//Active nodes of an interval is the intersection of leftEdge active nodes or open nodes & rightEdge active nodes
	//a.k.a the nodes that when leftEdge.time started were active or opened and kept active until rightEdge.time
	std::vector<std::string> activeNodes = (edges.at(leftEdge).activeNodes || edges.at(leftEdge).openingNodes) && edges.at(rightEdge).activeNodes;
	return activeNodes;
}





//uses attractions of this cluster instance, to generate other clusters
std::vector<Cluster> Divider::intervals(std::vector<TA*> attractions) {

	std::vector<Cluster> twClusters;
	std::map<std::string, int> profits;
	int attractionsTotalProfit = 0;
	int attractionsSize = attractions.size();

	for(const auto &a: attractions){
		profits[a->id] = a->profit;
		attractionsTotalProfit += a->profit;
	}
	

	std::vector<TimeEdge> timeEdges;

	while (attractions.size() != 0) {
		timeEdges = getTimeEdges(attractions);
		std::tuple<int, int> bestInterval = getBestInterval(timeEdges, attractions.size(), profits, attractionsSize);

		if (std::get<0>(bestInterval) != -1 && std::get<1>(bestInterval) != -1) {
			//get active nodes of interval
			std::vector < std::string > activeNodes = getActiveNodes(timeEdges, std::get<0>(bestInterval), std::get<1>(bestInterval));

			uint16_t activeNodesSize = activeNodes.size();
			if (activeNodesSize == 0) {
				std::cout << "No active nodes found" << std::endl;
				break;
			}

			//remove clustered attractions
			uint16_t found = 0;
			Cluster twCluster;
			for (std::vector<TA*>::iterator it = attractions.begin(); it != attractions.end();) {
				if (std::find(activeNodes.begin(), activeNodes.end(), (*it)->id) != activeNodes.end()) {
					twCluster.nodes.push_back(*it);
					it = attractions.erase(it);
					if (++found == activeNodesSize) {
						break;
					}
				}
				else {
					++it;
				}
			}
			twCluster.openTime = timeEdges.at(std::get<0>(bestInterval)).time;
			twCluster.closeTime = timeEdges.at(std::get<1>(bestInterval)).time;
 			twClusters.push_back(twCluster);

		}
	}


	return twClusters;
}

double Divider::calcScore(int intervalProfit, int totalProfit, double intervalDuration, double totalDuration) {
	if (totalProfit == 0 || totalDuration == 0) {
		std::cout << "calcScore: invalid totalProfit or totalDuration inputs" << std::endl;
		return 0;
	}

	double n = intervalProfit / (double)totalProfit;
	double l = intervalDuration / (double)totalDuration;
	return std::abs((n / l) - 1);
}

double Divider::calcScore2(int attractions, int totalAttractions, double intervalDuration, double totalDuration) {
	if (totalAttractions == 0 || totalDuration == 0) {
		std::cout << "calcScore: invalid totalNodesNum or totalDuration inputs" << std::endl;
		return 0;
	}

	std::cout << "attractions: " << attractions << ", totalAttractions: " << totalAttractions << ", intervalDuration: " << intervalDuration << ", totalDuration" << totalDuration << std::endl;

	double n = attractions / (double)totalAttractions;
	double l = intervalDuration / (double)totalDuration;
	return std::abs((n / l) - 1);
}

int collectProfit(std::vector<std::string> attractionIds, std::map<std::string, int>& profits) {
	int totalProfit = 0;
	for (const auto& id : attractionIds)
		totalProfit += profits[id];
	return totalProfit;
}



std::tuple<int, int> Divider::getBestInterval(std::vector<TimeEdge>& timeEdges, int nodesSize, std::map<std::string, int>& profits,int attractionsSize) {
	std::vector<TimeEdge> timeEdgesCopy = timeEdges;
	//double l; //time of interval / total time
	//double n; //active nodes of interval / total nodes
	double bestScore = DBL_MAX; //closer to 1 the better
	double tempScore;
	unsigned leftEdge, rightEdge;
	//double minDiff = DBL_MAX; //lesser the better
	std::tuple<int, int> bestInterval = std::make_tuple(0, 1);
	std::tuple<int, int> tempInterval = std::make_tuple(0, 1);
	//int maxActiveNodes = 0;
	size_t timeEdgesSize = timeEdges.size();
	double totalDuration = timeEdges.back().time - timeEdges.front().time;
	bool foundBetterInterval;
	std::vector<std::string> activeNodes;

	std::cout << "TimeEdges: (" << timeEdges.size() << ")" << std::endl;

	for (auto& edge : timeEdges) {
		std::cout << edge.id << "\t";
	}
	std::cout << std::endl;
	for (auto& edge : timeEdges) {
		std::cout << edge.time << "\t";
	}
	std::cout << std::endl;

	//we search first the small intervals trying to find the one with most active nodes
	for (size_t i = 1; i != timeEdgesSize; ++i) {
		leftEdge = i - 1;
		rightEdge = i;
		activeNodes = getActiveNodes(timeEdges, leftEdge, rightEdge);
		//int intervalProfit = collectProfit(activeNodes, profits);
		double intervalDuration = timeEdges.at(rightEdge).time - timeEdges.at(leftEdge).time;
		//double score = calcScore(intervalProfit, attractionsProfit, intervalDuration, totalDuration);
		double score = calcScore2(activeNodes.size(), attractionsSize, intervalDuration, totalDuration);
		std::cout << "score: " << score << std::endl;
		if (score < bestScore) {
			bestScore = score;
			std::get<0>(bestInterval) = i - 1;
			std::get<1>(bestInterval) = i;
		}
	}

	std::cout << "bestScore: " << bestScore << std::endl;


	leftEdge = std::get<0>(bestInterval);
	rightEdge = std::get<1>(bestInterval);

	std::cout << "leftEdge: " << leftEdge << std::endl;
	std::cout << "rightEdge:" << rightEdge << std::endl;

	std::cout << std::endl;

	while (true) {
		foundBetterInterval = false;
		
		//try decreasing left edge
		if (leftEdge != 0) {
			leftEdge--;
			activeNodes = getActiveNodes(timeEdges, leftEdge, rightEdge);
			//int intervalProfit = collectProfit(activeNodes, profits);
			/*tempScore = calcScore(intervalProfit, attractionsProfit, timeEdges.at(rightEdge).time - timeEdges.at(leftEdge).time, totalDuration);*/
			tempScore = calcScore2(activeNodes.size(), attractionsSize, timeEdges.at(rightEdge).time - timeEdges.at(leftEdge).time, totalDuration);
			std::cout << "left tempScore: " << tempScore << std::endl;
			if (tempScore < bestScore) {
				std::cout << "moving to left improves solution" << std::endl;
				bestScore = tempScore;
				std::get<0>(tempInterval) = leftEdge;
				std::get<1>(tempInterval) = rightEdge;
				foundBetterInterval = true;
			}
			else {
				std::cout << "moving to left doesn't improve solution" << std::endl;
			}
			leftEdge++;
		}

		//try increasing right edge
		if (rightEdge != timeEdgesSize - 1) {
			rightEdge++;
			activeNodes = getActiveNodes(timeEdges, leftEdge, rightEdge);
			//int intervalProfit = collectProfit(activeNodes, profits);
			//tempScore = calcScore(intervalProfit, attractionsProfit, timeEdges.at(rightEdge).time - timeEdges.at(leftEdge).time, totalDuration);
			tempScore = calcScore2(activeNodes.size(), attractionsSize, timeEdges.at(rightEdge).time - timeEdges.at(leftEdge).time, totalDuration);
			if (tempScore < bestScore) {
				std::cout << "mmoving to right improves solution" << std::endl;
				bestScore = tempScore;
				std::get<0>(tempInterval) = leftEdge;
				std::get<1>(tempInterval) = rightEdge;
				foundBetterInterval = true;
			}
			else {
				std::cout << "moving to right doesn't improve the solution" << std::endl;
			}
			rightEdge--;
		}

		//try decreasing left edge &+ increasing right edge
		if (leftEdge != 0 && rightEdge != timeEdgesSize - 1) {
			leftEdge--;
			rightEdge++;
			activeNodes = getActiveNodes(timeEdges, leftEdge, rightEdge);
			//int intervalProfit = collectProfit(activeNodes, profits);
			//tempScore = calcScore(intervalProfit, attractionsProfit, timeEdges.at(rightEdge).time - timeEdges.at(leftEdge).time, totalDuration);
			tempScore = calcScore2(activeNodes.size(), attractionsSize, timeEdges.at(rightEdge).time - timeEdges.at(leftEdge).time, totalDuration);
			if (tempScore < bestScore) {
				std::cout << "moving to left and right improves the solution" << std::endl;
				bestScore = tempScore;
				std::get<0>(tempInterval) = leftEdge;
				std::get<1>(tempInterval) = rightEdge;
				foundBetterInterval = true;
			}
			else {
				std::cout << "moving to left and right doesn't immprove the solution" << std::endl;
			}
			leftEdge++;
			rightEdge--;
		}

		if (!foundBetterInterval) {
			std::cout << "Didn't find any better interval so will break from the loop" << std::endl;
			break;
		}
		else {
			leftEdge = std::get<0>(tempInterval);
			rightEdge = std::get<1>(tempInterval);
		}
	}

	std::get<0>(bestInterval) = leftEdge;
	std::get<1>(bestInterval) = rightEdge;

	return bestInterval;
}

bool Divider::compareInterval(TimeWindowEvent te1, TimeWindowEvent te2) {
	return (te1.time < te2.time);
}


std::vector<TimeEdge> Divider::getTimeEdges(const std::vector<TA*>& attractions) {

	std::vector<TimeWindowEvent> twEvents;
	int index = 0;

	for (const TA* ta : attractions)
	{
		TimeWindowEvent t1(index++, ta->id, ta->timeWindow.openTime, timeWindowEventType::open);
		twEvents.push_back(t1);

		TimeWindowEvent t2(index++, ta->id, ta->timeWindow.closeTime, timeWindowEventType::close);
		twEvents.push_back(t2);
	}

	std::sort(twEvents.begin(), twEvents.end(), compareInterval);

	std::cout << "Printing time window events" << std::endl;
	for (auto& e : twEvents) {
		e.print();
	}

	//construct time edges
	size_t twEventsSize = twEvents.size();
	std::vector<TimeEdge> v;
	std::string currNodeId;
	index = 0; //reset index
	TimeEdge temp = TimeEdge();

	for (std::vector<TimeWindowEvent>::size_type i = 0; i < twEventsSize; ++i) {

		currNodeId = twEvents.at(i).nodeId;
		twEvents.at(i).type == timeWindowEventType::open ? temp.openingNodes.push_back(currNodeId) : temp.closingNodes.push_back(currNodeId);
		temp.events.push_back(twEvents.at(i));

		//if we are not in the last time window event
		if (i != twEventsSize - 1) {
			if (twEvents.at(i).time != twEvents.at(i + 1).time) {
				temp.id = index++;
				temp.time = twEvents.at(i).time;
				v.push_back(temp);
				temp = TimeEdge();
			}
		}
		else {
			temp.id = index;
			temp.time = twEvents.at(i).time;
			v.push_back(temp);
		}
	}


	//sort each opening and closing nodes for each time edge
	for (auto& edge : v) {
		std::sort(edge.openingNodes.begin(), edge.openingNodes.end());
		std::sort(edge.closingNodes.begin(), edge.closingNodes.end());
	}

	size_t timeEdgesSize = v.size();
	//set active nodes of each interval
	//<activeNodes_i = openingNodes_i-1 + activeNodes_i-1 - closingNodes_i-1>
	for (size_t i = 1; i != timeEdgesSize; ++i) {
		v.at(i).activeNodes = v.at(i - 1).activeNodes + v.at(i - 1).openingNodes;
		v.at(i).activeNodes = v.at(i).activeNodes - v.at(i - 1).closingNodes;
	}

	return v;
}