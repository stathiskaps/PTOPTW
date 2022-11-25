
#include "Controller.h"

Controller::Controller() {
	mOpenTime = 0;
	mCloseTime = 0;
	mWalks = 0;
	std::cout << "Constructed default controller component" << std::endl;
}

Controller::Controller(int openTime, int closeTime, std::vector<TA*> attractions, std::vector<Point> points, int walks) {
	mOpenTime = openTime;
	mCloseTime = closeTime;
	mAttractions = attractions;
	mPoints = points;
	mWalks = walks;
	std::cout << "Constructed controller component" << std::endl;
}

Controller::~Controller() {
	std::cout << "Desconstructed controller component" << std::endl;
}

std::tuple<std::vector<std::vector<double>>, double> calcTravelTimesMatrix(std::vector<Point>& points)
{
	int pointsSize = points.size();
	std::vector<std::vector<double>> ttMatrix;
	double totalTravelTime = 0;
	double meanTravelTime;
	double val;
	for (int i = 0; i < pointsSize; ++i) {
		std::vector<double> vec;
		for (int j = 0; j < pointsSize; j++) {
			val = GetEuclideanDistance(points.at(i).pos.lat, points.at(i).pos.lon, points.at(j).pos.lat, points.at(j).pos.lon);
			vec.push_back(val);
			totalTravelTime += val;
		}
		ttMatrix.push_back(vec);
	}

	meanTravelTime = totalTravelTime / ((double)pointsSize * pointsSize);

	return std::make_tuple(ttMatrix, meanTravelTime);
}

void Controller::exec() {
	
	// TaskManager* taskManager = TaskManager::GetInstance();
	// taskManager->Start();

	// Task task = {SEND_NODES, mAttractions, "Step 1: Parsing nodes"};
	// taskManager->QueueTask(task);


	std::this_thread::sleep_for(std::chrono::milliseconds(3000));

	TA* depot = mAttractions.at(0);
	mAttractions.erase(mAttractions.begin());

	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

	Divider divider(mAttractions, mWalks);
	std::vector<std::vector<Cluster>> clusters = divider.exec();
	
	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
	std::cout << "Divider duration = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;

	begin = std::chrono::steady_clock::now();

	//Add cluster centroid to points
	size_t index = mPoints.size();

	for (auto& d_cluster : clusters) {
		for (auto& tw_cluster : d_cluster) {
			tw_cluster.setCentroid(index++);
			mPoints.push_back(tw_cluster.centroid);
		}
	}

	//calc travel times of all points
	std::tuple<std::vector<std::vector<double>>, double> tuple = calcTravelTimesMatrix(mPoints); //TODO: delete pointer
	std::vector<std::vector<double>> ttMatrix = std::get<0>(tuple);

	for (auto& vec : ttMatrix) {
		for (auto& j : vec) {
			std::cout << j << " ";
		}
		std::cout << std::endl;
	}

	end = std::chrono::steady_clock::now();
	std::cout << "Processor duration = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;

}


