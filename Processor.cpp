#include "Processor.h"

Processor::Processor() {
	std::cout << "Constructed processor component" << std::endl;
}

Processor::Processor(std::vector<std::vector<Cluster>> clusters, std::vector<std::vector<double>> ttMatrix) {
	m_clusters = clusters;
	m_ttMatrix = ttMatrix;
	std::cout << "Constructed processor component" << std::endl;
}

Processor::~Processor() {
	std::cout << "Deconstructed processor component" << std::endl;
}

std::vector<std::vector<Solution>> Processor::exec() {
	std::vector<std::vector<Solution>> allSolutions;
	for (auto& d_cluster : m_clusters){
		std::vector<Solution> dSolutions;
		for (auto& tw_cluster : d_cluster) {
			TA* cluster_depot = new TA(tw_cluster.centroid); //todo: delete pointer
			/*OP op(tw_cluster.nodes, m_ttMatrix, cluster_depot, tw_cluster.openTime, tw_cluster.closeTime);
			Solution sol = op.solve();
			dSolutions.push_back(sol);
			delete(cluster_depot);*/
		}
		allSolutions.push_back(dSolutions);
	}
	return allSolutions;
}