#include "Unifier.h"

Unifier::Unifier() {
	std::cout << "Constructed unifier component" << std::endl;
}

Unifier::Unifier(std::vector<std::vector<Solution>> solutions, int tMax, std::vector<std::vector<double>> ttMatrix, TA* depot) {
	m_solutions = solutions;
	m_tMax = tMax;
	m_ttMatrix = ttMatrix;
	m_depot = depot;
	std::cout << "Constructed unifier component" << std::endl;
}

Unifier::~Unifier() {
	std::cout << "Deconstructed unifier component" << std::endl;
}

std::vector<Solution> Unifier::exec() {
	std::vector<Solution> solutions;
	
	for (auto &dSolution : m_solutions) {
		std::vector<TA*> touristAttractions;
		for (auto &twSolution : dSolution) {
			TA* ta = twSolution.toTA();
			if (ta != nullptr) {
				touristAttractions.push_back(ta);
			}
			
		}

		for (auto& ta : touristAttractions) {
			Sight* s;
			Route* r;

			if ((s = dynamic_cast<Sight*>(ta))) {
				std::cout << "TA " << ta->id << " is a sight" << std::endl;
				s->print();
			}
			else if ((r = dynamic_cast<Route*>(ta))) {
				std::cout << "TA " << ta->id << " is a route" << std::endl;
				r->print();
			}
		}

		OPTW optw(touristAttractions, m_ttMatrix, m_depot, OPEN_DAY_TIME, CLOSE_DAY_TIME);
		Solution sol = optw.solve();
		bool valid = optw.validate();
		std::string msg = valid ? "yes" : "no";
		std::cout << "valid solution? " << msg << std::endl;
		sol.print();

		solutions.push_back(sol);

		//use ilstop to find a solution about these node-routes
		//OPTW ilstop(routes, m_tMax, m_ttMatrix, m_depot, dSol);
		//Solution sol = ilstop.exec();
		//solutions.push_back(sol);
	}
	return solutions;
}