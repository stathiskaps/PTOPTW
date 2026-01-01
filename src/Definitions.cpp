#include "Definitions.h"

bool contains(std::list<std::string>& list, std::string val) {
	return (std::find(list.begin(), list.end(), val) != list.end());
}

bool contains(const std::vector<std::string>& v, std::string val) {
	return (std::find(v.begin(), v.end(), val) != v.end());
}

double GetEuclideanDistance(int x1, int y1, int x2, int y2) {
	double eu_dist;
	double nearest;
	eu_dist = pow(x2 - x1, 2) + pow(y2 - y1, 2);
	eu_dist = sqrt(eu_dist);
	nearest = (double)roundf(eu_dist * 100) / 100;
	return nearest;
}

namespace id {
	std::string generate()
	{
		const size_t length = 3;
		auto randchar = []() -> char
		{
			const char charset[] =
				"0123456789"
				"ABCDEFGHIJKLMNOPQRSTUVWXYZ"
				"abcdefghijklmnopqrstuvwxyz";
			const size_t max_index = (sizeof(charset) - 1);
			return charset[rand() % max_index];
		};
		std::string str(length, 0);
		std::generate_n(str.begin(), length, randchar);

		return str;
	}
}