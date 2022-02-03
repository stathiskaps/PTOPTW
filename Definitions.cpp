#include "Definitions.h"

bool contains(std::list<std::string>& list, std::string val) {
	return (std::find(list.begin(), list.end(), val) != list.end());
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