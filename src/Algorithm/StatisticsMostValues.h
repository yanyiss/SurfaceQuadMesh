#pragma once

#include <set>
#include <queue>

template <typename num_type, bool is_maximum>
class StatisticsMostValues
{
public:
	StatisticsMostValues(int size) : max_size(size) {};
	~StatisticsMostValues() = default;

	int size() { return queue.size(); }
	int capacity() { return max_size; }

	void update(int id, num_type x);

	void print(const char* title);
	std::vector<int> get_ids();
	std::vector<std::pair<int, num_type>> get_data();

	void get_ids(std::set<int>& id);

private:
	using datum_type = std::pair<num_type, int>;
	using compare_type = std::conditional_t<is_maximum, std::greater<datum_type>, std::less<datum_type>>;

	std::priority_queue<datum_type, std::vector<datum_type>, compare_type> queue;

	const int max_size;
};
