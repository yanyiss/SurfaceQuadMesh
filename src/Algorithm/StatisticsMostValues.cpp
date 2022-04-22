#include "StatisticsMostValues.h"

#include <iostream>

template <typename num_type, bool is_maximum>
void StatisticsMostValues<num_type, is_maximum>::update(int id, num_type x)
{
	queue.emplace(x, id);
	while (queue.size() > max_size) queue.pop();
}

template <typename num_type, bool is_maximum>
void StatisticsMostValues<num_type, is_maximum>::print(const char* title)
{
	std::cout << title << std::endl;
	while (!queue.empty())
	{
		std::cout << queue.top().second << " : " << queue.top().first << std::endl;
		queue.pop();
	}
}

template <typename num_type, bool is_maximum>
std::vector<int> StatisticsMostValues<num_type, is_maximum>::get_ids()
{
	std::vector<int> res;
	res.reserve(queue.size());

	while (!queue.empty())
	{
		res.emplace_back(queue.top().second);
		queue.pop();
	}

	return res;
}

template <typename num_type, bool is_maximum>
std::vector<std::pair<int, num_type>> StatisticsMostValues<num_type, is_maximum>::get_data()
{
	std::vector<std::pair<int, num_type>> res;
	res.reserve(queue.size());

	while (!queue.empty())
	{
		res.emplace_back(queue.top().second, queue.top().first);
		queue.pop();
	}

	return res;
}

template <typename num_type, bool is_maximum>
void StatisticsMostValues<num_type, is_maximum>::get_ids(std::set<int>& id)
{
	while (!queue.empty())
	{
		id.insert(queue.top().second);
		queue.pop();
	}
}

template StatisticsMostValues<double, true>;
template StatisticsMostValues<double, false>;