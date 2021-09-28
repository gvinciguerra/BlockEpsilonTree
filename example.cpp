#include <iostream>
#include <random>
#include <vector>
#include "BlockEpsilonTree.hpp"

int main() {
    std::mt19937 engine;
    std::geometric_distribution<uint32_t> distribution(0.85);
    std::vector<uint32_t> data(1000000);
    for (auto i = 1; i < data.size(); ++i)
        data[i] = data[i - 1] + distribution(engine) + 1;

    BlockEpsilonTree bet(data, 2);

    std::cout << "Bits per integer:      " << 8. * bet.size_in_bytes() / data.size() << std::endl
              << "Average depth:         " << bet.get_metadata()["average_depth"] << std::endl
              << "# of elements <= 500:  " << bet.rank(500) << std::endl
              << "10th smallest element: " << bet.select(10) << std::endl;

    return 0;
}
