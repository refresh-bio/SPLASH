#ifndef _HELMERT_DECOMPOSITION_H
#define _HELMERT_DECOMPOSITION_H

#include <vector>
#include <numeric>

//mkokot_TODO: consider Suppose for a fixed anchor, we observe >k=5 targets in M cell types with more than K cells (let’s take K=10).
class HelmertDecomposition {
	std::vector<double> pi;
	std::vector<double> P;
public:
	//n is a vector that for each cell type contains the number of columns in contignency table that are of this cell type
	//cell type id is this vector index
	HelmertDecomposition(const std::vector<uint32_t>& n) : pi(n.size()), P(n.size()) {
		auto sum = std::accumulate(n.begin(), n.end(), 0ul);

		pi[0] = (double)n[0] / sum;
		P[0] = pi[0];
		//mkokot_TODO: consider The order of the cells can be arbitrary, or could be sorted based on the number of cells in that from highest to lowest
		for (size_t i = 1; i < pi.size(); ++i) {
			pi[i] = (double)n[i] / sum;
			P[i] = P[i - 1] + pi[i]; 
		}
	}

	void get_row(uint32_t i, std::vector<double>& row) {
		row.resize(P.size());
		if (i == 0) {
			for (size_t j = 0; j < row.size(); ++j)
				row[j] = sqrt(pi[j]);
			return;
		}

		double precomputed = pi[i] / (P[i - 1] * P[i]);

		for (size_t j = 0; j < i; ++j)
			row[j] = sqrt(pi[j] * precomputed);
		
		row[i] = -sqrt(P[i - 1] / P[i]);
		for (size_t j = i + 1; j < row.size(); ++j)
			row[j] = 0;
	}
};

#endif // !_HELMERT_DECOMPOSITION_H
