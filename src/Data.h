#ifndef _DATA_H_
#define _DATA_H_

#include <stdio.h>

#include <vector>

class Data {
   private:
	int bin_capacity;
	int n_items;
	std::vector<int> weights;

   public:
	Data() {
		weights = {2, 1, 3, 3, 5};
		bin_capacity = 7;
		n_items = weights.size();
	}
	void readData(char* filePath);

	int getNItems();

	int getBinCapacity();

	int getItemWeight(unsigned int item);
};

#endif
