#ifndef SPACE_GRID_H_INCLUDED
#define SPACE_GRID_H_INCLUDED

#include <vector>
#include <map>
class SpaceGrid{
	public:
		void resize(int x, int y, int z, const particleGrid& ps);
		void refillGrid(const particleGrid& ps);

		void pushParticle(const particle& p);
		void getNeighbors(const particleGrid& ps, const particle& p, particleGrid& des);
		

		glm::vec3 positionToGridIndex(glm::vec3 p);
		int gridIndexToVecIndex(glm::vec3 index);
		int positionToVecIndex(glm::vec3 p);
		bool IfWithinBoundry(glm::vec3 gridIndex);
		
		glm::vec3 dim;

		std::vector<std::vector<int>> GridData;


};


#endif