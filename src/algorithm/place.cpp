/*!
 *  * Write **your own** code.
 *   * !!! Cheating will be strictly not accepted. !!!
 *    * If cheating is detected by the similarity check program and TA determine that you have cheated,
 *     * then you will get F grade or zero point for this term project.
 *      * You can use external libraries for only pure math libraries; i.e) fft, sparse matrix, etc
 *       * If you want to use external library, then please check whether it is okay by contact to TA.
 *        * */

#include "Circuit.h"
#include "matrixSolver.h"

namespace Placer {

// class for managing the Instances
// assigns each Instance an easily accessible index
/**
 * @class InstanceManager
 * @brief Manages the mapping between instance names and their corresponding indices.
 *
 * This class provides a bidirectional mapping between instance names and their
 * indices in a vector. It is designed to allow for efficient lookup of an instance's
 * index given its name, and vice versa. This can be particularly useful in scenarios
 * where instances need to be quickly accessed by name, or ordered by their index.
 *
 * @note The class assumes that the vector of instance pointers is correctly populated
 * and that each instance has a unique name. The constructor initializes the internal
 * data structures based on the provided vector, and the mapping is assumed to remain
 * valid for the lifetime of the `InstanceManager` object.
 *
 * Usage example:
 * @code
 *     std::vector<Instance*> instances = { ... };
 *     InstanceManager manager(instances);
 *     int index = manager.getIndexByName("instance_name");
 *     std::string name = manager.getNameByIndex(index);
 * @endcode
 */
class InstanceManager {
	public:
		/**
		 * @brief Constructs an `InstanceManager` and initializes the mapping from the provided vector.
		 * @param instance_pointers A vector of pointers to `Instance` objects to be managed.
		 */
		InstanceManager(const std::vector<Instance*>& instance_pointers)
			: instance_pointers_(instance_pointers) {
			
			// Reserving space for the vector to optimize performance
			index_to_instance_name_.reserve(instance_pointers.size());
			
			// Populate the map and vector
			for (int i = 0; i < instance_pointers.size(); ++i) {
				if (instance_pointers[i] != nullptr) { // Ensure the pointer is valid
					std::string instance_name = instance_pointers[i]->getName();
					
					// Map from name to index
					instance_name_to_index_[instance_name] = i;
					
					// Map from index to name
					index_to_instance_name_.push_back(instance_name);
				} else {
					// Handle nullptr if necessary, e.g., by inserting a placeholder or throwing an error
					index_to_instance_name_.push_back(""); // Placeholder for invalid pointer
				}
			}
		}
		/**
		 * @brief Retrieves the index of an instance given its name.
		 * @param name The name of the instance to find.
		 * @return The index of the instance, or -1 if the name is not found.
		 */		
		int getIndexByName(const std::string& name) const {
			auto it = instance_name_to_index_.find(name);
			if (it != instance_name_to_index_.end()) {
				return it->second;
			}
			// Handle the case where the name is not found
			return -1;
		}

		/**
		 * @brief Retrieves the name of an instance given its index.
		 * @param index The index of the instance to find.
		 * @return The name of the instance, or an empty string if the index is out of bounds.
		 */
		std::string getNameByIndex(int index) const {
			if (index >= 0 && index < index_to_instance_name_.size()) {
				return index_to_instance_name_[index];
			}
			// Handle the case where the index is out of bounds
			return "";
		}
	private:
		std::vector<Instance*> instance_pointers_; ///< The original vector of instance pointers.
		std::unordered_map<std::string, int> instance_name_to_index_; ///< Maps instance names to their indices.
		std::vector<std::string> index_to_instance_name_; ///< Maps indices to instance names.
};


class GridManager {
public:
    static const int gridSize = 20;
	int num_instances;
	int die_width;	// max column
	int die_height;	// max row
	int region_area;
	int region_row_width;
	int region_col_height;
	std::array<std::array<std::vector<int>, gridSize>, gridSize> grid;
	std::array<std::array<int, gridSize>, gridSize> gridCellArea;
	std::vector<int> cellAreaIndexStore; // stores cell area for each instance index 
public:
    GridManager(const std::vector<Instance*>& instance_pointers_, Die* die_)
	{
		num_instances = instance_pointers_.size();
		die_width = die_->getWidth();	// max column
		die_height = die_->getHeight();	// max row
		region_area = die_width * die_height / (gridSize * gridSize);
		region_row_width = die_width / gridSize;
		region_col_height = die_height / gridSize;

		cellAreaIndexStore.assign(num_instances, 0);
		for (int i = 0; i < num_instances; i++)
		{
			std::pair<int,int> instance_coordinate = instance_pointers_[i]->getCoordinate();
			int gridX = std::min(gridSize - 1, std::max(0, static_cast<int>(instance_coordinate.first / region_row_width)));
			int gridY = std::min(gridSize - 1, std::max(0, static_cast<int>(instance_coordinate.second / region_col_height)));
			grid[gridY][gridX].push_back(i);
			cellAreaIndexStore[i] = instance_pointers_[i]->getArea();
		}
		calculateAreaInRegion(0, 0, gridSize - 1, gridSize - 1);

    }

    // Method to calculate the area of cells in a region
    int calculateAreaInRegion(int row1, int col1, int row2, int col2) {
        int areaSum = 0;
        for (int i = row1; i <= row2; ++i) {
            for (int j = col1; j <= col2; ++j) {
				for (int k = 0; k < grid[i][j].size(); k++)
				{
					areaSum += cellAreaIndexStore[grid[i][j][k]];
				}
                
            }
        }
        return areaSum;
    }

    // Method to move a cell and update area sums
    void moveCell(int cellIndex, int fromRow, int fromCol, int toRow, int toCol) {
        // Move cell from grid region to new grid region
		auto it = std::find(grid[fromRow][fromCol].begin(), grid[fromRow][fromCol].end(), cellIndex);

		if (it != grid[fromRow][fromCol].end())
		{
			grid[fromRow][fromCol].erase(it);
			grid[toRow][toCol].push_back(cellIndex);
			updateAreaSums(fromRow, fromCol, -cellAreaIndexStore[cellIndex]);
			updateAreaSums(toRow, toCol, cellAreaIndexStore[cellIndex]);
		}
		else
		{
			//std::cout << "Cell not found in grid region" << std::endl;
		}

    }

    std::pair<int, int> getRegionCenter(int i, int j) {
        double regionWidth = die_width / gridSize;
        double regionHeight = die_height / gridSize;

        double centerX = (i * regionWidth) + (regionWidth / 2);
        double centerY = (j * regionHeight) + (regionHeight / 2);

		// change to int
		int intcenterX = std::min(die_width, std::max(0, static_cast<int>(centerX)));
		int intcenterY = std::min(die_height, std::max(0, static_cast<int>(centerY)));
        return {intcenterX, intcenterY};
    }

	// Calculate the maximum overshoot area for the grid
	int calculateMaxOvershootArea() {
		int maxOvershootArea = 0;
		for (int i = 0; i < gridSize; ++i) {
			for (int j = 0; j < gridSize; ++j) {
				int overshootArea = calculateOvershootArea(i, j);
				if (overshootArea > maxOvershootArea) {
					maxOvershootArea = overshootArea;
				}
			}
		}
		return maxOvershootArea;
	}

	// Calculate median overshoot area for the grid
	int calculateMedianOvershootArea() {
		std::vector<int> overshootAreas;
		for (int i = 0; i < gridSize; ++i) {
			for (int j = 0; j < gridSize; ++j) {
				overshootAreas.push_back(calculateOvershootArea(i, j));
			}
		}
		std::sort(overshootAreas.begin(), overshootAreas.end());
		return overshootAreas[overshootAreas.size() / 2];
	}

	// Calculate the overshoot area for a region
	int calculateOvershootArea(int i, int j)
	{
		int regionArea = die_width * die_height / (gridSize * gridSize);
		int regionWidth = die_width / gridSize;
		int regionHeight = die_height / gridSize;
		int regionAreaSum = calculateAreaInRegion(i, j, i, j);
		int overshootArea = regionAreaSum - regionArea;
		if (overshootArea < 0)
		{
			overshootArea = 0;
		}
		return overshootArea;
	}

	// Calculate the overshoot density for a region
	double calculateOvershootDensity(int i, int j)
	{
		int regionArea = die_width * die_height / (gridSize * gridSize);

		int regionAreaSum = calculateAreaInRegion(i, j, i, j);
		// cast to double to avoid integer division
		double overshootDensity = static_cast<double>(regionAreaSum ) / regionArea;
		if (overshootDensity < 0)
		{
			overshootDensity = 1.0;
		}
		return overshootDensity;
	}

	// Calculate the average overshoot area for all regions
	double calculateAverageOvershootArea() {
		double averageOvershootArea = 0.0;
		for (int i = 0; i < gridSize; ++i) {
			for (int j = 0; j < gridSize; ++j) {
				averageOvershootArea += calculateOvershootArea(i, j);
			}
		}
		return averageOvershootArea / (gridSize * gridSize);
	}
	// Proportion of regions with overshoot area > 0
	double calculateProportionOvershootArea() {
		double proportionOvershootArea = 0.0;
		for (int i = 0; i < gridSize; ++i) {
			for (int j = 0; j < gridSize; ++j) {
				if (calculateOvershootArea(i, j) > 0) {
					proportionOvershootArea += 1.0;
				}
			}
		}
		return proportionOvershootArea / (gridSize * gridSize);
	}
	// Average overshoot density
	double calculateAverageOvershootDensity() {
		double averageOvershootDensity = 0.0;
		for (int i = 0; i < gridSize; ++i) {
			for (int j = 0; j < gridSize; ++j) {
				averageOvershootDensity += calculateOvershootDensity(i, j);
			}
		}
		return averageOvershootDensity / (gridSize * gridSize);
	}
	// Find the region given instance index
	std::pair<int, int> findRegion(int instanceIndex) {
		for (int i = 0; i < gridSize; ++i) {
			for (int j = 0; j < gridSize; ++j) {
				auto it = std::find(grid[i][j].begin(), grid[i][j].end(), instanceIndex);
				if (it != grid[i][j].end())
				{
					return { i, j };
				}
			}
		}
		return { -1, -1 };
	}
	// Find region given location
	std::pair<int, int> findRegion(std::pair<int, int> location) {
		int regionWidth = die_width / gridSize;
		int regionHeight = die_height / gridSize;
		int gridX = std::min(gridSize - 1, std::max(0, static_cast<int>(location.first / regionWidth)));
		int gridY = std::min(gridSize - 1, std::max(0, static_cast<int>(location.second / regionHeight)));
		return { gridX, gridY };
	}

private:
    // Helper method to update area sums
    void updateAreaSums(int x, int y, int valueChange) {
        gridCellArea[x][y] += valueChange;
    }
};

void Circuit::quadraticPlacement(const string &output_path_name, const string &defName) {
	// matrix solve example.
	// You should refer below function when you get x respect to Ax = b
	// Below function is implemented in src/algorithm/math/matrixSolver.cpp
	// solve_example();  // You should erase this line.
	//this->placeExample();

	// Make hyper graph of circuit that treats a net as a node
	// From hyper graph form clique graph
	cout << "flag1" << endl;
	// Initiate a helper class for managing the instances and their indices
	InstanceManager instance_manager(instance_pointers_);

	// Adjacency list of clique graph
	/**
	 * @brief A data structure representing an adjacency list for a weighted graph.
	 *
	 * This data structure is composed of a vector, where each element is a list of pairs:
	 * - std::vector<...>: The outer structure, a dynamic array with elements that can be accessed via an index.
	 *   Each index corresponds to the instance number, which is directly related to the position of the instance
	 *   within the 'instance_pointers_' vector.
	 *
	 * - std::list<std::pair<int, int>>: Each element of the vector is a list, where each list corresponds to
	 *   all the other nodes that are connected to the node represented by the vector's index. A std::list is
	 *   chosen here for efficient insertions and deletions as the graph changes.
	 *
	 * - std::pair<int, int>: Inside each list, we have pairs. Each pair represents an edge in the graph,
	 *   where 'first' is an integer representing the index of the connected node and 'second' is an integer
	 *   representing the weight of the edge (e.g., the length of the net connecting two instances/cells).
	 *
	 * Example:
	 *   If instance 0 (from 'instance_pointers_' vector) is connected to instance 2 with a net of weight 5, then
	 *   the list at index 0 of the adjacency vector will contain a pair where the 'first' element of the pair is 2
	 *   (the index of the instance it's connected to) and the 'second' element is 5 (the weight of the connection).
	 *
	 * This structure is especially efficient for sparse graphs, where each node is not connected to many other nodes,
	 * thus saving space by only storing the connections that do exist. The instance number in the adjacency list
	 * corresponds to the index number of the instance in the 'instance_pointers_' vector, maintaining a consistent
	 * reference system across the data structure.
	 */




	int num_instances = instance_pointers_.size();

	const int INSTANCE_THRESHOLD = 1000; // Threshold for parallelization
	const int MAX_THREADS = (num_instances > 200000) ? 16 : 32; // Maximum number of threads, reduce number of threads when instance is high to limit memory use


	omp_set_num_threads(MAX_THREADS); // Set max threads for large data size


	std::vector<std::list<std::pair<int, int>>> adjacency_list(instance_pointers_.size());
	// Create data structure to store the weight of pad wire for each instance
	std::vector<int> pad_wire_weights(instance_pointers_.size(), 0);
	// Create b vector
	// bx vector:
	// If gate i connects to a pad at (xi, yi) with a wire with weight wi,
	// then set bx[i] = wi * xi
	std::valarray<double> bx_vector(num_instances);
	std::valarray<double> by_vector(num_instances);

	std::cout << "flag2" << endl;
	int test_idx = 0;

	// Initiate the adjacency list with the number of instances in the circuit
	// Also initiate the weight of pad wire for each instance
	// Checking structure:
	// instance -> pin -> net -> pin -> instance
	// Conditionally use parallelism for the outer loop (instances)
	if (num_instances > INSTANCE_THRESHOLD)
	{
		#pragma omp parallel
		{
			// privatization of variables for each thread, coalesced later
			std::vector<std::list<std::pair<int, int>>> local_adjacency_list(num_instances);
			std::vector<int> local_pad_wire_weights(num_instances, 0);
			std::valarray<double> local_bx_vector(0.0, num_instances);
			std::valarray<double> local_by_vector(0.0, num_instances);

				#pragma omp for nowait
				for (int idx = 0; idx < num_instances; idx++)
				{
					Instance *instance = instance_pointers_[idx];
					int instance_index = idx;
					// check if instance is nullptr
					if (!instance)
					{
						continue;
					}
					
					// Print out instance_index every 10000 instances to check on progress
					if (idx % 10000 == 0)
					{
						std::cout << "instance_index: " << idx << std::endl;
					}

					// Temporary map to store connections and their cumulative weights
					std::unordered_map<int, int> connection_weights;
					for (Pin *pin : instance->getPins())
					{
						// check if pin is nullptr
						if (!pin)
						{
							continue;
						}
						//cout << "pin " << endl;
						//cout << pin->getPinName() << endl;
						// Get the net connected to the pin
						Net *net = pin->getNet();
						// If the net is not connected to anything, skip it
						if (net == nullptr)
						{
							continue;
						}
						int current_net_weight = net->getWeight();

						// Iterate through the pins connected to the net
						for (Pin *connected_pin: net ->getConnectedPins())
						{
							if (!connected_pin)
							{
								continue;
							}
							// Check if the connected_pin is connected to an instance
							if(connected_pin-> getInstance())
							{
								//cout << "connected_pin " << connected_pin->getPinName() << endl;
								// Get the instance connected to the pin
								Instance *connected_instance = connected_pin->getInstance();
								// Skip self-pointing edges or edges to nullptr
								if (!connected_instance || connected_instance == instance)
								{
									continue;
								}
								// Get the index of the connected instance
								int connected_instance_index = instance_manager.getIndexByName(connected_instance->getName());
								
								// Accumulate weights for duplicate connections
								connection_weights[connected_instance_index] += current_net_weight;
								/*
								// Add the connection to the adjacency list
								adjacency_list[instance_manager.getIndexByName(instance->getName())].push_back(std::make_pair(connected_instance_index, current_net_weight));
								*/
							}
							// Check if the connected_pin is a block pin (pad)
							if (connected_pin->isBlockPin())
							{
								// Accumulate net weight for block pins
								local_pad_wire_weights[instance_index] += current_net_weight;

								std::pair<int, int> pin_coordinate = connected_pin->getCoordinate();
								// Accumulate net weight * x/y coordinate for block pins
								local_bx_vector[instance_index] += current_net_weight * pin_coordinate.first;
								local_by_vector[instance_index] += current_net_weight * pin_coordinate.second;							}
						}
					}
					// Add the accumulated connections to the adjacency list
					for (const std::pair<const int, int>& connection : connection_weights)
					{
						local_adjacency_list[instance_index].push_back(connection);
					}

				}
			

			// Critical section for merging results
			#pragma omp critical
			{
				for (int i = 0; i < num_instances; i++)
				{
					adjacency_list[i].insert(adjacency_list[i].end(), local_adjacency_list[i].begin(), local_adjacency_list[i].end());
					pad_wire_weights[i] += local_pad_wire_weights[i];
					bx_vector[i] += local_bx_vector[i];
					by_vector[i] += local_by_vector[i];

				}
			}


		}

	}
	else	// Serial execution for a small number of instances
	{
		for (Instance *instance : instance_pointers_)
		{
			// Iterate through pins of instance
			// From the pin get the corresponding net
			// From the net get the connected instance and weight
			// Use the information to update the adjacency list

			// check if instance is nullptr
			if (!instance)
			{
				continue;
			}
			// int instance_index = instance_manager.getIndexByName(instance->getName());
			int instance_index = test_idx;

			// Print out instance_index every 10000 instances to check on progress
			if (instance_index % 10000 == 0)
			{
				std::cout << "instance_index: " << instance_index << std::endl;
			}

			//std::cout << instance->getName() << endl;
			// std::cout << "instance_index: " << instance_index << " Test Index" << test_idx << endl;
			// assert (instance_index == test_idx);
			// Temporary map to store connections and their cumulative weights
			std::unordered_map<int, int> connection_weights;
			for (Pin *pin : instance->getPins())
			{
				// check if pin is nullptr
				if (!pin)
				{
					continue;
				}
				//cout << "pin " << endl;
				//cout << pin->getPinName() << endl;
				// Get the net connected to the pin
				Net *net = pin->getNet();
				// If the net is not connected to anything, skip it
				if (net == nullptr)
				{
					continue;
				}
				int current_net_weight = net->getWeight();

				// Iterate through the pins connected to the net
				for (Pin *connected_pin: net ->getConnectedPins())
				{
					if (!connected_pin)
					{
						continue;
					}
					// Check if the connected_pin is connected to an instance
					if(connected_pin-> getInstance())
					{
						//cout << "connected_pin " << connected_pin->getPinName() << endl;
						// Get the instance connected to the pin
						Instance *connected_instance = connected_pin->getInstance();
						// Skip self-pointing edges or edges to nullptr
						if (!connected_instance || connected_instance == instance)
						{
							continue;
						}
						// Get the index of the connected instance
						int connected_instance_index = instance_manager.getIndexByName(connected_instance->getName());
						
						// Accumulate weights for duplicate connections
						connection_weights[connected_instance_index] += current_net_weight;
						/*
						// Add the connection to the adjacency list
						adjacency_list[instance_manager.getIndexByName(instance->getName())].push_back(std::make_pair(connected_instance_index, current_net_weight));
						*/
					}
					// Check if the connected_pin is a block pin (pad)
					if (connected_pin->isBlockPin())
					{
						// Accumulate net weight for block pins
						pad_wire_weights[instance_index] += current_net_weight;

						std::pair<int, int> pin_coordinate = connected_pin->getCoordinate();
						// Accumulate net weight * x/y coordinate for block pins
						bx_vector[instance_index] += current_net_weight * pin_coordinate.first;
						by_vector[instance_index] += current_net_weight * pin_coordinate.second;
					}
				}
			}
			// Add the accumulated connections to the adjacency list
			for (const std::pair<const int, int>& connection : connection_weights)
			{
				adjacency_list[instance_index].push_back(connection);
			}
			test_idx++;
		}
	}

	std::cout << "flag3" << endl;
	//cout << endl;

	// Using the adjacency list and the pad wire weights, populate the A matrix
	// A matrix is where the ith diagonal is the sum of the weights of the connected pad nets that is connected to the ith instance and the sum of ith row of C matrix
	// For off diagonal elements, the value is the negative of the C matrix
	int non_zero_count = 0;
	// count non zero elements in adjacency list: non-diagonal elements
	for (const auto& instance_connections : adjacency_list)
	{
		non_zero_count += instance_connections.size();
	}
	// count non zero elements for (pad wire weights + row sum): diagonal elements
	// check if the sum of pad wire weights with row sum is zero
	for (int i = 0; i < pad_wire_weights.size(); i++)
	{
		// Sum of weights for the ith row in the adjacency list
		int row_sum = 0;

		// Calculate the sum of weights for the ith row in the adjacency list
		for (const auto& pair : adjacency_list[i])
		{
			row_sum += pair.second;
		}

		// Calculate the sum of the apd wire weight and the row sum
		int diagonal_sum = pad_wire_weights[i] + row_sum;
		// Check if the diagonal value is nonzero
		if (diagonal_sum != 0)
		{
			non_zero_count++;
		}
	}
	coo_matrix A_matrix;
	A_matrix.n = adjacency_list.size();
	A_matrix.nnz = non_zero_count;
	A_matrix.row.resize(non_zero_count);
	A_matrix.col.resize(non_zero_count);
	A_matrix.dat.resize(non_zero_count);

	std::cout << "flag4" << std::endl;

	std::valarray<int> A_row(non_zero_count);
	std::valarray<int> A_col(non_zero_count);
	std::valarray<double> A_dat(non_zero_count);
	// Populate the A matrix
	int idx = 0;
    for (int i = 0; i < adjacency_list.size(); i++) 
	{
        int sum_row_C = 0;

        // Populate off-diagonal elements (negative of C)
        for (const auto& pair : adjacency_list[i]) 
		{
			A_row[idx] = i;
			A_col[idx] = pair.first;
			A_dat[idx] = -static_cast<double>(pair.second);
            sum_row_C += pair.second;
            idx++;
        }

        // Add diagonal element
		A_row[idx] = i;
		A_col[idx] = i;
		A_dat[idx] = static_cast<double>(pad_wire_weights[i] + sum_row_C);

        idx++;
    }
	// Populate the A matrix
	A_matrix.row = A_row;
	A_matrix.col = A_col;
	A_matrix.dat = A_dat;
	// Change Parallelization parameters
	A_matrix.matrix_threshold = INSTANCE_THRESHOLD;
	A_matrix.matrix_max_num_threads = MAX_THREADS;
	// Compute preconditioner
	A_matrix.compute_preconditioner();

	std::cout << "flag5" << std::endl;

	//std::cout << "flag before matrix solve" << std::endl;
	// Solve Ax = b for bx and by
	std::valarray<double> x_vector(A_matrix.n);
	std::valarray<double> y_vector(A_matrix.n);
	// initialize x and y vector to some random value
	for (int i = 0; i < A_matrix.n; i++)
	{
		x_vector[i] = (double) random() / (double) RAND_MAX;
		y_vector[i] = (double) random() / (double) RAND_MAX;
	}

	std::cout << "flag6: before matrix solve" << std::endl;


	A_matrix.solve(bx_vector, x_vector);

	std::cout << "flag7: x solved" << std::endl;


	A_matrix.solve(by_vector, y_vector);

	std::cout << "flag8: y solved" << std::endl;


	// Set Coordinate (place) each instance/cell
	int die_width = die_->getWidth();
	int die_height = die_->getHeight();
	for (Instance *instance : instance_pointers_) 
	{
		double raw_x = x_vector[instance_manager.getIndexByName(instance->getName())];
		double raw_y = y_vector[instance_manager.getIndexByName(instance->getName())];
		// Make sure that the double raw_x, raw_y are int and fit within die
		int x = std::min(die_width, std::max(0, static_cast<int>(raw_x)));
		int y = std::min(die_height, std::max(0, static_cast<int>(raw_y)));
		instance->setCoordinate(x, y);
    }


/*
// Used to illustrate connection relationships between instances, pins, and nets.
// Does similar procedure to adjacency list but prints out the information instead of storing it.
// Starts from instance and does a sort of single layer tree traversal 
// to get to the connected instances.
for (Instance *instance : instance_pointers_) {
    if (!instance) continue; // Ensure the instance pointer is valid
    cout << "Instance: " << instance->getName() << endl;
    
    for (Pin *pin : instance->getPins()) {
        if (!pin) continue; // Ensure the pin pointer is valid
        cout << "\tPin: " << pin->getPinName() << endl;

        Net *net = pin->getNet();
        if (!net) continue; // Ensure the net pointer is valid
        cout << "\t\tNet: " << net->getName() << ", Weight: " << net->getWeight() << endl;

        // Iterate through the pins connected to the net
        for (Pin *connected_pin : net->getConnectedPins()) {
            if (!connected_pin) continue; // Ensure the connected_pin pointer is valid
            cout << "\t\t\tConnected Pin: " << connected_pin->getPinName() << endl;

            // Check if the connected_pin is connected to an instance
            Instance *connected_instance = connected_pin->getInstance();
            if (!connected_instance) continue; // Ensure the connected_instance pointer is valid

            cout << "\t\t\t\tConnected Instance: " << connected_instance->getName() << endl;
        }
    }
}


	cout << endl;


	// Used to print adjacency list
    for (size_t i = 0; i < adjacency_list.size(); ++i) 
	{
        std::cout << "Instance " << i << " is connected to:\n";
        for (const auto& edge : adjacency_list[i]) 
		{
            std::cout << "  Instance " << edge.first << " with weight " << edge.second << '\n';
        }
    }

	// From clique graph form a matrix(C, A, and b)
	// Solve Ax = b
	// Print the HPWL of placeExample

	cout <<"Access all Inances (cells) in the circuit" << endl;
	cout << "Number of instances: " << instance_pointers_.size() << endl;
	for (Instance *instance : instance_pointers_) {
		string cell_name = instance->getName();
		cout << cell_name << " ";
	}
	cout << endl << endl;
	// Do the same for nets
	cout << "Access all nets in the circuit" << endl;
	cout << "Number of nets: " << net_pointers_.size() << endl;
	for (Net *net : net_pointers_) {
		string net_name = net->getName();
		cout << net_name << " ";
	}
	cout << endl << endl;
	// Do the same for pins
	cout << "Access all pins in the circuit" << endl;
	cout << "Number of pins: " << pin_pointers_.size() << endl;
	for (Pin *pin : pin_pointers_) {
		string pin_name = pin->getPinName();
		cout << pin_name << " ";
	}
	cout << endl << endl;
	// Do the same for pads
	cout << "Access all pads in the circuit" << endl;
	cout << "Number of pads: " << pad_pointers_.size() << endl;
	for (Pin *pad : pad_pointers_) {
		string pad_name = pad->getPinName();
		cout << pad_name << " ";
	}
	cout << endl << endl;


*/

	// added as for large files memory error seems to happen when function exits
	// add a backup policy to save error before potential errors occur
  	string img_output_file_name = "qPlace_result_" + defName;
  	saveImg(output_path_name, img_output_file_name);


	saveImg("../output/qPlacer/", "FooImg");
	write("../output/qPlacer/Foo.def");
	
	write(output_path_name + defName);



	std::cout << "HPWL of placeExample: " << std::endl;
	std::cout << this->getHPWL() << std::endl;




 }
void Circuit::myPlacement() {
	// do random placement first
	this->placeExample();  // You should erase this line.

	std::cout << "After Initial Placement: " << std::endl;
	std::cout << "HPWL of placeExample: " << std::endl;
	std::cout << this->getHPWL() << std::endl;




	const int gridSize = 20;
	InstanceManager instance_manager(instance_pointers_);
	int num_instances = instance_pointers_.size();
	int die_width = die_->getWidth();	// max column
	int die_height = die_->getHeight();	// max row
	int region_area = die_width * die_height / (gridSize * gridSize);
	int region_row_width = die_width / gridSize;
	int region_col_height = die_height / gridSize;


	GridManager grid(instance_pointers_, die_);


	// print the average overshoot area, median overshoot area, and max overshoot area
	std::cout << "Average overshoot area: " << std::endl;
	std::cout << grid.calculateAverageOvershootArea() << std::endl;
	std::cout << "Median overshoot area: " << std::endl;
	std::cout << grid.calculateMedianOvershootArea() << std::endl;
	std::cout << "Max overshoot area: " << std::endl;
	std::cout << grid.calculateMaxOvershootArea() << std::endl;
	std::cout << "Proportion of regions with overshoot area > 0: " << std::endl;
	std::cout << grid.calculateProportionOvershootArea() << std::endl;
	std::cout << "Average overshoot density: " << std::endl;
	std::cout << grid.calculateAverageOvershootDensity() << std::endl;

	ulong initial_HPWL = this->getHPWL();
	double initial_Average_Overshoot = grid.calculateAverageOvershootArea();
	double initial_Median_Overshoot = grid.calculateMedianOvershootArea();
	double initial_Max_Overshoot = grid.calculateMaxOvershootArea();
	double initial_Proportion_Overshoot = grid.calculateProportionOvershootArea();
	double initial_Average_Overshoot_Density = grid.calculateAverageOvershootDensity();


	// Snake traverse the grid and move cell to next region when density is exceeded
    for (int i = 0; i < gridSize; ++i) 
	{
        // Determine direction: left to right for even rows, right to left for odd rows
        bool leftToRight = (i % 2 == 0);

        for (int j = 0; j < gridSize; ++j) 
		{
            int col = leftToRight ? j : gridSize - 1 - j;

			int densityThreshold = region_area;
            if (grid.calculateAreaInRegion(i, col, i, col) > densityThreshold) 
			{
				int overShootArea = grid.calculateAreaInRegion(i, col, i, col) - densityThreshold;
				int shiftArea = 0;
				std::vector<int> shiftCell;
				for (int k = 0; k < grid.grid[i][col].size(); k++)
				{
					int cellIndex = grid.grid[i][col][k];
					int cellArea = grid.cellAreaIndexStore[cellIndex];
					if (shiftArea + cellArea < overShootArea)
					{
						shiftArea += cellArea;
						shiftCell.push_back(cellIndex);
					}
					else
					{
						shiftArea += cellArea;
						shiftCell.push_back(cellIndex);
						break;
					}	
				}

				// Move the cell to the next region
                // Identify the next cell in the snaking path
                int nextRow = i + (col == gridSize - 1 ? 1 : 0);
                int nextCol = leftToRight ? (col + 1) % gridSize : (col - 1 + gridSize) % gridSize;
                if (nextRow < gridSize) 
				{
					for (int k = 0; k < shiftCell.size(); k++)
					{
						int cellIndex = shiftCell[k];
						grid.moveCell(cellIndex, i, col, nextRow, nextCol);
						std::pair<int, int> newCellCoordinate = grid.getRegionCenter(nextRow, nextCol);
						// to the new coordinate add a random value between -1 % of Width/Height and 1 % of Width/Height
						int randomX = rand() % (die_width / 100);
						int randomY = rand() % (die_height / 100);
						// check if the new coordinate is within the die
						int newX = std::min(die_width, std::max(0, static_cast<int>(newCellCoordinate.first + randomX)));
						int newY = std::min(die_height, std::max(0, static_cast<int>(newCellCoordinate.second + randomY)));
						instance_pointers_[cellIndex]->setCoordinate(newX, newY);
					}
                }
                // else: handle the case when the cell is at the end of the grid, give to left or right neigbor
				else
				{
					for (int k = 0; k < shiftCell.size(); k++)
					{
						if (k % 3 == 0)
						{
							int cellIndex = shiftCell[k];
							grid.moveCell(cellIndex, i, col, gridSize-1, gridSize-2);
							std::pair<int, int> newCellCoordinate = grid.getRegionCenter(nextRow, nextCol);
							// to the new coordinate add a random value between -1 % of Width/Height and 1 % of Width/Height
							int randomX = rand() % (die_width / 100);
							int randomY = rand() % (die_height / 100);
							// check if the new coordinate is within the die
							int newX = std::min(die_width, std::max(0, static_cast<int>(newCellCoordinate.first + randomX)));
							int newY = std::min(die_height, std::max(0, static_cast<int>(newCellCoordinate.second + randomY)));
							instance_pointers_[cellIndex]->setCoordinate(newX, newY);

						}
						else if (k % 3 == 1)
						{
							int cellIndex = shiftCell[k];
							grid.moveCell(cellIndex, i, col, gridSize-2, gridSize-1);
							std::pair<int, int> newCellCoordinate = grid.getRegionCenter(nextRow, nextCol);
							// to the new coordinate add a random value between -1 % of Width/Height and 1 % of Width/Height
							int randomX = rand() % (die_width / 100);
							int randomY = rand() % (die_height / 100);
							// check if the new coordinate is within the die
							int newX = std::min(die_width, std::max(0, static_cast<int>(newCellCoordinate.first + randomX)));
							int newY = std::min(die_height, std::max(0, static_cast<int>(newCellCoordinate.second + randomY)));
							instance_pointers_[cellIndex]->setCoordinate(newX, newY);

						}
						else
						{
							int cellIndex = shiftCell[k];
							grid.moveCell(cellIndex, i, col, gridSize-2, gridSize-2);
							std::pair<int, int> newCellCoordinate = grid.getRegionCenter(nextRow, nextCol);
							// to the new coordinate add a random value between -1 % of Width/Height and 1 % of Width/Height
							int randomX = rand() % (die_width / 100);
							int randomY = rand() % (die_height / 100);
							// check if the new coordinate is within the die
							int newX = std::min(die_width, std::max(0, static_cast<int>(newCellCoordinate.first + randomX)));
							int newY = std::min(die_height, std::max(0, static_cast<int>(newCellCoordinate.second + randomY)));
							instance_pointers_[cellIndex]->setCoordinate(newX, newY);

						}

					}
				}
            }
        }
    }

	std::cout << "After Top-Down Snake Traversal Grid Based spread out: " << std::endl;
	std::cout << "HPWL: " << std::endl;
	std::cout << this->getHPWL() << std::endl;
	// print overshoot areas
	std::cout << "Average overshoot area: " << std::endl;
	std::cout << grid.calculateAverageOvershootArea() << std::endl;
	std::cout << "Median overshoot area: " << std::endl;
	std::cout << grid.calculateMedianOvershootArea() << std::endl;
	std::cout << "Max overshoot area: " << std::endl;
	std::cout << grid.calculateMaxOvershootArea() << std::endl;
	std::cout << "Proportion of regions with overshoot area > 0: " << std::endl;
	std::cout << grid.calculateProportionOvershootArea() << std::endl;
	std::cout << "Average overshoot density: " << std::endl;
	std::cout << grid.calculateAverageOvershootDensity() << std::endl;

	ulong top_down_HPWL = this->getHPWL();
	double top_down_Average_Overshoot = grid.calculateAverageOvershootArea();
	double top_down_Median_Overshoot = grid.calculateMedianOvershootArea();
	double top_down_Max_Overshoot = grid.calculateMaxOvershootArea();
	double top_down_Proportion_Overshoot = grid.calculateProportionOvershootArea();
	double top_down_Average_Overshoot_Density = grid.calculateAverageOvershootDensity();

	

	// Do Bottom-Up Snake Traversal Grid Based spread out
	// Reverse snaking traverse the grid (bottom-up)
	for (int i = gridSize - 1; i >= 0; --i) 
	{
		bool leftToRight = ((gridSize - 1 - i) % 2 == 0);
		
		for (int j = 0; j < gridSize; ++j) 
		{
			int col = leftToRight ? j : gridSize - 1 - j;

			int densityThreshold = region_area;
			if (grid.calculateAreaInRegion(i, col, i, col) > densityThreshold)
			{
				int overShootArea = grid.calculateAreaInRegion(i, col, i, col) - densityThreshold;
				int shiftArea = 0;
				std::vector<int> shiftCell;
				for (int k = 0; k < grid.grid[i][col].size(); k++)
				{
					int cellIndex = grid.grid[i][col][k];
					int cellArea = grid.cellAreaIndexStore[cellIndex];
					if (shiftArea + cellArea < overShootArea)
					{
						shiftArea += cellArea;
						shiftCell.push_back(cellIndex);
					}
					else
					{
						shiftArea += cellArea;
						shiftCell.push_back(cellIndex);
						break;
					}
				}

				// Move the cell to the next region
				// Identify the next cell in the snaking path
				int nextRow = i - (leftToRight ? (col == gridSize - 1 ? 1 : 0) : (col == 0 ? 1 : 0));
				int nextCol = leftToRight ? (col + 1) % gridSize : (col - 1 + gridSize) % gridSize;
				if (nextRow >= 0)
				{
					for (int k = 0; k < shiftCell.size(); k++)
					{
						int cellIndex = shiftCell[k];
						grid.moveCell(cellIndex, i, col, nextRow, nextCol);
						std::pair<int, int> newCellCoordinate = grid.getRegionCenter(nextRow, nextCol);
						// to the new coordinate add a random value between -1 % of Width/Height and 1 % of Width/Height
						int randomX = rand() % (die_width / 100);
						int randomY = rand() % (die_height / 100);
						// check if the new coordinate is within the die
						int newX = std::min(die_width, std::max(0, static_cast<int>(newCellCoordinate.first + randomX)));
						int newY = std::min(die_height, std::max(0, static_cast<int>(newCellCoordinate.second + randomY)));
						instance_pointers_[cellIndex]->setCoordinate(newX, newY);
					}
				}
				else
				{
					for (int k = 0; k < shiftCell.size(); k++)
					{
						if (k % 3 == 0)
						{
							int cellIndex = shiftCell[k];
							grid.moveCell(cellIndex, i, col, 0, 1);
							std::pair<int, int> newCellCoordinate = grid.getRegionCenter(nextRow, nextCol);
							// to the new coordinate add a random value between -1 % of Width/Height and 1 % of Width/Height
							int randomX = rand() % (die_width / 100);
							int randomY = rand() % (die_height / 100);
							// check if the new coordinate is within the die
							int newX = std::min(die_width, std::max(0, static_cast<int>(newCellCoordinate.first + randomX)));
							int newY = std::min(die_height, std::max(0, static_cast<int>(newCellCoordinate.second + randomY)));
							instance_pointers_[cellIndex]->setCoordinate(newX, newY);

						}
						else if (k % 3 == 1)
						{
							int cellIndex = shiftCell[k];
							grid.moveCell(cellIndex, i, col, 1, 0);
							std::pair<int, int> newCellCoordinate = grid.getRegionCenter(nextRow, nextCol);
							// to the new coordinate add a random value between -1 % of Width/Height and 1 % of Width/Height
							int randomX = rand() % (die_width / 100);
							int randomY = rand() % (die_height / 100);
							// check if the new coordinate is within the die
							int newX = std::min(die_width, std::max(0, static_cast<int>(newCellCoordinate.first + randomX)));
							int newY = std::min(die_height, std::max(0, static_cast<int>(newCellCoordinate.second + randomY)));
							instance_pointers_[cellIndex]->setCoordinate(newX, newY);

						}
						else
						{
							int cellIndex = shiftCell[k];
							grid.moveCell(cellIndex, i, col, 1, 1);
							std::pair<int, int> newCellCoordinate = grid.getRegionCenter(nextRow, nextCol);
							// to the new coordinate add a random value between -1 % of Width/Height and 1 % of Width/Height
							int randomX = rand() % (die_width / 100);
							int randomY = rand() % (die_height / 100);
							// check if the new coordinate is within the die
							int newX = std::min(die_width, std::max(0, static_cast<int>(newCellCoordinate.first + randomX)));
							int newY = std::min(die_height, std::max(0, static_cast<int>(newCellCoordinate.second + randomY)));
							instance_pointers_[cellIndex]->setCoordinate(newX, newY);

						}


					}
					
				}
			}

		}
	}

	std::cout << "After Bottom-Up Snake Traversal Grid Based spread out: " << std::endl;
	std::cout << "HPWL: " << std::endl;
	std::cout << this->getHPWL() << std::endl;
	// print overshoot areas
	std::cout << "Average overshoot area: " << std::endl;
	std::cout << grid.calculateAverageOvershootArea() << std::endl;
	std::cout << "Median overshoot area: " << std::endl;
	std::cout << grid.calculateMedianOvershootArea() << std::endl;
	std::cout << "Max overshoot area: " << std::endl;
	std::cout << grid.calculateMaxOvershootArea() << std::endl;
	std::cout << "Proportion of regions with overshoot area > 0: " << std::endl;
	std::cout << grid.calculateProportionOvershootArea() << std::endl;
	std::cout << "Average overshoot density: " << std::endl;

	ulong bottom_up_HPWL = this->getHPWL();
	double bottom_up_Average_Overshoot = grid.calculateAverageOvershootArea();
	double bottom_up_Median_Overshoot = grid.calculateMedianOvershootArea();
	double bottom_up_Max_Overshoot = grid.calculateMaxOvershootArea();
	double bottom_up_Proportion_Overshoot = grid.calculateProportionOvershootArea();
	double bottom_up_Average_Overshoot_Density = grid.calculateAverageOvershootDensity();




	



	// Implement Simulated Annealing
	// Initialize temperature
	double temperature = 1000;
	// Initialize cooling rate
	double coolingRate = 0.03;
	// Initialize current solution
	double currentSolution = this->getHPWL();

	double current_Cost = 0.8 * bottom_up_Average_Overshoot_Density + 0.2 * (static_cast<double>(bottom_up_HPWL) / initial_HPWL);

	while (temperature > 1)
	{
		// Randomly select a cell
		int randomCellIndex = rand() % num_instances;
		int current_i = grid.findRegion(randomCellIndex).first;
		int current_j = grid.findRegion(randomCellIndex).second;
		std::pair<int, int> randomCellCoordinate = instance_pointers_[randomCellIndex]->getCoordinate();

		// Random move increment
		int random_x_increment = rand() % (die_width / 200);
		int random_y_increment = rand() % (die_height / 200);

		int new_x = std::min(die_width, std::max(0, static_cast<int>(randomCellCoordinate.first + random_x_increment)));
		int new_y = std::min(die_height, std::max(0, static_cast<int>(randomCellCoordinate.second + random_y_increment)));
		// make sure the new coordinate is within the die
		std::pair<int, int> newCellCoordinate(new_x, new_y);

		int new_i = grid.findRegion(newCellCoordinate).first;
		int new_j = grid.findRegion(newCellCoordinate).second;
		int action = rand() % 5;
		if (action == 0)
		{
			int anotherCellIndex = rand() % num_instances;
			if (randomCellIndex == anotherCellIndex)
			{
				continue;
			}
		}
		else
		{


			grid.moveCell(randomCellIndex, current_i, current_j, new_i, new_j);

			instance_pointers_[randomCellIndex]->setCoordinate(new_x, new_y);



		}
		double new_Cost = 0.8 * grid.calculateAverageOvershootDensity() + 0.2 * (static_cast<double>(this->getHPWL()) / initial_HPWL);

		double acceptanceProbability = exp((current_Cost - new_Cost) / temperature);
		if (new_Cost < current_Cost)
		{
			current_Cost = new_Cost;
		}
		else if (acceptanceProbability > static_cast<double>(rand()) / static_cast<double>(RAND_MAX))
		{
			current_Cost = new_Cost;
		}
		else // revert changes
		{
			grid.moveCell(randomCellIndex, new_i, new_j, current_i, current_j);

			instance_pointers_[randomCellIndex]->setCoordinate(randomCellCoordinate.first, randomCellCoordinate.second);
		}
		temperature = temperature * (1 - coolingRate);
	}

	std::cout << "After Simulated Annealing: " << std::endl;
	std::cout << "HPWL: " << std::endl;
	std::cout << this->getHPWL() << std::endl;
	// print overshoot areas
	std::cout << "Average overshoot area: " << std::endl;
	std::cout << grid.calculateAverageOvershootArea() << std::endl;
	std::cout << "Median overshoot area: " << std::endl;
	std::cout << grid.calculateMedianOvershootArea() << std::endl;
	std::cout << "Max overshoot area: " << std::endl;
	std::cout << grid.calculateMaxOvershootArea() << std::endl;
	std::cout << "Proportion of regions with overshoot area > 0: " << std::endl;
	std::cout << grid.calculateProportionOvershootArea() << std::endl;
	std::cout << "Average overshoot density: " << std::endl;
	std::cout << grid.calculateAverageOvershootDensity() << std::endl;


	ulong simulated_annealing_HPWL = this->getHPWL();
	double simulated_annealing_Average_Overshoot = grid.calculateAverageOvershootArea();
	double simulated_annealing_Median_Overshoot = grid.calculateMedianOvershootArea();
	double simulated_annealing_Max_Overshoot = grid.calculateMaxOvershootArea();
	double simulated_annealing_Proportion_Overshoot = grid.calculateProportionOvershootArea();
	double simulated_annealing_Average_Overshoot_Density = grid.calculateAverageOvershootDensity();


	// compare the HPWL and overshoot area before and after the grid based spread out
	std::cout << "Initial HPWL: " << initial_HPWL << std::endl;
	std::cout << "Top-Down HPWL: " << top_down_HPWL << std::endl;
	std::cout << "Bottom-Up HPWL: " << bottom_up_HPWL << std::endl;
	std::cout << "Simulated Annealing HPWL: " << simulated_annealing_HPWL << std::endl;
	std::cout << "Initial Average Overshoot Area: " << initial_Average_Overshoot << std::endl;
	std::cout << "Top-Down Average Overshoot Area: " << top_down_Average_Overshoot << std::endl;
	std::cout << "Bottom-Up Average Overshoot Area: " << bottom_up_Average_Overshoot << std::endl;
	std::cout << "Simulated Annealing Average Overshoot Area: " << simulated_annealing_Average_Overshoot << std::endl;
	std::cout << "Initial Median Overshoot Area: " << initial_Median_Overshoot << std::endl;
	std::cout << "Top-Down Median Overshoot Area: " << top_down_Median_Overshoot << std::endl;
	std::cout << "Bottom-Up Median Overshoot Area: " << bottom_up_Median_Overshoot << std::endl;
	std::cout << "Simulated Annealing Median Overshoot Area: " << simulated_annealing_Median_Overshoot << std::endl;
	std::cout << "Initial Max Overshoot Area: " << initial_Max_Overshoot << std::endl;
	std::cout << "Top-Down Max Overshoot Area: " << top_down_Max_Overshoot << std::endl;
	std::cout << "Bottom-Up Max Overshoot Area: " << bottom_up_Max_Overshoot << std::endl;
	std::cout << "Simulated Annealing Max Overshoot Area: " << simulated_annealing_Max_Overshoot << std::endl;
	std::cout << "Initial Proportion Overshoot Area: " << initial_Proportion_Overshoot << std::endl;
	std::cout << "Top-Down Proportion Overshoot Area: " << top_down_Proportion_Overshoot << std::endl;
	std::cout << "Bottom-Up Proportion Overshoot Area: " << bottom_up_Proportion_Overshoot << std::endl;
	std::cout << "Simulated Annealing Proportion Overshoot Area: " << simulated_annealing_Proportion_Overshoot << std::endl;
	std::cout << "Initial Average Overshoot Density: " << initial_Average_Overshoot_Density << std::endl;
	std::cout << "Top-Down Average Overshoot Density: " << top_down_Average_Overshoot_Density << std::endl;
	std::cout << "Bottom-Up Average Overshoot Density: " << bottom_up_Average_Overshoot_Density << std::endl;
	std::cout << "Simulated Annealing Average Overshoot Density: " << simulated_annealing_Average_Overshoot_Density << std::endl;



	






}
}
