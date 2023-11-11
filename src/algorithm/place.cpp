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

void Circuit::quadraticPlacement() {
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
	std::vector<std::list<std::pair<int, int>>> adjacency_list(instance_pointers_.size());
	// Create data structure to store the weight of pad wire for each instance
	std::vector<int> pad_wire_weights(instance_pointers_.size(), 0);
	std::cout << "flag2" << endl;
	int test_idx = 0;

	// Initiate the adjacency list with the number of instances in the circuit
	// Also initiate the weight of pad wire for each instance
	// Checking structure:
	// instance -> pin -> net -> pin -> instance
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

	std::cout << "flag5" << std::endl;


	// Create b vector
	// bx vector:
	// If gate i connects to a pad at (xi, yi) with a wire with weight wi,
	// then set bx[i] = wi * xi
	std::valarray<double> bx_vector(A_matrix.n);
	std::valarray<double> by_vector(A_matrix.n);

	for (Instance *instance: instance_pointers_)
	{
		if (!instance)
		{
			continue;
		}
		int instance_index = instance_manager.getIndexByName(instance->getName());
		
		for (Pin *pin : instance->getPins())
		{
			if (!pin)
			{
				continue;
			}
			Net *net = pin->getNet();
			if (!net)
			{
				continue;
			}
			int current_net_weight = net->getWeight();
			for (Pin *connected_pin: net->getConnectedPins())
			{
				// Check if the connected_pin is a block pin
				if (connected_pin && connected_pin->isBlockPin())
				{

					std::pair<int, int> pin_coordinate = connected_pin->getCoordinate();
					// Accumulate net weight * x/y coordinate for block pins
					bx_vector[instance_index] += current_net_weight * pin_coordinate.first;
					by_vector[instance_index] += current_net_weight * pin_coordinate.second;
				}
			}
		}
	}
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
	std::cout << "HPWL of placeExample: " << std::endl;
	std::cout << this->getHPWL() << std::endl;



 }
void Circuit::myPlacement() {
	this->placeExample();  // You should erase this line.
}
}
