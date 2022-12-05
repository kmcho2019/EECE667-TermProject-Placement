/*!
 * Write **your own** code.
 * !!! Cheating will be strictly not accepted. !!!
 * If cheating is detected by the similarity check program and TA determine that you have cheated,
 * then you will get F grade or zero point for this term project.
 * You can use external libraries for only pure math libraries; i.e) fft, sparse matrix, etc
 * If you want to use external library, then please check whether it is okay by contact to TA.
 * */

#include "Circuit.h"
#include "matrixSolver.h"
#include <ctime>
#include <set>

namespace Placer {
void Circuit::quadraticPlacement() {
  // matrix solve example.
  // You should refer below function when you get x respect to Ax = b
  // Below function is implemented in src/algorithm/math/matrixSolver.cpp
  
  // Build Metrix A
  coo_matrix A;

  // Give ID to instances
  unordered_map<string, int> instMap;
  int instCnt = 0;
  for (auto &inst : instance_pointers_) {
    instMap.insert(make_pair(inst->getName(), instCnt++));
  }

  vector<int> vector_row_idx;
  vector<int> vector_col_idx;
  vector<double> vector_data;
  vector<double> vector_bx(instCnt);
  vector<double> vector_by(instCnt);
  
  // int cnt = 0;
  // for(auto &net : net_pointers_) {
  //   if(net->getConnectedPins().size() > 50) {
  //     cout << net->getName() << endl;
  //     cnt ++;
  //   }
  // }
  // cout << "CNT : " << cnt << endl;
  // Adjacent Matrix & Get A elements
  // cout << "Inst CNT " << instCnt <<endl;
  for (auto &inst : instance_pointers_) {
    // cout << "Inst Name : "<<inst->getName() <<endl;
    int inst1Num = instMap.find(inst->getName())->second;
    double diagA = 0;
    double bx = 0;
    double by = 0;
    // For all pins in the instance
    for(auto &pin : inst->getPins()) {
      // cout << "In inst pin : "<<pin->getCoordinate().first << " : "<< pin->getCoordinate().second<<endl; 
      if(pin->getPinName() == "CK") continue;
      if(pin->getNet()) {
        Net *net = pin->getNet();
        if (net->getSignalType() != "POWER" && net->getSignalType() != "GROUND" && net->getSignalType() != "CLOCK" && net->getSignalType() != "RESET") {
          int netWeight = net->getWeight();
          if(net->getConnectedPins().size() > 10) continue;
          // For all connected pins
          for(auto &connectedPin : net->getConnectedPins()) {
            // When connected pin is in Instance
            if(connectedPin->isInstancePin()){
              int inst2Num = instMap.find(connectedPin->getInstance()->getName())->second;
              if(inst1Num != inst2Num) {
                vector_row_idx.push_back(inst1Num);
                vector_col_idx.push_back(inst2Num);
                vector_data.push_back(-1.0 * netWeight);

                bx += netWeight * (connectedPin->getCoordinate().first - pin->getCoordinate().first);
                by += netWeight * (connectedPin->getCoordinate().second - pin->getCoordinate().second);
                diagA += (double)netWeight;
                // cout << "inst pin : " << netWeight << ", " << inst->getName() << " and " << connectedPin->getInstance()->getName() << endl;
              }
            }
            // When connected pin is in blockPin
            else if(connectedPin->isBlockPin()) {
              // // cout << "Block : " << connectedPin->getCoordinate().first << connectedPin->getCoordinate().second <<endl;
              bx += netWeight * (connectedPin->getCoordinate().first - pin->getCoordinate().first);
              by += netWeight * (connectedPin->getCoordinate().second - pin->getCoordinate().second);
              diagA += (double)netWeight;
              // cout << "block pin : " << netWeight << ", " << inst->getName() << " and " << connectedPin->getCoordinate().first << " : "<< connectedPin->getCoordinate().second <<endl;
            }
          }
        }
      }
    }
    // Update diagonal value and matrix b
    vector_row_idx.push_back(inst1Num);
    vector_col_idx.push_back(inst1Num);
    vector_data.push_back(diagA);
    vector_bx[inst1Num] = (double)bx;
    vector_by[inst1Num] = (double)by;
  }
  
  // Update A and b for matrixSolver
  A.n = instCnt;
  A.nnz = vector_data.size();
  A.row.resize(vector_row_idx.size());
  A.col.resize(vector_col_idx.size());
  A.dat.resize(vector_data.size());

  A.row = valarray<int>(vector_row_idx.data(), A.nnz);
  A.col = valarray<int>(vector_col_idx.data(), A.nnz);
  A.dat = valarray<double>(vector_data.data(), A.nnz);
  valarray<double> x(0.0, A.n);
  valarray<double> y(0.0, A.n);
  valarray<double> bx(vector_bx.data(), A.n);
  valarray<double> by(vector_by.data(), A.n);
  cout << "A constructed" <<endl;
  // Solve Ax = b
  A.solve(bx, x);
  A.solve(by, y);
  cout << "A solved" <<endl;

  // Place with solutions
  for (auto &inst : instance_pointers_) {
    int inst1Num = instMap.find(inst->getName())->second;
    inst->setCoordinate(x[inst1Num], y[inst1Num]);
  }
}
template <typename T,typename U>
T maxClamp(const T x, const U max) {
  if(x > max) return max;
  return x;
}

template <typename T,typename U>
T minClamp(const T x, const U min) {
  if(x < min) return min;
  return x;
}

template <typename T,typename U>
T clamp(const T x, const U min, const U max) {
  if(x < min) return min;
  else if(x > max) return max;
  return x;
}

template <typename T,typename U>                                                   
std::pair<T,U> operator+(const std::pair<T,U> & l,const std::pair<T,U> & r) {   
  return {l.first+r.first,l.second+r.second};                                    
} 

template <typename T,typename U>                                                   
std::pair<T,T> operator*(const U &mul, const std::pair<T,T> & r) {   
  return {mul * r.first, mul * r.second};   
} 

class Bin {
public:
  Bin(){}
  pair<double, double> densityForce[4];
  double cell_area;
  double stayForce;
  void reset() {
    densityForce[0] = make_pair(0.0, 0.0);
    densityForce[1] = make_pair(0.0, 0.0);
    densityForce[2] = make_pair(0.0, 0.0);
    densityForce[3] = make_pair(0.0, 0.0);
    stayForce = -100.0;
    cell_area = 0.0;
  }
};

void Circuit::myPlacement() {
  // Do qplace
  // this->quadraticPlacement(); 
  clock_t start, end;
  start = clock();
  // Place with solutions
  long maxWidth = 0, maxHeight = 0;
  for(Instance *instance : instance_pointers_) {
    if(maxWidth < instance->getWidth()) {
      maxWidth = instance->getWidth();
    }
    if(maxHeight < instance->getHeight()) {
      maxHeight = instance->getHeight();
    }
  }
  uint die_width = die_->getWidth();
  uint die_height = die_->getHeight();
  // cout << die_width << " " <<  die_height << " " << (long long)die_width * die_height<<endl;
  // long long cellA = 0;
  // for (auto &inst : instance_pointers_) {
  //   cellA += (long long) inst->getArea();
  // }
  // cout << cellA << endl;
  int INITCNT = 3;
  for (auto &inst : instance_pointers_) {
    inst->setCoordinate(int(die_width/2), int(die_height/2));
  }

  for(int i = 0; i< INITCNT; i++) {
    for (auto &inst : instance_pointers_) {
      double hpwlForceX = 0.0;
      double hpwlForceY = 0.0;
      
      pair<int, int> instCoordinate = inst->getCoordinate();
      // For all pins in the instance
      for(auto &pin : inst->getPins()) {
        if(pin->getPinName() == "CK") continue;
        if(pin->getNet()) {
          Net *net = pin->getNet();
          // if (net->getSignalType() != "POWER" && net->getSignalType() != "GROUND" && net->getSignalType() != "CLOCK" && net->getSignalType() != "RESET") {
            int netWeight = net->getWeight();
            // For all connected pins
            pair<int, int> pinCoordinate = instCoordinate + pin->getCoordinate();
            for(auto &connectedPin : net->getConnectedPins()) {
              // When connected pin is in Instance
                // When connected pin is in Instance
                if(connectedPin->isInstancePin() && i>0){
                  Instance *connInst = connectedPin->getInstance();
                  pair<int, int> connPinCoordinate = connInst->getCoordinate() + connectedPin->getCoordinate();
                  hpwlForceX += (double)netWeight * (pinCoordinate.first - connPinCoordinate.first);
                  hpwlForceY += (double)netWeight * (pinCoordinate.second - connPinCoordinate.second);                
                }
                else if(connectedPin->isBlockPin()) {
                  pair<int, int> connPinCoordinate = connectedPin->getCoordinate();
                  hpwlForceX += (double)netWeight * (pinCoordinate.first - connPinCoordinate.first);
                  hpwlForceY += (double)netWeight * (pinCoordinate.second - connPinCoordinate.second);                  
                }
            }
          // }
        }
      }
      inst->hpwlForce = make_pair(hpwlForceX, hpwlForceY);
    }
    for (auto &inst : instance_pointers_) {
      pair<int, int> current = inst->getCoordinate();
      pair<double, double> hpwlForce = inst->hpwlForce;
      hpwlForce.first = -hpwlForce.first * 0.000001 * die_width;
      hpwlForce.second = -hpwlForce.second * 0.000001 * die_height;
      int newX = current.first + int(hpwlForce.first);
      if(newX > die_width - inst->getWidth()) newX = die_width - inst->getWidth();
      int newY = current.second + int(hpwlForce.second);
      if(newY > die_height - inst->getHeight()) newY = die_height - inst->getHeight();
      inst->setCoordinate(newX, newY);
    }
    ulong HPWL = 0;
    for (Net *net : net_pointers_) {
      HPWL += net->getHPWL();
    }

    cout << "INITIAL HPWL " << HPWL <<endl;
  }

  // int number_of_grid_X = floor(die_width / maxWidth) - 1;
  // int number_of_grid_Y = number_of_grid_X;
  int number_of_grid_X = 40;
  int number_of_grid_Y = 40;  
  bool fitWidth = true;
  bool fitHeight = true;
  cout << number_of_grid_X << " by " << number_of_grid_Y <<endl;

  int normal_bin_width = static_cast<int>(die_width / number_of_grid_X);
  int normal_bin_height = static_cast<int>(die_height / number_of_grid_Y); 
  if(die_width % normal_bin_width != 0) {
    fitWidth = false;
    number_of_grid_X += 1;
  }
  if(die_height % normal_bin_height != 0) {
    fitHeight = false;
    number_of_grid_Y += 1;
  }
  // Give ID to instances
  unordered_map<string, int> instMap;
  int instCnt = 0;
  for (auto &inst : instance_pointers_) {
    instMap.insert(make_pair(inst->getName(), instCnt++));
  }

  // Bin Construction
  vector<vector<Bin> > bins2D(number_of_grid_X + 4, vector<Bin> (number_of_grid_Y + 4));
  
  // Iterate until 
  bool condition = true;
  double coeff = 0.0;
  int iter = 0;
  ulong prevHPWL = 0;
  double prevDen = 0.0;

  while(condition) {
    //Bin reset
    for(int i = 0; i <= number_of_grid_X+3; i++) {
      for (int j = 0; j <= number_of_grid_Y+3; j++) {
        bins2D[i][j].reset();
      }
    }
    // Bin & inst update
    for(auto &inst : instance_pointers_) {
      int instID = instMap.find(inst->getName())->second;
      int position_x = inst->getCoordinate().first;
      int position_y = inst->getCoordinate().second;
      uint instWidth = inst->getWidth();
      uint instHeight = inst->getHeight();
      int bin_coordinate_x = floor(position_x / normal_bin_width) + 2;
      int bin_coordinate_y = floor(position_y / normal_bin_height) + 2;
      inst->binCoordinate = make_pair(bin_coordinate_x, bin_coordinate_y);
      // cout << position_x << " == " << instWidth << " == " << bin_coordinate_x << " == " <<normal_bin_width<<endl;
      double _overflowX = (double)position_x + (double)instWidth - (double)bin_coordinate_x * (double)normal_bin_width;
      double overflowX = (double)_overflowX / (double)instWidth;
      double _overflowY = (double)position_y + (double)instHeight - (double)bin_coordinate_y * (double)normal_bin_height;
      double overflowY = (double)_overflowY / (double)instHeight;
      // cout << overflowX << " and "<< overflowY <<endl;
      // if(overflowX > 0 || overflowY > 0) {
      //   cout << overflowX << " and "<< overflowY <<endl;
      // }
      if(overflowX > 0) {
        if(overflowY > 0) {
          inst->binType = 3;
          bins2D[bin_coordinate_x][bin_coordinate_y].cell_area += (double)inst->getArea() * (1 - overflowX) * (1 - overflowY);
          bins2D[bin_coordinate_x][bin_coordinate_y + 1].cell_area += (double)inst->getArea() * (1 - overflowX) * overflowY;        
          bins2D[bin_coordinate_x + 1][bin_coordinate_y].cell_area += (double)inst->getArea() * overflowX * (1 - overflowY);          
          bins2D[bin_coordinate_x + 1][bin_coordinate_y+1].cell_area += (double)inst->getArea() * overflowX * overflowY;     
        }
        else {
          inst->binType = 1;
          bins2D[bin_coordinate_x][bin_coordinate_y].cell_area += (double)inst->getArea() * (1 - overflowX);
          bins2D[bin_coordinate_x + 1][bin_coordinate_y].cell_area += (double)inst->getArea() * overflowX;
        }
      }
      else {
        if(overflowY > 0) {
          inst->binType = 2;
          bins2D[bin_coordinate_x][bin_coordinate_y].cell_area += (double)inst->getArea() * (1 - overflowY);
          bins2D[bin_coordinate_x][bin_coordinate_y + 1].cell_area += (double)inst->getArea() * overflowY;          
        }
        else {
          inst->binType = 0;
          bins2D[bin_coordinate_x][bin_coordinate_y].cell_area += (double)inst->getArea();
        }        
      }
    }

    // cout << "Cell area" <<endl;
    // Stay Force. Stay Force < 0 : have to move
    for (int i = 2; i <= number_of_grid_X+1; i++) {
      for (int j = 2; j <= number_of_grid_Y+1; j++) {
        long bin_width = normal_bin_width;
        if(!fitWidth && i == number_of_grid_X) {
          bin_width = die_width - bin_width * (number_of_grid_X - 1);
        }      
        long bin_height = normal_bin_height;
        if(!fitHeight && j == number_of_grid_Y) {
          bin_height = die_height - bin_height * (number_of_grid_Y - 1);
        }
        bins2D[i][j].stayForce = 1.0 - bins2D[i][j].cell_area / (double)(bin_width * bin_height);
      }
    }

    int left = number_of_grid_X/2 + 1;
    int down = number_of_grid_Y/2 + 1;
    int mini = 0, minj = 0;

    double minstayForce = 1.0, sumNegstay = 0.0, fd = 0.0, fl = 0.0;
    
    for (int i = 2; i <= number_of_grid_X+1; i++) {
      for (int j = 2 ; j <= number_of_grid_Y+1; j++) {
        if(minstayForce > bins2D[i][j].stayForce) {
          mini = i;
          minj = j;
          minstayForce = bins2D[i][j].stayForce;
        }
        if(bins2D[i][j].stayForce < 0) {
          sumNegstay += bins2D[i][j].stayForce;
          // if(i < left) {
          //   fl += bins2D[i][j].stayForce;
          // }
          // else if(i > number_of_grid_X - left) {
          //   fl -= bins2D[i][j].stayForce;
          // }
          // if(j < down) {
          //   fd += bins2D[i][j].stayForce;
          // }
          // else if(j>number_of_grid_Y - down){
          //   fd -= bins2D[i][j].stayForce;
          // }
        }
        if(iter % 10 == 0) cout << bins2D[i][j].stayForce << " ";
        // cout << bins2D[i][j].stayForce << " ";
      }
      if(iter % 10 == 0) cout << endl;
      // cout << endl;
    }
    if(minstayForce > 0) condition = false;
    //Global
    double globalForceX = 0.0, globalForceY = 0.0;
    double globalThreshold = 0.1;
    // globalForceX = (mini < left) ? 1.5 : -1.5;
    // globalForceY = (minj < down) ? 1.5 : -1.5;
    if(fl < sumNegstay * globalThreshold) globalForceX = 1.5;
    else if(fl > -sumNegstay * globalThreshold) globalForceX = -1.5;
    if(fd < sumNegstay * globalThreshold) globalForceY = 1.5;
    else if(fd > -sumNegstay * globalThreshold) globalForceY = -1.5;

    // Move Force
    for(auto &inst : instance_pointers_) {
      int bin_coordinate_x = inst->binCoordinate.first;
      int bin_coordinate_y = inst->binCoordinate.second;
      int binType = inst->binType;

      if(bins2D[bin_coordinate_x][bin_coordinate_y].densityForce[binType].first != 0 && bins2D[bin_coordinate_x][bin_coordinate_y].densityForce[binType].second != 0) {
        inst->densityForce = bins2D[bin_coordinate_x][bin_coordinate_y].densityForce[binType];
        continue;
      }

      double densityForceX = 0;
      double densityForceY = 0;

      if(binType == 0) {
        densityForceX -= bins2D[bin_coordinate_x - 1][bin_coordinate_y - 1].stayForce;
        densityForceX -= bins2D[bin_coordinate_x - 1][bin_coordinate_y].stayForce;
        densityForceX -= bins2D[bin_coordinate_x - 1][bin_coordinate_y + 1].stayForce;
        densityForceX += bins2D[bin_coordinate_x + 1][bin_coordinate_y - 1].stayForce;
        densityForceX += bins2D[bin_coordinate_x + 1][bin_coordinate_y].stayForce;
        densityForceX += bins2D[bin_coordinate_x + 1][bin_coordinate_y + 1].stayForce;
        densityForceY -= bins2D[bin_coordinate_x - 1][bin_coordinate_y - 1].stayForce;
        densityForceY -= bins2D[bin_coordinate_x][bin_coordinate_y - 1].stayForce;
        densityForceY -= bins2D[bin_coordinate_x + 1][bin_coordinate_y - 1].stayForce;
        densityForceY += bins2D[bin_coordinate_x - 1][bin_coordinate_y + 1].stayForce;
        densityForceY += bins2D[bin_coordinate_x][bin_coordinate_y + 1].stayForce;
        densityForceY += bins2D[bin_coordinate_x + 1][bin_coordinate_y + 1].stayForce;
      }
      else if(binType == 1) {
        densityForceX -= bins2D[bin_coordinate_x][bin_coordinate_y - 1].stayForce;
        densityForceX -= bins2D[bin_coordinate_x][bin_coordinate_y].stayForce;
        densityForceX -= bins2D[bin_coordinate_x][bin_coordinate_y + 1].stayForce;
        densityForceX += bins2D[bin_coordinate_x + 1][bin_coordinate_y - 1].stayForce;
        densityForceX += bins2D[bin_coordinate_x + 1][bin_coordinate_y].stayForce;
        densityForceX += bins2D[bin_coordinate_x + 1][bin_coordinate_y + 1].stayForce;
        densityForceY -= bins2D[bin_coordinate_x][bin_coordinate_y - 1].stayForce;
        densityForceY -= bins2D[bin_coordinate_x + 1][bin_coordinate_y - 1].stayForce;
        densityForceY += bins2D[bin_coordinate_x][bin_coordinate_y + 1].stayForce;
        densityForceY += bins2D[bin_coordinate_x + 1][bin_coordinate_y + 1].stayForce;
      }
      else if(binType == 2) {
        densityForceX -= bins2D[bin_coordinate_x][bin_coordinate_y].stayForce;
        densityForceX -= bins2D[bin_coordinate_x][bin_coordinate_y + 1].stayForce;
        densityForceX += bins2D[bin_coordinate_x + 1][bin_coordinate_y].stayForce;
        densityForceX += bins2D[bin_coordinate_x + 1][bin_coordinate_y + 1].stayForce;
        densityForceY -= bins2D[bin_coordinate_x - 1][bin_coordinate_y].stayForce;
        densityForceY -= bins2D[bin_coordinate_x][bin_coordinate_y].stayForce;
        densityForceY -= bins2D[bin_coordinate_x + 1][bin_coordinate_y].stayForce;
        densityForceY += bins2D[bin_coordinate_x - 1][bin_coordinate_y + 1].stayForce;
        densityForceY += bins2D[bin_coordinate_x][bin_coordinate_y + 1].stayForce;
        densityForceY += bins2D[bin_coordinate_x + 1][bin_coordinate_y + 1].stayForce;
      }
      else {
        densityForceX -= bins2D[bin_coordinate_x][bin_coordinate_y].stayForce;
        densityForceX += bins2D[bin_coordinate_x + 1][bin_coordinate_y].stayForce;
        densityForceY -= bins2D[bin_coordinate_x][bin_coordinate_y].stayForce;
        densityForceY += bins2D[bin_coordinate_x][bin_coordinate_y + 1].stayForce;
      }
      inst->densityForce = make_pair(densityForceX, densityForceY);
      bins2D[bin_coordinate_x][bin_coordinate_y].densityForce[binType] = make_pair(densityForceX, densityForceY);
    }

    // cout << "Spreading"<<endl;
    for(auto &inst : instance_pointers_) {
      int bin_coordinate_x = inst->binCoordinate.first;
      int bin_coordinate_y = inst->binCoordinate.second;
      if(bins2D[bin_coordinate_x][bin_coordinate_y].stayForce < 0) {
        pair<int, int> binCoord = make_pair(-bin_coordinate_x * normal_bin_width + inst->getWidth()/2, -bin_coordinate_y * normal_bin_height  + inst->getHeight()/2 );
        pair<int, int> inbinCoord = inst->getCoordinate() + binCoord;
        double outerCoeff = 0.5 * maxClamp((-1) * bins2D[bin_coordinate_x][bin_coordinate_y].stayForce, 1.0);

        double dx = inst->densityForce.first;
        double dy = inst->densityForce.second;
        double adx = abs(dx);
        double ady = abs(dy);
        double ratio0 = 1 / max(adx, ady);
        ratio0 = ratio0 * (1.29289 - (adx + ady) * ratio0 * 0.29289);
        inst->densityForce = make_pair(outerCoeff * dx * ratio0, outerCoeff * dy * ratio0);

        double fx = (double)(inbinCoord.first - normal_bin_width/2);
        if(fx < 0 && bins2D[bin_coordinate_x-1][bin_coordinate_y].stayForce <0) fx = minClamp(fx, -0.1);  
        if(fx > 0 && bins2D[bin_coordinate_x+1][bin_coordinate_y].stayForce <0) fx = maxClamp(fx, 0.1);
        double fy = (double)(inbinCoord.second - normal_bin_height/2);
        if(fy < 0 && bins2D[bin_coordinate_x][bin_coordinate_y-1].stayForce <0) fy = minClamp(fy, -0.1);  
        if(fy > 0 && bins2D[bin_coordinate_x][bin_coordinate_y+1].stayForce <0) fy = maxClamp(fy, 0.1);
        double ax = abs(fx);
        double ay = abs(fy);
        double ratio = 1 / max(ax, ay);
        ratio = ratio * (1.29289 - (ax + ay) * ratio * 0.29289);
        inst->spreadForce = make_pair(outerCoeff * fx * ratio, outerCoeff * fy * ratio);
      }
      else {
        inst->spreadForce = make_pair(0.0, 0.0);
      }
    }
    // cout << "Density Fin"<<endl;
    // HPWL Force
    for(auto &inst : instance_pointers_) {
      double hpwlForceX = 0.0;
      double hpwlForceY = 0.0;
      
      pair<int, int> instCoordinate = inst->getCoordinate();
      // For all pins in the instance
      for(auto &pin : inst->getPins()) {
        if(pin->getPinName() == "CK") continue;
        if(pin->getNet()) {
          Net *net = pin->getNet();
          if (net->getSignalType() != "POWER" && net->getSignalType() != "GROUND" && net->getSignalType() != "CLOCK" && net->getSignalType() != "RESET") {
            double netWeight = (double)net->getWeight()/net->getConnectedPins().size();
            // For all connected pins
            pair<int, int> pinCoordinate = instCoordinate + pin->getCoordinate();
            for(auto &connectedPin : net->getConnectedPins()) {
              // When connected pin is in Instance
              if(connectedPin->isInstancePin()){
                Instance *connInst = connectedPin->getInstance();
                if(bins2D[inst->binCoordinate.first][inst->binCoordinate.second].stayForce < 0 && (connInst->binCoordinate.first == inst->binCoordinate.first && connInst->binCoordinate.second == inst->binCoordinate.second)) continue;
                pair<int, int> connPinCoordinate = connInst->getCoordinate() + connectedPin->getCoordinate();
                hpwlForceX += netWeight * (pinCoordinate.first - connPinCoordinate.first);
                hpwlForceY += netWeight * (pinCoordinate.second - connPinCoordinate.second);                
              }
              else if(connectedPin->isBlockPin()) {
                pair<int, int> connPinCoordinate = connectedPin->getCoordinate();
                hpwlForceX += netWeight * (pinCoordinate.first - connPinCoordinate.first);
                hpwlForceY += netWeight * (pinCoordinate.second - connPinCoordinate.second);                  
              }
            }
          }
        }
      }
      inst->hpwlForce = make_pair(hpwlForceX, hpwlForceY);
    }
    
    for (auto &inst : instance_pointers_) {
      pair<double, double> hpwlForce = inst->hpwlForce;
      inst->hpwlForce = make_pair(hpwlForce.first/maxhpwlForceX, hpwlForce.second/maxhpwlForceY);
    }
    // cout << "HPWL Fin"<<endl;

    // Update coeff
    ulong HPWL = 0;
    for (Net *net : net_pointers_) {
      HPWL += net->getHPWL();
    }
    
    double lambda = (prevDen < sumNegstay) ? 0.01 : 0.03;
    double sfcoeff = 1.0;
    prevDen = sumNegstay;
    prevHPWL = HPWL;

    for (auto &inst : instance_pointers_) {
      pair<int, int> current = inst->getCoordinate();
      // cout << densityForce.first << ", "<<densityForce.second << " \t ";
      // if(inst->binCoordinate.first == mini && inst->binCoordinate.second == minj) cout << " SF "<< sfcoeff *inst->spreadForce.first << ", "<< sfcoeff *inst->spreadForce.second << " HPWL " << -lambda * inst->hpwlForce.first << ", " << -lambda * inst->hpwlForce.second << " GF " << globalForceX << ", " << globalForceY <<endl;
      pair<double, double> force = sfcoeff * inst->spreadForce + (-lambda) * inst->hpwlForce;
      force = force + make_pair(globalForceX, globalForceY);

      int newX = clamp(current.first + normal_bin_width * int(force.first), 0, (int)(die_width - inst->getWidth()));
      int newY = clamp(current.second + normal_bin_height * int(force.second), 0, (int)(die_height - inst->getHeight()));
      inst->setCoordinate(newX, newY);
    }
    end = clock();
    double tresult = (double)(end-start)/CLOCKS_PER_SEC;
    tresult = (double) (end-start)/CLOCKS_PER_SEC;
    cout << "iter " << iter++ << " HPWL : " << HPWL << " min : " << minstayForce << " sum : " << sumNegstay << "\tTIME : " << tresult << endl;

    // cout<< "Moved"<<endl;
  }
}
}
