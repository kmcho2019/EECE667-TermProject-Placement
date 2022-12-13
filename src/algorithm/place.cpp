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


namespace Placer {
void Circuit::quadraticPlacement() {
  // matrix solve example.
  // You should refer below function when you get x respect to Ax = b
  // Below function is implemented in src/algorithm/math/matrixSolver.cpp
  
  // Build Metrix A
  coo_matrix A;

  // Give ID to instances
  
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
  //   if(net->getConnectedPins().size() > 40) {
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
template <typename T>                                                   
std::pair<T,T> operator*(const double &mul, const std::pair<T,T> & r) {   
  return {mul * r.first, mul * r.second};   
} 

static double fastExp(double a) {
  a = 1.0 + a / 1024.0;
  a *= a;
  a *= a;
  a *= a;
  a *= a;
  a *= a;
  a *= a;
  a *= a;
  a *= a;
  a *= a;
  a *= a;
  return a;
}

double get_w_x(double x) {
  return 2.0 * 3.141592 / (double)number_of_grid_X * x;
}

double Net::getHPWL_WA(double waCoeffX, double waCoeffY) {
  sumPosX = 0.0;
  sumNegX = 0.0;
  sumCPosX = 0.0;
  sumCNegX = 0.0;
  sumPosY = 0.0;
  sumNegY = 0.0;
  sumCPosY = 0.0;
  sumCNegY = 0.0;

  vector<Pin*> pins;
  for (dbITerm* db_i_term: db_net_->getITerms()) {
    pins.push_back(data_mapping_->pin_map_i[db_i_term]);
  }
  for(dbBTerm* db_b_term: db_net_->getBTerms()){
    pins.push_back(data_mapping_->pin_map_b[db_b_term]);
  }

  for(auto &pin : this->getConnectedPins()) {
    pair<int, int> pinCoordinate = pin->getCoordinate();

    double wPinX = pinCoordinate.first * waCoeffX;
    double wPinY = pinCoordinate.second * waCoeffY;

    sumCPosX += pinCoordinate.first * fastExp(wPinX);
    sumCPosY += pinCoordinate.second * fastExp(wPinY);
    sumCNegX += pinCoordinate.first * fastExp(-wPinX);
    sumCNegY += pinCoordinate.second * fastExp(-wPinY);
    sumPosX += fastExp(wPinX);
    sumPosY += fastExp(wPinY);
    sumNegX += fastExp(-wPinX);
    sumNegY += fastExp(-wPinY);
  }
  return sumCPosX/sumPosX - sumCNegX/sumNegX + sumCPosY/sumPosY - sumCNegY/sumNegY;
}

void Net::calcHPWL_gradWA(double waCoeffX, double waCoeffY) {
  sumPosX = 0.0;
  sumNegX = 0.0;
  sumCPosX = 0.0;
  sumCNegX = 0.0;
  sumPosY = 0.0;
  sumNegY = 0.0;
  sumCPosY = 0.0;
  sumCNegY = 0.0;
  
  // cout << "calcHPWL_gradWA INIT" <<endl;
  for(auto &pin : this->getConnectedPins()) {
    pair<int, int> pinCoordinate = pin->getCoordinate();
    
    double wPinX = pinCoordinate.first * waCoeffX;
    double wPinY = pinCoordinate.second * waCoeffY;

    sumCPosX += pinCoordinate.first * fastExp(wPinX);
    sumCPosY += pinCoordinate.second * fastExp(wPinY);
    sumCNegX += pinCoordinate.first * fastExp(-wPinX);
    sumCNegY += pinCoordinate.second * fastExp(-wPinY);
    sumPosX += fastExp(wPinX);
    sumPosY += fastExp(wPinY);
    sumNegX += fastExp(-wPinX);
    sumNegY += fastExp(-wPinY);
  }

  for(auto &pin : this->getConnectedPins()) {
    pair<int, int> pinCoordinate = pin->getCoordinate();
    double wPinX = pinCoordinate.first * waCoeffX;
    double wPinY = pinCoordinate.second * waCoeffY;

    double ePoxX = fastExp(wPinX);
    double eNegX = fastExp(-wPinX);
    double ePoxY = fastExp(wPinY);
    double eNegY = fastExp(-wPinY);

    pin->gradWAX = ePoxX / sumPosX * (1 + pinCoordinate.first * waCoeffX - sumCPosX/sumPosX * waCoeffX) - eNegX / sumNegX * (1 - pinCoordinate.first * waCoeffX + sumCNegX/sumNegX * waCoeffX);
    pin->gradWAY = ePoxY / sumPosY * (1 + pinCoordinate.first * waCoeffY - sumCPosY/sumPosY * waCoeffY) - eNegY / sumNegY * (1 - pinCoordinate.first * waCoeffY + sumCNegY/sumNegY * waCoeffY);
  }
}

double calcAlphaHat(vector<double> &vX, vector<double> &vY, vector<double> &prev_vX, vector<double> &prev_vY, vector<double> &gradX, vector<double> &gradY, vector<double> &prev_gradX, vector<double> &prev_gradY) {
  double sum = 0.0;
  for(int i = 0; i < vX.size(); i++) {
    double diffX = (vX[i] - prev_vX[i]);
    double diffY = (vY[i] - prev_vY[i]);
    sum += diffX * diffX;
    sum += diffY * diffY;
  }
  double normV = sqrt(sum);

  sum = 0.0;
  for(int i = 0; i < gradX.size(); i++) {
    double diffX = (gradX[i] - prev_gradX[i]);
    double diffY = (gradY[i] - prev_gradY[i]);
    
    sum += diffX * diffX;
    sum += diffY * diffY;
  }
  double normgradV = sqrt(sum);

  double result = 0.0;
  if(normgradV != 0.0) result = normV / normgradV;

  return result;
}

double Circuit::initLambda() {
  uint die_width = die_->getWidth();
  uint die_height = die_->getHeight();
  double normal_bin_width = (double)die_width / number_of_grid_X;
  double normal_bin_height = (double)die_height / number_of_grid_Y; 

  double gamma_ = 80.0 * normal_bin_width;
  double gamma = 1.0/gamma_;

  vector<vector<Bin> > bins2D(number_of_grid_X, vector<Bin> (number_of_grid_Y));  
  vector<vector<double> > a(number_of_grid_X, vector<double> (number_of_grid_Y));

  for(auto &net : net_pointers_) {
    for(auto &pin : net->getConnectedPins()) {
      pin->gradWAX = 0;
      pin->gradWAY = 0;
    }
  }

  // HPWL WA
  double HPWL_WA = 0.0;
  for(auto &net : net_pointers_) {
    net->calcHPWL_gradWA(gamma, gamma);
  }
  
  // Bin & inst update
  for(auto &inst : instance_pointers_) {
    int instID = instMap.find(inst->getName())->second;
    int position_x = inst->getCoordinate().first;
    int position_y = inst->getCoordinate().second;
    uint instWidth = inst->getWidth();
    uint instHeight = inst->getHeight();
    int bin_coordinate_x = floor(position_x / normal_bin_width);
    int bin_coordinate_y = floor(position_y / normal_bin_height);
    inst->binCoordinate = make_pair(bin_coordinate_x, bin_coordinate_y);
    double _overflowX = (double)position_x + (double)instWidth - (double)(bin_coordinate_x + 1) * (double)normal_bin_width;
    double overflowX = (double)_overflowX / (double)instWidth;
    double _overflowY = (double)position_y + (double)instHeight - (double)(bin_coordinate_y + 1) * (double)normal_bin_height;
    double overflowY = (double)_overflowY / (double)instHeight;

    if(overflowX > 0) {
      if(overflowY > 0) {
        inst->binType = 3;
        inst->overflowX = overflowX;
        inst->overflowY = overflowY;
        bins2D[bin_coordinate_x][bin_coordinate_y].cell_area += (double)inst->getArea() * (1 - overflowX) * (1 - overflowY);
        bins2D[bin_coordinate_x][bin_coordinate_y + 1].cell_area += (double)inst->getArea() * (1 - overflowX) * overflowY;        
        bins2D[bin_coordinate_x + 1][bin_coordinate_y].cell_area += (double)inst->getArea() * overflowX * (1 - overflowY);          
        bins2D[bin_coordinate_x + 1][bin_coordinate_y+1].cell_area += (double)inst->getArea() * overflowX * overflowY;     
      }
      else {
        inst->binType = 1;
        inst->overflowX = overflowX;
        inst->overflowY = 0.0;
        bins2D[bin_coordinate_x][bin_coordinate_y].cell_area += (double)inst->getArea() * (1 - overflowX);
        bins2D[bin_coordinate_x + 1][bin_coordinate_y].cell_area += (double)inst->getArea() * overflowX;
      }
    }
    else {
      if(overflowY > 0) {
        inst->binType = 2;
        inst->overflowX = 0.0;
        inst->overflowY = overflowY;
        bins2D[bin_coordinate_x][bin_coordinate_y].cell_area += (double)inst->getArea() * (1 - overflowY);
        bins2D[bin_coordinate_x][bin_coordinate_y + 1].cell_area += (double)inst->getArea() * overflowY;          
      }
      else {
        inst->binType = 0;
        inst->overflowX = 0.0;
        inst->overflowY = 0.0;
        bins2D[bin_coordinate_x][bin_coordinate_y].cell_area += (double)inst->getArea();
      }        
    }
  }
  for (int i = 0; i < number_of_grid_X; i++) {
    for (int j = 0; j < number_of_grid_Y; j++) {
      double density = bins2D[i][j].cell_area / (normal_bin_width * normal_bin_height);
      bins2D[i][j].density = density;
    }
  }
  // Calc a
  for (int i = 0; i < number_of_grid_X; i++) {
    double w_u = get_w_x(i);

  }
  for (int i = 0; i < number_of_grid_X; i++) {
    double w_u = get_w_x(i);
    for (int j = 0; j < number_of_grid_Y; j++) {
      double sum = 0.0;
      double w_v = get_w_x(j); 

      for (int x = 0; x < number_of_grid_X; x++) {
        for (int y = 0; y < number_of_grid_Y; y++) {
          sum += bins2D[x][y].density * cos(w_u * x) * cos(w_v * y);
        }
      }
      a[i][j] = 1.0 / number_of_grid_X / number_of_grid_Y * sum;
    }
  }

  double sumWA = 0.0, sumDensity = 0.0;
  for(auto &inst : instance_pointers_) {
    int instNum = instMap.find(inst->getName())->second;

    double gradInstX = 0.0, gradInstY = 0.0;
    for(auto &pin : inst->getPins()) {
      gradInstX += pin->gradWAX;
      gradInstY += pin->gradWAY;
    }
    sumWA += abs(gradInstX);
    sumWA += abs(gradInstY);

    pair<int, int> coordinate = inst->getCoordinate();
    int x = coordinate.first;
    int y = coordinate.second;

    double sum_x = 0.0;
    double sum_y = 0.0;

    for (int u = 0; u < number_of_grid_X; u++) {
      double w_u = get_w_x(u);
      double w_u_2 = w_u * w_u;
      for (int v = 0; v < number_of_grid_Y; v++) {
        double w_v = get_w_x(v); 
        double w_v_2 = w_v * w_v;
        double coeff;
        if(w_u_2 == 0 && w_v_2 == 0) coeff = 0.0;
        else coeff = a[u][v] / (w_u_2 + w_v_2);
        sum_x += coeff * w_u * sin(w_u * x) * cos(w_v * y);
        sum_y += coeff * w_v * cos(w_u * x) * sin(w_v * y);
      }
    }
    sumDensity += abs(inst->getArea() * sum_x);
    sumDensity += abs(inst->getArea() * sum_y);
  }
  return sumWA/sumDensity;
}

void Circuit::calcGradient(vector<double> &gradX, vector<double> &gradY, double lambda) {
  // Bin Construction
  uint die_width = die_->getWidth();
  uint die_height = die_->getHeight();
  double normal_bin_width = (double)die_width / number_of_grid_X;
  double normal_bin_height = (double)die_height / number_of_grid_Y; 
  
  vector<vector<Bin> > bins2D(number_of_grid_X, vector<Bin> (number_of_grid_Y));
  vector<vector<double> > a(number_of_grid_X, vector<double> (number_of_grid_Y));

  for(auto &net : net_pointers_) {
    for(auto &pin : net->getConnectedPins()) {
      pin->gradWAX = 0;
      pin->gradWAY = 0;
    }
  }
  for(int i = 0; i< gradX.size();i++) {
    gradX[i] = 0;
    gradY[i] = 0;
  }

  // Bin & inst update
  for(auto &inst : instance_pointers_) {
    int instID = instMap.find(inst->getName())->second;
    int position_x = inst->getCoordinate().first;
    int position_y = inst->getCoordinate().second;
    uint instWidth = inst->getWidth();
    uint instHeight = inst->getHeight();
    int bin_coordinate_x = floor(position_x / normal_bin_width);
    int bin_coordinate_y = floor(position_y / normal_bin_height);
    inst->binCoordinate = make_pair(bin_coordinate_x, bin_coordinate_y);
    double _overflowX = (double)position_x + (double)instWidth - (double)(bin_coordinate_x + 1) * (double)normal_bin_width;
    double overflowX = (double)_overflowX / (double)instWidth;
    double _overflowY = (double)position_y + (double)instHeight - (double)(bin_coordinate_y + 1) * (double)normal_bin_height;
    double overflowY = (double)_overflowY / (double)instHeight;
    
    if(overflowX > 0) {
      if(overflowY > 0) {
        inst->binType = 3;
        inst->overflowX = overflowX;
        inst->overflowY = overflowY;
        bins2D[bin_coordinate_x][bin_coordinate_y].cell_area += (double)inst->getArea() * (1 - overflowX) * (1 - overflowY);
        bins2D[bin_coordinate_x][bin_coordinate_y + 1].cell_area += (double)inst->getArea() * (1 - overflowX) * overflowY;        
        bins2D[bin_coordinate_x + 1][bin_coordinate_y].cell_area += (double)inst->getArea() * overflowX * (1 - overflowY);          
        bins2D[bin_coordinate_x + 1][bin_coordinate_y+1].cell_area += (double)inst->getArea() * overflowX * overflowY;     
      }
      else {
        inst->binType = 1;
        inst->overflowX = overflowX;
        inst->overflowY = 0.0;
        bins2D[bin_coordinate_x][bin_coordinate_y].cell_area += (double)inst->getArea() * (1 - overflowX);
        bins2D[bin_coordinate_x + 1][bin_coordinate_y].cell_area += (double)inst->getArea() * overflowX;
      }
    }
    else {
      if(overflowY > 0) {
        inst->binType = 2;
        inst->overflowX = 0.0;
        inst->overflowY = overflowY;
        bins2D[bin_coordinate_x][bin_coordinate_y].cell_area += (double)inst->getArea() * (1 - overflowY);
        bins2D[bin_coordinate_x][bin_coordinate_y + 1].cell_area += (double)inst->getArea() * overflowY;          
      }
      else {
        inst->binType = 0;
        inst->overflowX = 0.0;
        inst->overflowY = 0.0;
        bins2D[bin_coordinate_x][bin_coordinate_y].cell_area += (double)inst->getArea();
      }        
    }
  }
  for (int i = 0; i < number_of_grid_X; i++) {
    for (int j = 0; j < number_of_grid_Y; j++) {
      double density = bins2D[i][j].cell_area / (normal_bin_width * normal_bin_height);
      bins2D[i][j].density = density;
    }
  }
  // HPWL WA
  double maxDensity = 0.0;
  for (int x = 0; x < number_of_grid_X; x++) {
    for (int y = 0; y < number_of_grid_Y; y++) {
      if(maxDensity < bins2D[x][y].density) maxDensity = bins2D[x][y].density;
    }
  }
  double tau = max(maxDensity - 1.2 , 0.0) * (double)normal_bin_width * (double)normal_bin_height / (double)total_cell_area;
  double gamma_ = 8.0 * normal_bin_width * pow(10, 20.0/9.0 * tau - 11.0/9.0);
  double gamma = 1.0/gamma_;
  double HPWL_WA = 0.0;
  for(auto &net : net_pointers_) {
    // HPWL_WA += net->getHPWL_WA();
    if (net->getSignalType() != "POWER" && net->getSignalType() != "GROUND" && net->getSignalType() != "CLOCK" && net->getSignalType() != "RESET") {
      net->calcHPWL_gradWA(gamma, gamma);
    }
  }
  
  // Calc a
  for (int i = 0; i < number_of_grid_X; i++) {
    double w_u = get_w_x(i); 
    for (int j = 0; j < number_of_grid_Y; j++) {
      double sum = 0.0;
      double w_v = get_w_x(j); 
      for (int x = 0; x < number_of_grid_X; x++) {
        for (int y = 0; y < number_of_grid_Y; y++) {
          sum += bins2D[x][y].density * cos(w_u * x) * cos(w_v * y);
        }
      }
      a[i][j] = 1.0 / number_of_grid_X / number_of_grid_Y * sum; 
    }
  }

  for (int x = 0; x < number_of_grid_X; x++) {
    for (int y = 0; y < number_of_grid_Y; y++) {
      
    }
  }
  // cout << "Final" << endl;
  // cout << "GRAD "<<endl;

  for(auto &inst : instance_pointers_) {
    int instNum = instMap.find(inst->getName())->second;

    double gradInstX = 0.0, gradInstY = 0.0;
    for(auto &pin : inst->getPins()) {
      gradInstX += pin->gradWAX;
      gradInstY += pin->gradWAY;
    }
    pair<int, int> coordinate = inst->getCoordinate();
    int x = coordinate.first;
    int y = coordinate.second;

    double sum_x = 0.0;
    double sum_y = 0.0;

    for (int u = 0; u < number_of_grid_X; u++) {
      double w_u = get_w_x(u);
      double w_u_2 = w_u * w_u;
      for (int v = 0; v < number_of_grid_Y; v++) {
        double w_v = get_w_x(v); 
        double w_v_2 = w_v * w_v;
        double coeff;
        if(w_u_2 == 0 && w_v_2 == 0) coeff = 0.0;
        else coeff = a[u][v] / (w_u_2 + w_v_2);
        sum_x += coeff * w_u * sin(w_u * x) * cos(w_v * y);
        sum_y += coeff * w_v * cos(w_u * x) * sin(w_v * y);
      }
    }

    gradX[instNum] += gradInstX + lambda * (inst->getArea() * sum_x);
    gradY[instNum] += gradInstY + lambda * (inst->getArea() * sum_y);
  }
}

bool Circuit::densityCheck(double normal_bin_width, double normal_bin_height) {
  // Bin & inst update
  vector<vector<Bin> > bins2D(number_of_grid_X, vector<Bin> (number_of_grid_Y));

  for(auto &inst : instance_pointers_) {
    int instID = instMap.find(inst->getName())->second;
    int position_x = inst->getCoordinate().first;
    int position_y = inst->getCoordinate().second;
    uint instWidth = inst->getWidth();
    uint instHeight = inst->getHeight();
    int bin_coordinate_x = floor(position_x / normal_bin_width);
    int bin_coordinate_y = floor(position_y / normal_bin_height);
    inst->binCoordinate = make_pair(bin_coordinate_x, bin_coordinate_y);
    // cout << position_x << " == " << instWidth << " == " << bin_coordinate_x << " == " <<normal_bin_width<<endl;
    double _overflowX = (double)position_x + (double)instWidth - (double)(bin_coordinate_x + 1) * (double)normal_bin_width;
    double overflowX = (double)_overflowX / (double)instWidth;
    double _overflowY = (double)position_y + (double)instHeight - (double)(bin_coordinate_y + 1) * (double)normal_bin_height;
    double overflowY = (double)_overflowY / (double)instHeight;
    // cout << overflowX << " and "<< overflowY <<endl;
    // if(overflowX > 0 || overflowY > 0) {
    //   cout << overflowX << " and "<< overflowY <<endl;
    // }
    if(overflowX > 0) {
      if(overflowY > 0) {
        inst->binType = 3;
        inst->overflowX = overflowX;
        inst->overflowY = overflowY;
        bins2D[bin_coordinate_x][bin_coordinate_y].cell_area += (double)inst->getArea() * (1 - overflowX) * (1 - overflowY);
        bins2D[bin_coordinate_x][bin_coordinate_y + 1].cell_area += (double)inst->getArea() * (1 - overflowX) * overflowY;        
        bins2D[bin_coordinate_x + 1][bin_coordinate_y].cell_area += (double)inst->getArea() * overflowX * (1 - overflowY);          
        bins2D[bin_coordinate_x + 1][bin_coordinate_y+1].cell_area += (double)inst->getArea() * overflowX * overflowY;     
      }
      else {
        inst->binType = 1;
        inst->overflowX = overflowX;
        inst->overflowY = 0.0;
        bins2D[bin_coordinate_x][bin_coordinate_y].cell_area += (double)inst->getArea() * (1 - overflowX);
        bins2D[bin_coordinate_x + 1][bin_coordinate_y].cell_area += (double)inst->getArea() * overflowX;
      }
    }
    else {
      if(overflowY > 0) {
        inst->binType = 2;
        inst->overflowX = 0.0;
        inst->overflowY = overflowY;
        bins2D[bin_coordinate_x][bin_coordinate_y].cell_area += (double)inst->getArea() * (1 - overflowY);
        bins2D[bin_coordinate_x][bin_coordinate_y + 1].cell_area += (double)inst->getArea() * overflowY;          
      }
      else {
        inst->binType = 0;
        inst->overflowX = 0.0;
        inst->overflowY = 0.0;
        bins2D[bin_coordinate_x][bin_coordinate_y].cell_area += (double)inst->getArea();
      }        
    }
  }
  // Calc a
  bool results = false;
  double worstDensity = 0.0;
  double sumDensity = 0.0;
  for (int i = 0; i < number_of_grid_X; i++) {
    for (int j = 0; j < number_of_grid_Y; j++) {
      double density = bins2D[i][j].cell_area / ((double)normal_bin_width * normal_bin_height);
      if(density >= 1.2) {
        if(worstDensity < density) worstDensity = density;
        sumDensity += density;
        results = true;
      }
    }
  }
  cout << " WORST " << worstDensity << " Neg sum : " << sumDensity <<endl;
  return results;
}

void Circuit::myPlacement() {
  // Do 
  clock_t start, end;
  start = clock();
  uint die_width = die_->getWidth();
  uint die_height = die_->getHeight();
  double normal_bin_width = (double)die_width / number_of_grid_X;
  double normal_bin_height = (double)die_height / number_of_grid_Y; 
  double dieArea = (double)die_width * die_height;

  for (auto &inst : instance_pointers_) {
    inst->setCoordinate(int(die_width/2 - inst->getWidth()/2), int(die_height/2  - inst->getHeight()/2));
  }
  // Filler insertion
  // double totalAreaFC = 1.2 * (die_->getArea() - total_cell_area) - total_cell_area;
  // double areaFC = total_cell_area / instance_pointers_.size();
  // int cntFC = (int)(totalAreaFC/areaFC);
  // for(int i = 0;i<cntFC; i++) {
  //   Instance *filler = new Instance();
  //   filler.
  //   instance_pointers_.push_back(&instance);
  //   data_mapping_.inst_map[instance.getDbInst()] = &instance;
  // }
  // Give ID to instances
  int instCnt = 0;
  for (auto &inst : instance_pointers_) {
    instMap.insert(make_pair(inst->getName(), instCnt++));
  }
  
  // Iterate until 
  bool condition = true;

  int iter = 0;
  double alpha_max = 0.044 * normal_bin_width;
  double lambda_0 = this->initLambda();
  double prev_a = 1.0;
  double prev_alpha = 0.044 * normal_bin_width;
  cout << "Init Lambda "<< lambda_0 << endl;

  double prev_lambda = lambda_0;
  double mew_0 = 1.1;

  long long prevHPWL = 0, HPWL = 0;
  prevHPWL = 0;
  for (Net *net : net_pointers_) {
    prevHPWL += (long long) net->getHPWL();
  }
  HPWL = prevHPWL;

  vector<double> prev_uX(instance_pointers_.size());
  vector<double> prev_uY(instance_pointers_.size());
  vector<double> prev_vX(instance_pointers_.size());
  vector<double> prev_vY(instance_pointers_.size());

  vector<double> uX(instance_pointers_.size());
  vector<double> uY(instance_pointers_.size());
  vector<double> vX(instance_pointers_.size());
  vector<double> vY(instance_pointers_.size());

  vector<double> prev_gradX(instance_pointers_.size());
  vector<double> prev_gradY(instance_pointers_.size());

  vector<double> gradX(instance_pointers_.size());
  vector<double> gradY(instance_pointers_.size());

  for(auto &inst : instance_pointers_) {
    int instNum = instMap.find(inst->getName())->second;
    pair<int, int> coordinate = inst->getCoordinate();
    vX[instNum] = coordinate.first;
    uX[instNum] = coordinate.first;
    vY[instNum] = coordinate.second;
    uY[instNum] = coordinate.second;
  }

  // cout << "start"<<endl;
  while(condition) {
    // Update coeff
    double mew = 1.1;
    double diff_ = 1.0 - (double)(HPWL - prevHPWL)/3.5e5;
    // cout << HPWL << " " << prevHPWL<< " " << diff_<< endl;
    if(diff_ <= -3.1) mew = 0.75;
    else if(diff_ >= 1) mew = 1.1;
    else {
      mew = clamp(pow(mew_0, diff_), 0.75, 1.1);
    }
    double lambda = mew * prev_lambda;
    // cout << lambda << " lambda " <<prev_lambda << endl;
    
    prevHPWL = HPWL;
    prev_lambda = lambda;

    calcGradient(gradX, gradY, lambda);
    
    // Initial grad
    if(iter == 0) {           
      for(auto &inst : instance_pointers_) {
        int instNum = instMap.find(inst->getName())->second;
        double H = 1.0 / (inst->getPins().size() + lambda * inst->getArea());
        prev_gradX[instNum] = H * gradX[instNum];
        prev_gradY[instNum] = H * gradY[instNum];
      }
    }

    // Back Tracking
    double sum = 0.0;
    double alpha_hat = calcAlphaHat(vX, vY, prev_vX, prev_vY, gradX, gradY, prev_gradX, prev_gradY);
    
    double epsilon = 0.95;
    double alpha_k_max = max(alpha_max, 2 * prev_alpha);
    double next_ak;

    vector<double> uX_hat(instance_pointers_.size()), uY_hat(instance_pointers_.size());
    vector<double> vX_hat(instance_pointers_.size()), vY_hat(instance_pointers_.size());

    for(int i = 0; i < vX_hat.size(); i++) {
      vX_hat[i] = vX[i] + (-alpha_hat) * gradX[i];
      vY_hat[i] = vY[i] + (-alpha_hat) * gradY[i];
    }
    placeMap(vX_hat, vY_hat);

    vector<double> next_gradX(instance_pointers_.size());
    vector<double> next_gradY(instance_pointers_.size());
    calcGradient(next_gradX, next_gradY, lambda);
    next_ak = calcAlphaHat(vX_hat, vY_hat, vX, vY, next_gradX, next_gradY, gradX, gradY);
    int wi = 0;
    while(alpha_hat > epsilon * next_ak && wi < 100) {
      wi++;
      alpha_hat = next_ak;

      for(int i = 0; i < vX_hat.size(); i++) {
        vX_hat[i] = vX[i] + (-alpha_hat) * gradX[i];
        vY_hat[i] = vY[i] + (-alpha_hat) * gradY[i];
      }
      placeMap(vX_hat, vY_hat);

      vector<double> next_gradX(instance_pointers_.size());
      vector<double> next_gradY(instance_pointers_.size());
      calcGradient(next_gradX, next_gradY, lambda);

      next_ak = calcAlphaHat(vX_hat, vY_hat, vX, vY, next_gradX, next_gradY, gradX, gradY);
      if(alpha_hat < 0.01 * alpha_k_max || alpha_hat > alpha_k_max) break;
      // cout << "WHILE " << alpha_hat << " : " << epsilon * next_ak << endl;
    }
    alpha_hat = clamp(alpha_hat, 0.01*alpha_k_max, alpha_k_max);

    // Nestrov
    for(int i = 0; i < uX.size(); i++) {
      uX_hat[i] = vX[i] + (-alpha_hat) * gradX[i];
      uY_hat[i] = vY[i] + (-alpha_hat) * gradY[i];
    } 
    double next_a = (1 + sqrt(4 * prev_a * prev_a + 1)) / 2.0;
    for(int i = 0; i < vX.size(); i++) {
      vX_hat[i] = uX_hat[i] + (prev_a - 1) * (uX_hat[i] - uX[i]) / next_a;
      vY_hat[i] = uY_hat[i] + (prev_a - 1) * (uY_hat[i] - uY[i]) / next_a;
    } 
    prev_alpha = alpha_hat;
    prev_a = next_a;
    for(int i = 0; i < vX.size(); i++) {
      prev_uX[i] = uX[i];
      prev_uY[i] = uY[i];
      prev_vX[i] = vX[i];
      prev_vY[i] = vY[i];
      uX[i] = uX_hat[i];
      uY[i] = uY_hat[i];
      vX[i] = vX_hat[i];
      vY[i] = vY_hat[i];
      prev_gradX[i] = gradX[i];
      prev_gradY[i] = gradY[i];
    }
    placeMap(vX_hat, vY_hat);

    end = clock();
    double tresult = (double)(end-start)/CLOCKS_PER_SEC;
    tresult = (double) (end-start)/CLOCKS_PER_SEC;
    // Update coeff
    HPWL = 0;
    for (Net *net : net_pointers_) {
      HPWL += (long long) net->getHPWL();
    }
    condition = densityCheck(normal_bin_width, normal_bin_height);

    cout << "iter " << iter++ << " HPWL : " << HPWL << "\tTIME : " << tresult << endl;

    string img_file_name = "result" + to_string(iter);
    saveImg(img_file_name);
  }
}
}
