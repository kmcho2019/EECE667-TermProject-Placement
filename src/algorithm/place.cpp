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
#define number_of_grid_X 32
#define number_of_grid_Y 32


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
  vector<float> vector_data;
  vector<float> vector_bx(instCnt);
  vector<float> vector_by(instCnt);
  
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
    float diagA = 0;
    float bx = 0;
    float by = 0;
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
                diagA += (float)netWeight;
                // cout << "inst pin : " << netWeight << ", " << inst->getName() << " and " << connectedPin->getInstance()->getName() << endl;
              }
            }
            // When connected pin is in blockPin
            else if(connectedPin->isBlockPin()) {
              // // cout << "Block : " << connectedPin->getCoordinate().first << connectedPin->getCoordinate().second <<endl;
              bx += netWeight * (connectedPin->getCoordinate().first - pin->getCoordinate().first);
              by += netWeight * (connectedPin->getCoordinate().second - pin->getCoordinate().second);
              diagA += (float)netWeight;
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
  //   vector_bx[inst1Num] = (float)bx;
  //   vector_by[inst1Num] = (float)by;
  // }
  
  // // Update A and b for matrixSolver
  // A.n = instCnt;
  // A.nnz = vector_data.size();
  // A.row.resize(vector_row_idx.size());
  // A.col.resize(vector_col_idx.size());
  // A.dat.resize(vector_data.size());

  // A.row = valarray<int>(vector_row_idx.data(), A.nnz);
  // A.col = valarray<int>(vector_col_idx.data(), A.nnz);
  // A.dat = valarray<float>(vector_data.data(), A.nnz);
  // valarray<float> x(0.0, A.n);
  // valarray<float> y(0.0, A.n);
  // valarray<float> bx(vector_bx.data(), A.n);
  // valarray<float> by(vector_by.data(), A.n);
  // cout << "A constructed" <<endl;
  // // Solve Ax = b
  // A.solve(bx, x);
  // A.solve(by, y);
  // cout << "A solved" <<endl;

  // // Place with solutions
  // for (auto &inst : instance_pointers_) {
  //   int inst1Num = instMap.find(inst->getName())->second;
  //   inst->setCoordinate(x[inst1Num], y[inst1Num]);
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
std::pair<T,T> operator*(const float &mul, const std::pair<T,T> & r) {   
  return {mul * r.first, mul * r.second};   
} 

static float fastExp(float a) {
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

float get_w_x(float x) {
  return 2.0 * 3.141592 / (float)number_of_grid_X * x;
}

float Net::getHPWL_WA(float waCoeffX, float waCoeffY) {
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

  for(auto &pin : getConnectedPins()) {
    pair<int, int> pinCoordinate = pin->getCoordinate();

    float wPinX = pinCoordinate.first * waCoeffX;
    float wPinY = pinCoordinate.second * waCoeffY;

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

void Net::calcHPWL_gradWA(float waCoeffX, float waCoeffY) {
  int maxX = 0, maxY = 0;
  int minX = 2100000000, minY = 2100000000;
  for(auto &pin : this->getConnectedPins()) {
    pair<int, int> pinCoordinate = pin->getCoordinate();
    if(maxX < pinCoordinate.first) maxX = pinCoordinate.first;
    if(maxY < pinCoordinate.second) maxY = pinCoordinate.second;
    if(minX > pinCoordinate.first) minX = pinCoordinate.first;
    if(minY > pinCoordinate.second) minY = pinCoordinate.second;
  }

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
    float wMaxPinX = (pinCoordinate.first - maxX) * waCoeffX;
    float wMaxPinY = (pinCoordinate.second - maxY) * waCoeffY;
    float wMinPinX = -(pinCoordinate.first - minX) * waCoeffX;
    float wMinPinY = -(pinCoordinate.second - minY) * waCoeffY;

    sumCPosX += pinCoordinate.first * fastExp(wMaxPinX);
    sumCPosY += pinCoordinate.second * fastExp(wMaxPinY);
    sumCNegX += pinCoordinate.first * fastExp(wMinPinX);
    sumCNegY += pinCoordinate.second * fastExp(wMinPinY);
    sumPosX += fastExp(wMaxPinX);
    sumPosY += fastExp(wMaxPinY);
    sumNegX += fastExp(wMinPinX);
    sumNegY += fastExp(wMinPinY);
  }
  
  for(auto &pin : this->getConnectedPins()) {
    pair<int, int> pinCoordinate = pin->getCoordinate();
    float wMaxPinX = (pinCoordinate.first - maxX) * waCoeffX;
    float wMaxPinY = (pinCoordinate.second - maxY) * waCoeffY;
    float wMinPinX = -(pinCoordinate.first - minX) * waCoeffX;
    float wMinPinY = -(pinCoordinate.second - minY) * waCoeffY;

    float ePoxX = fastExp(wMaxPinX);
    float ePoxY = fastExp(wMaxPinY);
    float eNegX = fastExp(wMinPinX);
    float eNegY = fastExp(wMinPinY);

    pin->gradWAX = ((1.0 + pinCoordinate.first * waCoeffX) * sumPosX - waCoeffX * sumCPosX) * ePoxX / sumPosX / sumPosX - ((1.0 - pinCoordinate.first * waCoeffX) * sumNegX + waCoeffX * sumCNegX) * eNegX / sumNegX / sumNegX; 
    pin->gradWAY = ((1.0 + pinCoordinate.first * waCoeffY) * sumPosY - waCoeffY * sumCPosY) * ePoxY / sumPosY / sumPosY - ((1.0 - pinCoordinate.first * waCoeffY) * sumNegY + waCoeffY * sumCNegY) * eNegY / sumNegY / sumNegY; 
  }
}

float calcAlphaHat(vector<float> &vX, vector<float> &vY, vector<float> &prev_vX, vector<float> &prev_vY, vector<float> &gradX, vector<float> &gradY, vector<float> &prev_gradX, vector<float> &prev_gradY) {
  float sum = 0.0;
  for(int i = 0; i < vX.size(); i++) {
    float diffX = (vX[i] - prev_vX[i]);
    float diffY = (vY[i] - prev_vY[i]);
    sum += diffX * diffX;
    sum += diffY * diffY;
  }
  float normV = sqrt(sum);

  sum = 0.0;
  for(int i = 0; i < gradX.size(); i++) {
    float diffX = (gradX[i] - prev_gradX[i]);
    float diffY = (gradY[i] - prev_gradY[i]);
    
    sum += diffX * diffX;
    sum += diffY * diffY;
  }
  float normgradV = sqrt(sum);

  float result = 0.0;
  if(normgradV != 0.0) result = normV / normgradV;

  return result;
}

float Circuit::initLambda() {
  uint die_width = die_->getWidth();
  uint die_height = die_->getHeight();
  float normal_bin_width = (float)die_width / number_of_grid_X;
  float normal_bin_height = (float)die_height / number_of_grid_Y; 

  float gamma_ = 80.0 * normal_bin_width;
  float gamma = 1.0/gamma_;

  vector<vector<float> > a(number_of_grid_X, vector<float> (number_of_grid_Y));

  for(auto &net : net_pointers_) {
    for(auto &pin : net->getConnectedPins()) {
      pin->gradWAX = 0;
      pin->gradWAY = 0;
    }
  }

  // HPWL WA
  float HPWL_WA = 0.0;
  for(auto &net : net_pointers_) {
    net->calcHPWL_gradWA(gamma, gamma);
  }

  for (int i = 0; i <= number_of_grid_X; ++i) {
    vector<Bin> bins1D;
    for (int j = 0; j <= number_of_grid_Y; ++j) {
      pair<int, int> lower_left{i * normal_bin_width, j * normal_bin_height};
      pair<int, int> upper_right{(i + 1) * normal_bin_width, (j + 1) * normal_bin_height};
      bins1D.emplace_back(floor(normal_bin_width) * floor(normal_bin_height), lower_left, upper_right);
    }
    bins2D.push_back(bins1D);
  }
  // get utility in each bins
  for (Instance *instance : instance_pointers_) {
    pair<int, int> instance_lower_left = instance->getCoordinate();
    pair<int, int> instance_upper_right{
        instance_lower_left.first + instance->getWidth(),
        instance_lower_left.second + instance->getHeight()
    };
    for (int i = 0; i <= number_of_grid_X; ++i) {
      for (int j = 0; j <= number_of_grid_Y; ++j) {
        bins2D[i][j].stdArea = 0;
        bins2D[i][j].fillerArea = 0;
      }
    }
    int left_idx = static_cast<int>(instance_lower_left.first / normal_bin_width);
    int right_idx = static_cast<int>(instance_upper_right.first / normal_bin_width);
    int lower_idx = static_cast<int>(instance_lower_left.second / normal_bin_width);
    int upper_idx = static_cast<int>(instance_upper_right.second / normal_bin_width);
    for (int j = left_idx; j <= right_idx; ++j) {
      for (int k = lower_idx; k <= upper_idx; ++k) {
        bins2D.at(j).at(k).getOverlapArea(instance);
      }
    }
  }
  for (int i = 0; i < number_of_grid_X; i++) {
    for (int j = 0; j < number_of_grid_Y; j++) {
      float density = (bins2D[i][j].stdArea + bins2D[i][j].fillerArea) / (normal_bin_width * normal_bin_height);
      bins2D[i][j].density = density;
    }
  }
  // For FFT
  wx.resize(number_of_grid_X);
  wx_sq.resize(number_of_grid_X);
  wy.resize(number_of_grid_X);
  wy_sq.resize(number_of_grid_X);
  for (int i = 0; i < number_of_grid_X; i++) {
    float temp = get_w_x(i);
    wx[i] = temp;
    wx_sq[i] = temp * temp;
    wy[i] = temp;
    wy_sq[i] = temp * temp;
  }
  float** binDensity = new float *[number_of_grid_X];
  float** electricPotential = new float *[number_of_grid_X];
  float** electricForceX = new float *[number_of_grid_X];
  float** electricForceY = new float *[number_of_grid_X];
  cosTable.resize(number_of_grid_X * 3 / 2, 0);
  workArea_.resize(round(sqrt(number_of_grid_X)) + 2, 0);

  for (int i = 0; i < number_of_grid_X; i++) {
    binDensity[i] = new float[number_of_grid_Y];
    electricPotential[i] = new float[number_of_grid_Y];
    electricForceX[i] = new float[number_of_grid_Y];
    electricForceY[i] = new float[number_of_grid_Y];

    for (int j = 0; j < number_of_grid_Y; j++) {
      binDensity[i][j] = bins2D[i][j].density;
      electricPotential[i][j] = bins2D[i][j].electricPotential;
      electricForceX[i][j] = bins2D[i][j].electricField_x;
      electricForceY[i][j] = bins2D[i][j].electricField_y;
    }
  }
  ddct2d(number_of_grid_X,number_of_grid_Y, -1, binDensity, NULL, workArea_.data(), cosTable.data());

  for (int i = 0; i < number_of_grid_X; ++i) {
    binDensity[i][0] *= 0.5;
  }
  for (int i = 0; i < number_of_grid_Y; ++i) {
    binDensity[0][i] *= 0.5;
  }
  for (int i = 0; i < number_of_grid_X; ++i) {
    for (int j = 0; j < number_of_grid_Y; ++j) {
      binDensity[i][j] *= 4.0 / number_of_grid_X / number_of_grid_Y;
    }
  }
  for (int i = 0; i < number_of_grid_X; i++) {
    float wx_ = wx[i];
    float wx2 = wx_sq[i];

    for (int j = 0; j < number_of_grid_Y; j++) {
      float wy_ = wy[j];
      float wy2 = wy_sq[j];

      float density = binDensity[i][j];
      float phi = 0;
      float electroX = 0, electroY = 0;

      if (i == 0 && j == 0) {
        phi = electroX = electroY = 0.0f;
      } else {
        phi = density / (wx2 + wy2);
        electroX = phi * wx_;
        electroY = phi * wy_;
      }
      electricPotential[i][j] = phi;
      electricForceX[i][j] = electroX;
      electricForceY[i][j] = electroY;
    }
  }

  ddct2d(number_of_grid_X, number_of_grid_Y, 1, electricPotential, NULL, workArea_.data(), cosTable.data());
  ddsct2d(number_of_grid_X, number_of_grid_Y, 1, electricForceX, NULL, workArea_.data(), cosTable.data());
  ddcst2d(number_of_grid_X, number_of_grid_Y, 1, electricForceY, NULL, workArea_.data(), cosTable.data());

  for (int i = 0; i < number_of_grid_X; i++) {
    for (int j = 0; j < number_of_grid_Y; j++) {
      bins2D[i][j].density = binDensity[i][j];
      bins2D[i][j].electricPotential = electricPotential[i][j];
      bins2D[i][j].electricField_x = electricForceX[i][j];
      bins2D[i][j].electricField_y = electricForceY[i][j];
    }
    delete[] binDensity[i];
    delete[] electricPotential[i];
    delete[] electricForceX[i];
    delete[] electricForceY[i];
  }
  delete[] binDensity;
  delete[] electricPotential;
  delete[] electricForceX;
  delete[] electricForceY;

  float sumWA = 0.0, sumDensity = 0.0;
  for(auto &inst : instance_pointers_) {
    int instNum = instMap.find(inst->getName())->second;
    float gradInstX = 0.0, gradInstY = 0.0;

    for(auto &pin : inst->getPins()) {
      gradInstX += pin->gradWAX;
      gradInstY += pin->gradWAY;
    }
    sumWA += abs(gradInstX);
    sumWA += abs(gradInstY);

    pair<int, int> coordinate = inst->binCoordinate;
    sumDensity += abs((double)inst->getArea() * bins2D[coordinate.first][coordinate.second].electricField_x);
    sumDensity += abs((double)inst->getArea() * bins2D[coordinate.first][coordinate.second].electricField_y);
  }
  return sumWA/sumDensity;
}

void Circuit::calcGradient(vector<float> &gradX, vector<float> &gradY, float lambda) {
  // Bin Construction
  uint die_width = die_->getWidth();
  uint die_height = die_->getHeight();
  float normal_bin_width = (float)die_width / number_of_grid_X;
  float normal_bin_height = (float)die_height / number_of_grid_Y;
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
  for (int i = 0; i <= number_of_grid_X; ++i) {
    for (int j = 0; j <= number_of_grid_Y; ++j) {
      bins2D[i][j].stdArea = 0;
      bins2D[i][j].fillerArea = 0;
    }
  }
  // get utility in each bins
  for (Instance *instance : instance_pointers_) {
    pair<int, int> instance_lower_left = instance->getCoordinate();
    pair<int, int> instance_upper_right{
        instance_lower_left.first + instance->getWidth(),
        instance_lower_left.second + instance->getHeight()
    };

    int left_idx = static_cast<int>(instance_lower_left.first / normal_bin_width);
    int right_idx = static_cast<int>(instance_upper_right.first / normal_bin_width);
    int lower_idx = static_cast<int>(instance_lower_left.second / normal_bin_height);
    int upper_idx = static_cast<int>(instance_upper_right.second / normal_bin_height);
    for (int j = left_idx; j <= right_idx; ++j) {
      for (int k = lower_idx; k <= upper_idx; ++k) {
        bins2D.at(j).at(k).getOverlapArea(instance);
      }
    }
  }
  for (int i = 0; i < number_of_grid_X; i++) {
    for (int j = 0; j < number_of_grid_Y; j++) {
      float density = (bins2D[i][j].stdArea + bins2D[i][j].fillerArea);
      bins2D[i][j].density = density;
    }
  }
  vector<vector<float> > a(number_of_grid_X, vector<float> (number_of_grid_Y));
  cosTable.resize(number_of_grid_X * 3 / 2, 0);
  workArea_.resize(round(sqrt(number_of_grid_X)) + 2, 0);

  float** binDensity = new float *[number_of_grid_X];
  float** electricPotential = new float *[number_of_grid_X];
  float** electricForceX = new float *[number_of_grid_X];
  float** electricForceY = new float *[number_of_grid_X];

  for (int i = 0; i < number_of_grid_X; i++) {
    binDensity[i] = new float[number_of_grid_Y];
    electricPotential[i] = new float[number_of_grid_Y];
    electricForceX[i] = new float[number_of_grid_Y];
    electricForceY[i] = new float[number_of_grid_Y];

    for (int j = 0; j < number_of_grid_Y; j++) {
      binDensity[i][j] = bins2D[i][j].density;
//      electricPotential[i][j] = bins2D[i][j].electricPotential;
//      electricForceX[i][j] = bins2D[i][j].electricField_x;
//      electricForceY[i][j] = bins2D[i][j].electricField_y;
      electricPotential[i][j] = 0.0f;
      electricForceX[i][j] = 0.0f;
      electricForceY[i][j] = 0.0f;
    }
  }

//  ddct2d(number_of_grid_X,number_of_grid_Y, -1, binDensity, NULL, workArea_.data(), cosTable.data());
  ddct2d(number_of_grid_X, number_of_grid_Y, -1, binDensity, NULL, (int*) &workArea_[0], (float*) &cosTable[0]);

  for (int i = 0; i < number_of_grid_X; ++i) {
    binDensity[i][0] *= 0.5;
  }
  for (int i = 0; i < number_of_grid_Y; ++i) {
    binDensity[0][i] *= 0.5;
  }
  for (int i = 0; i < number_of_grid_X; ++i) {
    for (int j = 0; j < number_of_grid_Y; ++j) {
      binDensity[i][j] *= 4.0 / number_of_grid_X / number_of_grid_Y;
    }
  }
  for (int i = 0; i < number_of_grid_X; i++) {
    float wx_ = wx[i];
    float wx2 = wx_sq[i];

    for (int j = 0; j < number_of_grid_Y; j++) {
      float wy_ = wy[j];
      float wy2 = wy_sq[j];

      float density = binDensity[i][j];
      float phi = 0;
      float electroX = 0, electroY = 0;

      if (i == 0 && j == 0) {
        phi = electroX = electroY = 0.0f;
      } else {
        phi = density / (wx2 + wy2);
        electroX = phi * wx_;
        electroY = phi * wy_;
      }
      electricPotential[i][j] = phi;
      electricForceX[i][j] = electroX;
      electricForceY[i][j] = electroY;
    }
  }

  ddct2d(number_of_grid_X, number_of_grid_Y, 1, electricPotential, NULL, (int*) &workArea_[0], (float*) &cosTable[0]);
  ddsct2d(number_of_grid_X, number_of_grid_Y, 1, electricForceX, NULL, (int*) &workArea_[0], (float*) &cosTable[0]);
  ddcst2d(number_of_grid_X, number_of_grid_Y, 1, electricForceY, NULL, (int*) &workArea_[0], (float*) &cosTable[0]);

  for (int i = 0; i < number_of_grid_X; i++) {
    for (int j = 0; j < number_of_grid_Y; j++) {
      bins2D[i][j].density = binDensity[i][j];
      bins2D[i][j].electricPotential = electricPotential[i][j];
      bins2D[i][j].electricField_x = electricForceX[i][j];
      bins2D[i][j].electricField_y = electricForceY[i][j];
    }
    delete[] binDensity[i];
    delete[] electricPotential[i];
    delete[] electricForceX[i];
    delete[] electricForceY[i];
  }
  delete[] binDensity;
  delete[] electricPotential;
  delete[] electricForceX;
  delete[] electricForceY;

  // HPWL WA
  float maxDensity = 0.0;
  for (int x = 0; x < number_of_grid_X; x++) {
    for (int y = 0; y < number_of_grid_Y; y++) {
      if(maxDensity < bins2D[x][y].density) maxDensity = bins2D[x][y].density;
    }
  }
  float tau = max(maxDensity - 1.2 , 0.0) * (float)normal_bin_width * (float)normal_bin_height / (float)total_cell_area;
  float gamma_ = 8.0 * normal_bin_width * pow(10, 20.0/9.0 * tau - 11.0/9.0);
  float gamma = 1.0/gamma_;
  float HPWL_WA = 0.0;
  for(auto &net : net_pointers_) {
    if (net->getSignalType() != "POWER" && net->getSignalType() != "GROUND" && net->getSignalType() != "CLOCK" && net->getSignalType() != "RESET") {
      net->calcHPWL_gradWA(gamma, gamma);
    }
  }
  float sum = 0.0f;
  float absSum = 0.0f;
  for(auto &inst : instance_pointers_) {
    int instNum = instMap.find(inst->getName())->second;

    float gradInstX = 0.0, gradInstY = 0.0;
    for(auto &pin : inst->getPins()) {
      gradInstX += pin->gradWAX;
      gradInstY += pin->gradWAY;
    }
    pair<int, int> coordinate = inst->binCoordinate;
    float temp = lambda * ((double)inst->getArea() * bins2D[coordinate.first][coordinate.second].electricField_x);
    sum += temp;
    absSum += abs(temp);
    gradX[instNum] += gradInstX + lambda * ((double)inst->getArea() * bins2D[coordinate.first][coordinate.second].electricField_x);
    gradY[instNum] += gradInstY + lambda * ((double)inst->getArea() * bins2D[coordinate.first][coordinate.second].electricField_y);
  }
//  cout<<"END"<<endl;
}

bool Circuit::densityCheck(float normal_bin_width, float normal_bin_height) {
  int die_width = die_->getWidth();
  int die_height = die_->getHeight();

  // Bin & inst update
  vector<vector<Bin>> bins2D;
  for (int i = 0; i <= number_of_grid_X; ++i) {
    vector<Bin> bins1D;
    for (int j = 0; j <= number_of_grid_Y; ++j) {
      pair<int, int> lower_left{i * normal_bin_width, j * normal_bin_height};
      pair<int, int> upper_right{(i + 1) * normal_bin_width, (j + 1) * normal_bin_height};
      bins1D.emplace_back(static_cast<int>(die_width * die_height), lower_left, upper_right);
    }
    bins2D.push_back(bins1D);
  }
  // get utility in each bins
  for (Instance *instance : instance_pointers_) {
    pair<int, int> instance_lower_left = instance->getCoordinate();
    pair<int, int> instance_upper_right{
        instance_lower_left.first + instance->getWidth(),
        instance_lower_left.second + instance->getHeight()
    };

    int left_idx = static_cast<int>(instance_lower_left.first / normal_bin_width);
    int right_idx = static_cast<int>(instance_upper_right.first / normal_bin_width);
    int lower_idx = static_cast<int>(instance_lower_left.second / normal_bin_width);
    int upper_idx = static_cast<int>(instance_upper_right.second / normal_bin_width);
    for (int j = left_idx; j <= right_idx; ++j) {
      for (int k = lower_idx; k <= upper_idx; ++k) {
        bins2D.at(j).at(k).getOverlapArea(instance);
      }
    }
  }

  // Calc a
  bool results = false;
  float worstDensity = 0.0;
  float sumDensity = 0.0;
  for (int i = 0; i < number_of_grid_X; i++) {
    for (int j = 0; j < number_of_grid_Y; j++) {
      float density = bins2D[i][j].stdArea / ((float)normal_bin_width * normal_bin_height);
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
  float normal_bin_width = (float)die_width / number_of_grid_X;
  float normal_bin_height = (float)die_height / number_of_grid_Y; 
  float dieArea = (float)die_width * die_height;

  // Give ID to instances
  int instCnt = 0;
  vector<int> fillerID;
  fillerID.resize(cntFC);
  int fillerCnt = 0;
//  cout << "Start"<<endl;
  for (auto &inst : instance_pointers_) {
    // cout << inst->name_ <<endl;
    if(inst->isFiller) fillerID[fillerCnt++] = instCnt;
    instMap.insert(make_pair(inst->getName(), instCnt++));    
  }
//  cout << "giveID";
  for (auto &inst : instance_pointers_) {
    if(inst->isFiller) continue;
    inst->setCoordinate(int(die_width/2 - inst->getWidth()/2), int(die_height/2  - inst->getHeight()/2));
  }
  // Iterate until 
  bool condition = true;

  int iter = 0;
  float alpha_max = 0.044 * normal_bin_width;
  float lambda_0 = initLambda();
  float prev_a = 1.0;
  float prev_alpha = 0.044 * normal_bin_width;
  cout << "Init Lambda "<< lambda_0 << endl;

  float prev_lambda = lambda_0;
  float mew_0 = 1.1;

  long long prevHPWL = 0, HPWL = 0;
  prevHPWL = 0;
  for (Net *net : net_pointers_) {
    prevHPWL += (long long) net->getHPWL();
  }
  HPWL = prevHPWL;

  vector<float> prev_uX(instance_pointers_.size());
  vector<float> prev_uY(instance_pointers_.size());
  vector<float> prev_vX(instance_pointers_.size());
  vector<float> prev_vY(instance_pointers_.size());

  vector<float> uX(instance_pointers_.size());
  vector<float> uY(instance_pointers_.size());
  vector<float> vX(instance_pointers_.size());
  vector<float> vY(instance_pointers_.size());

  vector<float> prev_gradX(instance_pointers_.size());
  vector<float> prev_gradY(instance_pointers_.size());

  vector<float> gradX(instance_pointers_.size());
  vector<float> gradY(instance_pointers_.size());

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
    float mew = 1.1;
    float diff_ = 1.0 - (float)(HPWL - prevHPWL)/3.5e6;
    // cout << HPWL << " " << prevHPWL<< " " << diff_<< endl;
    if(diff_ <= -3.1) mew = 0.75;
    else if(diff_ >= 1) mew = 1.1;
    else {
      mew = clamp(pow(mew_0, diff_), 0.75, 1.1);
    }
    float lambda = mew * prev_lambda;
    
    prevHPWL = HPWL;
    prev_lambda = lambda;

    calcGradient(gradX, gradY, lambda);
    
    // Initial grad
    if(iter == 0) {           
      for(auto &inst : instance_pointers_) {
        int instNum = instMap.find(inst->getName())->second;
        float H = 1.0 / (inst->getPins().size() + lambda * inst->getArea());
        prev_gradX[instNum] = H * gradX[instNum];
        prev_gradY[instNum] = H * gradY[instNum];
      }
    }

    // Back Tracking
    float sum = 0.0;
    float alpha_hat = calcAlphaHat(vX, vY, prev_vX, prev_vY, gradX, gradY, prev_gradX, prev_gradY);
    
    float epsilon = 0.95f;
    float alpha_k_max = max(alpha_max, 2 * prev_alpha);
    float next_ak;

    vector<float> uX_hat(instance_pointers_.size()), uY_hat(instance_pointers_.size());
    vector<float> vX_hat(instance_pointers_.size()), vY_hat(instance_pointers_.size());

    for(int i = 0; i < vX_hat.size(); i++) {
      vX_hat[i] = vX[i] + (-alpha_hat) * gradX[i];
      vY_hat[i] = vY[i] + (-alpha_hat) * gradY[i];
    }
    placeMap(vX_hat, vY_hat);

    vector<float> next_gradX(instance_pointers_.size());
    vector<float> next_gradY(instance_pointers_.size());
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

      vector<float> next_gradX(instance_pointers_.size());
      vector<float> next_gradY(instance_pointers_.size());
      calcGradient(next_gradX, next_gradY, lambda);

      next_ak = calcAlphaHat(vX_hat, vY_hat, vX, vY, next_gradX, next_gradY, gradX, gradY);
      if(alpha_hat < 0.01 * alpha_k_max || alpha_hat > alpha_k_max) break;
      // cout << "WHILE " << alpha_hat << " : " << epsilon * next_ak << endl;
    }
    alpha_hat = clamp(alpha_hat, 0.01f*alpha_k_max, alpha_k_max);

    // Nestrov
    for(int i = 0; i < uX.size(); i++) {
      uX_hat[i] = vX[i] + (-alpha_hat) * gradX[i];
      uY_hat[i] = vY[i] + (-alpha_hat) * gradY[i];
    } 
    float next_a = (1 + sqrt(4 * prev_a * prev_a + 1)) / 2.0;
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
    float tresult = (float)(end-start)/CLOCKS_PER_SEC;
    tresult = (float) (end-start)/CLOCKS_PER_SEC;
    // Update coeff
    HPWL = 0;
    for (Net *net : net_pointers_) {
      HPWL += (long long) net->getHPWL();
    }

    condition = densityCheck(normal_bin_width, normal_bin_height);

    if(iter % 5 == 0) {
      cout << "iter " << iter << " HPWL : " << HPWL << "\tTIME : " << tresult << endl;
      string img_file_name = "result" + to_string(iter);
      saveImg(img_file_name);
    }
    iter++;
  }
//  for(auto &id : fillerID) {
//    instance_pointers_.erase(instance_pointers_.begin()+id);
//  }
}
}
