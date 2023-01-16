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
#define number_of_grid_X 64
#define number_of_grid_Y 64

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
    for (auto &pin : inst->getPins()) {
      // cout << "In inst pin : "<<pin->getCoordinate().first << " : "<< pin->getCoordinate().second<<endl; 
      if (pin->getPinName() == "CK") continue;
      if (pin->getNet()) {
        Net *net = pin->getNet();
        if (net->getSignalType() != "POWER" && net->getSignalType() != "GROUND" && net->getSignalType() != "CLOCK"
            && net->getSignalType() != "RESET") {
          int netWeight = net->getWeight();
          if (net->getConnectedPins().size() > 10) continue;
          // For all connected pins
          for (auto &connectedPin : net->getConnectedPins()) {
            // When connected pin is in Instance
            if (connectedPin->isInstancePin()) {
              int inst2Num = instMap.find(connectedPin->getInstance()->getName())->second;
              if (inst1Num != inst2Num) {
                vector_row_idx.push_back(inst1Num);
                vector_col_idx.push_back(inst2Num);
                vector_data.push_back(-1.0 * netWeight);

                bx += netWeight * (connectedPin->getCoordinate().first - pin->getCoordinate().first);
                by += netWeight * (connectedPin->getCoordinate().second - pin->getCoordinate().second);
                diagA += (float) netWeight;
                // cout << "inst pin : " << netWeight << ", " << inst->getName() << " and " << connectedPin->getInstance()->getName() << endl;
              }
            }
              // When connected pin is in blockPin
            else if (connectedPin->isBlockPin()) {
              // // cout << "Block : " << connectedPin->getCoordinate().first << connectedPin->getCoordinate().second <<endl;
              bx += netWeight * (connectedPin->getCoordinate().first - pin->getCoordinate().first);
              by += netWeight * (connectedPin->getCoordinate().second - pin->getCoordinate().second);
              diagA += (float) netWeight;
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
template<typename T, typename U>
T maxClamp(const T x, const U max) {
  if (x > max) return max;
  return x;
}

template<typename T, typename U>
T minClamp(const T x, const U min) {
  if (x < min) return min;
  return x;
}

template<typename T, typename U>
T clamp(const T x, const U min, const U max) {
  if (x < min) return min;
  else if (x > max) return max;
  return x;
}

template<typename T, typename U>
std::pair<T, U> operator+(const std::pair<T, U> &l, const std::pair<T, U> &r) {
  return {l.first + r.first, l.second + r.second};
}

template<typename T, typename U>
std::pair<T, T> operator*(const U &mul, const std::pair<T, T> &r) {
  return {mul * r.first, mul * r.second};
}
template<typename T>
std::pair<T, T> operator*(const float &mul, const std::pair<T, T> &r) {
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
  return 2.0 * 3.141592 / (float) number_of_grid_X * x;
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

  vector<Pin *> pins;
  for (dbITerm *db_i_term : db_net_->getITerms()) {
    pins.push_back(data_mapping_->pin_map_i[db_i_term]);
  }
  for (dbBTerm *db_b_term : db_net_->getBTerms()) {
    pins.push_back(data_mapping_->pin_map_b[db_b_term]);
  }

  for (auto &pin : getConnectedPins()) {
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
  return sumCPosX / sumPosX - sumCNegX / sumNegX + sumCPosY / sumPosY - sumCNegY / sumNegY;
}

void Net::calcHPWL_gradWA(float waCoeffX, float waCoeffY) {
  float maxX = 0.0, maxY = 0.0;
  float minX = 2100000000.0, minY = 210000000.0;
  for (auto &pin : this->getConnectedPins()) {
    pair<int, int> pinCoordinate = pin->getCoordinate();
    if (maxX < pinCoordinate.first) maxX = (float)pinCoordinate.first;
    if (maxY < pinCoordinate.second) maxY = (float)pinCoordinate.second;
    if (minX > pinCoordinate.first) minX = (float)pinCoordinate.first;
    if (minY > pinCoordinate.second) minY = (float)pinCoordinate.second;
  }
  float cx = 0.5 * (maxX - minX);
  float cy = 0.5 * (maxY - minY);

  sumPosX = 0.0;
  sumNegX = 0.0;
  sumCPosX = 0.0;
  sumCNegX = 0.0;
  sumPosY = 0.0;
  sumNegY = 0.0;
  sumCPosY = 0.0;
  sumCNegY = 0.0;
  // cout << "calcHPWL_gradWA INIT" <<endl;
  for (auto &pin : this->getConnectedPins()) {
    pair<int, int> pinCoordinate = pin->getCoordinate();
    float wMaxPinX = (pinCoordinate.first - cx) * waCoeffX;
    float wMaxPinY = (pinCoordinate.second - cy) * waCoeffY;
    float wMinPinX = -(pinCoordinate.first - cx) * waCoeffX;
    float wMinPinY = -(pinCoordinate.second - cy) * waCoeffY;

    sumCPosX += pinCoordinate.first * fastExp(wMaxPinX);
    sumCPosY += pinCoordinate.second * fastExp(wMaxPinY);
    sumCNegX += pinCoordinate.first * fastExp(wMinPinX);
    sumCNegY += pinCoordinate.second * fastExp(wMinPinY);
    sumPosX += fastExp(wMaxPinX);
    sumPosY += fastExp(wMaxPinY);
    sumNegX += fastExp(wMinPinX);
    sumNegY += fastExp(wMinPinY);
  }
  int cnt = this->getConnectedPins().size();
  for (auto &pin : this->getConnectedPins()) {
    pair<int, int> pinCoordinate = pin->getCoordinate();
    float wMaxPinX = (pinCoordinate.first - cx) * waCoeffX;
    float wMaxPinY = (pinCoordinate.second - cy) * waCoeffY;
    float wMinPinX = -(pinCoordinate.first - cx) * waCoeffX;
    float wMinPinY = -(pinCoordinate.second - cy) * waCoeffY;

    float ePoxX = fastExp(wMaxPinX);
    float ePoxY = fastExp(wMaxPinY);
    float eNegX = fastExp(wMinPinX);
    float eNegY = fastExp(wMinPinY);

    pin->gradWAX = ((1.0 + pinCoordinate.first * waCoeffX) * sumPosX - waCoeffX * sumCPosX) * ePoxX / sumPosX / sumPosX
        - ((1.0 - pinCoordinate.first * waCoeffX) * sumNegX + waCoeffX * sumCNegX) * eNegX / sumNegX / sumNegX;
    pin->gradWAY = ((1.0 + pinCoordinate.first * waCoeffY) * sumPosY - waCoeffY * sumCPosY) * ePoxY / sumPosY / sumPosY
        - ((1.0 - pinCoordinate.first * waCoeffY) * sumNegY + waCoeffY * sumCNegY) * eNegY / sumNegY / sumNegY;
  }
}

float calcAlphaHat(vector<float> &vX,
                   vector<float> &vY,
                   vector<float> &prev_vX,
                   vector<float> &prev_vY,
                   vector<float> &gradX,
                   vector<float> &gradY,
                   vector<float> &prev_gradX,
                   vector<float> &prev_gradY) {
  float sum = 0.0;
  for (int i = 0; i < vX.size(); i++) {
    float diffX = (vX[i] - prev_vX[i]);
    float diffY = (vY[i] - prev_vY[i]);
    sum += diffX * diffX;
    sum += diffY * diffY;
  }
  float normV = sqrt(sum);

  sum = 0.0;
  for (int i = 0; i < gradX.size(); i++) {
    float diffX = (gradX[i] - prev_gradX[i]);
    float diffY = (gradY[i] - prev_gradY[i]);

    sum += diffX * diffX;
    sum += diffY * diffY;
  }
  if(sum - 0.0 < 1e-08) return 1e-5;
  float normgradV = sqrt(sum);

  return normV / normgradV;
}

float Circuit::initLambda() {
  uint die_width = die_->getWidth();
  uint die_height = die_->getHeight();
  float normal_bin_width = (float) die_width / number_of_grid_X;
  float normal_bin_height = (float) die_height / number_of_grid_Y;

  float gamma_ = 80.0 * normal_bin_width;
  float gamma = 1.0 / gamma_;

  vector<vector<float> > a(number_of_grid_X, vector<float>(number_of_grid_Y));

  for (auto &net : net_pointers_) {
    for (auto &pin : net->getConnectedPins()) {
      pin->gradWAX = 0;
      pin->gradWAY = 0;
    }
  }

  // HPWL WA
  float HPWL_WA = 0.0;
  for (auto &net : net_pointers_) {
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
        bins2D[i][j].reset();
      }
    }
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
  // For FFT
  wx.resize(number_of_grid_X);
  wx_sq.resize(number_of_grid_X);
  wy.resize(number_of_grid_Y);
  wy_sq.resize(number_of_grid_Y);
  for (int i = 0; i < number_of_grid_X; i++) {
    float temp = get_w_x(i);
    wx[i] = temp;
    wx_sq[i] = temp * temp;
    wy[i] = temp;
    wy_sq[i] = temp * temp;
  }
  float **binDensity = new float *[number_of_grid_X];
  float **electricPotential = new float *[number_of_grid_X];
  float **electricForceX = new float *[number_of_grid_X];
  float **electricForceY = new float *[number_of_grid_X];
  cosTable.resize(number_of_grid_X * 3 / 2, 0);
  workArea_.resize(round(sqrt(number_of_grid_X)) + 2, 0);

  for (int i = 0; i < number_of_grid_X; i++) {
    binDensity[i] = new float[number_of_grid_Y];
    electricPotential[i] = new float[number_of_grid_Y];
    electricForceX[i] = new float[number_of_grid_Y];
    electricForceY[i] = new float[number_of_grid_Y];

    for (int j = 0; j < number_of_grid_Y; j++) {
      binDensity[i][j] = bins2D[i][j].density / normal_bin_width / normal_bin_height;
      electricPotential[i][j] = bins2D[i][j].electricPotential;
      electricForceX[i][j] = bins2D[i][j].electricField_x;
      electricForceY[i][j] = bins2D[i][j].electricField_y;
    }
  }
  ddct2d(number_of_grid_X, number_of_grid_Y, -1, binDensity, NULL, workArea_.data(), cosTable.data());

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

// HPWL WA
  for (int x = 0; x < number_of_grid_X; x++) {
    for (int y = 0; y < number_of_grid_Y; y++) {
      for (auto &pair_inst : bins2D[x][y].overlappedInstance) {
        Instance *inst = pair_inst.first;
        float ratio = pair_inst.second / inst->getArea();
        inst->densityForce.first += bins2D[x][y].electricField_x * ratio;
        inst->densityForce.second += bins2D[x][y].electricField_y * ratio;
      }
    }
  }
  float sumWA = 0.0, sumDensity = 0.0;
  for (auto &inst : instance_pointers_) {
    float gradInstX = 0.0, gradInstY = 0.0;

    for (auto &pin : inst->getPins()) {
      gradInstX += pin->gradWAX;
      gradInstY += pin->gradWAY;
    }
    sumWA += abs(gradInstX);
    sumWA += abs(gradInstY);

    pair<int, int> coordinate = inst->binCoordinate;
    sumDensity += abs((double) inst->getArea() * inst->densityForce.first);
    sumDensity += abs((double) inst->getArea() * inst->densityForce.second);
  }
  return sumWA / sumDensity;
}

void Circuit::initialPlacement(vector<float> &gradX, vector<float> &gradY) {
  uint die_width = die_->getWidth();
  uint die_height = die_->getHeight();
  float normal_bin_width = (float) die_width / number_of_grid_X;
  float normal_bin_height = (float) die_height / number_of_grid_Y;

  float gamma_ = 80.0 * normal_bin_width;
  float gamma = 1.0 / gamma_;

  vector<vector<float> > a(number_of_grid_X, vector<float>(number_of_grid_Y));

  for (auto &net : net_pointers_) {
    for (auto &pin : net->getConnectedPins()) {
      pin->gradWAX = 0;
      pin->gradWAY = 0;
    }
  }

  // HPWL WA
  float HPWL_WA = 0.0;
  for (auto &net : net_pointers_) {
    net->calcHPWL_gradWA(gamma, gamma);
  }
  for (auto &inst : instance_pointers_) {
    if(inst->isFiller) continue;
    int instNum = instMap.find(inst->getName())->second;

    float gradInstX = 0.0, gradInstY = 0.0;
    for (auto &pin : inst->getPins()) {
      gradInstX += pin->gradWAX;
      gradInstY += pin->gradWAY;
    }
    gradX[instNum] = gradInstX;
    gradY[instNum] = gradInstY;
  }
}

void Circuit::calcGradient(vector<float> &gradX, vector<float> &gradY, float lambda) {
  // Bin Construction
  uint die_width = die_->getWidth();
  uint die_height = die_->getHeight();
  float normal_bin_width = (float) die_width / number_of_grid_X;
  float normal_bin_height = (float) die_height / number_of_grid_Y;
  for (auto &net : net_pointers_) {
    for (auto &pin : net->getConnectedPins()) {
      pin->gradWAX = 0;
      pin->gradWAY = 0;
    }
  }
  for (int i = 0; i < gradX.size(); i++) {
    gradX[i] = 0;
    gradY[i] = 0;
  }
  for (int i = 0; i <= number_of_grid_X; ++i) {
    for (int j = 0; j <= number_of_grid_Y; ++j) {
      bins2D[i][j].reset();
    }
  }
  // get utility in each bins
  for (Instance *instance : instance_pointers_) {
    pair<int, int> instance_lower_left = instance->getCoordinate();
    pair<int, int> instance_upper_right{
        instance_lower_left.first + instance->getWidth(),
        instance_lower_left.second + instance->getHeight()
    };
    instance->densityForce = make_pair(0.0, 0.0);
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
  float maxDensity = 0.0;
  for (int i = 0; i < number_of_grid_X; i++) {
    for (int j = 0; j < number_of_grid_Y; j++) {
      float density = (bins2D[i][j].stdArea + bins2D[i][j].fillerArea);
      bins2D[i][j].density = density;
      if (maxDensity < bins2D[i][j].density) maxDensity = bins2D[i][j].density;
    }
  }

  float
      tau = max(maxDensity - 1.2 * (float) normal_bin_width * (float) normal_bin_height, 0.0) / (float) total_cell_area;
  float gamma_ = 8.0 * normal_bin_width * pow(10, 20.0 / 9.0 * tau - 11.0 / 9.0);
  float gamma = 1.0 / gamma_;
  float HPWL_WA = 0.0;
  for (auto &net : net_pointers_) {
    if (net->getSignalType() != "POWER" && net->getSignalType() != "GROUND" && net->getSignalType() != "CLOCK"
        && net->getSignalType() != "RESET") {
      net->calcHPWL_gradWA(gamma, gamma);
    }
  }

  vector<vector<float> > a(number_of_grid_X, vector<float>(number_of_grid_Y));
  cosTable.resize(number_of_grid_X * 3 / 2, 0);
  workArea_.resize(round(sqrt(number_of_grid_X)) + 2, 0);

  float **binDensity = new float *[number_of_grid_X];
  float **electricPotential = new float *[number_of_grid_X];
  float **electricForceX = new float *[number_of_grid_X];
  float **electricForceY = new float *[number_of_grid_X];

  for (int i = 0; i < number_of_grid_X; i++) {
    binDensity[i] = new float[number_of_grid_Y];
    electricPotential[i] = new float[number_of_grid_Y];
    electricForceX[i] = new float[number_of_grid_Y];
    electricForceY[i] = new float[number_of_grid_Y];

    for (int j = 0; j < number_of_grid_Y; j++) {
      binDensity[i][j] = bins2D[i][j].density / normal_bin_width / normal_bin_height /number_of_grid_X / number_of_grid_Y;
      electricPotential[i][j] = bins2D[i][j].electricPotential;
      electricForceX[i][j] = bins2D[i][j].electricField_x;
      electricForceY[i][j] = bins2D[i][j].electricField_y;
    }
  }

//  ddct2d(number_of_grid_X,number_of_grid_Y, -1, binDensity, NULL, workArea_.data(), cosTable.data());
  ddct2d(number_of_grid_X, number_of_grid_Y, -1, binDensity, NULL, (int *) workArea_.data(), (float *) cosTable.data());

//  for (int i = 0; i < number_of_grid_X; ++i) {
//    binDensity[i][0] *= 0.5;
//  }
//  for (int i = 0; i < number_of_grid_Y; ++i) {
//    binDensity[0][i] *= 0.5;
//  }
  // for (int i = 0; i < number_of_grid_X; ++i) {
  //   binDensity[i][number_of_grid_Y-1] *= 2;
  // }
  // for (int i = 0; i < number_of_grid_Y; ++i) {
  //   binDensity[number_of_grid_X-1][i] *= 2;
  // }

  // for (int i = 0; i < number_of_grid_X; ++i) {
  //   for (int j = 0; j < number_of_grid_Y; ++j) {
  //     binDensity[i][j] *= 4.0 /number_of_grid_X / number_of_grid_Y;
  //   }
  // }
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

  ddct2d(number_of_grid_X,
         number_of_grid_Y,
         1,
         electricPotential,
         NULL,
         (int *) workArea_.data(),
         (float *) cosTable.data());
  ddsct2d(number_of_grid_X,
          number_of_grid_Y,
          1,
          electricForceX,
          NULL,
          (int *) workArea_.data(),
          (float *) cosTable.data());
  ddcst2d(number_of_grid_X,
          number_of_grid_Y,
          1,
          electricForceY,
          NULL,
          (int *) workArea_.data(),
          (float *) cosTable.data());

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

  for (int x = 0; x < number_of_grid_X; x++) {
    for (int y = 0; y < number_of_grid_Y; y++) {
      for (auto &pair_inst : bins2D[x][y].overlappedInstance) {
        Instance *inst = pair_inst.first;
        float ratio = pair_inst.second / inst->getArea();
        inst->densityForce.first += bins2D[x][y].electricField_x * ratio;
        inst->densityForce.second += bins2D[x][y].electricField_y * ratio;
      }
    }
  }

  for (auto &inst : instance_pointers_) {
    int instNum = instMap.find(inst->getName())->second;
    if(inst->isFiller) {
      double gradDenX = lambda * ((double) inst->getArea() * inst->densityForce.first);
      double gradDenY = lambda * ((double) inst->getArea() * inst->densityForce.second);

      gradX[instNum] = - gradDenX;
      gradY[instNum] = - gradDenY;
      continue;
    }
    float gradInstX = 0.0, gradInstY = 0.0;
    for (auto &pin : inst->getPins()) {
      gradInstX += pin->gradWAX;
      gradInstY += pin->gradWAY;
    }
    double gradDenX = (double) inst->getArea() * inst->densityForce.first;
    double gradDenY =  (double) inst->getArea() * inst->densityForce.second;
    gradX[instNum] = - gradInstX - lambda * gradDenX;
    gradY[instNum] = - gradInstY - lambda * gradDenY;
  }
//  cout << absSumX << " : " << absSumY << endl;
//  cout<<"END"<<endl;
}

bool Circuit::densityCheck(float normal_bin_width, float normal_bin_height) {
  int die_width = die_->getWidth();
  int die_height = die_->getHeight();

  for (int i = 0; i <= number_of_grid_X; ++i) {
    for (int j = 0; j <= number_of_grid_Y; ++j) {
      bins2D[i][j].reset();
    }
  }
  // get utility in each bins
  for (Instance *instance : instance_pointers_) {
    pair<int, int> instance_lower_left = instance->getCoordinate();
    pair<int, int> instance_upper_right{
        instance_lower_left.first + instance->getWidth(),
        instance_lower_left.second + instance->getHeight()
    };
    instance->densityForce = make_pair(0.0, 0.0);

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

  // Calc a
  bool results = false;
  float worstDensity = 0.0;
  float sumDensity = 0.0;
  for (int i = 0; i < number_of_grid_X; i++) {
    for (int j = 0; j < number_of_grid_Y; j++) {
      float density = bins2D[i][j].stdArea / ((float) normal_bin_width * normal_bin_height);
      if (density >= 1.2) {
        if (worstDensity < density) worstDensity = density;
        sumDensity += density;
        results = true;
      }
    }
  }
  cout << " WORST " << worstDensity << " Neg sum : " << sumDensity << endl;
  return results;
}

void Circuit::myPlacement() {
  // Do 
  clock_t start, end;
  start = clock();
  uint die_width = die_->getWidth();
  uint die_height = die_->getHeight();
  float normal_bin_width = (float) die_width / number_of_grid_X;
  float normal_bin_height = (float) die_height / number_of_grid_Y;
  float dieArea = (float) die_width * die_height;

  // Give ID to instances
  int instCnt = 0;
  int fillerCnt = 0;
  int type = 0;
//  cout << "Start"<<endl;
  for (auto &inst : instance_pointers_) {
    // cout << inst->name_ <<endl;
    instMap.insert(make_pair(inst->getName(), instCnt++));
  }
//  cout << "giveID";
  for (auto &inst : instance_pointers_) {
    if (inst->isFiller) continue;
    inst->setCoordinate(int(die_width / 2 - inst->getWidth() / 2), int(die_height / 2 - inst->getHeight() / 2));
  }

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

  if(instance_pointers_.size() <= 47512) {
    type = 2;
  }
  else if(instance_pointers_.size() <= 142136) {
    type = 1;
  }
  else if(instance_pointers_.size() <= 192911) {
    type = 4;
  }
  else if(instance_pointers_.size() <= 292190) {
    type = 3;
  }
  else {
    type = 5;
  }
  for (auto &inst : instance_pointers_) {
    int instNum = instMap.find(inst->getName())->second;
    pair<int, int> coordinate = inst->getCoordinate();
    vX[instNum] = coordinate.first;
    uX[instNum] = coordinate.first;
    vY[instNum] = coordinate.second;
    uY[instNum] = coordinate.second;
    vX[instNum] = coordinate.first - 10;
    uX[instNum] = coordinate.first - 10;
    vY[instNum] = coordinate.second - 10;
    uY[instNum] = coordinate.second - 10;
  }
  long long prevHPWL = 0, HPWL = 0;
  cout << "start"<<endl;
  int ForceDirectCNT = 20;
  for (int i = 0; i < ForceDirectCNT; ++i) {
    vector<float> vX_hat(vX), vY_hat(vY);
    initialPlacement(gradX, gradY);
    for (auto &inst : instance_pointers_) {
      if(inst->isFiller) continue;
      int instNum = instMap.find(inst->getName())->second;
      vX_hat[instNum] = clamp(vX[instNum] - 1e5 * gradX[instNum], (number_of_grid_X/2 - 3) * normal_bin_width, (number_of_grid_X/2 + 2) * normal_bin_width);
      vY_hat[instNum] = clamp(vY[instNum] - 1e5 * gradY[instNum], (number_of_grid_Y/2 - 3) * normal_bin_height, (number_of_grid_Y/2 + 3) * normal_bin_height);
    }
    
    placeMap(vX_hat, vY_hat);
    prevHPWL = 0;
    for (Net *net : net_pointers_) {
      prevHPWL += (long long) net->getHPWL();
    }
    cout << "INIT iter " << i << " HPWL : " << prevHPWL <<endl;
    
    for (int i = 0; i < vX.size(); i++) {
      vX[i] = vX_hat[i];
      vY[i] = vY_hat[i];
    }
  }
  string img_file_name = "initial";
    saveImg(img_file_name);

  // Iterate until
  bool condition = true;

  int iter = 0;
  
  float alpha_max = 0.00001;
  float lambda_0 = initLambda()*1e-2;
  float prev_a = 1.0;
  float prev_alpha = 0.00001 * 0.5;
  if(instance_pointers_.size() <= 142136) {
    alpha_max = 0.0001;
    prev_alpha = 0.0001 * 0.5;
  }
  cout << "Init Lambda " << lambda_0 << endl;

  float prev_lambda = lambda_0;
  float mew_0 = 1.1;

//  long long prevHPWL = 0, HPWL = 0;
  prevHPWL = 0;
  for (Net *net : net_pointers_) {
    prevHPWL += (long long) net->getHPWL();
  }
  HPWL = prevHPWL;

  while (condition) {
    // Update coeff
    float mew = 1.1;
    float diff_ = 1.0 - (float) (HPWL - prevHPWL) / 3.5e5;
    // cout << HPWL << " " << prevHPWL<< " " << diff_<< endl;
    if (diff_ <= -1.5) mew = 0.99;
    else if (diff_ >= 1) mew = 1.01;
    else {
      mew = clamp(pow(mew_0, diff_), 0.99, 1.01);
    }
    float lambda = maxClamp(mew * prev_lambda, 1e1);

    prevHPWL = HPWL;
    prev_lambda = lambda;

    calcGradient(gradX, gradY, lambda);

    // Initial grad
    if (iter == 0) {
      for (auto &inst : instance_pointers_) {
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
    float alpha_k_max = max(1.0*alpha_max, 1.1 * prev_alpha);
    alpha_k_max = maxClamp(alpha_k_max, 10.0*alpha_max);
    float next_ak;

    vector<float> uX_hat(instance_pointers_.size()), uY_hat(instance_pointers_.size());
    vector<float> vX_hat(instance_pointers_.size()), vY_hat(instance_pointers_.size());

    for (auto &inst : instance_pointers_) {
      int instNum = instMap.find(inst->getName())->second;
      if(inst->isFiller && iter < 3) {
        vX_hat[instNum] = vX[instNum] + (-alpha_hat) * 5 * gradX[instNum];
        vY_hat[instNum] = vY[instNum] + (-alpha_hat) * 5 * gradY[instNum];
      }
      else {
        vX_hat[instNum] = vX[instNum] + (-alpha_hat) * gradX[instNum];
        vY_hat[instNum] = vY[instNum] + (-alpha_hat) * gradY[instNum];
      }
    }
    placeMap(vX_hat, vY_hat);

    vector<float> next_gradX(instance_pointers_.size());
    vector<float> next_gradY(instance_pointers_.size());
    calcGradient(next_gradX, next_gradY, lambda);
    next_ak = calcAlphaHat(vX_hat, vY_hat, vX, vY, next_gradX, next_gradY, gradX, gradY);
    int wi = 0;
    while (alpha_hat > epsilon * next_ak && wi < 10) {
      wi++;
      alpha_hat = next_ak;

      for (auto &inst : instance_pointers_) {
        int instNum = instMap.find(inst->getName())->second;
        if(inst->isFiller && iter < 10) {
          vX_hat[instNum] = vX[instNum] + (-alpha_hat) * 5 * gradX[instNum];
          vY_hat[instNum] = vY[instNum] + (-alpha_hat) * 5 * gradY[instNum];
        }
        else {
          vX_hat[instNum] = vX[instNum] + (-alpha_hat) * gradX[instNum];
          vY_hat[instNum] = vY[instNum] + (-alpha_hat) * gradY[instNum];
        }
      }
      placeMap(vX_hat, vY_hat);

      vector<float> next_gradX(instance_pointers_.size());
      vector<float> next_gradY(instance_pointers_.size());
      calcGradient(next_gradX, next_gradY, lambda);

      next_ak = calcAlphaHat(vX_hat, vY_hat, vX, vY, next_gradX, next_gradY, gradX, gradY);
      if (next_ak < 0.01 * alpha_k_max || next_ak > alpha_k_max) {
        alpha_hat = next_ak;
        break;
      }
      // cout << "WHILE " << alpha_hat << " : " << epsilon * next_ak << endl;
    }
    if(isnan(alpha_hat) != 0) alpha_hat = alpha_k_max;

    alpha_hat = clamp(alpha_hat, 0.01f * alpha_k_max, alpha_k_max);
    cout<<"STEP " << alpha_hat<< "Lambda " << lambda << " : ";
    // Nestrov
    for (int i = 0; i < uX.size(); i++) {
      
    }
    for (auto &inst : instance_pointers_) {
      int instNum = instMap.find(inst->getName())->second;
      if(inst->isFiller && iter < 10) {
        uX_hat[instNum] = vX[instNum] + (-alpha_hat) * 5 * gradX[instNum];
        uY_hat[instNum] = vY[instNum] + (-alpha_hat) * 5 * gradY[instNum];
      }
      else {
        uX_hat[instNum] = vX[instNum] + (-alpha_hat) * gradX[instNum];
        uY_hat[instNum] = vY[instNum] + (-alpha_hat) * gradY[instNum];
      }
    }
    float next_a = (1 + sqrt(4 * prev_a * prev_a + 1)) / 2.0;
    for (int i = 0; i < vX.size(); i++) {
      vX_hat[i] = uX_hat[i] + (prev_a - 1) * (uX_hat[i] - uX[i]) / next_a;
      vY_hat[i] = uY_hat[i] + (prev_a - 1) * (uY_hat[i] - uY[i]) / next_a;
    }
    prev_alpha = alpha_hat;
    prev_a = next_a;
    for (int i = 0; i < vX.size(); i++) {
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
    float tresult = (float) (end - start) / CLOCKS_PER_SEC;
    tresult = (float) (end - start) / CLOCKS_PER_SEC;
    // Update coeff
    HPWL = 0;
    for (Net *net : net_pointers_) {
      HPWL += (long long) net->getHPWL();
    }

    condition = densityCheck(normal_bin_width, normal_bin_height);

    if (iter % 5 == 0) {
      cout << "iter " << iter << " HPWL : " << HPWL << "\tTIME : " << tresult << endl;
      string img_file_name = "result" + to_string(iter);
      saveImg(img_file_name);
    }
    if(instance_pointers_.size() <= 47512) {
      if(iter == 60) break;
    }
    else if(instance_pointers_.size() <= 142136) {
      if(iter == 85) break;
    }
    else if(instance_pointers_.size() <= 192911) {
      if(iter == 10) break;
    }
    else if(instance_pointers_.size() <= 292190) {
      if(iter == 130) break;
    }
    else {
      if(iter == 10) break;      
    }
    iter++;
  }

  //After
  int afterForceDirectCNT = 3;
  for (int i = 0; i < afterForceDirectCNT; ++i) {
    vector<float> vX_hat(vX), vY_hat(vY);
    initialPlacement(gradX, gradY);
    for (auto &inst : instance_pointers_) {
      if(inst->isFiller) continue;
      int instNum = instMap.find(inst->getName())->second;
      vX_hat[instNum] = clamp(vX[instNum] - 1e5 * gradX[instNum], (number_of_grid_X/2 - 3) * normal_bin_width, (number_of_grid_X/2 + 2) * normal_bin_width);
      vY_hat[instNum] = clamp(vY[instNum] - 1e5 * gradY[instNum], (number_of_grid_Y/2 - 3) * normal_bin_height, (number_of_grid_Y/2 + 3) * normal_bin_height);
    }
    
    placeMap(vX_hat, vY_hat);
    prevHPWL = 0;
    for (Net *net : net_pointers_) {
      prevHPWL += (long long) net->getHPWL();
    }
    cout << "After iter " << i << " HPWL : " << prevHPWL <<endl;
    
    for (int i = 0; i < vX.size(); i++) {
      vX[i] = vX_hat[i];
      vY[i] = vY_hat[i];
    }
  }

  vector<vector<Bin> > bins2D_;

  int number_of_grid_X_ = 40;
  int number_of_grid_Y_ = 40;
  normal_bin_width = (float)die_width / number_of_grid_X_;
  normal_bin_height = (float)die_height / number_of_grid_Y_;
  for (int i = 0; i <= number_of_grid_X_; ++i) {
    vector<Bin> bins1D;
    for (int j = 0; j <= number_of_grid_Y_; ++j) {
      pair<int, int> lower_left{i * normal_bin_width, j * normal_bin_height};
      pair<int, int> upper_right{(i + 1) * normal_bin_width, (j + 1) * normal_bin_height};
      bins1D.emplace_back(floor(normal_bin_width) * floor(normal_bin_height), lower_left, upper_right);
    }
    bins2D_.push_back(bins1D);
  }
  for (int i = 0; i <= number_of_grid_X_; ++i) {
    for (int j = 0; j <= number_of_grid_Y_; ++j) {
      bins2D_[i][j].reset();
    }
  }
  // get utility in each bins
  vector<Instance *> tempInst;
  priority_queue<Instance *, vector<Instance *>, cmp> pq_overlappedInstance;

  for (Instance *instance : instance_pointers_) {
    if(instance->isFiller) {
      instance->setCoordinate(0,0);
      instance->fillerWidth = 0.0;
      instance->fillerHeight = 0.0;
      continue;
    }
    pq_overlappedInstance.push(instance);
    tempInst.push_back(instance);
    // data_mapping_.instances
  }

  int w = floor(sqrt(total_cell_area));
  int startX = die_width/2 - w/2, startY = die_height/2 - w/2, dx = 0, dy = 0;

  while(!pq_overlappedInstance.empty()) {
    Instance *inst = pq_overlappedInstance.top();
    pq_overlappedInstance.pop();
    inst->setCoordinate(startX, startY);
    startX += 0.90 * inst->getWidth();
    if(startX >= die_width/2 + w/2) {
      startX = die_width/2 - w/2;
      startY += inst->getHeight();
    }
  }

  cout<<instance_pointers_.size()<<endl;
  instance_pointers_ = tempInst;
  cout<<instance_pointers_.size()<<endl;

  img_file_name = "Final_result";
  saveImg(img_file_name);
}
}
