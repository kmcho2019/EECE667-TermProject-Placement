///////////////////////////////////////////////////////////////////////////////
// Creator: Minjae Kim of CSDL, POSTECH
// Email:   kmj0824@postech.ac.kr
// GitHub:  ApeachM
//
// BSD 3-Clause License
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice, this
//   list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
//
// * Neither the name of the copyright holder nor the names of its
//   contributors may be used to endorse or promote products derived from
//   this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
///////////////////////////////////////////////////////////////////////////////

#ifndef PLACER_INCLUDE_DATASTRUCTURES_CIRCUIT_H_
#define PLACER_INCLUDE_DATASTRUCTURES_CIRCUIT_H_
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <random>
#include <queue>
#include "Parser.h"
#include "Instance.h"
#include "Net.h"
#include "Pin.h"
#include "Die.h"

namespace Placer {
using namespace odb;

struct cmp {
bool operator()(Instance* &a, Instance* &b) {
  pair<int, int> aC = a->getCoordinate();
  pair<int, int> bC = b->getCoordinate();
  if (aC.second == bC.second) {
      return aC.first < bC.first;
  }
  else {
      return aC.second < bC.second;
  }
}
};

class Bin {
 public:
  Bin() = default;
  explicit Bin(int area, pair<int, int> lower_left, pair<int, int> upper_right)
      : bin_area_(area), lower_left_(std::move(lower_left)), upper_right_(std::move(upper_right)) {}
  pair<int, int> lower_left_;
  pair<int, int> upper_right_;
  float bin_area_{0};
  float stdArea{0};
  float fillerArea{0};
  float density = 0.0;
  float electricPotential = 0.0;
  float electricField_x = 0.0;
  float electricField_y = 0.0;
  vector<pair<Instance *, float> > overlappedInstance;

//  void getOverlapArea(Instance *instance) {
//    if (bin_area_ == 0) {
//      assert(0);
//    } else {
//      pair<int, int> instance_lower_left = instance->getCoordinate();
//      pair<int, int> instance_upper_right
//          {instance_lower_left.first + instance->getWidth(), instance_lower_left.second + instance->getHeight()};
//      if (instance_upper_right.first <= lower_left_.first) {
//        return;
//      } else if (instance_upper_right.second <= lower_left_.second) {
//        return;
//      } else if (instance_lower_left.first >= upper_right_.first) {
//        return;
//      } else if (instance_lower_left.second >= upper_right_.second) {
//        return;
//      } else {
//        float dx = (float)instance_upper_right.first - (float)lower_left_.first;
//        float dy = (float)instance_upper_right.second - (float)lower_left_.second;
//        if(instance->isFiller) fillerArea += dx * dy;
//        else stdArea += dx * dy;
//      }
//    }
//  }

  void reset() {
    stdArea = 0;
    fillerArea = 0;
    density = 0.0;
    electricPotential = 0.0;
    electricField_x = 0.0;
    electricField_y = 0.0;
    overlappedInstance.clear();    
  }

  void getOverlapArea(Instance *instance) {    
    if (bin_area_ == 0) {
      assert(0);
    } else {
      int bin_width = upper_right_.first - lower_left_.first;
      
      pair<int, int> instance_lower_left = instance->getCoordinate();
      pair<int, int> instance_upper_right;
      float densityX = 1.0;

      if (instance->getWidth() >= bin_width) {
        instance_upper_right = make_pair(instance_lower_left.first + instance->getWidth(),
                                         instance_lower_left.second + instance->getHeight());
      } else {
        instance_upper_right =
            make_pair(instance_lower_left.first + bin_width, instance_lower_left.second + instance->getHeight());
        density = (float) instance->getWidth() / bin_width;
      }
      int rectLx = max(this->lower_left_.first, instance_lower_left.first);
      int rectLy = max(this->lower_left_.second, instance_lower_left.second);
      int rectUx = min(this->upper_right_.first, instance_upper_right.first);
      int rectUy = min(this->upper_right_.second, instance_upper_right.second);

      int overlapWidth = rectUx - rectLx;
      int overlapHeight = rectUy - rectLy;

      if (overlapWidth < 0 || overlapHeight < 0) {
        return;
      }
      float cell_area_ = (float)(rectUx - rectLx) * (rectUy - rectLy);
      
      float overlapArea = cell_area_ * density;
      overlappedInstance.push_back(make_pair(instance, overlapArea));
      if (instance->isFiller) fillerArea += overlapArea;
      else {
        stdArea += overlapArea;
      }
    }
  }
};

class Circuit {
 protected:
  Parser parser_;
  data_storage data_storage_;
  data_mapping data_mapping_;

  std::vector<Instance *> instance_pointers_;
  std::vector<Net *> net_pointers_;
  std::vector<Pin *> pin_pointers_;  // This vector includes instance pin pointers and pad pin pointers
  std::vector<Pin *> pad_pointers_;
  Die *die_ = nullptr;
  void init();

 public:
  Circuit() = default;
  ~Circuit() = default;
  void parse(const string &lef_name, const string &def_name);
  void write(const string &out_file_name);
  void quadraticPlacement();
  void myPlacement();
  void calcGradient(std::vector<float> &gradX, std::vector<float> &gradY);


  /// \brief
  /// get unit of micro
  /// \details
  /// the coordinate in this circuit is `return value`/1um.
  /// \example
  /// if the return value is 100, then
  /// (20000, 30000) means coordinate (200um, 300um)
  int getUnitOfMicro() const;

  /// \brief
  /// saveImg the circuit image.
  /// \details
  /// It saves the picture for the cells, pads, and nets in the circuit,
  /// in the output/images/file_name.png
  void saveImg(const string &output_path_name, const string &file_name);

  /// \brief
  /// return the HPWL of the total circuit
  ulong getHPWL();

  /// \brief
  /// Analyze the bench metrics and print them
  void analyzeBench();

  // etc
  int cntFC = 0;
  float total_cell_area = 0.0f;
  float areaFC = 0.0f, totalAreaFC = 0.0f;
  unordered_map<string, int> instMap;
  vector<float> wx;
  vector<float> wx_sq;
  vector<float> wy;
  vector<float> wy_sq;
  vector<float> cosTable;
  vector<int> workArea_;
  vector<vector<Bin> > bins2D;

  void howToUse();
  void placeExample();
  void dbTutorial() const;
  void initialPlacement(vector<float> &gradX, vector<float> &gradY);
  void calcGradient(vector<float> &gradX, vector<float> &gradY, float lambda);
  void placeMap(vector<float> &vX, vector<float> &vY);
  bool densityCheck(float normal_bin_width, float normal_bin_height);
  float initLambda();

  // FOR FFT
  void cdft(int n, int isgn, float *a, int *ip, float *w);
  void rdft(int n, int isgn, float *a, int *ip, float *w);
  void ddct(int n, int isgn, float *a, int *ip, float *w);
  void ddst(int n, int isgn, float *a, int *ip, float *w);
  void dfct(int n, float *a, float *t, int *ip, float *w);
  void dfst(int n, float *a, float *t, int *ip, float *w);

  void makewt(int nw, int *ip, float *w);
  void makeipt(int nw, int *ip);
  void makect(int nc, int *ip, float *c);

  void cftfsub(int n, float *a, int *ip, int nw, float *w);
  void cftbsub(int n, float *a, int *ip, int nw, float *w);
  void bitrv2(int n, int *ip, float *a);
  void bitrv2conj(int n, int *ip, float *a);
  void bitrv216(float *a);
  void bitrv216neg(float *a);
  void bitrv208(float *a);
  void bitrv208neg(float *a);
  void cftf1st(int n, float *a, float *w);
  void cftb1st(int n, float *a, float *w);

  void cftrec4(int n, float *a, int nw, float *w);
  int cfttree(int n, int j, int k, float *a, int nw, float *w);
  void cftleaf(int n, int isplt, float *a, int nw, float *w);
  void cftmdl1(int n, float *a, float *w);
  void cftmdl2(int n, float *a, float *w);
  void cftfx41(int n, float *a, int nw, float *w);
  void cftf161(float *a, float *w);
  void cftf162(float *a, float *w);
  void cftf081(float *a, float *w);
  void cftf082(float *a, float *w);
  void cftf040(float *a);
  void cftb040(float *a);
  void cftx020(float *a);
  void rftfsub(int n, float *a, int nc, float *c);
  void rftbsub(int n, float *a, int nc, float *c);
  void dctsub(int n, float *a, int nc, float *c);
  void dstsub(int n, float *a, int nc, float *c);

  // 2D fftsg
  void cdft2d(int n1, int n2, int isgn, float **a, float *t, int *ip, float *w);
  void rdft2d(int n1, int n2, int isgn, float **a, float *t, int *ip, float *w);
  void rdft2dsort(int n1, int n2, int isgn, float **a);
  void ddcst2d(int n1, int n2, int isgn, float **a, float *t, int *ip, float *w);
  void ddsct2d(int n1, int n2, int isgn, float **a, float *t, int *ip, float *w);
  void ddct2d(int n1, int n2, int isgn, float **a, float *t, int *ip, float *w);
  void ddst2d(int n1, int n2, int isgn, float **a, float *t, int *ip, float *w);
  void cdft2d_sub(int n1, int n2, int isgn, float **a, float *t, int *ip, float *w);
  void rdft2d_sub(int n1, int isgn, float **a);
  void ddxt2d_sub(int n1, int n2, int ics, int isgn, float **a, float *t, int *ip, float *w);
};

} // Placer

#endif //PLACER_INCLUDE_DATASTRUCTURES_CIRCUIT_H_
