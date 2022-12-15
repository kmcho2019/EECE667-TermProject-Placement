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
#include "Circuit.h"

using namespace std;

namespace Placer {

void Circuit::parse(const string &lef_name, const string &def_name) {
  parser_.readLef(lef_name);
  parser_.readDef(def_name);
  this->init();
}
void Circuit::init() {
  // You don't need to understand/see this code.
  // by TA
  dbBlock *block = parser_.db_database_->getChip()->getBlock();
  dbSet<dbInst> db_instances = block->getInsts();
  dbSet<dbNet> db_nets = block->getNets();
  auto lib = parser_.db_database_->findLib("nangate");

  /// Die setting
  data_storage_.die.setDbBlock(block);
  die_ = &data_storage_.die;
  total_cell_area = 0;
  /*!
   * @brief
   * Instance setting
   *
   * @details
   * 1. It makes real instance data and store in \c data_storage.instances. \n
   * 2. Then it makes pointer set for \c Circuit class, \n
   * 3. And also makes mapping from \c db_instance to instance pointer. \n\n
   * */
  // Filler insertion
  for (auto *db_inst: db_instances) {
    Instance instance(db_inst);
    total_cell_area += instance.getArea();
  }
  totalAreaFC = 1.2 * (die_->getArea() - total_cell_area) - total_cell_area;
  areaFC = total_cell_area / (float)db_instances.size();
  cntFC = floor(totalAreaFC/areaFC);

  data_storage_.instances.reserve(db_instances.size() + cntFC);  // real data for instance
  instance_pointers_.reserve(db_instances.size() + cntFC);  // pointer data for instances
  data_storage_.nets.reserve(db_nets.size());
  net_pointers_.reserve(db_nets.size());

  // 1. make real data for instances
  for (odb::dbInst *db_inst : db_instances) {
    Instance instance(db_inst);
    instance.setDataStorage(&data_storage_);
    instance.setDataMapping(&data_mapping_);
    if (instance.isPlaced()){
      auto coordinate = instance.getCoordinate();
      instance.setCoordinate(coordinate.first, coordinate.second);
    }
    data_storage_.instances.push_back(instance);
  }

  mt19937 genX(random_device{}());
  mt19937 genY(random_device{}());
  uniform_real_distribution<float> disX(0, die_->getWidth());
  uniform_real_distribution<float> disY(0, die_->getHeight());
  // for(auto n : lib->getMasters()) {
  //   cout << n->getName() <<endl;
  // }
  for(int i = 0;i<cntFC; i++) {
    Instance filler;
    string name =  "_filler" + to_string(i);
    // cout<<filler.name_<<endl;
    filler.isFiller = true;
    filler.name_ = name.c_str();
    filler.fillerWidth = sqrt(areaFC);
    filler.fillerHeight = sqrt(areaFC);
    filler.setCoordinate(floor(disX(genX)), floor(disY(genY)));
    // filler.setDataStorage(&data_storage_);
    // filler.setDataMapping(&data_mapping_);
    data_storage_.instances.push_back(filler);
  }

  // 2-3. make pointer set and map from db_instance to instance pointer
  for (auto &instance : data_storage_.instances) {
    instance_pointers_.push_back(&instance);
    if(!instance.isFiller) {
      data_mapping_.inst_map[instance.getDbInst()] = &instance;
    }
  }


  /*!
   * @brief
   * Pin setting
   *
   * @details
   * Same with above way
   */
  // 1. make real data
  // 1-1. Instance terminals
  for (auto instance : instance_pointers_) {
    if(instance->isFiller) continue;
    for (dbITerm *db_i_term : instance->getDbInst()->getITerms()) {
      Pin pin(db_i_term);
      pin.setDataStorage(&data_storage_);
      pin.setDataMapping(&data_mapping_);
      data_storage_.pins.push_back(pin);
    }
  }
  // 1-2. Block terminals
  for(dbBTerm* db_b_term: block->getBTerms()){
    Pin pin(db_b_term);
    pin.setDataStorage(&data_storage_);
    pin.setDataMapping(&data_mapping_);
    data_storage_.pins.push_back(pin);
  }

  // 2. make pointer set and map from db_pin to pin pointer
  for (auto & pin : data_storage_.pins) {
    Pin *pin_pointer = &pin;
    pin_pointers_.push_back(pin_pointer);
    if (pin_pointer->isInstancePin()) {
      data_mapping_.pin_map_i[pin_pointer->getDbITerm()] = pin_pointer;
    } else if (pin_pointer->isBlockPin()) {
      data_mapping_.pin_map_b[pin_pointer->getDbBTerm()] = pin_pointer;
      pad_pointers_.push_back(pin_pointer);
    }
  }

  /*!
   * @brief
   * Net setting
   *
   * @details
   * Same with above way
   */
  // 1. make real data
  for (odb::dbNet *db_net : db_nets) {
    Net net(db_net);
    net.setDataStorage(&data_storage_);
    net.setDataMapping(&data_mapping_);
    data_storage_.nets.push_back(net);
  }
  // 2. make pointer set and map from db_net to net pointer
  for (auto &net : data_storage_.nets) {
    net_pointers_.push_back(&net);
    data_mapping_.net_map[net.getDbNet()] = &net;
  }

//  for(auto &inst : instance_pointers_) {
//    cout<<inst->name_<<endl;
//  }
}

void Circuit::write(const string &out_file_name) {
  parser_.writeDef(out_file_name);
}
int Circuit::getUnitOfMicro() const {
  return parser_.db_database_->getTech()->getDbUnitsPerMicron();
}
ulong Circuit::getHPWL() {
  ulong HPWL = 0;
  for (Net *net : net_pointers_) {
    HPWL += net->getHPWL();
  }
  return HPWL;
}
void Circuit::analyzeBench() {
  cout << "================== bench analyze ==================" << endl;
  cout << "Instance #: " << instance_pointers_.size() << endl;
  cout << "Net #: " << net_pointers_.size() << endl;
  cout << "IO pad #: " << pad_pointers_.size() << endl;
  uint64 die_area = die_->getArea();

  cout << scientific << endl;
  cout << "Total Cell Area: " << static_cast<float>(total_cell_area) << endl;
  cout << "Die Area: " << static_cast<float>(die_area) << endl;
  cout << "===================================================" << endl;
}

void Circuit::placeMap(vector<float> &vX, vector<float> &vY) {
  uint die_width = die_->getWidth();
  uint die_height = die_->getHeight();
  for (auto &inst : instance_pointers_) {
    int instID = instMap.find(inst->getName())->second;
    int newX = clamp(int(vX[instID]), 0, (int) die_width - (int) inst->getWidth());
    int newY = clamp(int(vY[instID]), 0, (int) die_height - (int) inst->getHeight());
    inst->setCoordinate(newX, newY);
  }
}

} // Placer