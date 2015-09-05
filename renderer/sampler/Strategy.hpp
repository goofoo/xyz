#ifndef __TGIR_SAMPLER_STRATEGY_HPP__
#define __TGIR_SAMPLER_STRATEGY_HPP__

namespace tgir {
struct Ptrr : public hi::thread {
  void run();
};

struct Ptid : public hi::thread {
  void run();
};

/*
  struct Ptgw : public hi::thread
  {
    void run();
  };
*/
struct Ptpf1 : public hi::thread {
  void run();
};

struct Ptpf2 : public hi::thread {
  void run();
};

}  // end of namespace tgir

#endif __TGIR_SAMPLER_STRATEGY_HPP__
