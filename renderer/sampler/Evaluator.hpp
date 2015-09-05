#ifndef __TGIR_SAMPLER_EVALUATOR_HPP__
#define __TGIR_SAMPLER_EVALUATOR_HPP__

namespace tgir {
/// <summary>
/// Bidirectional Path Tracing
/// </summary>
struct BPT : public hi::thread {
  void run();
};

/// <summary>
/// Large-Step Metropolis Light Transport
/// </summary>
struct SimplifiedMLT : public hi::thread {
  void run();
};

/// <summary>
/// Original Metropolis Light Transport
/// </summary>
struct OriginalMLT : public hi::thread {
  void run();
};

/// <summary>
/// Full Metropolis Light Transport
/// </summary>
struct FullMLT : public hi::thread {
  void run();
};

/// <summary>
/// Replica-Exchange Light Transport
/// </summary>
struct RELT : public hi::thread {
  void run();
};

/// <summary>
/// Replica-Exchange Light Transport
/// </summary>
struct FullRELT : public hi::thread {
  void run();
};

}  // end of namespace tgir

#endif __TGIR_SAMPLER_EVALUATOR_HPP__
