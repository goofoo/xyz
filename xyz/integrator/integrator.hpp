#ifndef XYZ_PATHTRACER_HPP_
#define XYZ_PATHTRACER_HPP_
namespace xyz {
struct PathTracer : public hi::thread {
  virtual void run();
};

struct ImportanceDrivenPathTracer : public hi::thread {
  virtual void run();
};

struct PathTracer_MetropolisLightTransport : public hi::thread {
  virtual void run();
};

struct ImportanceDrivenPathTracer_MetropolisLightTransport : public hi::thread {
  virtual void run();
};

struct PathTracerWithGoWithTheWinnersStrategy : public hi::thread {
  virtual void run();
};

struct ImportanceDrivenPathTracerWithGoWithTheWinnersStrategy
    : public hi::thread {
  virtual void run();
};

struct PathTracerWithGoWithTheWinnersStrategy_MetropolisLightTransport
    : public hi::thread {
  virtual void run();
};

struct
    ImportanceDrivenPathTracerWithGoWithTheWinnersStrategy_MetropolisLightTransport
    : public hi::thread {
  virtual void run();
};
}

namespace xyz {
struct BidirectionalPathTracer : public hi::thread {
  virtual void run();
};

struct ImportanceDrivenBidirectionalPathTracer : public hi::thread {
  virtual void run();
};

struct BidirectionalPathTracer_MetropolisLightTransport : public hi::thread {
  virtual void run();
};

struct ImportanceDrivenBidirectionalPathTracer_MetropolisLightTransport
    : public hi::thread {
  virtual void run();
};
}

namespace xyz {
struct ResampledPathTracer : public hi::thread {
  virtual void run();
};
}

namespace xyz {
struct MetropolisLightTransportTestbed : public hi::thread {
  virtual void run();
};

struct PathTracer_PopulationLightTransportTestbed : public hi::thread {
  virtual void run();
};

struct BidirectionalPathTracer_PopulationLightTransportTestbed
    : public hi::thread {
  virtual void run();
};
}

#endif
