#pragma once

namespace hi {
class runnable abstract  // interface
    {
 public:
  virtual ~runnable();
  virtual void run() = 0;
};

}  // end of namespace runnable
