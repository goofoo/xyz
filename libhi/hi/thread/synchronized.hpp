#pragma once

#include <hi/thread/monitor_object.hpp>

namespace hi {
class synchronized sealed {
 public:
  inline synchronized(hi::monitor_object &monitor) : monitor_(monitor) {
    monitor_.enter();
  }

  inline ~synchronized() { monitor_.leave(); }

 private:
  hi::monitor_object &monitor_;

  // ñ¢íËã`ÇÃä÷êî
  synchronized(synchronized const &);
  synchronized &operator=(synchronized const &);
};

}  // end of namespace hi
