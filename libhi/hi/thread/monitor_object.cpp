#include <hi/thread/monitor_object.hpp>
#include <hi/thread/thread.hpp>

namespace hi {
monitor_object::monitor_object() : locking_thread_(nullptr), locking_count_(0) {
  ::InitializeCriticalSection(&cs_);
}

monitor_object::~monitor_object() { ::DeleteCriticalSection(&cs_); }

void monitor_object::enter() {
  ::EnterCriticalSection(&cs_);
  locking_thread_ = ::GetCurrentThread();
  locking_count_++;
}

void monitor_object::leave() {
  if (0 == --locking_count_) {
    locking_thread_ = nullptr;
  }
  ::LeaveCriticalSection(&cs_);
}

bool monitor_object::has_lock() const {
  return ::GetCurrentThread() == locking_thread_;
}

}  // end of namespace hi
