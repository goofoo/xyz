#pragma once

#include <hi/lang.hpp>

namespace hi {
class thread;
class synchronized;

/// <summary>
/// ロックできるモニタ
/// </summary>
class monitor_object {
  friend class hi::synchronized;

 public:
  monitor_object();
  ~monitor_object();

 private:
  void enter();
  void leave();
  bool has_lock() const;

 private:
  CRITICAL_SECTION cs_;
  HANDLE locking_thread_;
  std::size_t locking_count_;
  std::vector<hi::thread *> wait_set_;
};
}  // end of namespace hi
