#include <hi/thread/thread.hpp>
#include <hi/thread/synchronized.hpp>

namespace {
unsigned int __stdcall thread_main(void *thread_ptr) {
  reinterpret_cast<hi::thread *>(thread_ptr)->run();
  ::_endthreadex(0);
  return 0;
}

}  // end of unnamed namespace

namespace hi {
thread::thread(hi::runnable *const target)
    : target_(target), hThread_(nullptr) {
  // スレッド生成
  unsigned int thread_id = unsigned int();
  hThread_ = reinterpret_cast<HANDLE>(::_beginthreadex(
      nullptr, 0, thread_main, this, CREATE_SUSPENDED, &thread_id));
}

thread::~thread() {
  // スレッド破棄
  ::CloseHandle(hThread_);
}

void thread::start() const {
  // スレッド開始
  ::ResumeThread(hThread_);
}

void thread::join(hi::dword time) const {
  // スレッドが終了するまで待つ
  ::WaitForSingleObject(hThread_, time);
}

void thread::run() {
  if (target_) {
    target_->run();
  }
}

void thread::sleep(hi::dword const time) { ::Sleep(time); }

void thread::yield() { ::SwitchToThread(); }

hi::dword thread::get_logial_cpu_count() {
  DWORD_PTR process_affinity_mask;  // プロセスで利用しているCPU
  DWORD_PTR system_affinity_mask;   // システムで利用できるCPU

  ::GetProcessAffinityMask(::GetCurrentProcess(), &process_affinity_mask,
                           &system_affinity_mask);

  hi::dword count = 0;
  for (int i = 0; i < 32; ++i) {
    if ((system_affinity_mask >> i) & 1) {
      ++count;
    }
  }
  return count;
}

bool thread::reserve_cpu(hi::dword const core) {
  return 0 != ::SetThreadAffinityMask(::GetCurrentThread(), 1 << core);
}

void thread::set_priority(int const priority) {
  ::SetThreadPriority(::GetCurrentThread(), priority);
}
}  // end of namespace hi
