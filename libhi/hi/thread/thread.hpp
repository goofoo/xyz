#pragma once

#include <hi/lang.hpp>
#include <hi/thread/runnable.hpp>
#include <hi/thread/monitor_object.hpp>

namespace hi {
class thread : public hi::runnable {
  HI_DISALLOW_COPY_AND_ASSIGN(thread);

  friend class hi::monitor_object;

 public:
  thread(hi::runnable *const kicker = nullptr);
  virtual ~thread();

  /// <summary>
  /// スレッドの処理を開始する.
  /// </summary>
  void start() const;

  /// <summary>
  /// スレッドが終了するまで待つ.
  /// </summary>
  void join(hi::dword time = ~0) const;

  /// <summary>
  /// スレッドが実行する処理を書く．
  /// </summary>
  virtual void run();

  /// <summary>
  /// スレッドの処理をいったん休止する.
  /// </summary>
  static void sleep(hi::dword time = 0);

  /// <summary>
  /// 別のスレッドに処理を切り替える.
  /// </summary>
  static void yield();

  /// <summary>
  /// 論理的なCPU数を返す.
  /// </summary>
  static hi::dword get_logial_cpu_count();

  /// <summary>
  /// スレッドをどのCPUで実行するか設定する.
  /// </summary>
  static bool reserve_cpu(hi::dword core = 0);

  /// <summary>
  /// スレッドの優先度を設定する.
  /// </summary>
  static void set_priority(int priority);

 private:
  // インスタンス変数
  hi::runnable *const target_;
  HANDLE hThread_;
};

//[TODO]mutex, semaphoreを実装する

}  // end of namespace hi
