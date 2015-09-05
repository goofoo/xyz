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
  /// �X���b�h�̏������J�n����.
  /// </summary>
  void start() const;

  /// <summary>
  /// �X���b�h���I������܂ő҂�.
  /// </summary>
  void join(hi::dword time = ~0) const;

  /// <summary>
  /// �X���b�h�����s���鏈���������D
  /// </summary>
  virtual void run();

  /// <summary>
  /// �X���b�h�̏�������������x�~����.
  /// </summary>
  static void sleep(hi::dword time = 0);

  /// <summary>
  /// �ʂ̃X���b�h�ɏ�����؂�ւ���.
  /// </summary>
  static void yield();

  /// <summary>
  /// �_���I��CPU����Ԃ�.
  /// </summary>
  static hi::dword get_logial_cpu_count();

  /// <summary>
  /// �X���b�h���ǂ�CPU�Ŏ��s���邩�ݒ肷��.
  /// </summary>
  static bool reserve_cpu(hi::dword core = 0);

  /// <summary>
  /// �X���b�h�̗D��x��ݒ肷��.
  /// </summary>
  static void set_priority(int priority);

 private:
  // �C���X�^���X�ϐ�
  hi::runnable *const target_;
  HANDLE hThread_;
};

//[TODO]mutex, semaphore����������

}  // end of namespace hi
