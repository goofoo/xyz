#pragma once

// �X���b�h�֌W�̃��[�e�B���e�B
namespace hi {
/// <summary>
/// Compare-and-Exchange.
/// </summary>
/// <param name="_val">�f�[�^�ւ̃|�C���^</param>
/// <param name="_old">�f�[�^�̌Â��l</param>
/// <param name="_new">�ݒ肷��V�����l</param>
/// <returns>����������true</returns>
/// <history>
/// </history>
inline bool compare_and_exchange(dword volatile *const _val, dword const _old,
                                 dword const _new) {
    __asm
    {
      xor edx,  edx                     ; ���������ǂ����̃t���O
      mov ecx, _new                     ; �ݒ肷��l
      mov ebx, _val                     ; �f�[�^�ւ̎Q��
      mov eax, _old                     ; �f�[�^�̌Â��l
      lock cmpxchg dword ptr[ebx], ecx  ; ��r���Č���(�A�g�~�b�N)
      jnz __compare_and_exchange_exit   ; ���s������W�����v
      inc edx                           ; �t���O�𗧂Ă�
    __compare_and_exchange_exit:
      mov eax, edx                      ; �߂�l��ݒ�
    }
}

/* ��x�������s����֐�
inline void thread_once()
{
  static DWORD volatile flag; // �ÓI�ϐ��̏����l��0

  if ( !hi::thread_cmpxchg_dword( &flag, 0, 1 ) )
  {
    return; // ���łɎ��s���Ă���
  }

  // doSomething
}
 */

// Sequence lock
// o Optimistic lock (�y�ϓI�ȃ��b�N)
// o �C�ӂ̃f�[�^ + counter
// o �ǂݍ��݃X���b�h�����Ȃ� lock-free
// o �������݃X���b�h�� lock �͕K�v
//   o counter �������Ȃ�J���A��Ȃ��L���
//
// �ǂݍ���
// 1. Read counter
// 2. Read data
// 3. Read counter
// 1.counter ����� 1.counter != 3.counter �Ȃ玸�s
// data ��j�����ă��g���C
//
// ��������
// 1. Counter �������Ȃ� CAS ���߂� +1
// 2. data ����������
// 3. Counter ������� +1

// Read Copy Update (RCU)
// o �Z���������X�g
// o �������݂̒x��������
// o �A�g�~�b�N���߂��s�v
// o data �̃R�s�[������āA��������ŁA�Ȃ�������

// Double-ended Queue (Deque)
// o N.Arora et al., "Thread scheduling for multiprogrammed multiprocessors,"
// SPAA 1998.
//   o OS �����̃^�X�N�L���[�̂��߂̍l����ꂽ deque
//   o �Б������L�X���b�h�p�A�����Б��͑��X���b�h�p
//   o push �� pop �� lock-free ���A�ʏ�̓A�g�~�b�N���߂��s�v
//   o Sun HotSpot VM �̕��� GC �Ȃǂŗ��p����Ă���
// Owner thread <-> |||||||||| -> Other threads

// ���̑�
// o Deque
//   o M.Micchel, "CAS-based lock-free algorithm for shared dequeues," EuroPar
//   2003.
// o Bidirectional linked list
//   o H.Sundell, "Lock-free and practical doubly linked list-based deques using
//   single-word compare-and-swap," OPODIS 2004
//     o NOBLE - alibrary of non-blocking synchronization protocols
//     http://www.cs.chalmers.se/~noble/

}  // end of namespace hi
