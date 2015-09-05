#pragma once

namespace hi {
/// <summary>
/// Compare-and-Exchange.
/// </summary>
/// <param name="_val">�f�[�^�ւ̃|�C���^</param>
/// <param name="_old">�f�[�^�̌Â��l</param>
/// <param name="_new">�ݒ肷��V�����l</param>
/// <returns>����������true</returns>
/// <history>
/// [kitaoka] 2007/04/13 Created
/// [kitaoka] 2007/06/10 Modify
/// </history>
inline bool dword_cmpxchg(dword volatile &_val, dword &_old, dword const _new) {
  UNREFERENCED_PARAMETER(_val);
  UNREFERENCED_PARAMETER(_old);
  UNREFERENCED_PARAMETER(_new);
  //__asm
  //{
  //  mov esi, _old                      ; �O�ɓǂݏo�������̌Â��l
  //  mov eax, dword ptr[esi]            ;
  //  mov ebx, _new                      ; �V�����ݒ肵�����l
  //  mov edi, _val                      ;
  //  lock cmpxchg dword ptr[edi], ebx   ; ��r���Č���(�A�g�~�b�N)
  //  jnz __exit_failed_dword_cmpxchg    ;
  //  mov eax, 1                         ; �����F�߂�l��ݒ�
  //  jmp __exit_succeed_dword_cmpxchg   ;
  //__exit_failed_dword_cmpxchg:         ;
  //  mov dword ptr[esi], eax            ; ���s�F�Â��l�Ƃ��Č��݂̒l��ݒ�
  //  xor eax, eax                       ;       �߂�l��ݒ�
  //__exit_succeed_dword_cmpxchg:        ;
  //}
}

inline bool qword_cmpxchg(qword volatile &_val, qword &_old,
                          qword const &_new) {
  UNREFERENCED_PARAMETER(_val);
  UNREFERENCED_PARAMETER(_old);
  UNREFERENCED_PARAMETER(_new);
  //__asm
  //{
  //  mov esi, _old                      ; �O�ɓǂݏo�������̌Â��l
  //  mov eax, dword ptr[esi  ]          ;
  //  mov edx, dword ptr[esi+4]          ;
  //  mov edi, _new                      ; �V�����ݒ肵�����l
  //  mov ebx, dword ptr[edi  ]          ;
  //  mov ecx, dword ptr[edi+4]          ;
  //  mov edi, _val                      ; NOTE: 64bit
  //  CPU�̏ꍇ�́A���̕ӂ̏������قȂ�
  //  lock cmpxchg8b qword ptr[edi]      ; ��r���Č���(�A�g�~�b�N) (64bit
  //  CPU�Ȃ�cmpxchgq)
  //  jnz __exit_failed_qword_cmpxchg    ;
  //  mov eax, 1                         ; �����F�߂�l��ݒ�
  //  jmp __exit_succeed_qword_cmpxchg   ;
  //__exit_failed_qword_cmpxchg:         ;
  //  mov dword ptr[esi  ], eax          ; ���s�F�Â��l�Ƃ��Č��݂̒l��ݒ�
  //  mov dword ptr[esi+4], edx          ;
  //  xor eax, eax                       ;       �߂�l��ݒ�
  //__exit_succeed_qword_cmpxchg:        ;
  //}
}

}  // end of namespace hi
