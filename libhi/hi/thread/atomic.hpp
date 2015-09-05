#pragma once

namespace hi {
/// <summary>
/// Compare-and-Exchange.
/// </summary>
/// <param name="_val">データへのポインタ</param>
/// <param name="_old">データの古い値</param>
/// <param name="_new">設定する新しい値</param>
/// <returns>成功したらtrue</returns>
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
  //  mov esi, _old                      ; 前に読み出した時の古い値
  //  mov eax, dword ptr[esi]            ;
  //  mov ebx, _new                      ; 新しく設定したい値
  //  mov edi, _val                      ;
  //  lock cmpxchg dword ptr[edi], ebx   ; 比較して交換(アトミック)
  //  jnz __exit_failed_dword_cmpxchg    ;
  //  mov eax, 1                         ; 成功：戻り値を設定
  //  jmp __exit_succeed_dword_cmpxchg   ;
  //__exit_failed_dword_cmpxchg:         ;
  //  mov dword ptr[esi], eax            ; 失敗：古い値として現在の値を設定
  //  xor eax, eax                       ;       戻り値を設定
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
  //  mov esi, _old                      ; 前に読み出した時の古い値
  //  mov eax, dword ptr[esi  ]          ;
  //  mov edx, dword ptr[esi+4]          ;
  //  mov edi, _new                      ; 新しく設定したい値
  //  mov ebx, dword ptr[edi  ]          ;
  //  mov ecx, dword ptr[edi+4]          ;
  //  mov edi, _val                      ; NOTE: 64bit
  //  CPUの場合は、この辺の処理が異なる
  //  lock cmpxchg8b qword ptr[edi]      ; 比較して交換(アトミック) (64bit
  //  CPUならcmpxchgq)
  //  jnz __exit_failed_qword_cmpxchg    ;
  //  mov eax, 1                         ; 成功：戻り値を設定
  //  jmp __exit_succeed_qword_cmpxchg   ;
  //__exit_failed_qword_cmpxchg:         ;
  //  mov dword ptr[esi  ], eax          ; 失敗：古い値として現在の値を設定
  //  mov dword ptr[esi+4], edx          ;
  //  xor eax, eax                       ;       戻り値を設定
  //__exit_succeed_qword_cmpxchg:        ;
  //}
}

}  // end of namespace hi
