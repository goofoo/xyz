#pragma once

// スレッド関係のユーティリティ
namespace hi {
/// <summary>
/// Compare-and-Exchange.
/// </summary>
/// <param name="_val">データへのポインタ</param>
/// <param name="_old">データの古い値</param>
/// <param name="_new">設定する新しい値</param>
/// <returns>成功したらtrue</returns>
/// <history>
/// </history>
inline bool compare_and_exchange(dword volatile *const _val, dword const _old,
                                 dword const _new) {
    __asm
    {
      xor edx,  edx                     ; 成功したどうかのフラグ
      mov ecx, _new                     ; 設定する値
      mov ebx, _val                     ; データへの参照
      mov eax, _old                     ; データの古い値
      lock cmpxchg dword ptr[ebx], ecx  ; 比較して交換(アトミック)
      jnz __compare_and_exchange_exit   ; 失敗したらジャンプ
      inc edx                           ; フラグを立てる
    __compare_and_exchange_exit:
      mov eax, edx                      ; 戻り値を設定
    }
}

/* 一度だけ実行する関数
inline void thread_once()
{
  static DWORD volatile flag; // 静的変数の初期値は0

  if ( !hi::thread_cmpxchg_dword( &flag, 0, 1 ) )
  {
    return; // すでに実行している
  }

  // doSomething
}
 */

// Sequence lock
// o Optimistic lock (楽観的なロック)
// o 任意のデータ + counter
// o 読み込みスレッドだけなら lock-free
// o 書き込みスレッドは lock は必要
//   o counter が偶数なら開放、奇数なら占有状態
//
// 読み込み
// 1. Read counter
// 2. Read data
// 3. Read counter
// 1.counter が奇数か 1.counter != 3.counter なら失敗
// data を破棄してリトライ
//
// 書き込み
// 1. Counter が偶数なら CAS 命令で +1
// 2. data を書き換え
// 3. Counter をさらに +1

// Read Copy Update (RCU)
// o 短い方向リスト
// o 書き込みの遅延を許す
// o アトミック命令が不要
// o data のコピーを作って、書き込んで、つなぎかえる

// Double-ended Queue (Deque)
// o N.Arora et al., "Thread scheduling for multiprogrammed multiprocessors,"
// SPAA 1998.
//   o OS 内部のタスクキューのための考えられた deque
//   o 片側が所有スレッド用、もう片側は他スレッド用
//   o push も pop も lock-free かつ、通常はアトミック命令も不要
//   o Sun HotSpot VM の並列 GC などで利用されている
// Owner thread <-> |||||||||| -> Other threads

// その他
// o Deque
//   o M.Micchel, "CAS-based lock-free algorithm for shared dequeues," EuroPar
//   2003.
// o Bidirectional linked list
//   o H.Sundell, "Lock-free and practical doubly linked list-based deques using
//   single-word compare-and-swap," OPODIS 2004
//     o NOBLE - alibrary of non-blocking synchronization protocols
//     http://www.cs.chalmers.se/~noble/

}  // end of namespace hi
