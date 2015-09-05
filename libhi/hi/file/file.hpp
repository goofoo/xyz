#pragma once

#include <hi/lang.hpp>

namespace hi {
/// <summary>
/// ディレクトリを作成します．
/// </summary>
bool mkdir(std::tstring dirname);

/// <summary>
/// ディレクトリ内のファイルを列挙します．
/// </summary>
bool enum_files(std::tstring const &dirname, std::tstring const &filter,
                std::vector<std::tstring> &files);

/// <summary>
/// 絶対パスから拡張子を除いたファイル名を取得します．
/// ただし，パスの最後から最初に見つかるピリオドまでを拡張子とします．
/// </summary>
void path2name(std::tstring const &path, std::tstring &name);
}
