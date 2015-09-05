#include <hi/file/file.hpp>
#include <hi/util.hpp>

namespace hi {
bool mkdir(std::tstring dirname) {
#define make_dir()                        \
  if (-1 == ::_tmkdir(dirname.c_str())) { \
    errno_t err;                          \
    ::_get_errno(&err);                   \
    if (err != EEXIST) {                  \
      return false;                       \
    }                                     \
  }

  for (std::size_t i = 0, size = dirname.length(); i < size; ++i) {
    if ((_TEXT('\\') == dirname[i]) || (_TEXT('/') == dirname[i])) {
      std::tchar_t tmp = dirname[i];
      dirname[i] = _TEXT('\0');
      make_dir();
      dirname[i] = tmp;
    }
  }
  make_dir();

#undef make_dir
  return true;
}

bool enum_files(std::tstring const &dirname, std::tstring const &filter,
                std::vector<std::tstring> &files) {
  files.clear();

  std::tstring query(dirname);  // 検索クエリ
  std::replace(query.begin(), query.end(), _TEXT('/'), _TEXT('\\'));
  if (_TEXT('\\') != query[query.size() - 1]) {
    query += _TEXT("\\");
  }
  query += _TEXT("*");

  // ファイルの列挙
  WIN32_FIND_DATA findFileData;
  HANDLE hFindFile = ::FindFirstFile(query.c_str(), &findFileData);

  if (INVALID_HANDLE_VALUE == hFindFile) {
    return false;
  }

  try {
    do {
      if ((findFileData.dwFileAttributes &
           FILE_ATTRIBUTE_DIRECTORY)  // ディレクトリは無視
          ||
          !hi::ends_with(findFileData.cFileName,
                         filter))  // 拡張子でフィルタリング
      {
        continue;
      }

      files.push_back(findFileData.cFileName);
    } while (::FindNextFile(hFindFile, &findFileData));
  } catch (...) {
    ::FindClose(hFindFile);
    throw;
  }
  ::FindClose(
      hFindFile);  //[TODO]ハンドル管理クラスを作って自動解放するようにする

  return true;
}

void path2name(std::tstring const &path, std::tstring &name) {
  name = path;

  std::replace(name.begin(), name.end(), _TEXT('/'), _TEXT('\\'));

  // 拡張子を削る
  {
    std::size_t const i = name.find_last_of(_TEXT('.'));
    if (std::tstring::npos != i) {
      name.resize(i);
    }
  }

  // パスを削る
  {
    std::size_t const i = name.find_last_of(_TEXT('\\'));
    if (std::tstring::npos != i) {
      name.erase(0, i + 1);
    }
  }
}
}
