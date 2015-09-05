#pragma once

#include <hi/lang.hpp>

namespace hi {
/// <summary>
/// �f�B���N�g�����쐬���܂��D
/// </summary>
bool mkdir(std::tstring dirname);

/// <summary>
/// �f�B���N�g�����̃t�@�C����񋓂��܂��D
/// </summary>
bool enum_files(std::tstring const &dirname, std::tstring const &filter,
                std::vector<std::tstring> &files);

/// <summary>
/// ��΃p�X����g���q���������t�@�C�������擾���܂��D
/// �������C�p�X�̍Ōォ��ŏ��Ɍ�����s���I�h�܂ł��g���q�Ƃ��܂��D
/// </summary>
void path2name(std::tstring const &path, std::tstring &name);
}
